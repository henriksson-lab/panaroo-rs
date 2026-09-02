#!/usr/bin/env bash
# Benchmark harness: panaroo-rs vs the patched Python reference, wall time and peak RSS.
#
#   tests/bench/run.sh [DATASET] [options] [-- EXTRA_PANAROO_ARGS...]
#
#   DATASET      ci (default) | smoke | tiny | scale, or a path to an input list
#
# Options:
#   -n N         timed repeats per configuration (default 3)
#   -T LIST      comma-separated thread counts (default 1,10)
#   -m MODE      --clean-mode, passed to BOTH sides (default strict)
#   -r DIR       reference tree (default tests/parity/build/reference)
#   -o DIR       results/scratch dir (default tests/bench/build)
#   -k           keep the per-run output directories (default: deleted after measuring)
#   -V           skip the output-equivalence verification
#   -P           python only (skip the Rust side)
#
# Requires the pinned parity env:  conda activate panaroo-parity
#
# Everything (both implementations, cd-hit, mafft) runs under THIS script's PATH, so both
# sides provably resolve the same external binaries; the resolved paths and versions are
# recorded in the env report and asserted at startup.
set -euo pipefail

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
parity_build="$repo/tests/parity/build"

ref="$parity_build/reference"
outroot="$here/build"
repeats=3
threadlist="1,10"
mode=strict
keep=0
verify=1
python_only=0

dataset="${1:-ci}"
if [[ $# -gt 0 && "$dataset" != -* ]]; then shift; else dataset=ci; fi
while getopts "n:T:m:r:o:kVP" opt; do
  case "$opt" in
    n) repeats="$OPTARG" ;;
    T) threadlist="$OPTARG" ;;
    m) mode="$OPTARG" ;;
    r) ref="$OPTARG" ;;
    o) outroot="$OPTARG" ;;
    k) keep=1 ;;
    V) verify=0 ;;
    P) python_only=1 ;;
    *) sed -n '2,20p' "$0" >&2; exit 2 ;;
  esac
done
shift $((OPTIND - 1))
[[ "${1:-}" == "--" ]] && shift
extra=("$@")

IFS=',' read -r -a threads <<< "$threadlist"

die() { echo "bench: $*" >&2; exit 1; }

# ---------------------------------------------------------------- inputs and prereqs ---
if [[ -f "$dataset" ]]; then
  input="$dataset"; name="$(basename "$(dirname "$dataset")")"
else
  input="$parity_build/data/$dataset/input.txt"; name="$dataset"
fi
[[ -f "$input" ]] || die "no input list at $input -- run tests/parity/data/fetch.sh $name"
grep -qv '^/' "$input" && die "input list $input has relative paths; regenerate with fetch.sh"
n_genomes="$(wc -l < "$input")"

[[ -x /usr/bin/time ]] || die "/usr/bin/time (GNU time) not found; install the 'time' package"
/usr/bin/time -v /bin/true 2>/dev/null >/dev/null || die "/usr/bin/time does not support -v (not GNU time?)"

[[ -d "$ref/panaroo" ]] || die "no reference tree at $ref -- run tests/parity/reference/apply.sh"
[[ -f "$ref/PARITY_REFERENCE" ]] || die "$ref is not a parity reference build (no PARITY_REFERENCE stamp)"

command -v cd-hit  >/dev/null || die "cd-hit not on PATH -- conda activate panaroo-parity"
command -v python  >/dev/null || die "python not on PATH -- conda activate panaroo-parity"
cdhit_bin="$(command -v cd-hit)"
mafft_bin="$(command -v mafft || echo 'NOT FOUND')"
python_bin="$(command -v python)"

# Capture tool versions up front. `cd-hit -h` exits non-zero and a `grep -m1` in the
# pipeline makes it die of SIGPIPE, either of which trips `set -o pipefail`; so each of
# these is guarded with `|| true` and consumes its input fully.
cdhit_ver="$( { cd-hit -h 2>&1 || true; } | awk '/CD-HIT version/{if(!v)v=$0} END{gsub(/^[ \t=]+|[ \t=]+$/,"",v); print v}' || true)"
mafft_ver="$( { mafft --version 2>&1 || true; } | awk 'NR==1' || true)"
python_ver="$(python -V 2>&1 || true)"

max_threads=0
for t in "${threads[@]}"; do (( t > max_threads )) && max_threads=$t; done
ncpu="$(nproc)"
# Physical vs logical matters for choosing -t: this host is 20 physical cores with
# hyperthreading, so -t 20 is one thread per physical core and -t 40 doubles up.
ncpu_phys="$(lscpu 2>/dev/null | awk -F: '/^Core\(s\) per socket/{gsub(/ /,"",$2);c=$2} /^Socket\(s\)/{gsub(/ /,"",$2);s=$2} END{if(c&&s) print c*s}' || true)"
threads_per_core="$(lscpu 2>/dev/null | awk -F: '/^Thread\(s\) per core/{gsub(/ /,"",$2); print $2}' || true)"
load_start="$(cut -d' ' -f1-3 /proc/loadavg)"
load_1min="${load_start%% *}"
# Nothing here can make the box idle; the least it can do is say so.
idle_note="load average at start was $load_start on $ncpu logical cores"
if awk -v l="${load_1min/,/.}" -v n="$ncpu" 'BEGIN{exit !(l > 0.25*n)}'; then
  idle_note="MACHINE WAS NOT IDLE -- $idle_note; absolute times and thread scaling are contaminated"
  echo "bench: WARNING $idle_note" >&2
else
  idle_note="machine looked idle enough -- $idle_note"
fi
(( max_threads > ncpu )) && echo "bench: WARNING requested -t $max_threads > $ncpu cores" >&2

# ------------------------------------------------------------------------- rust build ---
rust_ok=1
rsbin="$repo/target/release/panaroo"
if [[ "$python_only" == 1 ]]; then
  rust_ok=0
  rust_note="skipped (-P)"
else
  # The CLI binary is behind a Cargo feature that is OFF by default: without --features
  # cli the [[bin]] target is simply skipped and cargo exits 0 having built nothing. So
  # delete the binary first -- if the build silently skips it, the file is absent and we
  # fail loudly instead of timing a stale binary from some earlier build.
  echo "== cargo build --release --features cli =="
  rm -f "$rsbin"
  if cargo build --release --features cli --manifest-path "$repo/Cargo.toml" 2>&1 | tail -5; then
    if [[ -x "$rsbin" ]]; then
      rust_note="ok ($("$rsbin" --version 2>/dev/null | head -1 || echo 'no --version'))"
    else
      rust_ok=0
      rust_note="build reported success but $rsbin does not exist (missing --features cli?)"
    fi
  else
    rust_ok=0; rust_note="cargo build --release --features cli FAILED"
  fi
  [[ "$rust_ok" == 0 ]] && echo "bench: WARNING Rust side unavailable: $rust_note" >&2
fi

# --------------------------------------------------------------------------- scratch ---
if repo_rev="$(git -C "$repo" rev-parse --short HEAD 2>/dev/null)"; then
  git -C "$repo" diff --quiet 2>/dev/null && repo_rev="$repo_rev (clean)" || repo_rev="$repo_rev (DIRTY working tree)"
else
  repo_rev="n/a (no commit on HEAD)"
fi

stamp="$(date +%Y%m%d-%H%M%S)"
run_dir="$outroot/run-$stamp"
mkdir -p "$run_dir/logs" "$run_dir/out"
results="$run_dir/results.tsv"
envfile="$run_dir/env.txt"
ln -sfn "run-$stamp" "$outroot/latest"

printf '%s\n' "run_id	dataset	n_genomes	impl	threads	rep	phase	exit_status	wall_s	user_s	sys_s	cpu_s	time_maxrss_kb	tree_peak_rss_kb	tree_peak_pss_kb	tree_peak_nproc	clean_mode	extra_args	log	cmdline" > "$results"

# ------------------------------------------------------------------ environment report ---
{
  echo "panaroo-rs benchmark environment"
  echo "run_id                 $stamp"
  echo "date                   $(date -Is)"
  echo "host                   $(hostname)"
  echo "kernel                 $(uname -sr)"
  echo "cpu_model              $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//')"
  echo "cpu_cores_logical      $ncpu"
  echo "cpu_cores_physical     ${ncpu_phys:-unknown} (threads per core: ${threads_per_core:-unknown})"
  echo "mem_total              $(grep -m1 MemTotal /proc/meminfo | awk '{print $2" "$3}')"
  echo "loadavg_at_start       $load_start"
  echo "governor               $(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo unknown)"
  echo
  echo "repo_rev               $repo_rev"
  echo "rustc                  $(rustc --version 2>/dev/null || echo 'n/a')"
  echo "cargo_profile          release"
  echo "rust_binary            $rsbin ($rust_note)"
  echo
  echo "reference_tree         $ref"
  sed 's/^/reference_stamp       /' "$ref/PARITY_REFERENCE"
  echo "python                 $python_bin  $python_ver"
  echo "cd-hit                 $cdhit_bin  ${cdhit_ver:-version?}"
  echo "mafft                  $mafft_bin  ${mafft_ver:-version?}"
  echo "conda_env              ${CONDA_DEFAULT_ENV:-<none>}"
  echo
  echo "dataset                $name ($n_genomes genomes, list $input)"
  echo "clean_mode             $mode"
  echo "extra_args             ${extra[*]:-<none>}"
  echo "thread_counts          ${threads[*]}"
  echo "timed_repeats          $repeats (plus 1 discarded warm-up per configuration)"
  echo
  echo "Both sides are launched from this shell, so they inherit one PATH and therefore"
  echo "resolve the SAME cd-hit / mafft binaries listed above. With no -a/--alignment in"
  echo "extra_args, mafft is never invoked by either side."
} > "$envfile"
cat "$envfile"

# ---------------------------------------------------------------- python import cost ---
# Fixed interpreter+import overhead is genuinely part of Python's cost and stays inside
# the measured runs; it is measured separately here only so the reader can see how much
# of the gap is startup rather than compute.
py_import_s="NA"
if [[ -z "${BENCH_SKIP_IMPORT:-}" ]]; then
  best=""
  for _ in 1 2 3; do
    tf="$run_dir/logs/pyimport.time"
    ( cd "$ref" && exec env PYTHONHASHSEED=0 /usr/bin/time -v -o "$tf" python -c "import panaroo" ) >/dev/null 2>&1 || true
    s="$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' "$tf")"
    s="$(python3 -c "import sys;p=sys.argv[1].split(':');print(sum(float(x)*60**i for i,x in enumerate(reversed(p))))" "$s")"
    if [[ -z "$best" ]] || python3 -c "import sys;sys.exit(0 if float(sys.argv[1])<float(sys.argv[2]) else 1)" "$s" "$best"; then best="$s"; fi
  done
  py_import_s="$best"
fi
echo
echo "python interpreter startup + 'import panaroo': ${py_import_s}s (best of 3)"

# ---------------------------------------------------------------------- measurement ---
mapfile -t infiles < "$input"

# measure IMPL THREADS REP PHASE OUTDIR
#
# Runs one configuration under /usr/bin/time -v while a /proc sampler follows the whole
# process tree, and appends one row to the TSV.
measure() {
  local impl="$1" th="$2" rep="$3" phase="$4" out="$5"
  local tag="$impl-t$th-$phase$rep"
  local tf="$run_dir/logs/$tag.time"
  local sf="$run_dir/logs/$tag.sample.json"
  local lf="$run_dir/logs/$tag.log"
  local cmdline pid st

  rm -rf "$out"; mkdir -p "$out"

  if [[ "$impl" == python ]]; then
    cmdline="(cd $ref && PYTHONHASHSEED=0 $python_bin -m panaroo -i $input -o $out --clean-mode $mode -t $th ${extra[*]:-})"
    ( cd "$ref" && exec env PYTHONHASHSEED=0 /usr/bin/time -v -o "$tf" \
        python -m panaroo -i "$input" -o "$out" --clean-mode "$mode" -t "$th" ${extra[@]+"${extra[@]}"} ) >"$lf" 2>&1 &
  else
    cmdline="$rsbin -i <$n_genomes files> -o $out --clean-mode $mode -t $th ${extra[*]:-}"
    ( exec /usr/bin/time -v -o "$tf" \
        "$rsbin" -i "${infiles[@]}" -o "$out" --clean-mode "$mode" -t "$th" ${extra[@]+"${extra[@]}"} ) >"$lf" 2>&1 &
  fi
  pid=$!
  python3 "$here/rss_sample.py" --pid "$pid" --interval 0.1 --out "$sf" &
  local spid=$!
  set +e; wait "$pid"; st=$?; set -e
  wait "$spid" || true

  local wall=NA maxrss=NA tr=NA tp=NA tn=NA user=NA sys=NA cpu=NA
  if [[ -s "$tf" ]]; then
    wall="$(awk -F': ' '/Elapsed \(wall clock\)/{print $2}' "$tf")"
    wall="$(python3 -c "import sys;p=sys.argv[1].split(':');print('%.3f'%sum(float(x)*60**i for i,x in enumerate(reversed(p))))" "${wall:-0}")"
    maxrss="$(awk -F': ' '/Maximum resident set size/{print $2}' "$tf")"
    # User+System is CPU time summed over all threads/processes: at -t N it separates
    # single-core efficiency from how well each side actually keeps N cores busy.
    user="$(awk -F': ' '/User time \(seconds\)/{print $2}' "$tf")"
    sys="$(awk -F': ' '/System time \(seconds\)/{print $2}' "$tf")"
    cpu="$(awk -v u="${user:-0}" -v s="${sys:-0}" 'BEGIN{printf "%.2f", u+s}')"
  fi
  if [[ -s "$sf" ]]; then
    read -r tr tp tn < <(python3 -c "
import json,sys
d=json.load(open(sys.argv[1]))
print(d['peak_tree_rss_kb'], d['peak_tree_pss_kb'] if d['peak_tree_pss_kb'] is not None else 'NA', d['peak_tree_nproc'])
" "$sf")
  fi

  printf '%s\n' "$stamp	$name	$n_genomes	$impl	$th	$rep	$phase	$st	$wall	$user	$sys	$cpu	$maxrss	$tr	$tp	$tn	$mode	${extra[*]:-}	$lf	$cmdline" >> "$results"
  printf '  %-6s -t%-3s %-8s  wall %8ss  cpu %9ss  time-v maxRSS %8s kB  tree PSS %8s kB  (%s procs)  exit %s\n' \
    "$impl" "$th" "$phase$rep" "$wall" "$cpu" "$maxrss" "$tp" "$tn" "$st"
  [[ "$st" == 0 ]] || echo "     ^ FAILED; log: $lf" >&2
}

impls=(python)
[[ "$rust_ok" == 1 ]] && impls+=(rust)

echo
echo "== warm-up (discarded: fills the page cache and any import cache for both sides) =="
for th in "${threads[@]}"; do
  for impl in "${impls[@]}"; do
    measure "$impl" "$th" 0 warmup "$run_dir/out/$impl-t$th-warmup"
  done
done

# --------------------------------------------------------- output equivalence check ---
# A speed number for a run that computed something different is not a speed number.
verified="not checked"
if [[ "$verify" == 1 && "$rust_ok" == 1 ]]; then
  echo
  echo "== output equivalence (warm-up outputs, tests/parity/canonicalise.py) =="
  allok=1
  for th in "${threads[@]}"; do
    echo "-- -t $th"
    if python3 "$repo/tests/parity/canonicalise.py" \
         "$run_dir/out/python-t$th-warmup" "$run_dir/out/rust-t$th-warmup" \
         > "$run_dir/logs/verify-t$th.txt" 2>&1; then
      echo "   OK: identical (tier 0/1) -- see logs/verify-t$th.txt"
    else
      echo "   MISMATCH -- see $run_dir/logs/verify-t$th.txt"
      sed -n '/DIFFERS\|MISSING/p' "$run_dir/logs/verify-t$th.txt" | sed 's/^/     /'
      allok=0
    fi
  done
  if [[ "$allok" == 1 ]]; then verified="identical output (tier 0/1) at all thread counts"
  else verified="OUTPUTS DIFFER -- timings below compare runs that computed different things"; fi
fi

echo
echo "== timed repeats (interleaved python/rust so machine drift hits both equally) =="
for rep in $(seq 1 "$repeats"); do
  for th in "${threads[@]}"; do
    for impl in "${impls[@]}"; do
      measure "$impl" "$th" "$rep" timed "$run_dir/out/$impl-t$th-r$rep"
    done
  done
done

if [[ "$keep" == 0 ]]; then rm -rf "$run_dir/out"; fi

# ----------------------------------------------------------------------- summary ---
echo
echo "==================================== SUMMARY ===================================="
echo "dataset $name ($n_genomes genomes)   clean-mode $mode   extra args: ${extra[*]:-<none>}"
echo "repeats $repeats (+1 warm-up discarded)   host $(hostname), $ncpu logical cores"
echo "cpu   $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//')"
echo "load  start $load_start / end $(cut -d' ' -f1-3 /proc/loadavg) (1/5/15 min)"
echo "idle  $idle_note"
echo "rust binary:  $rsbin -- $rust_note"
echo "verification: $verified"
echo "python startup + import panaroo: ${py_import_s}s (inside every python measurement below)"
echo
python3 "$here/summarise.py" "$results" \
  --env "$envfile" --verified "$verified" --markdown "$run_dir/results.md"
cp -f "$run_dir/results.md" "$here/results.md" 2>/dev/null || true
cat <<'CAVEAT'

--------------------------------------------------------------------------------
HOW TO READ THE MEMORY COLUMNS

  time -v maxRSS   getrusage(RUSAGE_CHILDREN) "Maximum resident set size", exactly as
                   GNU time prints it. This is the peak of the SINGLE LARGEST process in
                   the tree -- it is NOT a sum. Both implementations spawn cd-hit (and
                   mafft when -a is given), and the Python side forks a multiprocessing
                   pool, so at -t 10 this number can badly understate the machine-level
                   footprint. It is kept because it is exact and trivially reproducible.

  tree peak PSS    Peak of (sum of Pss over the root process and every live descendant),
                   sampled every 100 ms. PSS splits shared pages between the processes
                   mapping them, so a forked pool is not double-counted. THIS is the
                   number to quote for "how much memory did the run need".

  tree peak RSS    Same sampling, summing VmRSS instead. Over-counts pages shared between
                   forked workers; useful only as an upper bound.

  Both sampled columns can MISS a spike shorter than the 100 ms interval, so they are
  lower bounds. procs = the largest number of live processes seen in one sample.

  Rust parallelism is threads inside one process; Python's is forked processes. That is a
  real difference in the thing being measured, not an artefact: at -t 10 the Rust maxRSS
  is a whole-run figure while the Python maxRSS is one worker's.
--------------------------------------------------------------------------------
CAVEAT
echo
echo "results  $results"
echo "table    $here/results.md   (README-ready markdown, with python/rust ratios)"
echo "env      $envfile"
echo "logs     $run_dir/logs"
