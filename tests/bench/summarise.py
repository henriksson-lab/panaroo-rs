#!/usr/bin/env python3
"""Render the bench summary from a results TSV.

    summarise.py RESULTS.tsv [--markdown OUT.md] [--env ENV.txt] [--verified TEXT]

Only rows with phase=timed and exit_status=0 are aggregated; warm-up rows and failures
are reported separately. Every metric is min / median / max across repeats, because with
3 repeats a single scheduling hiccup would dominate a mean.

--markdown writes a README-ready table with python/rust ratios. Ratios are computed from
the MEDIANS, and any cell whose two min-max ranges overlap is flagged, because such a
ratio is not distinguishable from noise.
"""
import argparse, csv, os
from statistics import median

FOOTNOTE = """
**How to read the memory columns.** `time -v maxRSS` is GNU `/usr/bin/time`'s "Maximum
resident set size", i.e. `getrusage(RUSAGE_CHILDREN)`. For a pipeline like this one it is
the peak of the **single largest process in the tree, never a sum**: both implementations
shell out to `cd-hit` (and to `mafft` when `-a/--alignment` is given), and the Python side
additionally forks a `multiprocessing` pool, so at higher thread counts this column can
badly understate what the machine actually had to hold. It is reported because it is
exact and trivially reproducible.

`tree peak PSS` is the peak of the summed *proportional* set size over the root process
and every live descendant, sampled from `/proc` every 100 ms. PSS charges each shared page
to the processes mapping it in proportion, so a forked worker pool is not double-counted.
**This is the column to quote for "how much memory did the run need".** Being sampled, it
is a lower bound: a spike shorter than 100 ms is missed.

One asymmetry is real rather than an artefact: Rust parallelises with threads inside one
process, Python with forked processes. At high `-t` the Rust `maxRSS` is essentially the
whole run's footprint while the Python `maxRSS` is one worker's share. Compare the **PSS**
columns across implementations; compare `maxRSS` only within one implementation.
"""


def fmt_secs(v):
    return f"{v:.2f}"


def fmt_mb(kb):
    return f"{kb / 1024.0:.0f}"


def stat3(values, fmt):
    return f"{fmt(min(values))} / {fmt(median(values))} / {fmt(max(values))}"


def nums(rows, key):
    out = []
    for r in rows:
        if r.get(key) not in (None, "", "NA"):
            out.append(float(r[key]))
    return out


def overlap(a, b):
    """True if the two min-max ranges intersect -- i.e. the ratio is within the noise."""
    if not a or not b:
        return False
    return min(a) <= max(b) and min(b) <= max(a)


def text_table(groups, order):
    med_wall = {}
    hdr = (f"{'impl':<10} {'thr':>4} {'n':>2}  {'wall s (min/med/max)':>26}  "
           f"{'cpu s user+sys':>26}  "
           f"{'time -v maxRSS MB':>26}  {'tree peak PSS MB':>26}  "
           f"{'tree peak RSS MB':>26}  {'procs':>5}")
    print(hdr)
    print("-" * len(hdr))
    for key in order:
        rs = groups[key]
        w = nums(rs, "wall_s")
        c = nums(rs, "cpu_s")
        m = nums(rs, "time_maxrss_kb")
        p = nums(rs, "tree_peak_pss_kb")
        t = nums(rs, "tree_peak_rss_kb")
        nproc = max(nums(rs, "tree_peak_nproc") or [0])
        med_wall[key] = median(w)
        print(f"{key[0]:<10} {key[1]:>4} {len(rs):>2}  {stat3(w, fmt_secs):>26}  "
              f"{(stat3(c, fmt_secs) if c else 'NA'):>26}  "
              f"{stat3(m, fmt_mb):>26}  {(stat3(p, fmt_mb) if p else 'NA'):>26}  "
              f"{(stat3(t, fmt_mb) if t else 'NA'):>26}  {int(nproc):>5}")
    return med_wall


def sig2(x):
    """Two significant figures -- ratios from few repeats do not deserve more."""
    if x != x or x in (float("inf"), float("-inf")):
        return "n/a"
    if x >= 10:
        return f"{x:.0f}"
    if x >= 1:
        return f"{x:.1f}"
    return f"{x:.2g}"


def markdown(groups, order, rows, envfile, verified, out):
    threads = sorted({k[1] for k in order})
    dataset = rows[0]["dataset"] if rows else "?"
    n_genomes = rows[0]["n_genomes"] if rows else "?"
    mode = rows[0]["clean_mode"] if rows else "?"
    extra = rows[0]["extra_args"] if rows else ""
    n_reps = max(len(v) for v in groups.values()) if groups else 0
    single = n_reps < 2
    warmup_ran = any(r["phase"] == "warmup" for r in rows)
    warmup_text = "preceded by one discarded warm-up" if warmup_ran else "with no discarded warm-up"

    env = {}
    if envfile and os.path.exists(envfile):
        for line in open(envfile):
            parts = line.rstrip("\n").split(None, 1)
            if len(parts) == 2:
                env.setdefault(parts[0], parts[1].strip())

    L = []
    L.append("# panaroo-rs vs Python Panaroo -- benchmark\n")
    L.append(f"Dataset `{dataset}` ({n_genomes} genomes), `--clean-mode {mode}`"
             + (f", extra args `{extra}`" if extra else ", no other non-default flags")
             + f", `-t {', '.join(str(t) for t in threads)}`.\n")

    if single:
        L.append("> **These are SINGLE-RUN measurements.** One timed run per "
                 f"configuration, {warmup_text}. **No run-to-run "
                 "variance was measured**, so there is no spread to report and the "
                 "ratios below are approximate: read them as rough magnitudes, not as "
                 "established figures. A ratio near 1 (say 0.9x-1.1x) is **not** "
                 "evidence of a real difference on this evidence. Re-run with "
                 "`tests/bench/run.sh <dataset> -n 5` for figures with a measured "
                 "spread.\n")
    else:
        L.append(f"{n_reps} timed repeats per configuration, {warmup_text}; "
                 "repeats interleaved python/rust so machine drift hits both "
                 "equally. Cells show median (min-max).\n")

    if verified:
        L.append(f"Output equivalence: **{verified}**\n")

    L.append("| machine | |")
    L.append("|---|---|")
    for k, label in (("cpu_model", "CPU"),
                     ("cpu_cores_physical", "physical cores"),
                     ("cpu_cores_logical", "logical cores"),
                     ("mem_total", "RAM"), ("kernel", "kernel"),
                     ("loadavg_at_start", "load average at start (1/5/15 min)"),
                     ("rustc", "rustc"), ("python", "python"),
                     ("cd-hit", "cd-hit"), ("mafft", "mafft"),
                     ("reference_tree", "python reference")):
        if k in env:
            L.append(f"| {label} | {env[k]} |")
    L.append("")
    if "cpu_cores_physical" in env:
        L.append(f"Thread counts are chosen against the **physical** core count: this host "
                 f"is one socket, {env['cpu_cores_physical'].split()[0]} physical cores "
                 f"with hyperthreading for {env.get('cpu_cores_logical', '?')} logical "
                 f"CPUs. `-t 20` is therefore **one thread per physical core** -- it "
                 f"saturates the machine without putting two threads on one core. A "
                 f"reader who assumes 40 cores would misread it.\n")

    L.append("| config | Python wall s | Rust wall s | speedup | Python CPU s | Rust CPU s "
             "| CPU ratio | Python maxRSS MB | Rust maxRSS MB | RSS ratio "
             "| Python tree PSS MB | Rust tree PSS MB | PSS ratio |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")

    flagged = []

    def cell(vals, fmt):
        if len(vals) < 2:
            return fmt(vals[0])
        return f"{fmt(median(vals))} ({fmt(min(vals))}-{fmt(max(vals))})"

    for th in threads:
        py = groups.get(("python", th))
        rs = groups.get(("rust", th))
        if not py or not rs:
            have = "python" if py else ("rust" if rs else "neither")
            L.append(f"| {dataset}, -t {th} | " + " | ".join([f"(only {have} measured)"] * 12) + " |")
            continue
        row = [f"{dataset}, -t {th}"]
        for key, fmt, label in (("wall_s", fmt_secs, "wall time"),
                                ("cpu_s", fmt_secs, "CPU time"),
                                ("time_maxrss_kb", fmt_mb, "maxRSS"),
                                ("tree_peak_pss_kb", fmt_mb, "tree PSS")):
            a, b = nums(py, key), nums(rs, key)
            if not a or not b:
                row += ["NA", "NA", "NA"]
                continue
            mark = ""
            if not single and overlap(a, b):
                mark = " \u2020"
                flagged.append(f"`-t {th}` {label}: the python and rust min-max ranges "
                               "overlap, so this ratio is not distinguishable from noise")
            ratio = median(a) / median(b) if median(b) else float("nan")
            row += [cell(a, fmt), cell(b, fmt), f"{sig2(ratio)}x{mark}"]
        L.append("| " + " | ".join(row) + " |")

    L.append("")
    L.append("All ratios are **python / rust**, computed from the medians and quoted to two "
             "significant figures (>1 means Rust is faster / uses less memory).")
    L.append("")
    L.append("`Python wall s` vs `Python CPU s` (user+system, summed over all threads and "
             "child processes) separates two different things: wall time at `-t N` mixes "
             "single-core efficiency with how well each implementation actually keeps N "
             "cores busy, while the CPU-time ratio isolates total work done. Note both "
             "sides spend part of that CPU time inside the same `cd-hit` binary.")
    if flagged:
        L.append("")
        for f in dict.fromkeys(flagged):
            L.append(f"\u2020 {f}.")
    L.append(FOOTNOTE.strip())

    if env.get("loadavg_at_start"):
        try:
            load1 = float(env["loadavg_at_start"].split()[0].replace(",", "."))
            ncpu = float(env.get("cpu_cores_logical", "0"))
            if ncpu and load1 > 0.25 * ncpu:
                L.append(f"\n**Warning: the machine was not idle.** Load average at the "
                         f"start of the run was {env['loadavg_at_start']} on "
                         f"{env['cpu_cores_logical']} logical cores. Absolute times, and "
                         "multi-thread scaling in particular, are contaminated by the "
                         "competing load; re-run on an idle box before quoting these "
                         "figures.")
        except ValueError:
            pass

    with open(out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("results")
    ap.add_argument("--markdown")
    ap.add_argument("--env")
    ap.add_argument("--verified", default="")
    args = ap.parse_args()

    with open(args.results) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    timed = [r for r in rows if r["phase"] == "timed" and r["exit_status"] == "0"]
    failed = [r for r in rows if r["exit_status"] != "0"]

    groups = {}
    for r in timed:
        groups.setdefault((r["impl"], int(r["threads"])), []).append(r)
    order = sorted(groups, key=lambda k: (k[1], k[0]))

    med_wall = text_table(groups, order)

    threads = sorted({k[1] for k in order})
    lines = []
    for th in threads:
        py, rs = ("python", th), ("rust", th)
        if py in med_wall and rs in med_wall and med_wall[rs] > 0:
            lines.append(f"  -t {th}: rust is {med_wall[py] / med_wall[rs]:.2f}x faster "
                         f"(median wall {med_wall[py]:.2f}s vs {med_wall[rs]:.2f}s)")
    if lines:
        print()
        print("speedup (median wall, python / rust):")
        print("\n".join(lines))

    for impl in ("python", "rust"):
        ks = [k for k in med_wall if k[0] == impl]
        if len(ks) > 1:
            base = min(ks, key=lambda k: k[1])
            print()
            print(f"thread scaling, {impl} (vs -t {base[1]}):")
            for k in sorted(ks, key=lambda k: k[1]):
                print(f"  -t {k[1]}: {med_wall[base] / med_wall[k]:.2f}x")

    if failed:
        print()
        print(f"!! {len(failed)} run(s) exited non-zero:")
        for r in failed:
            print(f"   {r['impl']} -t{r['threads']} {r['phase']} rep{r['rep']} "
                  f"exit={r['exit_status']} log={r['log']}")

    if args.markdown and timed:
        out = markdown(groups, order, rows, args.env, args.verified, args.markdown)
        print()
        print(f"markdown table written to {out}")


if __name__ == "__main__":
    main()
