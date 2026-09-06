# panaroo-rs vs Python Panaroo -- benchmark

Dataset `tiny` (4 genomes), `--clean-mode strict`, extra args `-a core`, `-t 1`.

> **These are SINGLE-RUN measurements.** One timed run per configuration, with no discarded warm-up. **No run-to-run variance was measured**, so there is no spread to report and the ratios below are approximate: read them as rough magnitudes, not as established figures. A ratio near 1 (say 0.9x-1.1x) is **not** evidence of a real difference on this evidence. Re-run with `tests/bench/run.sh <dataset> -n 5` for figures with a measured spread.

Output equivalence: **identical output (tier 0/1) at all thread counts**

| machine | |
|---|---|
| CPU | Intel(R) Xeon(R) Gold 6138 CPU @ 2.00GHz |
| physical cores | 20 (threads per core: 2) |
| logical cores | 40 |
| RAM | 196439188 kB |
| kernel | Linux 6.8.0-58-generic |
| load average at start (1/5/15 min) | 1.56 1.61 1.89 |
| rustc | rustc 1.98.0 (88d9e12ae 2026-08-18) |
| python | /home/mahogny/miniconda3/envs/panaroo-parity/bin/python  Python 3.11.16 |
| cd-hit | /home/mahogny/miniconda3/envs/panaroo-parity/bin/cd-hit  CD-HIT version 4.8.1 (built on Apr 24 2025) |
| mafft | /home/mahogny/miniconda3/envs/panaroo-parity/bin/mafft  v7.526 (2024/Apr/26) |
| python reference | /home/mahogny/github/claude/panaroo-rs/tests/parity/build/reference |

Thread counts are chosen against the **physical** core count: this host is one socket, 20 physical cores with hyperthreading for 40 logical CPUs. `-t 20` is therefore **one thread per physical core** -- it saturates the machine without putting two threads on one core. A reader who assumes 40 cores would misread it.

| config | Python wall s | Rust wall s | speedup | Python CPU s | Rust CPU s | CPU ratio | Python maxRSS MB | Rust maxRSS MB | RSS ratio | Python tree PSS MB | Rust tree PSS MB | PSS ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| tiny, -t 1 | 411.29 | 83.65 | 4.9x | 486.71 | 83.45 | 5.8x | 356 | 302 | 1.2x | 352 | 300 | 1.2x |

All ratios are **python / rust**, computed from the medians and quoted to two significant figures (>1 means Rust is faster / uses less memory).

## Embedded MAFFT profiling notes

Single-thread embedded MAFFT remains the main hot path for `tiny -a core`. After flattening
`local_align` scratch storage and the hot `profile_align_imp_multimtx` `h`/`ijp` arena, a
`perf record` sample
(`/tmp/panaroo-rs-tiny-core-after-profile-flat-hijp-20260906.perf.data`, 16,028 samples)
attributed ~24% of cycles to `mafft_align::local::local_align` and ~22% to
`mafft_align::profile::profile_align_imp_multimtx`; the largest individual CD-HIT symbols
were ~4-7%.

After the later flat-`cpmx1` and local `u8`/traceback-preallocation changes, the current
profile sample (`/tmp/panaroo-rs-tiny-core-after-local-u8-prealloc-20260906.perf.data`,
16,628 samples) attributed ~26% of cycles to `mafft_align::local::local_align` and ~22%
to `mafft_align::profile::profile_align_imp_multimtx`; CD-HIT remained secondary.

On the 451-file `tiny` core-alignment MAFFT corpus
(`/tmp/panaroo-rs-tiny-core-unaligned-inputs-20260906/unaligned_gene_sequences`), the
embedded `mafft-rs` path at that point ran in `66.714606s` with `451/451` outputs byte-identical
to `/tmp/panaroo-rs-tiny-core-profile-unchecked-20260906/aligned_gene_sequences`.
The corresponding full embedded Panaroo check wrote
`/tmp/panaroo-rs-tiny-core-local-scratch-20260906` in `80.67s` wall / `79.97s` user /
`0.67s` sys, with `core_gene_alignment.aln`, `core_gene_alignment_filtered.aln`,
`gene_presence_absence.csv`, and the whole `aligned_gene_sequences/` directory
byte-identical to the prior embedded baseline.

After also flattening the hot `profile_align_imp_multimtx` `h`/`ijp` arena, the same MAFFT
corpus ran in `64.823496s` with `451/451` byte-identical outputs. The corresponding full
embedded Panaroo check wrote `/tmp/panaroo-rs-tiny-core-profile-flat-hijp-20260906` in
`80.95s` wall / `80.40s` user / `0.52s` sys, with the same checked outputs
byte-identical to the prior embedded baseline.

Flattening `profile_align_imp_multimtx`'s remaining `cpmx1` sparse profile lowered the
same MAFFT corpus again to `63.812844s` with `451/451` byte-identical outputs. The
corresponding full embedded Panaroo check wrote
`/tmp/panaroo-rs-tiny-core-profile-flat-cpmx1-20260906` in `81.63s` wall / `80.86s` user /
`0.74s` sys, with the same checked outputs byte-identical to the prior embedded baseline.
A local scratch no-clear experiment stayed byte-identical but regressed the corpus to
`66.768811s`, so it was reverted.

Changing `local_align`'s pre-mapped `seq2` alphabet buffer from `usize` to `u8` and
preallocating traceback output vectors lowered the same MAFFT corpus to `62.763034s` with
`451/451` byte-identical outputs. The corresponding full embedded Panaroo check wrote
`/tmp/panaroo-rs-tiny-core-local-u8-prealloc-20260906` in `83.01s` wall / `82.21s` user /
`0.77s` sys, with the same checked outputs byte-identical to the prior embedded baseline;
single full-run wall times in this range have been noisy, so the direct MAFFT corpus is
the cleaner speed signal for this micro-change.

Splitting `local_align`'s valid-row DP loop into real sequence columns plus the final C
boundary column removed a per-cell `j < m` branch in the hot path. The same MAFFT corpus
then ran in `60.508935s` and `60.536705s` on two passes, both with `451/451`
byte-identical outputs. The corresponding full embedded Panaroo check wrote
`/tmp/panaroo-rs-tiny-core-local-split-boundary-20260906` in `81.88s` wall / `81.08s`
user / `0.77s` sys, with `core_gene_alignment.aln`,
`core_gene_alignment_filtered.aln`, `gene_presence_absence.csv`, and the whole
`aligned_gene_sequences/` directory byte-identical to the previous embedded baseline.
A refreshed profile of the same embedded run
(`/tmp/panaroo-rs-tiny-core-local-split-boundary-20260906.perf.data`, 81,417 samples)
still put `local_align` first at `25.34%`, followed by
`profile_align_imp_multimtx` at `22.01%`; the largest embedded CD-HIT symbols were
`diag_test_aapn` `6.53%`, `WordTable::count_words` `6.18%`, and
`simd_diag_aa_kernel::<16>` `6.03%`.

Pooling the two short dense `bcarr` buffers used by `profile_align_imp_multimtx`'s initial
row/column score setup into the existing per-thread `DpScratch` removed two heap
allocations per profile alignment call without changing accumulation order. The MAFFT
corpus then ran in `60.058220s` and `60.181905s`, both with `451/451` byte-identical
outputs. The corresponding full embedded Panaroo check wrote
`/tmp/panaroo-rs-tiny-core-profile-bcarr-pool-20260906` in `79.32s` wall / `78.57s`
user / `0.72s` sys, with `core_gene_alignment.aln`,
`core_gene_alignment_filtered.aln`, `gene_presence_absence.csv`, and the whole
`aligned_gene_sequences/` directory byte-identical to the previous embedded baseline.
A refreshed integrated perf sample after this change
(`/tmp/panaroo-rs-tiny-core-profile-bcarr-pool-refreshed-20260906.perf.data`, 84,620
samples) had the same overall shape: `local_align` `26.04%`,
`profile_align_imp_multimtx` `22.73%`, `diag_test_aapn` `6.41%`,
`WordTable::count_words` `5.95%`, `simd_diag_aa_kernel::<16>` `5.94%`, and
`local_band_align` `4.23%`.

A `local_align` valid-index fast path stayed byte-identical but regressed the corpus to
`66.652177s`, so it was reverted. Reusing pooled `scarr` for profile boundary
accumulators also stayed byte-identical but did not beat the retained result
(`62.838892s`), so it was reverted too.
Avoiding the full flat-`ijp` zero-fill in `local_align` when scratch capacity was already
large enough also stayed byte-identical on the corpus, but regressed to `62.985460s`, so
it was reverted.
Adding a score-only `local_align_score` and using it for DP-based `adjust_direction`
comparisons passed local score-equivalence tests, adjust-direction tests, and the pinned
Panaroo cluster fixture, but the default corpus uses k-mer adjust-direction instead of the
DP mode and the extra code regressed the corpus to `62.505789s`, so it was reverted.

Pre-mapping `seq1` into a scratch alphabet-code buffer also stayed byte-identical but
regressed the corpus to `65.827766s`, so it was reverted.

Two later profile-allocation experiments were also rejected. Reusing the existing `scarr`
scratch for the two dense profile boundary accumulators stayed byte-identical but ran the
corpus in `64.480212s`. Preallocating the profile traceback gaptables to `n + m` also
stayed byte-identical but ran the corpus in `64.375289s`. Both were reverted.

Changing `match_calc_row_into` to accumulate each sparse profile score in a local scalar
and store `output[j]` once preserved profile tests and the pinned Panaroo cluster fixture,
and the corpus stayed `451/451` byte-identical, but total time regressed to `61.178939s`;
the edit was reverted.

Moving the disabled-by-default Rust profile diagnostic dump blocks into cold noinline
helpers also preserved corpus byte parity, but regressed the same 451-file MAFFT corpus
to `60.669690s` against the retained `60.058220s` / `60.181905s` baseline, so it was
reverted.
Caching the same diagnostic environment gates with `OnceLock` preserved `451/451` corpus
byte parity but regressed the corpus to `63.752122s`, so that experiment was reverted too.

A CD-HIT protein SIMD counter-zeroing experiment (`span.div_ceil(32) * 32` instead of
`span + 64`) produced byte-identical direct `cd-hit` output and `.clstr` on the tiny
`combined_protein_CDS.fasta`, but direct timings were indistinguishable from baseline
(`0.61s` wall each), so it was not retained.

Larger direct CD-HIT checks on the CI Panaroo intermediates give a clearer baseline for
future CD-HIT work. On
`/tmp/panaroo-py-ci-reference-20260905/combined_protein_CDS.fasta` with Panaroo's protein
flags (`-T 1 -c 0.98 -s 0.98 -aL 0.0 -AL 99999999 -aS 0.0 -AS 99999999 -M 0 -d 999 -g 1
-n 2`), bundled C++ `cd-hit` took `7.24s` wall / `7.18s` user / `0.05s` sys / `34560`
maxRSS, while `cdhit-rs` took `6.11s` wall / `6.08s` user / `0.03s` sys / `32784`
maxRSS; output FASTA and `.clstr` were byte-identical. The matching Rust perf sample
(`/tmp/cdhit-rs-ci-protein-20260906.perf.data`, 6,106 samples) attributed `86.64%` of
cycles to `WordTable::count_words`, `5.73%` to `local_band_align`, and only `0.12%` to
`diag_test_aapn`.

On `/tmp/panaroo-py-ci-reference-20260905/combined_DNA_CDS.fasta` with Panaroo's
`cd-hit-est` flags (`-T 1 -c 0.99 -s 0.0 -aL 0.0 -AL 99999999 -aS 99999999 -AS 99999999
-r 1 -M 0 -d 999 -mask NX -n 7`), bundled C++ `cd-hit-est` took `24.05s` wall / `23.95s`
user / `0.08s` sys / `77440` maxRSS, while `cdhit-rs` took `22.69s` wall / `22.60s` user /
`0.08s` sys / `87596` maxRSS; output FASTA and `.clstr` were byte-identical.
Adding the same direct-count/materialisation strategy to the non-fragment EST path kept
the output byte-identical to both the previous Rust output and bundled C++ output (FASTA
`28ab97f90b150da694c54fd2b6bc704f49593851b95ac8478e184e141258e494`, `.clstr`
`000d4cc96fe7cbf3d42dd828cbd569a0b6b35b37cfe54a25cde7486763ddccbd`) and lowered the
same run to `17.31s` and `17.18s` wall. The direct EST variant preserves generic
`count_words(est=true)`'s leading negative-word skip for windows containing `N`; it only
changes the candidate count representation before materialisation.
Removing the now-dead post-materialisation `ic.count < required_aan` checks preserved both
protein and EST output hashes. It gave only a marginal EST sample (`17.14s` wall) and
regressed the protein direct run to `4.97s` wall against the retained `4.73s`-`4.88s`
range, so the mixed tradeoff was reverted.

A `WordTable::count_words` inner-loop prefetch for future `index_mapping[ic.index]` stayed
test-clean but regressed the CI protein direct run to `7.49s` wall, so it was reverted.
Refreshing the direct Rust protein profile after the current retained CD-HIT changes wrote
`/tmp/cdhit-rs-ci-protein-refreshed-20260906.perf.data` with 6,111 samples and again put
`WordTable::count_words` first at `86.01%`, followed by `local_band_align` `5.60%`,
`check_one_aa` `2.12%`, and `WorkingBuffer::encode_words` `1.14%`. Annotating
`count_words` showed the main remaining cost inside the long row-entry loop, especially
updating existing `look_counts.items[idm - 1].count`.
A `count_words` specialization for the common `word_encodes_no[k] == 1` case passed
release tests but regressed the same direct protein run to `6.28s` wall, so it was
reverted.

A retained CD-HIT change adds a non-fragment protein `count_words_direct` path:
`index_mapping[rep]` now holds the accumulated count directly during `check_one_aa`, and
`look_counts` is materialised once after the hot scan, keeping only candidates that reached
`required_aan`. This avoids the previous
`index_mapping -> look_counts.items[idm - 1].count` dependent random write and does the
low-count skip while clearing the direct counters; the fragment and EST paths still use the
original 1-based mapping. `cargo test --release` passed. On the
same CI protein intermediate, the direct run stayed byte-identical to the retained baseline
(FASTA `11e21947454d8e471f8a57b4865ba6790aef4c51e2a18bf7a2ac2aa7b19afbb2`, `.clstr`
`274c5cd9ceaeb4e63bdcf94733d40d93530c2805ae9e2c0e7994927ce5f0786b`) and improved from
`6.11s` wall / `6.08s` user to `4.73s`, `4.84s`, and `4.88s` wall across three retained
runs. The final step was only hoisting `ic.index as usize` once per direct-loop row entry,
after the materialisation filter had measured `5.01s` and `4.94s` wall. A fresh perf sample
before the final materialisation filter
(`/tmp/cdhit-rs-ci-protein-direct-count-20260906.perf.data`, 16,254 samples) now puts
`count_words_direct` at `84.56%`, followed by `local_band_align` `7.16%`, `check_one_aa`
`2.23%`, and `WorkingBuffer::encode_words` `1.05%`.
A refreshed perf sample after the materialisation filter
(`/tmp/cdhit-rs-ci-protein-direct-filter-20260906.perf.data`, 5,397 samples, profiling
overhead run `5.47s`) put `count_words_direct` at `86.35%`, `local_band_align` at
`7.13%`, and `check_one_aa` down at `0.13%`; further direct protein wins therefore need
to attack the core counter loop rather than candidate filtering.
Removing the now-redundant leading clear loop from `count_words_direct` preserved the
same CI protein output hashes, but direct timings were `5.20s` and `5.09s` wall against
the retained `5.04s` / `5.07s` repeat range, so it was reverted.
Increasing `count_words_direct`'s prefetch distances from the retained `24/6` to `48/12`
also preserved release tests but clearly regressed the same direct protein run to `5.85s`
and `5.88s` wall, so it was reverted.
Decreasing the direct-path distances to `12/3` also preserved the same CI protein output
hashes, but regressed to `5.40s` and `5.38s` wall, so it was reverted as well.
An explicit two-phase `count_words_direct` split that separated candidate-adding words
from late existing-only words also preserved the same output hashes, but regressed to
`5.23s` and `5.39s` wall, so it was reverted.
Writing only the candidate index on first hit, then avoiding same-slot copies during
materialisation, also preserved the same output hashes, but did not beat the retained
path available at the time (`5.04s` and `5.02s` wall), so it was reverted.
Hoisting the materialisation-loop `index as usize` cast separately also stayed test-clean
but did not improve the retained direct path (`4.82s` and `4.94s` wall), so it was
reverted.
Forcing the direct-loop index cast through `u32` also stayed byte-identical but was neutral
to slower (`4.83s` and `4.94s` wall), so it was reverted.
Manually unrolling the direct row-entry loop by two was much worse: the first benchmark
run exceeded `90s` before being terminated, versus the retained `~5s` path. It was
reverted without further repeats.
Replacing the direct-loop count clamp with `ic.count.min(j1)` also preserved byte-identical
outputs but regressed to `5.08s` and `5.00s` wall, so it was reverted.

A Panaroo no-align profile (`-a core --aligner none`, `-t 1`) of
`target/profiling/panaroo` wrote
`/tmp/panaroo-rs-tiny-noalign-profiling-20260906` and captured 34,151 samples in
`/tmp/panaroo-rs-tiny-noalign-profiling-20260906.perf.data`. The hot symbols were embedded
CD-HIT, not Panaroo graph/output code: `diag_test_aapn` `15.77%`,
`simd_diag_aa_kernel::<16>` `15.04%`, `WordTable::count_words` `15.00%`,
`local_band_align` `10.88%`, `diag_test_aapn_est` `8.28%`, and
`simd_diag_kernel::<16>` `7.90%`. Panaroo-native symbols were down in noise
(`prokka::translate` `0.20%`, `generate_network` `0.01%`), so after alignment work the
next useful single-thread speed target is still embedded CD-HIT.

Re-running the same tiny no-align/core Panaroo path after `count_words_direct` wrote
`/tmp/panaroo-rs-tiny-noalign-direct-count-20260906` in `33.83s` wall / `33.67s` user /
`0.15s` sys. `diff -qr` against `/tmp/panaroo-rs-tiny-noalign-profiling-20260906` differed
only in `alignment_resume_state.json`'s timestamp; the Panaroo-generated protein CD-HIT
outputs matched exactly (FASTA `9f23d6aeeaeca15d1892229e8d803f56d1e8847776a798658f078f59da0e47b2`,
`.clstr` `f70c66f47d5c97939160888a0db01faaa71a122ce4ca925b2fbd2fdf887367d7`).
After adding the direct-path materialisation filter, the same no-align/core integration
run wrote `/tmp/panaroo-rs-tiny-noalign-direct-filter-20260906` in `35.00s` wall /
`34.81s` user / `0.18s` sys. This was slower than the previous no-align sample, so the
standalone CD-HIT timings above remain the cleaner signal for this micro-change. The
integration outputs still matched exactly apart from `alignment_resume_state.json`'s
timestamp, and the generated protein CD-HIT hashes stayed unchanged (FASTA
`9f23d6aeeaeca15d1892229e8d803f56d1e8847776a798658f078f59da0e47b2`, `.clstr`
`f70c66f47d5c97939160888a0db01faaa71a122ce4ca925b2fbd2fdf887367d7`).
After adding the direct EST path, the same no-align/core integration run wrote
`/tmp/panaroo-rs-tiny-noalign-direct-est-20260906` in `34.70s` wall / `34.53s` user /
`0.17s` sys. `diff -qr` against the previous retained no-align run again differed only in
`alignment_resume_state.json`; the generated protein CD-HIT hashes stayed unchanged
(FASTA `9f23d6aeeaeca15d1892229e8d803f56d1e8847776a798658f078f59da0e47b2`, `.clstr`
`f70c66f47d5c97939160888a0db01faaa71a122ce4ca925b2fbd2fdf887367d7`), and
`combined_DNA_CDS.fasta` hashed to
`8527b1e9223c8db26956c00487f8ec7dc724e00fc79321828ee8cf443d27ccca`.
A refreshed no-align/core perf sample after the direct EST path
(`/tmp/panaroo-rs-tiny-noalign-direct-est-20260906.perf.data`, 32,936 samples) again
matched outputs apart from the resume timestamp. The main symbols were now protein
diagonal/filter work: `diag_test_aapn` `16.09%`, `simd_diag_aa_kernel::<10>` `15.49%`,
`local_band_align` `10.98%`, `diag_test_aapn_est` `8.96%`,
`simd_diag_kernel::<10>` `8.11%`, `count_words_direct_est` `7.46%`, and
`count_words_direct` `4.96%`.
Replacing the protein SIMD diagonal kernel's `cpx += acc * w` update with exact
special cases for the only protein complexity weights (`cpx += acc` for weight 1 and
`cpx += acc + acc` for weight 2) preserved the standalone CI protein CD-HIT output hashes
(FASTA `11e21947454d8e471f8a57b4865ba6790aef4c51e2a18bf7a2ac2aa7b19afbb2`, `.clstr`
`274c5cd9ceaeb4e63bdcf94733d40d93530c2805ae9e2c0e7994927ce5f0786b`) and the no-align
Panaroo outputs still differed from the retained direct-EST run only in
`alignment_resume_state.json`. The no-align/core integration samples improved to `33.28s`
and `33.13s` wall. A fresh perf sample
(`/tmp/panaroo-rs-tiny-noalign-aa-cpx-add-20260906.perf.data`, 32,808 samples) kept the
same broad profile shape: `diag_test_aapn` `16.88%`, `simd_diag_aa_kernel::<10>` `15.26%`,
`local_band_align` `11.06%`, `diag_test_aapn_est` `9.01%`, `simd_diag_kernel::<10>`
`8.12%`, `count_words_direct_est` `7.53%`, and `count_words_direct` `4.91%`. The isolated
CI protein CD-HIT timings for this change were noisy/slower (`5.10s` and `5.08s` wall),
so the retained signal is the no-align Panaroo workload where the protein diagonal kernel
is actually a first-order cost.
Splitting `local_band_align`'s traceback loop into 454 and non-454 macro-specialised paths
passed release tests but regressed the same no-align/core Panaroo workload to `34.30s` and
`34.35s` wall against the retained `33.13s`-`33.28s` range, so it was reverted.
Reducing SIMD diagonal accumulator zeroing from the retained slack lengths (`span + 64`
for protein and `span + 128` for nucleotide) to rounded store widths also passed the SIMD
count parity tests, but regressed the same no-align/core workload to `35.46s` and `35.64s`
wall, so it was reverted.
Hoisting `local_band_align`'s per-cell `j == len2` boundary calculation into a row-local
`end_j1` value also passed release tests and kept no-align outputs identical except for
`alignment_resume_state.json`, but it was neutral/noisy (`33.13s` then `33.37s` wall
against the retained `33.13s`-`33.28s` range), so it was reverted.
Splitting the protein SIMD query entries from packed `position | code << 32` into separate
`Vec<u32>` positions and `Vec<u16>` codes passed focused SIMD parity tests, but the no-align
tiny/core workload regressed to `35.05s` wall at
`/tmp/panaroo-rs-tiny-noalign-aa-split-entry-1-20260906`, so it was reverted.
Passing the protein SIMD query groups directly to the kernel, instead of temporarily moving
`cnt/cpx` out of `SimdDiag`, also passed focused SIMD parity tests, but timed at `33.44s`
and `33.38s` wall against the retained `33.13s`-`33.28s` range, so it was reverted.
Combining the protein SIMD weight-1 and weight-2 accumulators in registers and flushing
`cnt/cpx` once per span passed focused SIMD parity tests and preserved no-align/core outputs
except for `alignment_resume_state.json`. The tiny/core no-align integration samples improved
to `33.01s` and `33.17s` wall at
`/tmp/panaroo-rs-tiny-noalign-aa-combined-acc-{1,2}-20260906`, a small repeatable win over
the retained `33.13s`-`33.28s` range.
The direct CI protein CD-HIT output hashes were unchanged:
FASTA `11e21947454d8e471f8a57b4865ba6790aef4c51e2a18bf7a2ac2aa7b19afbb2`,
`.clstr` `274c5cd9ceaeb4e63bdcf94733d40d93530c2805ae9e2c0e7994927ce5f0786b`.
Removing the now-redundant protein SIMD `cnt/cpx` zero-fill and old-value loads, storing
the freshly computed combined accumulators directly, passed focused SIMD parity tests and
improved the same no-align/core run again to `32.72s` and `32.83s` wall at
`/tmp/panaroo-rs-tiny-noalign-aa-direct-store-{1,2}-20260906`.
Direct single-thread CD-HIT comparisons on the Panaroo CI inputs now show the Rust binary
ahead of the original C++ binaries with byte-identical outputs:
protein `cd-hit` `5.01s` wall vs C++ `7.02s` (FASTA
`11e21947454d8e471f8a57b4865ba6790aef4c51e2a18bf7a2ac2aa7b19afbb2`, `.clstr`
`274c5cd9ceaeb4e63bdcf94733d40d93530c2805ae9e2c0e7994927ce5f0786b`), and nucleotide
`cd-hit-est` `17.01s` wall vs C++ `23.74s` (FASTA
`28ab97f90b150da694c54fd2b6bc704f49593851b95ac8478e184e141258e494`, `.clstr`
`000d4cc96fe7cbf3d42dd828cbd569a0b6b35b37cfe54a25cde7486763ddccbd`).
A refreshed full embedded tiny/core run after the same direct-store change wrote
`/tmp/panaroo-rs-tiny-core-direct-store-full-20260906` in `78.93s` wall / `78.21s` user /
`0.68s` sys, slightly ahead of the retained `79.32s` full-run sample. Its
`core_gene_alignment.aln`, `core_gene_alignment_filtered.aln`, `gene_presence_absence.csv`,
and all `451` files under `aligned_gene_sequences/` were byte-identical to
`/tmp/panaroo-rs-tiny-core-profile-bcarr-pool-20260906`.
Adding a separate `u16` direct-count map for ordinary CDS-length protein/EST queries passed
release tests, but direct CD-HIT timings regressed to `5.10s` for protein and `17.55s` for
EST against the retained `5.01s` / `17.01s` samples, so it was reverted.
Hoisting the direct-count "can a new candidate still reach `required_aan`?" test from
`aan_no - k + 1 < min` to a precomputed cutoff also passed release tests, but regressed the
same direct timings to `5.10s` for protein and `17.55s` for EST, so it was reverted.

A full aligned tiny/core Panaroo run after the same CD-HIT change wrote
`/tmp/panaroo-rs-tiny-core-direct-count-20260906` in `83.85s` wall / `83.08s` user /
`0.75s` sys. This is slower than the prior `79.32s` full-run sample despite the standalone
CD-HIT win, so treat it as noisy full-pipeline timing, not a regression signal. The
meaningful outputs matched the retained `bcarr` run byte-for-byte: all `451` files under
`aligned_gene_sequences/`, the protein CD-HIT FASTA and `.clstr`, and the primary graph,
CSV, reference FASTA, and core-alignment outputs. The raw root-directory `diff -qr` is noisy
because the older tree contains duplicate root-level `.aln.fas` files that this current
run did not emit.

Removing the `black_box` barrier from the protein SIMD diagonal kernel passed release tests
and preserved the known CI protein output hashes, but direct timings were neutral/noisy
(`6.08s`, then `6.18s` wall against a retained `6.11s` baseline), so it was reverted.

`Python wall s` vs `Python CPU s` (user+system, summed over all threads and child processes) separates two different things: wall time at `-t N` mixes single-core efficiency with how well each implementation actually keeps N cores busy, while the CPU-time ratio isolates total work done. Note both sides spend part of that CPU time inside the same `cd-hit` binary.
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
