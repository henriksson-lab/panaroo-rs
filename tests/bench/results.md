# panaroo-rs vs Python Panaroo -- benchmark

Dataset `ci` (4 genomes), `--clean-mode strict`, extra args `--alignment core`, `-t 20`.

> **These are SINGLE-RUN measurements.** One timed run per configuration, preceded by one discarded warm-up. **No run-to-run variance was measured**, so there is no spread to report and the ratios below are approximate: read them as rough magnitudes, not as established figures. A ratio near 1 (say 0.9x-1.1x) is **not** evidence of a real difference on this evidence. Re-run with `tests/bench/run.sh <dataset> -n 5` for figures with a measured spread.

Output equivalence: **identical output (tier 0/1) at all thread counts**

| machine | |
|---|---|
| CPU | Intel(R) Xeon(R) Gold 6138 CPU @ 2.00GHz |
| physical cores | 20 (threads per core: 2) |
| logical cores | 40 |
| RAM | 196439188 kB |
| kernel | Linux 6.8.0-58-generic |
| load average at start (1/5/15 min) | 11.32 11.14 10.74 |
| rustc | rustc 1.98.0 (88d9e12ae 2026-08-18) |
| python | /home/mahogny/miniconda3/envs/panaroo-parity/bin/python  Python 3.11.16 |
| cd-hit | /home/mahogny/miniconda3/envs/panaroo-parity/bin/cd-hit  CD-HIT version 4.8.1 (built on Apr 24 2025) |
| mafft | /home/mahogny/miniconda3/envs/panaroo-parity/bin/mafft  v7.526 (2024/Apr/26) |
| python reference | /data/henriksson/github/claude/panaroo-rs/tests/parity/build/reference |

Thread counts are chosen against the **physical** core count: this host is one socket, 20 physical cores with hyperthreading for 40 logical CPUs. `-t 20` is therefore **one thread per physical core** -- it saturates the machine without putting two threads on one core. A reader who assumes 40 cores would misread it.

| config | Python wall s | Rust wall s | speedup | Python CPU s | Rust CPU s | CPU ratio | Python maxRSS MB | Rust maxRSS MB | RSS ratio | Python tree PSS MB | Rust tree PSS MB | PSS ratio |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| ci, -t 20 | 379.10 | 209.86 | 1.8x | 5047.27 | 4115.32 | 1.2x | 388 | 389 | 1x | 2283 | 962 | 2.4x |

All ratios are **python / rust**, computed from the medians and quoted to two significant figures (>1 means Rust is faster / uses less memory).

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


**Warning: the machine was not idle.** Load average at the start of the run was 11.32 11.14 10.74 on 40 logical cores. Absolute times, and multi-thread scaling in particular, are contaminated by the competing load; re-run on an idle box before quoting these figures.
