#!/usr/bin/env python3
"""Sample the total memory footprint of a process TREE at fixed intervals.

Why this exists: `/usr/bin/time -v` reports "Maximum resident set size" from
getrusage(RUSAGE_CHILDREN), which is the MAXIMUM OVER individual processes, never a
sum. Panaroo (both implementations) spawns cd-hit and mafft, and the Python side
additionally forks a multiprocessing pool -- so at -t 10 the true footprint of the run
can be ~10x what /usr/bin/time prints. This sampler sums the whole live tree instead.

Two numbers are produced:

  peak_tree_rss_kb  sum of VmRSS over root + all descendants alive at the same instant.
                    OVERSTATES the truth: pages shared between processes (libc, the
                    forked parent's copy-on-write heap) are counted once per process.

  peak_tree_pss_kb  same, but summing Pss from /proc/PID/smaps_rollup. PSS divides each
                    shared page by the number of processes mapping it, so a forked pool
                    is accounted correctly. This is the honest number, and it is the one
                    to quote. It can be unavailable (older kernels, permissions), in
                    which case it is reported as null.

Both are SAMPLED, at --interval seconds. A spike shorter than the interval is missed, so
these are lower bounds on the true peak; /usr/bin/time's value is exact (for the single
largest process) and is kept as the reproducible baseline column.

Usage: rss_sample.py --pid PID [--interval 0.1] [--out FILE] [--max-seconds N]
"""
import argparse, json, os, sys, time


def read_stat_ppid(pid):
    try:
        with open(f"/proc/{pid}/stat", "rb") as fh:
            data = fh.read()
    except (OSError, IOError):
        return None
    # comm may contain spaces/parens; ppid is the field after the final ')'
    try:
        rest = data[data.rindex(b")") + 2:].split()
        return int(rest[1])
    except (ValueError, IndexError):
        return None


def read_rss_kb(pid):
    """Resident set size in kB, from statm (cheap: 2 fields, no parsing of status)."""
    try:
        with open(f"/proc/{pid}/statm", "rb") as fh:
            fields = fh.read().split()
        return int(fields[1]) * (os.sysconf("SC_PAGE_SIZE") // 1024)
    except (OSError, IOError, ValueError, IndexError):
        return 0


def read_pss_kb(pid):
    """Proportional set size in kB from smaps_rollup, or None if unavailable."""
    try:
        with open(f"/proc/{pid}/smaps_rollup", "rb") as fh:
            for line in fh:
                if line.startswith(b"Pss:"):
                    return int(line.split()[1])
    except (OSError, IOError, ValueError, IndexError):
        return None
    return None


def all_pids():
    for name in os.listdir("/proc"):
        if name.isdigit():
            yield int(name)


def tree(root):
    """PIDs of root plus every descendant currently alive."""
    children = {}
    for pid in all_pids():
        ppid = read_stat_ppid(pid)
        if ppid is not None:
            children.setdefault(ppid, []).append(pid)
    seen, stack = [], [root]
    while stack:
        pid = stack.pop()
        seen.append(pid)
        stack.extend(children.get(pid, ()))
    return seen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pid", type=int, required=True)
    ap.add_argument("--interval", type=float, default=0.1)
    ap.add_argument("--out", default="-")
    ap.add_argument("--max-seconds", type=float, default=86400.0)
    args = ap.parse_args()

    peak_rss = 0
    peak_pss = 0
    pss_ok = False
    peak_nproc = 0
    n_samples = 0
    started = time.time()

    while True:
        if not os.path.exists(f"/proc/{args.pid}"):
            break
        if time.time() - started > args.max_seconds:
            break
        pids = tree(args.pid)
        rss = 0
        pss = 0
        for pid in pids:
            rss += read_rss_kb(pid)
            v = read_pss_kb(pid)
            if v is not None:
                pss += v
                pss_ok = True
        n_samples += 1
        peak_nproc = max(peak_nproc, len(pids))
        peak_rss = max(peak_rss, rss)
        peak_pss = max(peak_pss, pss)
        time.sleep(args.interval)

    result = {
        "peak_tree_rss_kb": peak_rss,
        "peak_tree_pss_kb": peak_pss if pss_ok else None,
        "peak_tree_nproc": peak_nproc,
        "samples": n_samples,
        "interval_s": args.interval,
    }
    text = json.dumps(result)
    if args.out == "-":
        sys.stdout.write(text + "\n")
    else:
        with open(args.out, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
