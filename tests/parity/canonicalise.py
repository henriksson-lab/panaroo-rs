#!/usr/bin/env python3
"""Tier-1 canonicalisation of Panaroo output (PORTING_PLAN.md 6.0).

Normalises away renderings that depend on Python `set` iteration order, so two runs that
produced the same pangenome compare equal even if they emitted its members in a different
sequence. It does NOT normalise away structural differences: node IDs, edges, membership,
counts and every numeric field are left untouched.

Handled:
  *.gml   `seqIDs` in two encodings -- a stringized list on one line (pre_filt_graph.gml,
          written while seqIDs is still a Python set) and a real GML list emitted as
          repeated keys (final_graph.gml, written after main converts it to a list).
          `geneIDs`/`genomeIDs` are ';'-joined strings.
  *.csv   ';'-joined tokens inside a cell; rows re-sorted, since row order depends on
          `entry_count`, which depends on component iteration order.

Usage:
    canonicalise.py FILE            -> canonical form on stdout
    canonicalise.py DIR_A DIR_B     -> compare two output dirs, exit 1 on any difference
"""
import csv, io, os, re, sys

_LIST_LINE = re.compile(r'^(\s*)(seqIDs)\s+"\[(.*)\]"\s*$')
_STR_LINE = re.compile(r'^(\s*)(geneIDs|genomeIDs)\s+"(.*)"\s*$')
_REPEATED = re.compile(r'^(\s*)(seqIDs)\s+"([^"]*)"\s*$')


def canon_gml(text):
    out, i, lines = [], 0, text.splitlines(keepends=True)
    while i < len(lines):
        line = lines[i]

        m = _LIST_LINE.match(line)
        if m:  # stringized Python set: seqIDs "['a','b']"
            ind, key, inner = m.groups()
            toks = sorted(t.strip().strip("'") for t in inner.split(",") if t.strip())
            out.append(f'{ind}{key} "[{",".join(chr(39)+t+chr(39) for t in toks)}]"\n')
            i += 1
            continue

        m = _STR_LINE.match(line)
        if m:  # ';'-joined string, optionally wrapped in the stringizer's quotes
            ind, key, val = m.groups()
            quoted = val.startswith("'") and val.endswith("'")
            body = val[1:-1] if quoted else val
            body = ";".join(sorted(body.split(";")))
            out.append(f'{ind}{key} "{chr(39)+body+chr(39) if quoted else body}"\n')
            i += 1
            continue

        m = _REPEATED.match(line)
        if m:  # real GML list: a run of repeated `seqIDs "..."` keys
            ind, key = m.group(1), m.group(2)
            vals, j = [], i
            while j < len(lines):
                mm = _REPEATED.match(lines[j])
                if not mm or mm.group(2) != key:
                    break
                vals.append(mm.group(3))
                j += 1
            for v in sorted(vals):
                out.append(f'{ind}{key} "{v}"\n')
            i = j
            continue

        out.append(line)
        i += 1
    return "".join(out)


def canon_csv(text):
    rows = list(csv.reader(io.StringIO(text)))
    if not rows:
        return text
    body = [[";".join(sorted(c.split(";"))) if ";" in c else c for c in r] for r in rows[1:]]
    body.sort()
    buf = io.StringIO()
    w = csv.writer(buf, lineterminator="\n")
    w.writerow(rows[0])
    w.writerows(body)
    return buf.getvalue()


def canon(path):
    text = open(path, newline="").read()
    if os.path.basename(path) == "alignment_resume_state.json":
        # `started_at` is a wall-clock timestamp -- it cannot match between two runs and is
        # not part of the parity contract. Everything else in the manifest is compared.
        import json
        d = json.loads(text)
        d.pop("started_at", None)
        return json.dumps(d, indent=2, sort_keys=True) + "\n"
    if path.endswith(".gml"):
        return canon_gml(text)
    if path.endswith(".csv"):
        return canon_csv(text)
    return text


def main():
    if len(sys.argv) == 2:
        sys.stdout.write(canon(sys.argv[1]))
        return 0
    if len(sys.argv) != 3:
        sys.exit(__doc__)

    a, b = sys.argv[1], sys.argv[2]

    def walk(root):
        out = set()
        for dirpath, _dirnames, filenames in os.walk(root):
            rel = os.path.relpath(dirpath, root)
            for f in filenames:
                out.add(f if rel == "." else os.path.join(rel, f))
        return out

    names = sorted(walk(a) | walk(b))
    bad = 0
    for n in names:
        pa, pb = os.path.join(a, n), os.path.join(b, n)
        if not (os.path.isfile(pa) and os.path.isfile(pb)):
            print(f"  {n:36s} MISSING on one side"); bad += 1; continue
        raw = open(pa, "rb").read() == open(pb, "rb").read()
        if raw:
            print(f"  {n:36s} identical (tier 0)")
            continue
        if canon(pa) == canon(pb):
            print(f"  {n:36s} identical after canonicalisation (tier 1)")
        else:
            print(f"  {n:36s} DIFFERS")
            bad += 1
    print()
    print("FAIL" if bad else "OK")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
