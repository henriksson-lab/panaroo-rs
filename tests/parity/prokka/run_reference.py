#!/usr/bin/env python3
"""Run only the reference Python's process_prokka_input, for the Phase 3 checkpoint.

    run_reference.py INPUT_LIST OUT_DIR [N_CPU]
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "build", "reference"))
from panaroo.prokka import process_prokka_input

input_list, out_dir = sys.argv[1], sys.argv[2]
n_cpu = int(sys.argv[3]) if len(sys.argv) > 3 else 1
files = [l.strip() for l in open(input_list) if l.strip()]
os.makedirs(out_dir, exist_ok=True)
if not out_dir.endswith("/"):
    out_dir += "/"
process_prokka_input(files, out_dir, False, True, n_cpu, 11)
