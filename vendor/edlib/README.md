# Vendored: edlib

Unmodified source of the **edlib** sequence alignment library.

- Upstream: <https://github.com/Martinsos/edlib>
- Version: **v1.2.7**
- Licence: MIT, Copyright (c) 2014 Martin Šošić — see `LICENSE` in this directory
- Files: `edlib.h`, `edlib.cpp`, taken verbatim from `edlib/include/` and `edlib/src/`

```
sha256  faf76e864ff653b91b91f56a09746fa70a4c7bd2180cd17c138cc98116a64944  edlib.h
sha256  a5f4a292131d1a31eb2c36359ffbe627b09a8e4927098cd62b95166a82502d48  edlib.cpp
```

## Why vendored rather than reimplemented

This is the only third-party library the port **links** rather than reimplements. Panaroo's
Python `edlib` package wraps this same C++ code, and `run_pw`, `search_dna` and
`align_dna_cdhit` feed edit distances straight into threshold comparisons that decide
whether two genes merge. A different aligner — or a different version of this one — would
change the pangenome. Reproducing its behaviour by reimplementation is not realistic and
would not be trustworthy.

The `edlib_rs` crate wraps the same library but pulls in cmake and libclang/bindgen at build
time; vendoring the two source files and compiling them with the `cc` crate pins the version
exactly and removes those build dependencies.

Built by `build.rs` at the repository root. Bindings are hand-written in
`src/support/edlib.rs` — see that file for the API surface Panaroo actually uses.
