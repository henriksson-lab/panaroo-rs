//! Compile the vendored edlib C++ library.
//!
//! See `vendor/edlib/README.md` for provenance (edlib v1.2.7, MIT). Bindings are
//! hand-written in `src/support/edlib.rs`; no bindgen, so no libclang at build time.
fn main() {
    println!("cargo:rerun-if-changed=vendor/edlib/edlib.cpp");
    println!("cargo:rerun-if-changed=vendor/edlib/edlib.h");
    cc::Build::new()
        .cpp(true)
        .file("vendor/edlib/edlib.cpp")
        .include("vendor/edlib")
        .flag_if_supported("-std=c++11")
        .opt_level(3)
        .warnings(false)
        .compile("edlib");
}
