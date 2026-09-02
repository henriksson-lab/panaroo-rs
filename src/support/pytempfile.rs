//! CPython's `tempfile.mkdtemp`.
//!
//! # Provenance
//!
//! No Panaroo counterpart — this is CPython standard library behaviour that Panaroo relies
//! on. Reimplemented from the documented semantics and the reference implementation in
//! CPython 3.11 `Lib/tempfile.py` (`_RandomNameSequence`, `mkdtemp`). No source was copied.
//! See `NOTICE.md`.
//!
//! Why it matters that this is not just `create_dir("tmp_panaroo")`: `__main__.py::main`
//! does
//!
//! ```python
//! temp_dir = os.path.join(tempfile.mkdtemp(dir=args.output_dir), "")
//! os.environ['TMPDIR'] = temp_dir
//! ```
//!
//! so the working directory has a **fresh random name on every run**. Two Panaroo runs
//! writing into the same output directory therefore do not collide, and a crashed run
//! leaves its scratch behind rather than being silently reused by the next one. A fixed
//! name changes both of those.
//!
//! # Name format (CPython `_RandomNameSequence`)
//!
//! `"tmp" + 8 characters` drawn uniformly from `[a-z0-9_]` (`characters` in tempfile.py).
//! CPython retries up to `TMP_MAX` (10000) times on collision; a name that already exists
//! is not an error, just another draw.
//!
//! The directory is created with mode `0o700`, as CPython does.

use std::io::ErrorKind;
use std::path::{Path, PathBuf};

/// CPython `tempfile._RandomNameSequence.characters`.
const CHARACTERS: &[u8] = b"abcdefghijklmnopqrstuvwxyz0123456789_";

/// CPython `tempfile.TMP_MAX`.
const TMP_MAX: u32 = 10_000;

/// One 8-character random suffix.
///
/// CPython seeds a `random.Random()` per process; the values are not reproducible across
/// runs and nothing downstream may depend on them, so any decent entropy source is a
/// faithful substitute. We use `RandomState`, which is seeded from the OS.
fn random_suffix() -> String {
    use std::hash::{BuildHasher, Hasher};
    let mut h = std::collections::hash_map::RandomState::new().build_hasher();
    h.write_usize(
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .map(|d| d.subsec_nanos() as usize)
            .unwrap_or(0),
    );
    let mut bits = h.finish();
    let mut s = String::with_capacity(8);
    for _ in 0..8 {
        s.push(CHARACTERS[(bits % CHARACTERS.len() as u64) as usize] as char);
        bits /= CHARACTERS.len() as u64;
    }
    s
}

/// `tempfile.mkdtemp(dir=dir)` — creates the directory and returns its path.
///
/// Unlike Python's, this returns the error rather than raising; callers panic with the
/// message Python would print.
pub fn mkdtemp(dir: &Path) -> std::io::Result<PathBuf> {
    for _ in 0..TMP_MAX {
        let candidate = dir.join(format!("tmp{}", random_suffix()));
        match std::fs::create_dir(&candidate) {
            Ok(()) => {
                #[cfg(unix)]
                {
                    use std::os::unix::fs::PermissionsExt;
                    // CPython creates the directory 0o700.
                    let _ = std::fs::set_permissions(
                        &candidate,
                        std::fs::Permissions::from_mode(0o700),
                    );
                }
                return Ok(candidate);
            }
            // Name collision: draw again, exactly as CPython does.
            Err(e) if e.kind() == ErrorKind::AlreadyExists => continue,
            Err(e) => return Err(e),
        }
    }
    Err(std::io::Error::new(
        ErrorKind::AlreadyExists,
        "No usable temporary directory name found",
    ))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn makes_a_fresh_directory_each_time() {
        let base = std::env::temp_dir().join(format!("panaroo_rs_mkdtemp_{}", std::process::id()));
        std::fs::create_dir_all(&base).unwrap();

        let a = mkdtemp(&base).unwrap();
        let b = mkdtemp(&base).unwrap();
        assert_ne!(a, b, "mkdtemp must not reuse a name");
        assert!(a.is_dir() && b.is_dir());

        for p in [&a, &b] {
            let name = p.file_name().unwrap().to_str().unwrap();
            assert!(name.starts_with("tmp"), "{name}");
            assert_eq!(name.len(), 11, "tmp + 8 chars, got {name}");
            assert!(
                name[3..].bytes().all(|c| CHARACTERS.contains(&c)),
                "unexpected characters in {name}"
            );
        }
        std::fs::remove_dir_all(&base).unwrap();
    }
}
