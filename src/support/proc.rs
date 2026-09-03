//! `subprocess` — the wrapped external programs stay external.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable behaviour
//! of [CPython](https://python.org/) 3.11 `subprocess` (Python Software Foundation License
//! 2.0), for the four call shapes Panaroo uses. No CPython source was copied. See
//! `NOTICE.md`.
//!
//! Command strings are assembled by the *call sites*, byte-identically to the Python,
//! because they are printed verbatim under `--verbose` and because cd-hit's clustering
//! depends on its exact flags. This module only runs them.
//!
//! Wrapped programs: `cd-hit`, `cd-hit-est`, `mafft`, `muscle`, `muscle-super5`, `prank`,
//! `clustalo`, `famsa`.
//!
//! `shell=True` means `/bin/sh -c <string>`, which is what the redirections Panaroo appends
//! (`> /dev/null`, `> out.aln`) rely on.
//!
//! # Portability
//!
//! Every call goes through [`shell`], which picks the platform's shell rather than
//! hardcoding `/bin/sh`. On Unix that is `/bin/sh -c`, byte-for-byte what CPython's
//! `shell=True` does, so behaviour there is unchanged. On Windows it is `cmd /C`, which is
//! also what CPython's `shell=True` does.
//!
//! One thing does NOT translate by itself: the call sites append the literal `> /dev/null`
//! (`cdhit.rs`), because that is what the Python appends. `cmd` has no `/dev/null` — the
//! equivalent sink is `NUL`. [`shell`] rewrites that one redirection on Windows. This is a
//! deliberate, documented divergence from the Python string: upstream Panaroo does not run
//! on Windows at all, so there is no reference behaviour to be faithful to, and the
//! alternative is a command that silently writes a file called `dev\null`.
use super::pyfmt::py_repr_bytes;
use std::process::{Command, Stdio};

/// The platform shell, prepared with `cmd` as its command string.
///
/// Unix: `/bin/sh -c <cmd>`, exactly what CPython's `shell=True` spawns.
/// Windows: `cmd /C <cmd>`, likewise — with `> /dev/null` rewritten to `> NUL`, since the
/// call sites emit the Python's literal redirection and `cmd` would otherwise create a
/// file named `dev\null`.
fn shell(cmd: &str) -> Command {
    #[cfg(windows)]
    {
        let translated = cmd.replace("> /dev/null", "> NUL");
        let mut c = Command::new("cmd");
        c.arg("/C").arg(translated);
        c
    }
    #[cfg(not(windows))]
    {
        let mut c = Command::new("/bin/sh");
        c.arg("-c").arg(cmd);
        c
    }
}

/// `subprocess.run(cmd, shell=True, check=True)` — panics on a non-zero exit, as
/// `CalledProcessError` would.
pub fn run_shell_check(cmd: &str) {
    let st = shell(cmd)
        .status()
        .unwrap_or_else(|e| panic!("failed to run {cmd:?}: {e}"));
    if !st.success() {
        panic!(
            "CalledProcessError: Command {cmd:?} returned non-zero exit status {}.",
            st.code().unwrap_or(-1)
        );
    }
}

/// `subprocess.run(cmd, shell=True)` — returns the exit status, no check.
pub fn run_shell(cmd: &str) -> i32 {
    shell(cmd)
        .status()
        .unwrap_or_else(|e| panic!("failed to run {cmd:?}: {e}"))
        .code()
        .unwrap_or(-1)
}

/// `subprocess.run(cmd, stdout=PIPE, stderr=PIPE, shell=True)`
pub struct CompletedProcess {
    pub returncode: i32,
    pub stdout: Vec<u8>,
    pub stderr: Vec<u8>,
}

pub fn run_shell_capture(cmd: &str) -> CompletedProcess {
    let o = shell(cmd)
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .unwrap_or_else(|e| panic!("failed to run {cmd:?}: {e}"));
    CompletedProcess {
        returncode: o.status.code().unwrap_or(-1),
        stdout: o.stdout,
        stderr: o.stderr,
    }
}

/// `str(subprocess.run(cmd, stdout=PIPE, shell=True))`
///
/// **This is not a convenience — the version probes genuinely parse the repr of the
/// `CompletedProcess` object rather than its stdout.** `cdhit::check_cdhit_version` does
///
/// ```python
/// p = str(subprocess.run(cdhit_exec + ' -h', stdout=subprocess.PIPE, shell=True))
/// find_ver = re.search(r'CD-HIT version \d+\.\d+', p)
/// ```
///
/// so the text searched is `CompletedProcess(args='...', returncode=N, stdout=b'...')`
/// with the stdout rendered through Python's `bytes` repr. Reproduced here so a version
/// string that happens to straddle an escape behaves identically.
///
/// Note only `stdout` appears: `stderr` is not captured, so it is absent from the repr —
/// which means a tool that prints its version banner to stderr will not be detected. That
/// is upstream behaviour, preserved.
pub fn run_shell_capture_repr(cmd: &str) -> String {
    let o = shell(cmd)
        .stdout(Stdio::piped())
        .output()
        .unwrap_or_else(|e| panic!("failed to run {cmd:?}: {e}"));
    format!(
        "CompletedProcess(args='{}', returncode={}, stdout={})",
        cmd,
        o.status.code().unwrap_or(-1),
        py_repr_bytes(&o.stdout)
    )
}

/// `str(subprocess.run(cmd, stdout=PIPE, stderr=PIPE, shell=True))`
///
/// As [`run_shell_capture_repr`], but with `stderr` captured too — so the repr carries a
/// `stderr=b'...'` field as well, and a version banner printed to stderr is visible.
/// `generate_alignments::check_aligner_install` needs this form; `cdhit::check_cdhit_version`
/// needs the stdout-only one.
pub fn run_shell_capture_repr_both(cmd: &str) -> String {
    let c = run_shell_capture(cmd);
    format!(
        "CompletedProcess(args='{}', returncode={}, stdout={}, stderr={})",
        cmd,
        c.returncode,
        py_repr_bytes(&c.stdout),
        py_repr_bytes(&c.stderr)
    )
}

/// `subprocess.run(argv)` — an argv list, no shell.
pub fn run_argv(argv: &[String]) -> CompletedProcess {
    let o = Command::new(&argv[0])
        .args(&argv[1..])
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .output()
        .unwrap_or_else(|e| panic!("failed to run {argv:?}: {e}"));
    CompletedProcess {
        returncode: o.status.code().unwrap_or(-1),
        stdout: o.stdout,
        stderr: o.stderr,
    }
}

/// `subprocess.Popen(argv, stdout=PIPE, stderr=PIPE).communicate()` — argv, no shell.
pub fn popen_argv(argv: &[String]) -> (Vec<u8>, Vec<u8>) {
    let c = run_argv(argv);
    (c.stdout, c.stderr)
}

/// `subprocess.Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE).communicate()`
pub fn popen_communicate(cmd: &str) -> (Vec<u8>, Vec<u8>) {
    let c = run_shell_capture(cmd);
    (c.stdout, c.stderr)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn capture_repr_has_the_shape_the_version_probes_regex_over() {
        // python: str(subprocess.run("echo hi", stdout=subprocess.PIPE, shell=True))
        //   == "CompletedProcess(args='echo hi', returncode=0, stdout=b'hi\\n')"
        //
        // `cmd`'s echo terminates with CRLF where `sh`'s uses LF, and it is the raw bytes
        // that go through `py_repr_bytes`. The shape under test -- the repr wrapper the
        // version probes regex over -- is the same on both.
        #[cfg(windows)]
        let expected = r"CompletedProcess(args='echo hi', returncode=0, stdout=b'hi\r\n')";
        #[cfg(not(windows))]
        let expected = r"CompletedProcess(args='echo hi', returncode=0, stdout=b'hi\n')";
        assert_eq!(run_shell_capture_repr("echo hi"), expected);
    }

    #[test]
    fn run_shell_reports_exit_status() {
        // `exit 0` rather than `true`: `true` is a Unix builtin that `cmd` does not have.
        assert_eq!(run_shell("exit 0"), 0);
        assert_eq!(run_shell("exit 3"), 3);
    }

    #[test]
    #[should_panic(expected = "CalledProcessError")]
    fn check_panics_like_called_process_error() {
        run_shell_check("exit 1");
    }
}
