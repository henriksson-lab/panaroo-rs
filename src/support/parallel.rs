//! `joblib.Parallel` / `joblib.delayed`.
//!
//! # Provenance
//!
//! No Python counterpart in Panaroo — infrastructure. Reproduces the observable behaviour
//! of [joblib](https://joblib.readthedocs.io/) 1.4.2 (BSD 3-Clause) for the one call shape
//! Panaroo uses. No joblib source was copied. See `NOTICE.md`.
//!
//! # What has to hold
//!
//! **Result order.** `Parallel(n_jobs=n)(delayed(f)(x) for x in xs)` returns results in the
//! order of `xs`, whatever order they complete in. Everything downstream zips these against
//! the inputs, so a reordered result list corrupts output.
//!
//! **Pool shape is *not* observable.** Several call sites create a fresh pool per batch or
//! per cd-hit cluster (PORTING_PLAN.md §9 items 2 and 9). That is a performance bug, but it
//! cannot change results, so this uses one pool throughout. It is the single place where
//! "same logic" is relaxed, and it is safe precisely because joblib's contract here is
//! order-preserving and the workers are pure.
//!
//! joblib's `prefer="threads"` sites are thread-based upstream and the rest are
//! process-based; since the translated workers are pure, threads serve for both.

use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::Mutex;

/// joblib's `n_jobs` -> a concrete worker count.
///
/// A **negative** `n_jobs` is meaningful in joblib and reachable here, because Panaroo
/// passes `--threads` straight through: `-1` means all CPUs, `-2` all but one, and so on
/// (`n_cpus + 1 + n_jobs`). `0` is an error. Treating the value as unsigned instead would
/// turn `-1` into a nonsense worker count.
pub fn n_jobs_to_workers(n_jobs: i64) -> usize {
    let cpus = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1) as i64;
    let w = if n_jobs < 0 {
        cpus + 1 + n_jobs
    } else if n_jobs == 0 {
        panic!("ValueError: n_jobs == 0 in Parallel has no meaning");
    } else {
        n_jobs
    };
    w.max(1) as usize
}

/// `Parallel(n_jobs=n_cpu)(delayed(f)(x) for x in items)`.
///
/// Results come back in `items` order.
pub fn parallel_map<T, R, F>(n_cpu: i64, items: Vec<T>, f: F) -> Vec<R>
where
    T: Send,
    R: Send,
    F: Fn(T) -> R + Sync + Send,
{
    let n = items.len();
    if n == 0 {
        return Vec::new();
    }
    let workers = n_jobs_to_workers(n_cpu).min(n);
    if workers == 1 {
        return items.into_iter().map(f).collect();
    }

    // Slots are filled by index, so completion order never reaches the caller.
    let slots: Vec<Mutex<Option<R>>> = (0..n).map(|_| Mutex::new(None)).collect();
    let inputs: Vec<Mutex<Option<T>>> = items.into_iter().map(|t| Mutex::new(Some(t))).collect();
    let next = AtomicUsize::new(0);

    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= n {
                    break;
                }
                let item = inputs[i]
                    .lock()
                    .unwrap()
                    .take()
                    .expect("each index taken once");
                let r = f(item);
                *slots[i].lock().unwrap() = Some(r);
            });
        }
    });

    slots
        .into_iter()
        .map(|m| m.into_inner().unwrap().expect("every slot filled"))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn results_come_back_in_input_order() {
        // Deliberately make later items finish first; the output must still be in order.
        let items: Vec<usize> = (0..64).collect();
        let out = parallel_map(8, items, |x| {
            std::thread::sleep(std::time::Duration::from_micros((64 - x) as u64 * 20));
            x * 2
        });
        assert_eq!(out, (0..64).map(|x| x * 2).collect::<Vec<_>>());
    }

    #[test]
    fn negative_n_jobs_follows_joblib() {
        // joblib: -1 => all CPUs, -2 => all but one, positive => itself
        let cpus = std::thread::available_parallelism()
            .map(|n| n.get())
            .unwrap_or(1);
        assert_eq!(n_jobs_to_workers(-1), cpus);
        assert_eq!(n_jobs_to_workers(-2), (cpus - 1).max(1));
        assert_eq!(n_jobs_to_workers(4), 4);
        assert_eq!(n_jobs_to_workers(1), 1);
    }

    #[test]
    #[should_panic(expected = "n_jobs == 0")]
    fn zero_n_jobs_is_an_error_like_joblib() {
        n_jobs_to_workers(0);
    }

    #[test]
    fn empty_and_single_threaded_paths() {
        let e: Vec<usize> = parallel_map(4, Vec::<usize>::new(), |x| x);
        assert!(e.is_empty());
        assert_eq!(parallel_map(1, vec![1, 2, 3], |x| x + 1), [2, 3, 4]);
    }
}
