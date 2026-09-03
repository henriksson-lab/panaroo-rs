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
use std::sync::{Condvar, Mutex};

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

/// `Parallel(n_jobs=n_cpu)` over a list that is **consumed in order as it completes**,
/// instead of being collected in full first.
///
/// # Why this exists rather than `parallel_map` + a loop
///
/// `prokka.py::process_prokka_input` (lines 307-317) does this:
///
/// ```python
/// job_list = list(enumerate(gff_list))
/// job_list = [job_list[i:i + n_cpu] for i in range(0, len(job_list), n_cpu)]
/// for job in tqdm(job_list, disable=quiet):
///     gene_sequence_list = Parallel(n_jobs=n_cpu)(
///         delayed(get_gene_sequences)(gff, gff_no, filter_seqs, trans_table)
///         for gff_no, gff in job)
///     for i, gene_seq in enumerate(gene_sequence_list):
///         output_files(gene_seq[0], gene_seq[1], ...)
/// ```
///
/// i.e. it chunks the job list into slices of `n_cpu`, runs a fresh pool per slice, and
/// writes that slice's results out **before starting the next slice**. Two things about
/// that shape are observable, and this helper preserves both:
///
///  1. **the sink sees results in input order** — inside a batch joblib returns them in
///     order, and the batches themselves are in order; and
///  2. **at most `n_cpu` results are alive at once** — a batch is written out and dropped
///     before the next is computed.
///
/// Property 2 is the one `parallel_map` loses. Collecting every result first and only then
/// looping over them holds all `N` results simultaneously, where CPython holds `n_cpu`. For
/// `process_prokka_input` each result is ~9-10 MB (one genome's `SeqRecord`s), so at
/// `N = 500` that is ~5 GB against CPython's ~80 MB — a regression introduced by the port,
/// not inherited from it.
///
/// What is deliberately **not** reproduced is the barrier between batches: joblib tears down
/// and respawns a pool per slice, which is PORTING_PLAN.md §9 item 9's complaint. Here one
/// pool runs throughout and a worker may start item `i + 1` while item `i` is still being
/// written. That cannot change results — the workers are pure and the sink order is fixed by
/// index — it only removes idle time at the batch boundary.
///
/// # How the bound is enforced
///
/// A worker will not *begin* item `i` until `i < emitted + capacity`, where `emitted` counts
/// items already handed to the sink and `capacity` is the worker count. So the only results
/// that can be alive are those with an index in `[emitted, emitted + capacity)`, plus the one
/// currently inside `sink`: at most `capacity + 1`.
///
/// This cannot deadlock. Indices are handed out in increasing order, so when the consumer is
/// waiting for index `k` every index `< k` has been emitted, i.e. `emitted == k`; the gate
/// for index `k` is then `k < k + capacity`, which holds for any `capacity >= 1`. The worker
/// holding the index the consumer needs therefore never blocks.
pub fn parallel_map_ordered_consume<T, R, F, S>(n_cpu: i64, items: Vec<T>, f: F, mut sink: S)
where
    T: Send,
    R: Send,
    F: Fn(T) -> R + Sync + Send,
    S: FnMut(usize, R),
{
    let n = items.len();
    if n == 0 {
        return;
    }
    let workers = n_jobs_to_workers(n_cpu).min(n);
    if workers == 1 {
        for (i, t) in items.into_iter().enumerate() {
            sink(i, f(t));
        }
        return;
    }

    // Results that have been computed but not yet consumed, plus the emit cursor. One lock
    // covers both so the backpressure gate and the hand-off cannot race.
    struct State<R> {
        // `done[i]` is `Some` once item `i` has been computed and not yet taken. At most
        // `capacity` of these are `Some` at any moment -- the `Vec` spine is `N` pointers,
        // which is not the memory this helper exists to bound.
        done: Vec<Option<R>>,
        // How many items have been handed to the sink.
        emitted: usize,
    }

    let capacity = workers;
    let inputs: Vec<Mutex<Option<T>>> = items.into_iter().map(|t| Mutex::new(Some(t))).collect();
    let next = AtomicUsize::new(0);
    let state = Mutex::new(State {
        done: (0..n).map(|_| None).collect::<Vec<Option<R>>>(),
        emitted: 0,
    });
    let cv = Condvar::new();

    std::thread::scope(|scope| {
        for _ in 0..workers {
            scope.spawn(|| loop {
                let i = next.fetch_add(1, Ordering::Relaxed);
                if i >= n {
                    break;
                }
                // Backpressure: do not compute a result the consumer is not ready for.
                {
                    let mut st = state.lock().unwrap();
                    while i >= st.emitted + capacity {
                        st = cv.wait(st).unwrap();
                    }
                }
                let item = inputs[i]
                    .lock()
                    .unwrap()
                    .take()
                    .expect("each index taken once");
                let r = f(item);
                state.lock().unwrap().done[i] = Some(r);
                cv.notify_all();
            });
        }

        // The consumer runs on this thread, so `sink` need not be `Send` and can hold the
        // output handles by mutable reference.
        for k in 0..n {
            let r = {
                let mut st = state.lock().unwrap();
                while st.done[k].is_none() {
                    st = cv.wait(st).unwrap();
                }
                let r = st.done[k].take().expect("just checked");
                st.emitted = k + 1;
                r
            };
            cv.notify_all();
            sink(k, r);
        }
    });
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

    #[test]
    fn ordered_consume_sinks_in_input_order() {
        // Same shape as the parallel_map test: later items finish first, and the sink must
        // still see 0, 1, 2, ... with the matching values.
        let items: Vec<usize> = (0..64).collect();
        let mut seen: Vec<(usize, usize)> = Vec::new();
        parallel_map_ordered_consume(
            8,
            items,
            |x| {
                std::thread::sleep(std::time::Duration::from_micros((64 - x) as u64 * 20));
                x * 2
            },
            |i, r| seen.push((i, r)),
        );
        assert_eq!(seen, (0..64).map(|x| (x, x * 2)).collect::<Vec<_>>());
    }

    #[test]
    fn ordered_consume_bounds_live_results() {
        // The point of the helper: results are not all alive at once. A result increments a
        // live counter when it is built and decrements it when dropped, so the peak is the
        // number simultaneously held. The gate allows indices in [emitted, emitted+capacity)
        // plus the one inside the sink, so the bound is capacity + 1.
        use std::sync::Arc;

        struct Live {
            live: Arc<AtomicUsize>,
        }
        impl Drop for Live {
            fn drop(&mut self) {
                self.live.fetch_sub(1, Ordering::SeqCst);
            }
        }

        let workers = n_jobs_to_workers(4);
        let live = Arc::new(AtomicUsize::new(0));
        let peak = Arc::new(Mutex::new(0usize));
        let items: Vec<usize> = (0..200).collect();

        let (l, p) = (live.clone(), peak.clone());
        let mut count = 0usize;
        parallel_map_ordered_consume(
            4,
            items,
            move |_x| {
                let now = l.fetch_add(1, Ordering::SeqCst) + 1;
                let mut pk = p.lock().unwrap();
                if now > *pk {
                    *pk = now;
                }
                drop(pk);
                Live { live: l.clone() }
            },
            |_i, r| {
                count += 1;
                drop(r);
            },
        );

        assert_eq!(count, 200);
        assert_eq!(live.load(Ordering::SeqCst), 0, "every result was dropped");
        let observed = *peak.lock().unwrap();
        assert!(observed > 0);
        assert!(
            observed <= workers + 1,
            "peak live results {observed} exceeded capacity {workers} + 1"
        );
    }
}
