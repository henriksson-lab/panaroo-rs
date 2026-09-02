//! `intbitset.intbitset` — a C-extension bitset of non-negative integers.
//!
//! Panaroo stores genome membership in these (`G.nodes[n]['members']`,
//! `G.edges[e]['members']`). Semantics used by the pipeline:
//!
//! ```text
//! intbitset([i])      construct from an iterable
//! a | b               union                a |= b
//! a & b               intersection
//! a.add(i)            a.discard(i)
//! i in a              len(a)
//! a.copy()
//! a.intersection(b)
//! for i in a          ascending order
//! ```
//!
//! # Provenance
//!
//! Reimplements the observable behaviour of [intbitset](https://github.com/inveniosoftware-contrib/intbitset)
//! 3.1.0 — **LGPL-3.0-or-later**. No intbitset source was consulted; the semantics
//! reproduced (what a bitset of non-negative integers does) are not specific to that
//! implementation, and this crate neither links against nor redistributes it. See
//! `NOTICE.md`.
//!
//! Iteration is ascending, which is not an implementation detail we are free to change:
//! `";".join(str(m) for m in members)` lands in the GML `genomeIDs` attribute.

#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct IntBitSet {
    words: Vec<u64>,
}

impl IntBitSet {
    /// `intbitset()`
    pub fn new() -> Self {
        IntBitSet { words: Vec::new() }
    }

    /// `intbitset(iterable)`
    pub fn from_iter_ints<I: IntoIterator<Item = usize>>(it: I) -> Self {
        let mut s = IntBitSet::new();
        for i in it {
            s.add(i);
        }
        s
    }

    fn ensure(&mut self, word: usize) {
        if self.words.len() <= word {
            self.words.resize(word + 1, 0);
        }
    }

    /// `s.add(i)`
    pub fn add(&mut self, i: usize) {
        let (w, b) = (i / 64, i % 64);
        self.ensure(w);
        self.words[w] |= 1u64 << b;
    }

    /// `s.discard(i)` — a no-op if absent.
    pub fn discard(&mut self, i: usize) {
        let (w, b) = (i / 64, i % 64);
        if w < self.words.len() {
            self.words[w] &= !(1u64 << b);
        }
    }

    /// `i in s`
    pub fn contains(&self, i: usize) -> bool {
        let (w, b) = (i / 64, i % 64);
        w < self.words.len() && (self.words[w] >> b) & 1 == 1
    }

    /// `len(s)`
    pub fn len(&self) -> usize {
        self.words.iter().map(|w| w.count_ones() as usize).sum()
    }

    pub fn is_empty(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }

    /// `a | b`
    pub fn union(&self, other: &Self) -> Self {
        let mut out = self.clone();
        out.union_update(other);
        out
    }

    /// `a |= b`
    pub fn union_update(&mut self, other: &Self) {
        if other.words.len() > self.words.len() {
            self.words.resize(other.words.len(), 0);
        }
        for (i, &w) in other.words.iter().enumerate() {
            self.words[i] |= w;
        }
    }

    /// `a & b` / `a.intersection(b)`
    pub fn intersection(&self, other: &Self) -> Self {
        let n = self.words.len().min(other.words.len());
        IntBitSet {
            words: (0..n).map(|i| self.words[i] & other.words[i]).collect(),
        }
    }

    /// `a.isdisjoint(b)`
    pub fn isdisjoint(&self, other: &Self) -> bool {
        let n = self.words.len().min(other.words.len());
        (0..n).all(|i| self.words[i] & other.words[i] == 0)
    }

    /// `s.copy()`
    pub fn copy(&self) -> Self {
        self.clone()
    }

    /// `for i in s` — ascending.
    pub fn iter(&self) -> impl Iterator<Item = usize> + '_ {
        self.words.iter().enumerate().flat_map(|(wi, &w)| {
            (0..64).filter_map(move |b| {
                if (w >> b) & 1 == 1 {
                    Some(wi * 64 + b)
                } else {
                    None
                }
            })
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn v(s: &IntBitSet) -> Vec<usize> {
        s.iter().collect()
    }

    #[test]
    fn iteration_is_ascending_regardless_of_insertion_order() {
        // python: list(intbitset([100, 3, 64, 0])) -> [0, 3, 64, 100]
        let s = IntBitSet::from_iter_ints([100, 3, 64, 0]);
        assert_eq!(v(&s), [0, 3, 64, 100]);
        assert_eq!(s.len(), 4);
    }

    #[test]
    fn set_algebra() {
        let a = IntBitSet::from_iter_ints([1, 2, 3, 70]);
        let b = IntBitSet::from_iter_ints([3, 70, 200]);
        assert_eq!(v(&a.union(&b)), [1, 2, 3, 70, 200]);
        assert_eq!(v(&a.intersection(&b)), [3, 70]);
        assert!(!a.isdisjoint(&b));
        assert!(a.isdisjoint(&IntBitSet::from_iter_ints([4, 5])));

        let mut c = a.clone();
        c.union_update(&b);
        assert_eq!(v(&c), [1, 2, 3, 70, 200]);
    }

    #[test]
    fn add_discard_contains() {
        let mut s = IntBitSet::new();
        assert!(s.is_empty());
        s.add(5);
        s.add(300);
        assert!(s.contains(5) && s.contains(300) && !s.contains(6));
        s.discard(5);
        s.discard(999); // discarding an absent element is a no-op, as in Python
        assert_eq!(v(&s), [300]);
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn copy_is_independent() {
        let a = IntBitSet::from_iter_ints([1, 2]);
        let mut b = a.copy();
        b.add(3);
        assert_eq!(v(&a), [1, 2]);
        assert_eq!(v(&b), [1, 2, 3]);
    }

    #[test]
    fn intersection_across_different_word_counts() {
        // The shorter operand bounds the result; no out-of-range indexing.
        let a = IntBitSet::from_iter_ints([1]);
        let b = IntBitSet::from_iter_ints([1, 500]);
        assert_eq!(v(&a.intersection(&b)), [1]);
        assert_eq!(v(&b.intersection(&a)), [1]);
    }
}
