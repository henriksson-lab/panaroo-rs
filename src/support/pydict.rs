//! CPython `dict`: insertion-ordered mapping.
//!
//! No Python counterpart — infrastructure. CPython has guaranteed insertion order since
//! 3.7 and preserves it across deletions (the entries array is compacted in order on
//! resize), so an `IndexMap` with `shift_remove` reproduces it exactly. Unlike
//! [`super::pyset`], no hash-table emulation is needed: only *order* is observable here,
//! and insertion order is a documented language guarantee rather than an implementation
//! accident.
//!
//! **Removal must go through [`PyDict::pop`], never a swap-remove.** A swap-remove silently
//! reorders, which would break byte parity on everything downstream of iteration order —
//! GML attribute order, `G.neighbors()` order feeding BFS, GFF3 attribute order.

use indexmap::IndexMap;
use std::hash::Hash;

#[derive(Debug, Clone)]
pub struct PyDict<K: Hash + Eq, V> {
    inner: IndexMap<K, V>,
}

impl<K: Hash + Eq, V> Default for PyDict<K, V> {
    fn default() -> Self {
        Self::new()
    }
}

impl<K: Hash + Eq, V> PyDict<K, V> {
    /// `dict()`
    pub fn new() -> Self {
        PyDict {
            inner: IndexMap::new(),
        }
    }

    /// `d[k] = v` — an existing key keeps its position, as in CPython.
    pub fn insert(&mut self, k: K, v: V) -> Option<V> {
        self.inner.insert(k, v)
    }

    /// `d[k]` / `d.get(k)`
    pub fn get(&self, k: &K) -> Option<&V> {
        self.inner.get(k)
    }

    pub fn get_mut(&mut self, k: &K) -> Option<&mut V> {
        self.inner.get_mut(k)
    }

    /// `k in d`
    pub fn contains_key(&self, k: &K) -> bool {
        self.inner.contains_key(k)
    }

    /// `del d[k]` — order-preserving, so later keys keep their relative order.
    pub fn pop(&mut self, k: &K) -> Option<V> {
        self.inner.shift_remove(k)
    }

    /// `d.setdefault(k, default)`
    pub fn setdefault(&mut self, k: K, default: V) -> &mut V {
        self.inner.entry(k).or_insert(default)
    }

    /// `defaultdict(T)[k]` — insert `T()` if absent, then return a mutable reference.
    /// Insertion position follows first access, as it does for a Python `defaultdict`.
    pub fn entry_or_default(&mut self, k: K) -> &mut V
    where
        V: Default,
    {
        self.inner.entry(k).or_default()
    }

    /// `len(d)`
    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }

    /// `d.clear()`
    pub fn clear(&mut self) {
        self.inner.clear();
    }

    /// `for k in d` — insertion order.
    pub fn keys(&self) -> impl Iterator<Item = &K> {
        self.inner.keys()
    }

    /// `d.values()` — insertion order.
    pub fn values(&self) -> impl Iterator<Item = &V> {
        self.inner.values()
    }

    pub fn values_mut(&mut self) -> impl Iterator<Item = &mut V> {
        self.inner.values_mut()
    }

    /// `d.items()` — insertion order.
    pub fn items(&self) -> impl Iterator<Item = (&K, &V)> {
        self.inner.iter()
    }
}

impl<K: Hash + Eq, V> FromIterator<(K, V)> for PyDict<K, V> {
    fn from_iter<I: IntoIterator<Item = (K, V)>>(iter: I) -> Self {
        PyDict {
            inner: IndexMap::from_iter(iter),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preserves_insertion_order() {
        let d: PyDict<&str, i32> = [("zz", 1), ("aa", 2), ("mm", 3)].into_iter().collect();
        assert_eq!(d.keys().copied().collect::<Vec<_>>(), ["zz", "aa", "mm"]);
    }

    #[test]
    fn reassignment_keeps_position() {
        // CPython: `d['zz'] = 9` does not move 'zz' to the end.
        let mut d: PyDict<&str, i32> = [("zz", 1), ("aa", 2)].into_iter().collect();
        d.insert("zz", 9);
        assert_eq!(d.keys().copied().collect::<Vec<_>>(), ["zz", "aa"]);
        assert_eq!(d.get(&"zz"), Some(&9));
    }

    #[test]
    fn deletion_preserves_order_of_the_rest() {
        let mut d: PyDict<&str, i32> = [("a", 1), ("b", 2), ("c", 3)].into_iter().collect();
        assert_eq!(d.pop(&"b"), Some(2));
        assert_eq!(d.keys().copied().collect::<Vec<_>>(), ["a", "c"]);
    }
}
