/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
//! **UNSTABLE / EXPERIMENTAL**
//!
//! Clubcard is a compact data-structure that solves the exact membership query problem.
//!
//! Given some universe of objects U, a subset R of U, and two hash functions defined on U (as
//! described below), Clubcard outputs a compact encoding of the function `f : U -> {0, 1}` defined
//! by `f(x) = 0 if and only if x ∈ R`.
//!
//! Clubcard is based on the Ribbon filters from
//! - <https://arxiv.org/pdf/2103.02515>, and
//! - <https://arxiv.org/pdf/2109.01892>
//!
//! And some ideas from Mike Hamburg's RWC 2022 talk
//! - <https://rwc.iacr.org/2022/program.php#abstract-talk-39>
//! - <https://youtu.be/Htms5rNy7B8?list=PLeeS-3Ml-rpovBDh6do693We_CP3KTnHU&t=2357>
//!
//! The construction will be described in detail in a forthcoming paper.
//!
//! At a high level, a clubcard is a pair of GF(2) matrices (X, Y). The two hash functions
//! (h and g) map elements of U to vectors in the domain of X and Y respectively.
//!
//! The matrix X is a solution to `A * X = 0` where the rows of A are obtained by hashing the
//! elements of R with h. The number of columns in X is ~ log(|U\R| / |R|).
//!
//! The matrix Y is a solution to `B * Y = C` where the rows of B are obtained by hashing, with g,
//! the elements u ∈ U for which h(u) * X = 0. The matrix Y has one column. The rows of C encode
//! membership in R.
//!
//! Given a partition P of the set U, this library will construct the matrices A and B such that
//! they are (essentially) block diagonal with blocks indexed by p in P. The blocks are sorted by
//! rank := log(|U_p \ R_p| / |R_p|). The key observation is that we do not need to encode the
//! columns of X beyond the rank'th position in each block, as we are free to set these columns
//! equal to 0.
//!
//! Clubcard was developed to replace the use of Bloom cascades in CRLite. In a preliminary
//! experiment using a real-world collection of 12.2M revoked certs and 789.2M non-revoked certs,
//! the currently-deployed Bloom cascade implementation of CRLite produces a 19.8MB filter in 293
//! seconds (on a Ryzen 3975WX with 64GB of RAM). This Clubcard implementation produces a 9.2MB
//! filter in 190 seconds.
//!
//#![warn(missing_docs)]

extern crate bincode;

#[cfg(test)]
use rand::distributions::Distribution;
use rand::{thread_rng, Rng};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::cmp::{max, min};
use std::collections::BTreeMap;
use std::fmt;

/// An Equation\<W\> is a representation of a GF(2) linear functional
///     a(x) = b + sum_i a_i x_i
/// where a_i is equal to zero except for i in a block of 64*W coefficients
/// starting at i=s. We say an Equation is /aligned/ if a_s = 1.
/// (Note: a_i above denotes the i-th bit, not the i'th 64-bit limb.)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Equation<const W: usize> {
    s: usize,    // the row number
    a: [u64; W], // the non-trivial columns
    b: u8,       // the constant term
                 // TODO? save some space by using one bit of s for b.
}

impl<const W: usize> Equation<W> {
    /// Construct the equation a(x) = 0.
    pub fn zero() -> Self {
        Equation {
            s: 0,
            a: [0u64; W],
            b: 0,
        }
    }

    /// Construct the equation a(x) = x_i
    #[cfg(test)]
    pub fn std(i: usize) -> Self {
        let mut a = [0u64; W];
        a[0] = 1;
        Equation { s: i, a, b: 0 }
    }

    /// Construct an random aligned equation using the given distribution for s.
    #[cfg(test)]
    pub fn rand(s_dist: &impl Distribution<usize>) -> Self {
        let mut rng = rand::thread_rng();
        let s = s_dist.sample(&mut rng);
        let mut a = [0u64; W];
        for a_i in a.iter_mut() {
            *a_i = rng.gen();
        }
        a[0] |= 1;
        Equation {
            s,
            a,
            b: rng.gen::<u8>() & 1,
        }
    }

    /// Is this a(x) = 1 or a(x) = 0?
    pub fn is_zero(&self) -> bool {
        // TODO: is_const? or maybe this gets the point across.
        self.a == [0u64; W]
    }

    /// Adds `other` into `self`, i.e. sets self.a ^= other.a and self.b ^= other.b and then aligns
    /// the result.
    pub fn add(&mut self, other: &Equation<W>) {
        assert!(self.s == other.s);
        // Add the equations in GF(2)
        for i in 0..W {
            self.a[i] ^= other.a[i];
        }
        self.b ^= other.b;
        // Exit early if this equation is now zero.
        if self.is_zero() {
            return;
        }
        // Shift until there is a non-zero bit in the lowest limb.
        while self.a[0] == 0 {
            self.a.rotate_left(1);
        }
        // Shift first non-zero bit to position 0.
        let k = self.a[0].trailing_zeros();
        if k == 0 {
            return;
        }
        for i in 0..W - 1 {
            self.a[i] >>= k;
            self.a[i] |= self.a[i + 1] << (64 - k);
        }
        self.a[W - 1] >>= k;
        // Update the starting position
        self.s += k as usize;
    }

    /// Computes a(z) = sum a_i z_i.
    pub fn eval(&self, z: &[u64]) -> u8 {
        // Compute a(z), noting that this only depends
        // on 64*W bits of z starting from position s.
        let limb = self.s / 64;
        let shift = self.s % 64;
        let mut r = 0;
        for i in limb..min(z.len(), limb + W) {
            let mut tmp = z[i] >> shift;
            if i + 1 < z.len() && shift != 0 {
                tmp |= z[i + 1] << (64 - shift);
            }
            r ^= tmp & self.a[i - limb];
        }
        (r.count_ones() & 1) as u8
    }
}

/// A struct that implements Filterable can be hashed into an equation or
/// stored in a secondary retrieval mechanism.
/// TODO: document choice of W.
pub trait Filterable<const W: usize> {
    /// A good implementation of as_equation will produce equations { s, a, b } such that
    ///     (1) s is uniform in {0, 1, ..., m-1},
    ///     (2) a satisfies the alignment requirement (a\[0\] & 1 == 1) but is otherwise uniformly random,
    ///     (3) b is zero if and only if the item should be in the filter, and b = 1 otherwise.
    fn as_equation(&self, m: usize) -> Equation<W>;

    /// The shard that this item belongs in.
    fn shard(&self) -> &[u8];

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8];

    /// Whether this item should be included in an exact filter.
    fn included(&self) -> bool;
}

/// Marker type for checking that, for example, only Exact ribbons are passed to functions such as
/// Clubcard::collect_exact_ribbons.
pub struct Exact;
pub type ExactRibbon<const W: usize, T> = Ribbon<W, T, Exact>;

/// Marker type for checking that, for example, only Approximate ribbons are passed to functions such as
/// Clubcard::collect_approximate_ribbons.
pub struct Approximate;
pub type ApproximateRibbon<const W: usize, T> = Ribbon<W, T, Approximate>;

/// A RibbonBuilder collects a set of items for insertion into a Ribbon. If the optional filter is
/// provided, then only items that are contained in the filter will be inserted.
pub struct RibbonBuilder<'a, const W: usize, T: Filterable<W>> {
    /// block id.
    pub id: Vec<u8>,
    /// items to be inserted.
    pub items: Vec<T>,
    /// filter for pruning insertions.
    pub filter: Option<&'a ShardedRibbonFilter<W, T, Approximate>>,
    /// size of the universe that contains self.items
    universe_size: usize,
}

impl<'a, const W: usize, T: Filterable<W>> fmt::Display for RibbonBuilder<'a, W, T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "RibbonBuilder({:?}, {}, {})",
            self.id,
            self.items.len(),
            self.filter.is_some()
        )
    }
}

impl<'a, const W: usize, T: Filterable<W>> RibbonBuilder<'a, W, T> {
    pub fn new(
        id: &impl AsRef<[u8]>,
        filter: Option<&'a ShardedRibbonFilter<W, T, Approximate>>,
    ) -> RibbonBuilder<'a, W, T> {
        RibbonBuilder {
            id: AsRef::<[u8]>::as_ref(id).to_vec(),
            items: vec![],
            filter,
            universe_size: 0,
        }
    }

    /// Queue `item` for insertion into the ribbon (if it is contained in the provided filter).
    pub fn insert(&mut self, item: T) {
        if let Some(filter) = self.filter {
            if filter.contains(&item) {
                self.items.push(item);
            }
        } else {
            self.items.push(item);
        }
    }

    pub fn set_universe_size(&mut self, universe_size: usize) {
        self.universe_size = universe_size;
    }
}

impl<'a, const W: usize, T: Filterable<W>> From<RibbonBuilder<'a, W, T>>
    for ApproximateRibbon<W, T>
{
    /// Denote the inserted set by R and the universe by U.
    /// The ribbon returned by ApproximateRibbon::from encodes a function f : U -> {0, 1} where
    /// f(x) = 0 if and only if x is in R union S where S is a (random) subset of U \ R of size
    /// ~|R|. In other words, the ribbon solves the approximate membership query problem with a
    /// false positive rate roughly 2^-r = |R| / (|U| - |R|).
    /// The size of this ribbon is proportional to r|R|.
    fn from(mut builder: RibbonBuilder<'a, W, T>) -> ApproximateRibbon<W, T> {
        assert!(builder.items.len() <= builder.universe_size);
        let mut out =
            ApproximateRibbon::new(&builder.id, builder.items.len(), builder.universe_size);
        for item in builder.items.drain(..) {
            out.insert(item);
        }
        // Insertions should not fail for a homogeneous system.
        assert!(out.errors.is_empty());
        out
    }
}

impl<'a, const W: usize, T: Filterable<W>> From<RibbonBuilder<'a, W, T>> for ExactRibbon<W, T> {
    /// Denote the inserted set by R and the universe by U.
    /// The ribbon returned by ExactRibbon::from encodes the function "f(x) = 0 iff x in R". The
    /// size of this ribbon is proportional to |U|. In the typical use case, the set U is the
    /// result of filtering a larger universe with a false positive rate of 2^-r. This allows for
    /// exact encoding of R-membership using a pair of filters of total size ~(r+2)|R|.
    fn from(mut builder: RibbonBuilder<'a, W, T>) -> ExactRibbon<W, T> {
        assert!(builder.universe_size == 0 || builder.universe_size == builder.items.len());
        let mut out = ExactRibbon::new(&builder.id, builder.items.len());
        // By inserting the included items first, we ensure that any errors that occur during
        // insertion are for excluded items.
        let mut excluded = vec![];
        for item in builder.items.drain(..) {
            if item.included() {
                out.insert(item);
            } else {
                excluded.push(item);
            }
        }
        for item in excluded.drain(..) {
            out.insert(item);
        }
        out
    }
}

/// A compact representation of a linear system AX = B
pub struct Ribbon<const W: usize, T: Filterable<W>, ApproxOrExact> {
    /// A block identifier. Used to build an index for partitioned filters.
    id: Vec<u8>,
    /// The overhead.
    epsilon: f64,
    /// Equal to (1+epsilon) * |R|
    m: usize,
    /// The rank is round(-log2(subset_size / (universe_size - subset_size)))
    rank: usize,
    /// A linear system in which each equation has s in {0, ..., m-1}
    rows: Vec<Equation<W>>,
    /// A (typically short) list of items that failed insertion
    errors: Vec<T>,
    /// Marker for whether this is an Approximate or an Exact filter.
    phantom: std::marker::PhantomData<ApproxOrExact>,
}

impl<const W: usize, T: Filterable<W>, ApproxOrExact> fmt::Display for Ribbon<W, T, ApproxOrExact> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "ribbon({:?}): m: {}, rows: {}, rank: {}, errors: {}, epsilon: {}, overhead {}",
            self.id,
            self.m,
            self.rows.len(),
            self.rank,
            self.errors.len(),
            self.epsilon,
            (self.rows.iter().filter(|eq| eq.is_zero()).count() as f64 / (self.rows.len() as f64))
        )
    }
}

impl<const W: usize, T: Filterable<W>> ApproximateRibbon<W, T> {
    /// Construct an empty ribbon to encode a set R of size `subset_size` in a universe U of size
    /// `universe_size`.
    pub fn new(id: &impl AsRef<[u8]>, subset_size: usize, universe_size: usize) -> Self {
        assert!(subset_size <= universe_size);

        // TODO: Tune epsilon as a function of the inputs. Numerical experiments?
        let epsilon = 0.033;
        let m = ((1.0 + epsilon) * (subset_size as f64)).floor() as usize;

        let rank = if subset_size == 0 || 2 * subset_size >= universe_size {
            0
        } else {
            (((universe_size - subset_size) as f64) / (subset_size as f64))
                .log2()
                .round() as usize
        };

        Ribbon {
            id: AsRef::<[u8]>::as_ref(id).to_vec(),
            rows: vec![Equation::zero(); m],
            m,
            epsilon,
            rank,
            errors: vec![],
            phantom: std::marker::PhantomData,
        }
    }
}

impl<const W: usize, T: Filterable<W>> ExactRibbon<W, T> {
    /// Construct an empty ribbon to encode a set R of size `subset_size` in a universe U of size
    /// `universe_size`.
    pub fn new(id: &impl AsRef<[u8]>, size: usize) -> Self {
        // TODO: Tune epsilon as a function of the inputs. Numerical experiments?
        let epsilon = 0.033;
        let m = ((1.0 + epsilon) * (size as f64)).floor() as usize;

        Ribbon {
            id: AsRef::<[u8]>::as_ref(id).to_vec(),
            rows: vec![Equation::zero(); m],
            m,
            epsilon,
            rank: 0,
            errors: vec![],
            phantom: std::marker::PhantomData,
        }
    }
}

impl<const W: usize, T: Filterable<W>, ApproxOrExact> Ribbon<W, T, ApproxOrExact> {
    /// Hash the item to an Equation and insert it into the system.
    pub fn insert(&mut self, item: T) -> bool {
        let eq = item.as_equation(self.m);
        assert!(eq.is_zero() || eq.a[0] & 1 == 1);
        let rv = self.insert_equation(eq);
        if !rv {
            self.errors.push(item)
        }
        rv
    }

    /// Insert an equation into the system using Algorithm 1 from <https://arxiv.org/pdf/2103.02515>
    pub fn insert_equation(&mut self, mut eq: Equation<W>) -> bool {
        loop {
            if eq.is_zero() {
                return eq.b == 0; /* redundant (b=0) or inconsistent (b!=0) */
            }
            if eq.s >= self.rows.len() {
                // TODO: could be smarter here
                self.rows.resize_with(eq.s + 1, Equation::zero);
            }
            let cur = &mut self.rows[eq.s];
            if cur.is_zero() {
                *cur = eq;
                return true; /* inserted */
            }
            eq.add(cur);
        }
    }

    /// Solve the system using back-substitution. If this is a block in a larger system, the `tail`
    /// argument should be set to the the solution vector for the block to the right of this one.
    pub fn solve(&self, tail: &[u64]) -> Vec<u64> {
        let mut z = vec![0u64; ((self.rows.len() + 63) / 64) + tail.len()];
        // insert tail into z starting at bit self.rows.len()
        let k = self.rows.len() / 64;
        let p = self.rows.len() % 64;
        if p == 0 {
            z[k..(tail.len() + k)].copy_from_slice(tail);
        } else {
            for i in 0..tail.len() {
                z[k + i] |= tail[i] << p;
                z[k + i + 1] = tail[i] >> (64 - p)
            }
        }

        // Solve by back substitution
        for i in (0..self.rows.len()).rev() {
            let limb = i / 64;
            let pos = i % 64;
            let z_i = if self.rows[i].is_zero() {
                // Row i has a zero in column i, so we're free to choose.
                // TODO: We want multiple calls to solve() to give a different
                // solutions (when the system is suitably under-determined),
                // but it might be nice if this was deterministic.
                thread_rng().gen::<u8>()
            } else {
                // Let z' be the vector we get by setting bit i of z to z'_i.
                // Since z_i is zero, and row i has a one in column i, we have
                // row_i(z') = z'_i ^ row_i(z).
                // We want row_i(z') = b, so we must choose
                // z'_i = row_i(z) ^ b.
                self.rows[i].eval(&z) ^ self.rows[i].b
            };
            z[limb] |= ((z_i & 1) as u64) << pos;
        }
        z
    }
}

/// Helper struct for building block systems. Don't construct this
/// directly, just do `ShardedRibbonFilter::from(vec![r1, r2, r3]);`
struct ShardedRibbonFilterBuilder<const W: usize, T: Filterable<W>, ApproxOrExact> {
    blocks: Vec<Ribbon<W, T, ApproxOrExact>>,
}

impl<const W: usize, T: Filterable<W>, ApproxOrExact>
    ShardedRibbonFilterBuilder<W, T, ApproxOrExact>
{
    /// Sort ribbons by descending rank (descending simplifies indexing).
    fn sort(&mut self) {
        self.blocks.sort_unstable_by(|a, b| b.rank.cmp(&a.rank));
    }

    /// Solve the (block) system.
    /// The blocks are sorted by descending rank. We need at least one solution (i.e. column
    /// vector) per block, but we need no more than i solutions for a block of rank i.
    /// Concretely, suppose the ranks are [4, 2, 1, 0]. Then our solution can look like
    ///      block 0: | | | |
    ///      block 1: | | 0 0
    ///      block 2: | 0 0 0
    ///      block 3: | 0 0 0
    /// Since we serialize the block identifiers, offsets, and ranks in the final filter, we
    /// don't need to encode the zeros.
    fn solve(&self) -> Vec<Vec<u64>> {
        let Some(first) = self.blocks.first() else {
            return vec![];
        };
        let mut sols = vec![];
        for i in 0..max(1, first.rank) {
            // Back substitution across blocks.
            let mut tail = vec![];
            for j in (0..self.blocks.len()).rev() {
                if max(1, self.blocks[j].rank) > i {
                    tail = self.blocks[j].solve(&tail);
                }
            }
            sols.push(tail);
        }
        sols
    }

    /// Do it! Sort the blocks, solve the system, and build an index.
    pub fn finalize(&mut self) -> ShardedRibbonFilter<W, T, ApproxOrExact> {
        self.sort();
        let solution = self.solve();
        // construct the index---a map from a block identifier to that
        // block's offset in the solution vector.
        let mut index = ShardedRibbonFilterIndex::new();
        let mut offset = 0;
        for block in &self.blocks {
            let errors = block
                .errors
                .iter()
                .filter_map(|x| {
                    (!x.included()).then(|| x.discriminant().to_vec())
                })
                .collect();
            index.insert(
                block.id.clone(),
                (
                    offset,
                    block.m,
                    block.rank,
                    errors,
                ),
            );
            offset += block.rows.len();
        }
        ShardedRibbonFilter {
            index,
            solution,
            phantom: std::marker::PhantomData,
            phantom2: std::marker::PhantomData,
        }
    }
}

/// Metadata
type ShardedRibbonFilterIndex = BTreeMap<
    /* block id */ Vec<u8>,
    (
        /* offset */ usize,
        /* m */ usize,
        /* rank */ usize,
        /* exclude errors */ Vec<Vec<u8>>,
    ),
>;

/// A solution to a ribbon system, along with metadata necessary for querying it.
pub struct ShardedRibbonFilter<const W: usize, T: Filterable<W>, ApproxOrExact> {
    index: ShardedRibbonFilterIndex,
    solution: Vec<Vec<u64>>,
    phantom: std::marker::PhantomData<T>,
    phantom2: std::marker::PhantomData<ApproxOrExact>,
}

impl<const W: usize, T: Filterable<W>, ApproxOrExact> fmt::Display
    for ShardedRibbonFilter<W, T, ApproxOrExact>
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "ShardedRibbonFilter({:?})", self.index)
    }
}

impl<const W: usize, T: Filterable<W>, ApproxOrExact> ShardedRibbonFilter<W, T, ApproxOrExact> {
    /// Check if this filter contains the given item in the given block.
    pub fn contains(&self, item: &T) -> bool {
        let Some((offset, m, rank, exclude_errors)) = self.index.get(item.shard())
        else {
            return false;
        };
        // Empty blocks do not contain anything,
        // despite having inner product 0 with everything.
        if *m == 0 {
            return false;
        }
        let mut eq = item.as_equation(*m);
        eq.s += *offset;
        // eq.b is irrelevant here. We don't know it when querying.
        for i in 0..max(1, *rank) {
            if eq.eval(&self.solution[i]) != 0 {
                return false;
            }
        }
        for error in exclude_errors {
            if error == item.discriminant() {
                return false;
            }
        }
        true
    }
}

impl<const W: usize, T: Filterable<W>, ApproxOrExact> From<Vec<Ribbon<W, T, ApproxOrExact>>>
    for ShardedRibbonFilter<W, T, ApproxOrExact>
{
    fn from(blocks: Vec<Ribbon<W, T, ApproxOrExact>>) -> ShardedRibbonFilter<W, T, ApproxOrExact> {
        ShardedRibbonFilterBuilder { blocks }.finalize()
    }
}

/// A pair of ribbon filters that, together, solve the exact membership query problem.
pub struct ClubcardBuilder<const W: usize, T: Filterable<W>> {
    /// An approximate membership query filter to whittle down the universe
    /// to a managable size.
    approx_filter: Option<ShardedRibbonFilter<W, T, Approximate>>,
    /// An exact membership query filter to confirm membership in R for items that
    /// pass through the approximate filter.
    exact_filter: Option<ShardedRibbonFilter<W, T, Exact>>,
}

impl<const W: usize, T: Filterable<W>> Default for ClubcardBuilder<W, T> {
    fn default() -> Self {
        ClubcardBuilder {
            approx_filter: None,
            exact_filter: None,
        }
    }
}

impl<const W: usize, T: Filterable<W>> ClubcardBuilder<W, T> {
    pub fn new() -> Self {
        ClubcardBuilder::default()
    }

    pub fn get_approx_builder(&self, block: &impl AsRef<[u8]>) -> RibbonBuilder<'static, W, T> {
        RibbonBuilder::new(block, None)
    }

    pub fn get_exact_builder<'a>(&'a self, block: &impl AsRef<[u8]>) -> RibbonBuilder<'a, W, T> {
        RibbonBuilder::new(block, self.approx_filter.as_ref())
    }

    pub fn collect_approx_ribbons(&mut self, ribbons: Vec<ApproximateRibbon<W, T>>) {
        self.approx_filter = Some(ShardedRibbonFilter::from(ribbons));
    }

    pub fn collect_exact_ribbons(&mut self, ribbons: Vec<Ribbon<W, T, Exact>>) {
        self.exact_filter = Some(ShardedRibbonFilter::from(ribbons));
    }

    pub fn build(self) -> Clubcard {
        let mut index: ClubcardIndex = BTreeMap::new();

        assert!(self.approx_filter.is_some());
        let approx_filter = self.approx_filter.unwrap();
        for (shard, (offset, m, rank, exclude_errors)) in approx_filter.index {
            let mut meta = ClubcardShardMeta::default();
            meta.approx_filter_offset = offset;
            meta.approx_filter_m = m;
            meta.approx_filter_rank = rank;
            assert!(exclude_errors.len() == 0);
            meta.exclude_errors = exclude_errors;
            index.insert(shard, meta);
        }

        assert!(self.exact_filter.is_some());
        let mut exact_filter = self.exact_filter.unwrap();
        for (shard, (offset, m, rank, exclude_errors)) in exact_filter.index {
            assert!(rank == 0);
            let meta = index.get_mut(&shard).unwrap();
            meta.exact_filter_offset = offset;
            meta.exact_filter_m = m;
            meta.exclude_errors.extend(exclude_errors);
        }

        assert!(exact_filter.solution.len() == 1);
        let exact_filter = exact_filter.solution.pop().unwrap();

        Clubcard {
            index,
            approx_filter: approx_filter.solution,
            exact_filter,
        }
    }
}

// TODO: We should have a way to say that a shard encodes |U\R| rather than |R|,
// for the case where |R| = (1-ϵ)|U|.
#[derive(Default, Serialize, Deserialize)]
struct ClubcardShardMeta {
    approx_filter_offset: usize,
    approx_filter_m: usize,
    approx_filter_rank: usize,
    exact_filter_offset: usize,
    exact_filter_m: usize,
    exclude_errors: Vec<Vec<u8>>,
}

type ClubcardIndex = BTreeMap</* block id */ Vec<u8>, ClubcardShardMeta>;

#[derive(Serialize, Deserialize)]
pub struct Clubcard {
    index: ClubcardIndex,
    approx_filter: Vec<Vec<u64>>,
    exact_filter: Vec<u64>,
}

impl<const W: usize, T: Filterable<W>> From<ClubcardBuilder<W, T>> for Clubcard {
    fn from(builder: ClubcardBuilder<W, T>) -> Clubcard {
        builder.build()
    }
}

impl Clubcard {
    pub fn to_bytes(&self) -> Vec<u8> {
        bincode::serialize(self).unwrap()
    }

    pub fn from_bytes(&self, bytes: &[u8]) -> Self {
        bincode::deserialize(bytes).unwrap()
    }

    pub fn contains<const W: usize>(&self, item: &impl Filterable<W>) -> bool {
        let Some(meta) = self.index.get(item.shard()) else {
            return false;
        };

        if meta.approx_filter_m == 0 {
            return false;
        }

        let mut approx_eq = item.as_equation(meta.approx_filter_m);
        approx_eq.s += meta.approx_filter_offset;

        for i in 0..max(1, meta.approx_filter_rank) {
            if approx_eq.eval(&self.approx_filter[i]) != 0 {
                return false;
            }
        }

        if meta.exact_filter_m == 0 {
            return false;
        }

        let mut exact_eq = item.as_equation(meta.exact_filter_m);
        exact_eq.s += meta.exact_filter_offset;

        if exact_eq.eval(&self.exact_filter) != 0 {
            return false;
        }

        for error in &meta.exclude_errors {
            if error == item.discriminant() {
                return false;
            }
        }

        true
    }
}

//XXX: move this to rust-create-cascade
//TODO: can we avoid copying members?
#[derive(Clone, Debug)]
pub struct CRLiteKey(
    /// issuer spki hash
    pub [u8; 32],
    /// serial number. TODO: Can we use an array here?
    pub Vec<u8>,
    /// revocation status
    pub bool,
);

impl Filterable<4> for CRLiteKey {
    fn as_equation(&self, m: usize) -> Equation<4> {
        let mut hasher = Sha256::new();
        hasher.update(self.0);
        hasher.update(&self.1);
        let mut a = [0u64; 4];
        for (i, x) in hasher
            .finalize()
            .as_slice()
            .chunks(std::mem::size_of::<u64>())
            .enumerate()
        {
            a[i] = u64::from_le_bytes(x.try_into().unwrap());
        }
        a[0] |= 1;
        // TODO better position selection
        let s = (a[0] % max(1, m) as u64) as usize;
        let b = if self.2 { 0 } else { 1 };
        Equation { s, a, b }
    }

    fn shard(&self) -> &[u8] {
        &self.0[..]
    }

    fn discriminant(&self) -> &[u8] {
        &self.1
    }

    fn included(&self) -> bool {
        self.2
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::distributions::Uniform;

    impl<const W: usize> Filterable<W> for Equation<W> {
        fn as_equation(&self, _m: usize) -> Equation<W> {
            self.clone()
        }

        fn shard(&self) -> &[u8] {
            &[]
        }

        fn discriminant(&self) -> &[u8] {
            unsafe { std::mem::transmute(&self.a[..]) }
        }

        fn included(&self) -> bool {
            self.b == 0
        }
    }

    #[test]
    fn test_equation_add() {
        let mut e1 = Equation {
            s: 127,
            a: [0b11],
            b: 1,
        };
        let e2 = Equation {
            s: 127,
            a: [0b01],
            b: 1,
        };
        e1.add(&e2);
        // test that shifting works
        assert!(e1.s == 128);
        assert!(e1.a[0] == 0b1);
        assert!(e1.b == 0);

        let mut e1 = Equation {
            s: 127,
            a: [0b11, 0b1110, 0b1, 0],
            b: 1,
        };
        let e2 = Equation {
            s: 127,
            a: [0b01, 0b0100, 0b0, 0],
            b: 1,
        };
        e1.add(&e2);
        // test that shifting works
        assert!(e1.s == 128);
        assert!(e1.a[0] == 0b1);
        // test that bits move between limbs
        assert!(e1.a[1] == (1 << 63) | 0b101);
        assert!(e1.a[2] == 0);
        assert!(e1.a[3] == 0);
        assert!(e1.b == 0);
    }

    #[test]
    fn test_equation_eval() {
        for s in 0..64 {
            let eq = Equation {
                s,
                a: [0xffffffffffffffff, 0, 0, 0],
                b: 0,
            };
            assert!(0 == eq.eval(&[]));
            for i in 0..64 {
                assert!(((i >= eq.s) as u8) == eq.eval(&[1 << i, 0]));
                assert!(((i < eq.s) as u8) == eq.eval(&[0, 1 << i]));
                assert!(0 == eq.eval(&[0, 0, 1 << i]));
            }
        }
    }

    #[test]
    fn test_solve_identity() {
        let n = 1024;
        let mut builder = RibbonBuilder::new(&[], None);
        for i in 0usize..n {
            let mut eq = Equation::<1>::std(i);
            eq.b = (i % 2) as u8;
            builder.insert(eq);
        }
        let ribbon = ExactRibbon::from(builder);
        let filter = ShardedRibbonFilter::from(vec![ribbon]);
        for i in 0usize..n {
            let eq = Equation::<1>::std(i);
            assert!(eq.eval(&filter.solution[0]) == (i % 2) as u8);
        }
    }

    #[test]
    fn test_solve_empty() {
        let builder = RibbonBuilder::<4, Equation<4>>::new(&"0", None);
        let ribbon = ApproximateRibbon::from(builder);
        let filter = ShardedRibbonFilter::from(vec![ribbon]);
        assert!(!filter.contains(&Equation::std(0)));
    }

    #[test]
    fn test_solve_random() {
        let n = 1024;
        const W: usize = 2;
        let mut r = Ribbon::<W, Equation<W>, Exact>::new(&"0", n);
        let mut s_dist = Uniform::new(0, r.m);
        let mut eqs = Vec::with_capacity(n);
        for _ in 0..n {
            let eq = Equation::<W>::rand(&mut s_dist);
            eqs.push(eq.clone());
            r.insert(eq);
        }
        let x = r.solve(&[]);
        for eq in &eqs {
            assert!(eq.eval(&x) == eq.b);
        }
    }

    #[test]
    fn test_clubcard() {
        let subset_sizes = [1 << 16, 1 << 15, 1 << 14, 1 << 13];
        let universe_size = 1 << 18;

        let mut clubcard_builder = ClubcardBuilder::new();
        let mut approx_builders = vec![];
        for (i, n) in subset_sizes.iter().enumerate() {
            let mut r = clubcard_builder.get_approx_builder(&[i as u8; 32]);
            for j in 0usize..*n {
                let eq = CRLiteKey([i as u8; 32], j.to_le_bytes().to_vec(), true);
                r.insert(eq);
            }
            r.set_universe_size(universe_size);
            approx_builders.push(r)
        }

        let approx_ribbons = approx_builders
            .drain(..)
            .map(ApproximateRibbon::from)
            .collect();

        println!("Approx ribbons:");
        for r in &approx_ribbons {
            println!("\t{}", r);
        }

        clubcard_builder.collect_approx_ribbons(approx_ribbons);

        let mut exact_builders = vec![];
        for (i, n) in subset_sizes.iter().enumerate() {
            let mut r = clubcard_builder.get_exact_builder(&[i as u8; 32]);
            for j in 0usize..universe_size {
                let item = CRLiteKey([i as u8; 32], j.to_le_bytes().to_vec(), j < *n);
                r.insert(item);
            }
            exact_builders.push(r)
        }

        let exact_ribbons = exact_builders.drain(..).map(ExactRibbon::from).collect();

        println!("Exact ribbons:");
        for r in &exact_ribbons {
            println!("\t{}", r);
        }

        clubcard_builder.collect_exact_ribbons(exact_ribbons);

        let clubcard = clubcard_builder.build();
        let size = 8 * clubcard.to_bytes().len();
        println!("Serialized clubcard size: {}kB", size / 8 / 1024);

        let sum_subset_sizes: usize = subset_sizes.iter().sum();
        let sum_universe_sizes: usize = subset_sizes.len() * universe_size;
        let min_size = (sum_subset_sizes as f64)
            * ((sum_universe_sizes as f64) / (sum_subset_sizes as f64)).log2()
            + 1.44 * ((sum_subset_sizes) as f64);
        println!("Size overhead {}", size as f64 / min_size);
        println!("Checking construction");
        println!(
            "\texpecting {} included, {} excluded",
            sum_subset_sizes,
            subset_sizes.len() * universe_size - sum_subset_sizes
        );

        let mut included = 0;
        let mut excluded = 0;
        for i in 0..subset_sizes.len() {
            for j in 0..universe_size {
                let item = CRLiteKey(
                    [i as u8; 32],
                    j.to_le_bytes().to_vec(),
                    /* unused */ false,
                );
                if clubcard.contains(&item) {
                    included += 1;
                } else {
                    excluded += 1;
                }
            }
        }
        println!("\tfound {} included, {} excluded", included, excluded);
        assert!(sum_subset_sizes == included);
        assert!(sum_universe_sizes - sum_subset_sizes == excluded);
    }
}
