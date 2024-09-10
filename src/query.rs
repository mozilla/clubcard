/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

use crate::clubcard::ClubcardShardMeta;
#[cfg(test)]
use rand::{distributions::Distribution, Rng};
use std::cmp::min;

/// An Equation\<W\> is a representation of a GF(2) linear functional
///     a(x) = b + sum_i a_i x_i
/// where a_i is equal to zero except for i in a block of 64*W coefficients
/// starting at i=s. We say an Equation is /aligned/ if a_s = 1.
/// (Note: a_i above denotes the i-th bit, not the i'th 64-bit limb.)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Equation<const W: usize> {
    /* TODO restrict visibility */
    pub s: usize,    // the row number
    pub a: [u64; W], // the non-trivial columns
    pub b: u8,       // the constant term
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
    fn included(&self) -> bool {
        false
    }
}

pub trait Queryable<const W: usize> {
    fn as_approx_query(&self, meta: &ClubcardShardMeta) -> Equation<W>;
    fn as_exact_query(&self, meta: &ClubcardShardMeta) -> Equation<W>;

    /// The shard that this item belongs in.
    fn shard(&self) -> &[u8];

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8];
}

impl<const W: usize, T: Filterable<W>> Queryable<W> for T {
    fn as_approx_query(&self, meta: &ClubcardShardMeta) -> Equation<W> {
        let mut approx_eq = self.as_equation(meta.approx_filter_m);
        approx_eq.s += meta.approx_filter_offset;
        approx_eq
    }

    fn as_exact_query(&self, meta: &ClubcardShardMeta) -> Equation<W> {
        let mut approx_eq = self.as_equation(meta.exact_filter_m);
        approx_eq.s += meta.exact_filter_offset;
        approx_eq
    }

    /// The shard that this item belongs in.
    fn shard(&self) -> &[u8] {
        Filterable::shard(self)
    }

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8] {
        Filterable::discriminant(self)
    }
}
