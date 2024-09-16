/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

use crate::clubcard::ClubcardIndexEntry;
use std::cmp::min;

/// An Equation\<W\> is a representation of a GF(2) linear functional
///     a(x) = b + sum_i a_i x_i
/// where a_i is equal to zero except for i in a block of 64*W coefficients
/// starting at i=s. We say an Equation is /aligned/ if a_s = 1.
/// (Note: a_i above denotes the i-th bit, not the i'th 64-bit limb.)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Equation<const W: usize> {
    pub(crate) s: usize,    // the row number
    pub(crate) a: [u64; W], // the non-trivial columns
    pub(crate) b: u8,       // the constant term
                            // TODO? save some space by using one bit of s for b.
}

impl<const W: usize> Equation<W> {
    /// Construct the aligned equation equivalent to Equation { s, a, b }
    pub fn new(s: usize, a: [u64; W], b: u8) -> Equation<W> {
        let mut eq = Equation { s: 0, a, b };
        eq.add(&Equation::zero());
        eq.s += s;
        eq
    }

    /// Construct the equation a(x) = 0.
    pub fn zero() -> Self {
        Equation {
            s: 0,
            a: [0u64; W],
            b: 0,
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

    /// The block that this item belongs in.
    fn block_id(&self) -> &[u8];

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8];

    /// Whether this item should be included in an exact filter.
    fn included(&self) -> bool {
        false
    }
}

/// Queryable is the same as Filterable except it does not have the `included()` method---one
/// needs a Filterable to construct a filter, but a Queryable suffices to query a filter.
/// A default implementation is provided for any Filterable.
pub trait Queryable<const W: usize> {
    /// Given the metadata describing a block of a clubcard computes h(self).
    fn as_approx_query(&self, meta: &ClubcardIndexEntry) -> Equation<W>;

    /// Given the metadata describing a block of a clubcard computes g(self).
    fn as_exact_query(&self, meta: &ClubcardIndexEntry) -> Equation<W>;

    /// The block that this item belongs in.
    fn block_id(&self) -> &[u8];

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8];
}

impl<const W: usize, T: Filterable<W>> Queryable<W> for T {
    fn as_approx_query(&self, meta: &ClubcardIndexEntry) -> Equation<W> {
        let mut approx_eq = self.as_equation(meta.approx_filter_m);
        approx_eq.s += meta.approx_filter_offset;
        approx_eq
    }

    fn as_exact_query(&self, meta: &ClubcardIndexEntry) -> Equation<W> {
        let mut approx_eq = self.as_equation(meta.exact_filter_m);
        approx_eq.s += meta.exact_filter_offset;
        approx_eq
    }

    /// The block that this item belongs in.
    fn block_id(&self) -> &[u8] {
        Filterable::block_id(self)
    }

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8] {
        Filterable::discriminant(self)
    }
}

#[cfg(test)]
mod tests {
    use crate::Equation;

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
}
