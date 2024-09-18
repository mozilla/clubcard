/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

use crate::clubcard::ClubcardIndexEntry;
use crate::equation::Equation;

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

pub trait AsQuery<const W: usize> {
    fn as_approx_query(&self, meta: &ClubcardIndexEntry) -> Equation<W>;

    fn as_exact_query(&self, meta: &ClubcardIndexEntry) -> Equation<W>;

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8];
}

impl<const W: usize, T: Filterable<W>> AsQuery<W> for T {
    fn as_approx_query(&self, meta: &ClubcardIndexEntry) -> Equation<W> {
        let mut approx_eq = self.as_equation(meta.approx_filter_m);
        approx_eq.s += meta.approx_filter_offset;
        approx_eq
    }

    fn as_exact_query(&self, meta: &ClubcardIndexEntry) -> Equation<W> {
        let mut exact_eq = self.as_equation(meta.exact_filter_m);
        exact_eq.s += meta.exact_filter_offset;
        exact_eq
    }

    fn discriminant(&self) -> &[u8] {
        Filterable::discriminant(self)
    }
}

pub trait Queryable<const W: usize>: AsQuery<W> {
    type UniverseMetadata;
    type PartitionMetadata;

    /// Use the PartitionMetadata associated with a Clubcard to identify the
    /// block of the Clubcard that this item belongs to.
    /// Returns None if this item does not belong to any block.
    fn block_id(&self, meta: &Self::PartitionMetadata) -> Option<&[u8]>;

    /// Use the UniverseMetadata associated with a Clubcard to determine
    /// whether or not this item is in the Clubcard's encoded universe.
    fn in_universe(&self, meta: &Self::UniverseMetadata) -> bool;
}
