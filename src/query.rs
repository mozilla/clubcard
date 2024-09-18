/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

use crate::clubcard::ClubcardIndexEntry;
use crate::equation::Equation;

pub trait AsQuery<const W: usize> {
    /// The block that this item belongs in.
    fn block_id(&self) -> &[u8];

    /// Hash this item to a homogeneous equation (s, a) such that
    ///     (1) s is uniform in {0, 1, ..., m-1},
    ///     (2) a satisfies the alignment requirement (a\[0\] & 1 == 1) but is otherwise uniformly random,
    fn as_query(&self, m: usize) -> Equation<W>;

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8];

    #[doc(hidden)]
    fn as_approx_query(&self, meta: &ClubcardIndexEntry) -> Equation<W> {
        let mut approx_eq = self.as_query(meta.approx_filter_m);
        approx_eq.s += meta.approx_filter_offset;
        approx_eq
    }

    #[doc(hidden)]
    fn as_exact_query(&self, meta: &ClubcardIndexEntry) -> Equation<W> {
        let mut exact_eq = self.as_query(meta.exact_filter_m);
        exact_eq.s += meta.exact_filter_offset;
        exact_eq
    }
}

/// A Filterable is an item that can be inserted into a RibbonBuilder.
pub trait Filterable<const W: usize>: AsQuery<W> {
    /// Whether this item should be included in an exact filter.
    fn included(&self) -> bool {
        false
    }
}

/// A Queryable is an item that can be passed to Clubcard::contains.
pub trait Queryable<const W: usize>: AsQuery<W> {
    type UniverseMetadata;
    type PartitionMetadata;

    /// Use the PartitionMetadata associated with a Clubcard to identify the
    /// block of the Clubcard that this item belongs to.
    /// Returns None if this item does not belong to any block.
    fn get_block_id(&self, meta: &Self::PartitionMetadata) -> Option<&[u8]>;

    /// Use the UniverseMetadata associated with a Clubcard to determine
    /// whether or not this item is in the Clubcard's encoded universe.
    fn in_universe(&self, meta: &Self::UniverseMetadata) -> bool;
}
