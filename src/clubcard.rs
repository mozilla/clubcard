/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

use crate::query::Queryable;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

#[derive(Debug)]
pub enum ClubcardError {
    Deserialize,
    Serialize,
    UnsupportedVersion,
}

/// Metadata needed to compute membership in a clubcard.
#[derive(Default, Serialize, Deserialize)]
pub struct ClubcardShardMeta {
    /// Description of the hash function h.
    pub approx_filter_m: usize,
    /// Description of the hash function g.
    pub exact_filter_m: usize,
    /// The number of columns in X.
    pub approx_filter_rank: usize,
    /// An offset t such that [0^t || h(u)] * X = h(u) * Xi, where i is the shard identifier.
    pub approx_filter_offset: usize,
    /// An offset t such that [0^t || g(u)] * Y = g(u) * Yi, where i is the shard identifier.
    pub exact_filter_offset: usize,
    /// Whether to invert the output of queries to this shard.
    pub inverted: bool,
    /// A list of elements of Ui \ Ri that are not correctly encoded by this shard.
    pub exceptions: Vec<Vec<u8>>,
}

/// Lookup table from shard identifiers to shard metadata.
pub type ClubcardIndex = BTreeMap</* block id */ Vec<u8>, ClubcardShardMeta>;

/// A queryable Clubcard
#[derive(Serialize, Deserialize)]
pub struct Clubcard<const W: usize, T: Queryable<W>> {
    /// Shard lookup table
    pub(crate) index: ClubcardIndex,
    /// The matrix X
    pub(crate) approx_filter: Vec<Vec<u64>>,
    /// The matrix Y
    pub(crate) exact_filter: Vec<u64>,
    pub(crate) phantom: std::marker::PhantomData<T>,
}

impl<const W: usize, T: Queryable<W>> fmt::Display for Clubcard<W, T> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let approx_size = 8 * self.approx_filter.iter().map(|x| x.len()).sum::<usize>();
        let exact_size = 8 * self.exact_filter.len();
        let exceptions = self
            .index
            .values()
            .map(|meta| meta.exceptions.len())
            .sum::<usize>();
        writeln!(
            f,
            "Clubcard of size {} ({} + {})",
            approx_size + exact_size,
            approx_size,
            exact_size
        )?;
        writeln!(f, "- exceptions: {}", exceptions)?;
        writeln!(
            f,
            "- Serialized size: {}",
            self.to_bytes().map_or(0, |x| x.len())
        )
    }
}

impl<const W: usize, T: Queryable<W>> Clubcard<W, T> {
    const SERIALIZATION_VERSION: u16 = 0xffff;

    /// Serialize this clubcard.
    pub fn to_bytes(&self) -> Result<Vec<u8>, ClubcardError> {
        let mut out = u16::to_le_bytes(Self::SERIALIZATION_VERSION).to_vec();
        bincode::serialize_into(&mut out, self).map_err(|_| ClubcardError::Serialize)?;
        Ok(out)
    }

    /// Deserialize a clubcard.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, ClubcardError> {
        let (version_bytes, rest) = bytes.split_at(std::mem::size_of::<u16>());
        let Ok(version_bytes) = version_bytes.try_into() else {
            return Err(ClubcardError::Deserialize);
        };
        let version = u16::from_le_bytes(version_bytes);
        if version != Self::SERIALIZATION_VERSION {
            return Err(ClubcardError::UnsupportedVersion);
        }
        bincode::deserialize(rest).map_err(|_| ClubcardError::Deserialize)
    }

    /// Check whether item is in the set encoded by this clubcard.
    pub fn contains(&self, item: &T) -> bool {
        let Some(meta) = self.index.get(item.shard()) else {
            return false;
        };

        let result = (|| {
            // All queries evaluate to 0 on an empty filter, but logically
            // such a filter does not include anything. So we handle it as a
            // special case.
            if meta.approx_filter_m == 0 {
                return false;
            }

            // Check if h(item) * X is 0
            let approx_query = item.as_approx_query(meta);
            for i in 0..meta.approx_filter_rank {
                if approx_query.eval(&self.approx_filter[i]) != 0 {
                    return false;
                }
            }

            // Check if g(item) * X is 0
            let exact_query = item.as_exact_query(meta);
            if exact_query.eval(&self.exact_filter) != 0 {
                return false;
            }

            for exception in &meta.exceptions {
                if exception == item.discriminant() {
                    return false;
                }
            }
            true
        })();
        result ^ meta.inverted
    }
}
