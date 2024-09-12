/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

use crate::query::Queryable;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fmt;

#[derive(Default, Serialize, Deserialize)]
pub struct ClubcardShardMeta {
    pub approx_filter_offset: usize,
    pub approx_filter_m: usize,
    pub approx_filter_rank: usize,
    pub exact_filter_offset: usize,
    pub exact_filter_m: usize,
    pub inverted: bool,
    pub exceptions: Vec<Vec<u8>>,
}

pub type ClubcardIndex = BTreeMap</* block id */ Vec<u8>, ClubcardShardMeta>;

#[derive(Serialize, Deserialize)]
pub struct Clubcard {
    pub index: ClubcardIndex,
    pub approx_filter: Vec<Vec<u64>>,
    pub exact_filter: Vec<u64>,
}

impl fmt::Display for Clubcard {
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
        writeln!(f, "- Serialized size: {}", self.to_bytes().len())
    }
}

impl Clubcard {
    pub fn to_bytes(&self) -> Vec<u8> {
        bincode::serialize(self).unwrap()
    }

    pub fn from_bytes(bytes: &[u8]) -> Self {
        bincode::deserialize(bytes).unwrap()
    }

    pub fn contains<const W: usize>(&self, item: &impl Queryable<W>) -> bool {
        let Some(meta) = self.index.get(item.shard()) else {
            return false;
        };

        let result = (|| {
            if meta.approx_filter_m == 0 {
                return false;
            }

            let approx_eq = item.as_approx_query(meta);
            for i in 0..meta.approx_filter_rank {
                if approx_eq.eval(&self.approx_filter[i]) != 0 {
                    return false;
                }
            }

            if meta.exact_filter_m == 0 {
                return false;
            }

            let exact_eq = item.as_exact_query(meta);
            if exact_eq.eval(&self.exact_filter) != 0 {
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
