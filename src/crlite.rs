use crate::{AsQuery, ClubcardIndexEntry, Equation, Filterable, Queryable};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::cmp::max;
use std::collections::HashMap;

#[cfg(feature = "builder")]
use base64::Engine;
#[cfg(feature = "builder")]
use std::io::Read;

type LogId = [u8; 32];
type TimestampInterval = (u64, u64);

#[derive(Serialize, Deserialize)]
pub struct CRLiteCoverage(pub HashMap<LogId, TimestampInterval>);

#[cfg(feature = "builder")]
impl CRLiteCoverage {
    pub fn from_mozilla_ct_logs_json<T>(reader: T) -> Self
    where
        T: Read,
    {
        #[allow(non_snake_case)]
        #[derive(Deserialize)]
        struct MozillaCtLogsJson {
            LogID: String,
            MaxTimestamp: u64,
            MinTimestamp: u64,
        }

        let mut coverage = HashMap::new();
        let json_entries: Vec<MozillaCtLogsJson> = match serde_json::from_reader(reader) {
            Ok(json_entries) => json_entries,
            _ => return CRLiteCoverage(Default::default()),
        };
        for entry in json_entries {
            let mut log_id = [0u8; 32];
            match base64::prelude::BASE64_STANDARD.decode(&entry.LogID) {
                Ok(bytes) if bytes.len() == 32 => log_id.copy_from_slice(&bytes),
                _ => continue,
            };
            coverage.insert(log_id, (entry.MinTimestamp, entry.MaxTimestamp));
        }
        CRLiteCoverage(coverage)
    }
}

#[derive(Clone, Debug)]
pub struct CRLiteBuilderItem {
    /// issuer spki hash
    pub issuer: [u8; 32],
    /// serial number. TODO: smallvec?
    pub serial: Vec<u8>,
    /// revocation status
    pub revoked: bool,
}

impl CRLiteBuilderItem {
    pub fn revoked(issuer: [u8; 32], serial: Vec<u8>) -> Self {
        Self {
            issuer,
            serial,
            revoked: true,
        }
    }

    pub fn not_revoked(issuer: [u8; 32], serial: Vec<u8>) -> Self {
        Self {
            issuer,
            serial,
            revoked: false,
        }
    }
}

impl Filterable<4> for CRLiteBuilderItem {
    fn as_equation(&self, m: usize) -> Equation<4> {
        let mut eq = CRLiteQuery::from(self).as_equation(m);
        eq.b = if self.revoked { 0 } else { 1 };
        eq
    }

    fn block_id(&self) -> &[u8] {
        self.issuer.as_ref()
    }

    fn discriminant(&self) -> &[u8] {
        &self.serial
    }

    fn included(&self) -> bool {
        self.revoked
    }
}

#[derive(Clone, Debug)]
pub struct CRLiteQuery<'a> {
    pub issuer: &'a [u8; 32],
    pub serial: &'a [u8],
    pub log_timestamps: Option<&'a [([u8; 32], u64)]>,
}

impl<'a> CRLiteQuery<'a> {
    fn as_equation(&self, m: usize) -> Equation<4> {
        let mut digest = [0u8; 32];
        let mut hasher = Sha256::new();
        hasher.update(self.issuer);
        hasher.update(&self.serial);
        hasher.finalize_into((&mut digest).into());

        let mut a = [0u64; 4];
        for (i, x) in digest
            .chunks_exact(8) // TODO: use array_chunks::<8>() when stable
            .map(|x| TryInto::<[u8; 8]>::try_into(x).unwrap())
            .map(u64::from_le_bytes)
            .enumerate()
        {
            a[i] = x;
        }
        a[0] |= 1;
        let s = (a[3] as usize) % max(1, m);
        Equation { s, a, b: 0 }
    }
}

impl<'a> From<&'a CRLiteBuilderItem> for CRLiteQuery<'a> {
    fn from(item: &'a CRLiteBuilderItem) -> Self {
        Self {
            issuer: &item.issuer,
            serial: &item.serial,
            log_timestamps: None,
        }
    }
}

impl<'a> AsQuery<4> for CRLiteQuery<'a> {
    fn as_approx_query(&self, meta: &ClubcardIndexEntry) -> Equation<4> {
        let mut approx_eq = self.as_equation(meta.approx_filter_m);
        approx_eq.s += meta.approx_filter_offset;
        approx_eq
    }

    fn as_exact_query(&self, meta: &ClubcardIndexEntry) -> Equation<4> {
        let mut exact_eq = self.as_equation(meta.exact_filter_m);
        exact_eq.s += meta.exact_filter_offset;
        exact_eq
    }

    fn discriminant(&self) -> &[u8] {
        &self.serial
    }
}

impl<'a> Queryable<4> for CRLiteQuery<'a> {
    type UniverseMetadata = CRLiteCoverage;

    // The set of CRLiteKeys is partitioned by issuer, and each
    // CRLiteKey knows its issuer. So there's no need for additional
    // partition metadata.
    type PartitionMetadata = ();

    fn block_id(&self, _meta: &Self::PartitionMetadata) -> Option<&[u8]> {
        Some(self.issuer.as_ref())
    }

    fn in_universe(&self, universe: &Self::UniverseMetadata) -> bool {
        let Some(log_timestamps) = self.log_timestamps else {
            return false;
        };
        for (log_id, timestamp) in log_timestamps {
            if let Some((low, high)) = universe.0.get(log_id) {
                if low <= timestamp && timestamp <= high {
                    return true;
                }
            }
        }
        false
    }
}
