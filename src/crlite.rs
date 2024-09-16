use crate::{AsQuery, ClubcardIndexEntry, Equation, Filterable, Queryable};
use sha2::{Digest, Sha256};
use std::cmp::max;

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

    /// A unique identifier for this item. If this item cannot be inserted into the linear system,
    /// then we will store its `included()` status in a secondary retrieval mechanism keyed by
    /// `discriminant()`.
    fn discriminant(&self) -> &[u8] {
        &self.serial
    }
}

impl<'a> Queryable<4> for CRLiteQuery<'a> {
    // TODO: Replace this with HashMap<ct log id, (min timestamp, max timestamp)>
    type UniverseMetadata = ();

    // The set of CRLiteKeys is partitioned by issuer, and each
    // CRLiteKey knows its issuer. So there's no need for additional
    // partition metadata.
    type PartitionMetadata = ();

    fn block_id(&self, _meta: &Self::PartitionMetadata) -> Option<&[u8]> {
        Some(self.issuer.as_ref())
    }

    fn in_universe(&self, _universe: &Self::UniverseMetadata) -> bool {
        true
    }
}
