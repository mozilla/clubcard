use crate::{Equation, Filterable};
use sha2::{Digest, Sha256};
use std::cmp::max;

#[derive(Clone, Debug)]
pub struct CRLiteKey {
    /// issuer spki hash
    pub issuer: [u8; 32],
    /// serial number. TODO: smallvec?
    pub serial: Vec<u8>,
    /// revocation status
    pub revoked: bool,
}

impl CRLiteKey {
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

    pub fn query(issuer: [u8; 32], serial: Vec<u8>) -> Self {
        Self {
            issuer,
            serial,
            revoked: false, /* unused */
        }
    }
}

impl Filterable<4> for CRLiteKey {
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
        let b = if self.revoked { 0 } else { 1 };
        Equation { s, a, b }
    }

    fn shard(&self) -> &[u8] {
        self.issuer.as_ref()
    }

    fn discriminant(&self) -> &[u8] {
        &self.serial
    }

    fn included(&self) -> bool {
        self.revoked
    }
}
