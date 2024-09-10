use crate::{Equation, Filterable};
use sha2::{Digest, Sha256};
use std::cmp::max;

#[derive(Clone, Debug)]
pub struct CRLiteKey(
    /// issuer spki hash
    pub [u8; 32],
    /// serial number. TODO: smallvec?
    pub Vec<u8>,
    /// revocation status
    pub bool,
);

impl Filterable<4> for CRLiteKey {
    fn as_equation(&self, m: usize) -> Equation<4> {
        let mut digest = [0u8; 32];
        let mut hasher = Sha256::new();
        hasher.update(self.0);
        hasher.update(&self.1);
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
        let b = if self.2 { 0 } else { 1 };
        Equation { s, a, b }
    }

    fn shard(&self) -> &[u8] {
        self.0.as_ref()
    }

    fn discriminant(&self) -> &[u8] {
        &self.1
    }

    fn included(&self) -> bool {
        self.2
    }
}
