/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

use crate::Equation;
use crate::crlite::{CRLiteCoverage, CRLiteQuery};
use serde::Deserialize;
use std::collections::HashMap;

use crate::Filterable;
use base64::Engine;
use std::io::Read;

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

impl<'a> From<&'a CRLiteBuilderItem> for CRLiteQuery<'a> {
    fn from(item: &'a CRLiteBuilderItem) -> Self {
        Self {
            issuer: &item.issuer,
            serial: &item.serial,
            log_timestamps: None,
        }
    }
}


#[cfg(test)]
mod tests {
    use crate::builder::*;
    use crate::crlite::*;
    use crate::*;
    use std::collections::HashMap;

    #[test]
    fn test_crlite_clubcard() {
        let subset_sizes = [1 << 17, 1 << 16, 1 << 15, 1 << 14, 1 << 13];
        let universe_size = 1 << 18;

        let mut clubcard_builder = ClubcardBuilder::new();
        let mut approx_builders = vec![];
        for (i, n) in subset_sizes.iter().enumerate() {
            let mut r = clubcard_builder.get_approx_builder(&[i as u8; 32]);
            for j in 0usize..*n {
                let eq = CRLiteBuilderItem::revoked([i as u8; 32], j.to_le_bytes().to_vec());
                r.insert(eq);
            }
            r.set_universe_size(universe_size);
            approx_builders.push(r)
        }

        let approx_ribbons = approx_builders
            .drain(..)
            .map(ApproximateRibbon::from)
            .collect();

        println!("Approx ribbons:");
        for r in &approx_ribbons {
            println!("\t{}", r);
        }

        clubcard_builder.collect_approx_ribbons(approx_ribbons);

        let mut exact_builders = vec![];
        for (i, n) in subset_sizes.iter().enumerate() {
            let mut r = clubcard_builder.get_exact_builder(&[i as u8; 32]);
            for j in 0usize..universe_size {
                let item = if j < *n {
                    CRLiteBuilderItem::revoked([i as u8; 32], j.to_le_bytes().to_vec())
                } else {
                    CRLiteBuilderItem::not_revoked([i as u8; 32], j.to_le_bytes().to_vec())
                };
                r.insert(item);
            }
            exact_builders.push(r)
        }

        let exact_ribbons = exact_builders.drain(..).map(ExactRibbon::from).collect();

        println!("Exact ribbons:");
        for r in &exact_ribbons {
            println!("\t{}", r);
        }

        clubcard_builder.collect_exact_ribbons(exact_ribbons);

        let mut log_coverage = HashMap::new();
        log_coverage.insert([0u8; 32], (0u64, u64::MAX));

        let clubcard =
            clubcard_builder.build::<CRLiteQuery>(CRLiteCoverage(log_coverage), Default::default());
        let size = 8 * clubcard
            .to_bytes()
            .expect("serialization should succeed")
            .len();
        println!("Serialized clubcard size: {}kB", size / 8 / 1024);

        let sum_subset_sizes: usize = subset_sizes.iter().sum();
        let sum_universe_sizes: usize = subset_sizes.len() * universe_size;
        let min_size = (sum_subset_sizes as f64)
            * ((sum_universe_sizes as f64) / (sum_subset_sizes as f64)).log2()
            + 1.44 * ((sum_subset_sizes) as f64);
        println!("Size overhead {}", size as f64 / min_size);
        println!("Checking construction");
        println!(
            "\texpecting {} included, {} excluded",
            sum_subset_sizes,
            subset_sizes.len() * universe_size - sum_subset_sizes
        );

        let mut included = 0;
        let mut excluded = 0;
        for i in 0..subset_sizes.len() {
            let issuer = [i as u8; 32];
            for j in 0..universe_size {
                let serial = j.to_le_bytes();
                let item = CRLiteQuery {
                    issuer: &issuer,
                    serial: &serial,
                    log_timestamps: None,
                };
                if clubcard.unchecked_contains(&item) {
                    included += 1;
                } else {
                    excluded += 1;
                }
            }
        }
        println!("\tfound {} included, {} excluded", included, excluded);
        assert!(sum_subset_sizes == included);
        assert!(sum_universe_sizes - sum_subset_sizes == excluded);

        // Test that querying a serial from a never-before-seen issuer results in a non-member return.
        let issuer = [subset_sizes.len() as u8; 32];
        let serial = 0usize.to_le_bytes();
        let item = CRLiteQuery {
            issuer: &issuer,
            serial: &serial,
            log_timestamps: None,
        };
        assert!(!clubcard.unchecked_contains(&item));

        assert!(subset_sizes.len() > 0 && subset_sizes[0] > 0 && subset_sizes[0] < universe_size);
        let issuer = [0u8; 32];
        let revoked_serial = 0usize.to_le_bytes();
        let nonrevoked_serial = (universe_size - 1).to_le_bytes();

        // Test that calling contains() a without a timestamp results in a NotInUniverse return
        let item = CRLiteQuery {
            issuer: &issuer,
            serial: &revoked_serial,
            log_timestamps: None,
        };
        assert!(matches!(
            clubcard.contains(&item),
            Membership::NotInUniverse
        ));

        // Test that calling contains() without a timestamp in a covered interval results in a
        // Member return.
        let timestamps = [([0u8; 32], 100)];
        let item = CRLiteQuery {
            issuer: &issuer,
            serial: &revoked_serial,
            log_timestamps: Some(&timestamps),
        };
        assert!(matches!(clubcard.contains(&item), Membership::Member));

        // Test that calling contains() without a timestamp in a covered interval results in a
        // Member return.
        let timestamps = [([0u8; 32], 100)];
        let item = CRLiteQuery {
            issuer: &issuer,
            serial: &nonrevoked_serial,
            log_timestamps: Some(&timestamps),
        };
        assert!(matches!(clubcard.contains(&item), Membership::Nonmember));

        // Test that calling contains() without a timestamp in a covered interval results in a
        // Member return.
        let timestamps = [([1u8; 32], 100)];
        let item = CRLiteQuery {
            issuer: &issuer,
            serial: &revoked_serial,
            log_timestamps: Some(&timestamps),
        };
        assert!(matches!(
            clubcard.contains(&item),
            Membership::NotInUniverse
        ));
    }
}
