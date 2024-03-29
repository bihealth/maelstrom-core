use std::collections::{HashMap, HashSet};

use rust_htslib::{bcf, bcf::Read};

use crate::common::stats::Stats;

/// Name of the autosomes.
static NAMES_AUTOSOMES: &[&str] = &[
    "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17",
    "18", "19", "20", "21", "22", "chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8",
    "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
    "chr19", "chr20", "chr21", "chr22", // The following are for testing only.
    "one", "two",
];

/// Name of the gonomosomes.
static NAMES_ALLOSOMES: &[&str] = &["X", "Y", "chrX", "chrY"];

/// Store information about read depth.
pub struct MedianReadDepthInfo {
    /// Median read depth by chromosome.
    pub by_chrom: HashMap<String, f64>,
    /// Median read depth on autosomes.
    pub on_autosomes: f64,
}

/// Load DoC from file and compute median.
pub fn load_doc_median(path: &str) -> Result<MedianReadDepthInfo, anyhow::Error> {
    let autosomes: HashSet<String> = NAMES_AUTOSOMES.iter().map(|&s| s.to_string()).collect();
    let all_chroms: HashSet<String> = (&[NAMES_ALLOSOMES, NAMES_AUTOSOMES])
        .concat()
        .iter()
        .map(|s| s.to_string())
        .collect();
    let mut reader = bcf::Reader::from_path(path)?;
    let mut rcvs_by_chrom: Vec<Vec<f64>> = vec![vec![]; reader.header().contig_count() as usize];

    let mut record = reader.empty_record();
    while let Some(result) = reader.read(&mut record) {
        result?;
        rcvs_by_chrom[record.rid().unwrap() as usize]
            .push(record.format(b"RCV").float()?[0][0].into());
    }

    let mut by_chrom: HashMap<String, f64> = HashMap::new();
    let mut rcvs_autosomes: Vec<f64> = Vec::new();
    for rid in 0..(reader.header().contig_count()) {
        let chrom = String::from_utf8(reader.header().rid2name(rid)?.to_vec())?;
        if all_chroms.contains(&chrom) {
            let median = if rcvs_by_chrom.is_empty() {
                0.0
            } else if !rcvs_by_chrom[rid as usize].is_empty() {
                (&rcvs_by_chrom[rid as usize]).median()
            } else {
                0.0
            };
            by_chrom.insert(chrom.clone(), median);
        }
        if autosomes.contains(&chrom) {
            rcvs_autosomes.append(&mut rcvs_by_chrom[rid as usize].to_vec());
        }
    }

    let on_autosomes = if rcvs_autosomes.is_empty() {
        0.0
    } else {
        (&rcvs_autosomes).median()
    };

    Ok(MedianReadDepthInfo {
        by_chrom,
        on_autosomes,
    })
}
