/// Command to collect depth of coverage information from BAM file.
pub mod aggregation;
pub mod reference;

use bio_types::genome::{AbstractInterval, Interval};
use clap::Args as ClapArgs;
use console::{Emoji, Term};
use indicatif::{ProgressBar, ProgressStyle};
use itertools::sorted;
use rust_htslib::{bam, bam::Read as BamRead, bcf, bcf::Read as BcfRead};
use separator::Separatable;
use tempfile::tempdir;

use crate::common::bam::{build_chroms_bam, samples_from_file};
use crate::common::bcf::guess_bcf_format;
use crate::common::doc::{load_doc_median, MedianReadDepthInfo};
use crate::common::Args as CommonArgs;
use crate::common::{prefix_lines, ArgInterval};
use aggregation::{BamRecordAggregator, CoverageAggregator, FragmentsAggregator};
use reference::ReferenceStats;

/// Enum for selecting sequencing technology
#[derive(PartialEq, Eq, Debug, Clone, clap::ValueEnum, enum_utils::FromStr)]
pub enum CountKind {
    /// Count fragments starting in bin
    FragmentCount,
    /// Consider base-wise coverage
    Coverage,
}

#[derive(ClapArgs, Debug, Clone)]
pub struct DocConfig {
    /// The coverage kind to count
    #[clap(short = 'k', long = "count-kind", value_enum, default_value_t = CountKind::Coverage)]
    count_kind: CountKind,
    /// The window length
    #[clap(long = "window-length", default_value_t = 1_000)]
    window_length: usize,
    /// Minimal MAPQ when collecting depth of coverage information
    #[clap(long = "min-mapq")]
    min_mapq: Option<u8>,
    /// Minimal unclipped bases when collecting depth of coverage information
    #[clap(long = "min-unclipped")]
    min_unclipped: Option<f32>,
}

/// Command line options
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Optional list of regions to call.
    #[clap(short = 'R', long = "regions")]
    regions: Option<Vec<ArgInterval>>,
    /// Path to reference FASTA file.
    #[clap(short = 'r', long = "reference")]
    path_reference: Option<String>,
    /// Path to input BAM file.
    #[clap(short = 'i', long = "in", required = true)]
    path_input: String,
    /// Path to output VCF/BCF file.
    #[clap(short = 'o', long = "out", required = true)]
    path_output: String,
    /// Force overwriting output file.
    #[clap(short = 'f', long = "force", action, default_value_t = false)]
    force_overwrite: bool,
    #[clap(flatten)]
    doc_config: DocConfig,
}

/// Build header for the coverage output BCF file.
fn build_header(samples: &[String], contigs: &[Interval]) -> bcf::Header {
    let mut header = bcf::Header::new();

    // Put overall meta information into the BCF header.
    let now = chrono::Utc::now();
    if !cfg!(test) {
        header.push_record(format!("##fileDate={}", now.format("%Y%m%d")).as_bytes());
    } else {
        header.push_record(b"##fileDate=20200828");
    }

    // Add samples to BCF header.
    for sample in samples {
        header.push_sample(sample.as_bytes());
    }

    // Put contig information into BCF header.
    contigs.iter().for_each(|contig| {
        header.push_record(
            format!(
                "##contig=<ID={},length={}>",
                contig.contig(),
                contig.range().end
            )
            .as_bytes(),
        );
    });

    // Push the relevant header records.
    // TODO: later decide about commented-out lines
    let lines = vec![
        // Define ALT column <WINDOW>/<TARGET>
        "##ALT=<ID=WINDOW,Description=\"Record describes a window for read or coverage \
         counting\">",
        // INFO fields describing the window
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Window end\">",
        "##INFO=<ID=MAPQ,Number=1,Type=Float,Description=\"Mean MAPQ value across samples \
         for approximating mapability\">",
        "##INFO=<ID=GC,Number=1,Type=Float,Description=\"Reference GC fraction, if reference \
         FASTA file was given\">",
        "##INFO=<ID=GAP,Number=0,Type=Flag,Description=\"Window overlaps with N in \
         reference (gap)\">",
        // Generic FORMAT fields
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
        "##FORMAT=<ID=MQ,Number=1,Type=Float,Description=\"Mean read MAPQ from region\">",
        // The meaning of coverage differs between the counting approaches (fragments vs. coverage).
        "##FORMAT=<ID=RCV,Number=1,Type=Float,Description=\"Raw coverage value\">",
        "##FORMAT=<ID=RCVSD,Number=1,Type=Float,Description=\"Raw coverage standard deviation\">",
        "##FORMAT=<ID=CV,Number=1,Type=Float,Description=\"Normalized coverage value\">",
        "##FORMAT=<ID=CVSD,Number=1,Type=Float,Description=\"Normalized coverage standard deviation\">",
    ];
    for line in lines {
        header.push_record(line.as_bytes());
    }

    header
}

/// Build bcf::Writer with appropriate header.
fn build_bcf_writer(
    path: &str,
    samples: &[String],
    contigs: &[Interval],
) -> Result<bcf::Writer, anyhow::Error> {
    let guessed = guess_bcf_format(path);

    let header = build_header(samples, contigs);
    Ok(bcf::Writer::from_path(
        &path,
        &header,
        guessed.uncompressed,
        guessed.format,
    )?)
}

/// Process one region.
fn process_region(
    term: &Term,
    _common_args: &CommonArgs,
    args: &Args,
    contig: &Interval,
    bcf_writer: &mut bcf::Writer,
) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "Processing contig {}:{}-{}",
        contig.contig(),
        (contig.range().start + 1).separated_string(),
        contig.range().end.separated_string()
    ))?;

    let window_length = args.doc_config.window_length;
    let ref_stats = if let Some(path) = &args.path_reference {
        Some(ReferenceStats::from_path(
            path,
            contig.contig(),
            window_length,
        )?)
    } else {
        None
    };

    let mut aggregator: Box<dyn BamRecordAggregator> = match args.doc_config.count_kind {
        CountKind::FragmentCount => Box::new(FragmentsAggregator::new(
            args.doc_config.clone(),
            contig.clone(),
        )),
        CountKind::Coverage => Box::new(CoverageAggregator::new(
            args.doc_config.clone(),
            contig.clone(),
        )),
    };

    // TODO: 2 blocks -> function
    // Jump to region with BAM reader.
    let mut bam_reader = bam::IndexedReader::from_path(&args.path_input)?;
    let tid: u32 = bam_reader.header().tid(contig.contig().as_bytes()).unwrap();
    bam_reader.fetch((tid, contig.range().start, contig.range().end))?;

    let progress_bar = {
        let prog_bar =
            ProgressBar::new(((contig.range().end - contig.range().start) / 1_000) as u64);
        prog_bar.set_style(
            ProgressStyle::default_bar()
                .template(
                    "scanning {msg:.green.bold} [{elapsed_precise}] [{wide_bar:.cyan/blue}] \
            {pos:>7}/{len:7} Kbp {elapsed}/{eta}",
                )
                .unwrap()
                .progress_chars("=>-"),
        );
        prog_bar.set_message(String::from(contig.contig()));
        Some(prog_bar)
    };

    // Main loop for region: pass all BAM records in region through aggregator.
    term.write_line("Computing coverage...")?;
    aggregator.put_fetched_records(&mut bam_reader, &|pos| {
        if let Some(prog_bar) = &progress_bar {
            if pos >= 0 {
                prog_bar.set_position((pos / 1_000) as u64);
            }
        }
    })?;
    term.write_line(&format!(
        "Processed {}, skipped {} records ({:.2}% were processed)",
        aggregator.num_processed().separated_string(),
        aggregator.num_skipped().separated_string(),
        100.0 * (aggregator.num_processed() - aggregator.num_skipped()) as f64
            / aggregator.num_processed() as f64
    ))?;

    // TODO: full block -> function
    // Create the BCF records for this region.
    term.write_line("Writing BCF with coverage information...")?;
    for region_id in 0..aggregator.num_regions() {
        let stats = aggregator.get_stats(region_id);
        let mut record = bcf_writer.empty_record();
        let rid = bcf_writer.header().name2rid(contig.contig().as_bytes())?;

        // Columns: CHROM, POS, ID, REF, ALT, (FILTER)
        let pos = stats.start;
        let window_end = stats.end;
        let alleles_v = vec![Vec::from("N"), Vec::from("<WINDOW>")];
        let alleles = alleles_v
            .iter()
            .map(|x| x.as_slice())
            .collect::<Vec<&[u8]>>();

        record.set_rid(Some(rid));
        record.set_pos(pos as i64);
        record.set_id(format!("{}:{}-{}", &contig.contig(), pos + 1, window_end).as_bytes())?;
        record.set_alleles(&alleles)?;

        // Columns: INFO
        record.push_info_integer(b"END", &[window_end as i32])?;
        if let Some(ref_stats) = ref_stats.as_ref() {
            let bucket = pos / window_length;
            let gc = ref_stats.gc_content[bucket];
            if !gc.is_nan() {
                record.push_info_float(b"GC", &[gc])?;
            }
            if ref_stats.has_gap[bucket] {
                record.push_info_flag(b"GAP")?;
            }
        }

        // Columns: FORMAT/GT
        record.push_format_integer(b"GT", &[0, 0])?;

        // Columns: FORMAT/CV etc.
        record.push_format_float(b"RCV", &[stats.cov])?;
        if let Some(cov_sd) = stats.cov_sd {
            record.push_format_float(b"RCVSD", &[cov_sd])?;
        } else {
            record.push_format_float(b"RCVSD", &[0.0_f32])?;
        }

        record.push_format_float(b"MQ", &[stats.mean_mapq])?;
        record.push_info_float(b"MAPQ", &[stats.mean_mapq])?;

        bcf_writer.write(&record)?;
    }

    Ok(())
}

fn perform_final_write(
    path_in: &str,
    path_out: &str,
    doc_median_info: &MedianReadDepthInfo,
) -> Result<(), anyhow::Error> {
    let mut reader = bcf::Reader::from_path(&path_in)?;
    let mut header = bcf::Header::from_template(reader.header());
    let sample = std::str::from_utf8(reader.header().samples()[0])?;
    // NB: we need to prefix the underscore because htslib does not digits in front of keys
    let by_contig = sorted(
        doc_median_info
            .by_chrom
            .iter()
            .map(|(k, v)| format!("_{}={}", k, v)),
    )
    .collect::<Vec<String>>();
    header.push_record(
        format!(
            "##median-coverage=<ID={},autosomes={},{}>",
            &sample,
            doc_median_info.on_autosomes,
            by_contig.join(",")
        )
        .as_bytes(),
    );

    let guessed = guess_bcf_format(path_out);
    let mut writer =
        bcf::Writer::from_path(&path_out, &header, guessed.uncompressed, guessed.format)?;
    let mut record = reader.empty_record();
    while let Some(result) = reader.read(&mut record) {
        result?;
        writer.translate(&mut record);

        let rcvs = record.format(b"RCV").float()?;
        let rcvsds = record.format(b"RCVSD").float()?;

        let cv = &[if doc_median_info.on_autosomes > 1e-6 {
            rcvs[0][0] / doc_median_info.on_autosomes as f32
        } else {
            0.0_f32
        }];
        let cvsd = &[if doc_median_info.on_autosomes > 1e-6 {
            rcvsds[0][0] / doc_median_info.on_autosomes as f32
        } else {
            0.0_f32
        }];

        record.push_format_float(b"CV", cv)?;
        record.push_format_float(b"CVSD", cvsd)?;

        writer.write(&record)?;
    }

    Ok(())
}

fn perform_collection(
    term: &Term,
    common_args: &CommonArgs,
    args: &Args,
    _regions: &Option<Vec<Interval>>,
) -> Result<(), anyhow::Error> {
    // Create output file writer and kick off processing.  This is done in its own block such
    // that the file is definitely closed when building the index below.
    let contigs = {
        let bam_reader = bam::IndexedReader::from_path(&args.path_input)?;
        build_chroms_bam(bam_reader.header(), None)?
    };

    let regions = if let Some(regions) = &args.regions {
        regions
            .iter()
            .map(|itv| itv.to_interval())
            .collect::<Vec<Interval>>()
    } else {
        contigs.clone()
    };

    let samples = samples_from_file(&args.path_input)?;

    // Write to temporary directory.
    term.write_line("Scan BAM file for coverage information; write results to temporary file.")?;
    let tmp_dir = tempdir()?;
    let tmp_path = tmp_dir.path();
    let tmp_out = tmp_path.join("tmp.bcf").to_str().unwrap().to_string();
    {
        let mut writer = build_bcf_writer(&tmp_out, &samples, &contigs)?;
        for region in &regions {
            process_region(term, common_args, args, region, &mut writer)?;
        }
    }

    term.write_line("Done scanning BAM. Will now compute per-contig coverage medians.")?;
    let doc_median_info = load_doc_median(&tmp_out)?;

    term.write_line("Done computing per-contig coverage medians. Building final coverage file.")?;
    perform_final_write(&tmp_out, &args.path_output, &doc_median_info)?;

    // Close temporary directory to we can handle any errors here.
    tmp_dir.close()?;

    Ok(())
}

pub fn run(term: &Term, common_args: &CommonArgs, args: &Args) -> Result<(), anyhow::Error> {
    term.write_line(&format!(
        "{}bam-collect-doc -- collect depth of coverage from BAM",
        Emoji("🧬 ", "")
    ))?;
    term.write_line(&format!("{}configuration:", Emoji("🔧 ", "")))?;
    term.write_line(&prefix_lines(
        "   ",
        &format!("common = {:#?}", &common_args),
    ))?;
    let regions: Option<Vec<Interval>> = args
        .regions
        .as_ref()
        .map(|regions| regions.iter().map(|itv| itv.to_interval()).collect());
    term.write_line(&prefix_lines("   ", &format!("sequencing = {:#?}", &args)))?;
    perform_collection(term, common_args, args, &regions)
}

#[cfg(test)]
mod tests {
    use crate::cmd_bam_collect_doc::DocConfig;

    use super::CountKind;
    use anyhow;
    use console::Term;
    use std::fs;
    use tempdir::TempDir;

    use pretty_assertions::assert_eq;

    /// Helper that runs `perform_collection()` and compares the result.
    fn _perform_collection_and_test(
        tmp_dir: &TempDir,
        path_input: &str,
        path_expected: &str,
        count_kind: CountKind,
    ) -> Result<(), anyhow::Error> {
        let path_output = String::from(tmp_dir.path().join("out.vcf").to_str().unwrap());
        let common_args = super::CommonArgs {
            verbose: clap_verbosity_flag::Verbosity::new(0, 0),
        };
        let args = super::Args {
            regions: None,
            path_reference: None,
            path_input: String::from(path_input),
            path_output: path_output.clone(),
            force_overwrite: false,
            doc_config: DocConfig {
                count_kind: count_kind,
                window_length: 100,
                min_mapq: None,
                min_unclipped: None,
            },
        };
        let term = Term::stderr();

        super::perform_collection(&term, &common_args, &args, &None)?;

        assert_eq!(
            fs::read_to_string(path_expected)?,
            fs::read_to_string(path_output)?
        );

        Ok(())
    }

    #[test]
    fn test_perform_collection_examples_fragments() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/bam-collect-doc-1/ex.sorted.bam",
            "./src/tests/data/bam-collect-doc-1/ex.expected.fragments.vcf",
            CountKind::FragmentCount,
        )?;
        Ok(())
    }

    #[test]
    fn test_perform_collection_examples_coverage() -> Result<(), anyhow::Error> {
        let tmp_dir = TempDir::new("tests")?;
        _perform_collection_and_test(
            &tmp_dir,
            "./src/tests/data/bam-collect-doc-1/ex.sorted.bam",
            "./src/tests/data/bam-collect-doc-1/ex.expected.coverage.vcf",
            CountKind::Coverage,
        )?;
        Ok(())
    }
}
