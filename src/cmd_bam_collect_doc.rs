/// Command to collect depth of coverage information from BAM file.
use bio_types::genome::Interval;
use clap::Args as ClapArgs;
use console::{Emoji, Term};

use crate::common::Args as CommonArgs;
use crate::common::{prefix_lines, ArgInterval};

/// Enum for selecting sequencing technology
#[derive(PartialEq, Eq, Debug, Clone, clap::ValueEnum, enum_utils::FromStr)]
pub enum CountKind {
    /// Count fragments starting in bin
    FragmentCount,
    /// Consider base-wise coverage
    Coverage,
}

/// Command line options
#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Optional list of regions to call.
    #[clap(short = 'r', long = "regions")]
    regions: Option<Vec<ArgInterval>>,
    /// Path to input file.
    #[clap(short = 'i', long = "in", required = true)]
    path_input: String,
    /// Path to output file.
    #[clap(short = 'o', long = "out", required = true)]
    path_output: String,
    /// Force overwriting output file.
    #[clap(short = 'f', long = "force", action, default_value_t = false)]
    force_overwrite: bool,
    /// The coverage kind to count
    #[clap(short = 'k', long = "count-kind", value_enum, default_value_t = CountKind::Coverage)]
    count_kind: CountKind,
}

fn perform_collection(
    _term: &Term,
    _common_args: &CommonArgs,
    _args: &Args,
    _regions: &Option<Vec<Interval>>,
) -> Result<(), anyhow::Error> {
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
        .map(|regions| regions.iter().map(|itv| itv.interval()).collect());
    term.write_line(&prefix_lines("   ", &format!("sequencing = {:#?}", &args)))?;
    perform_collection(term, common_args, args, &regions)
}

#[cfg(test)]
mod tests {
    use super::CountKind;
    use anyhow;
    use console::Term;
    use file_diff::diff;
    use tempdir::TempDir;

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
            path_input: String::from(path_input),
            path_output: path_output.clone(),
            force_overwrite: false,
            count_kind: count_kind,
        };
        let term = Term::stderr();

        super::perform_collection(&term, &common_args, &args, &None)?;

        assert!(!diff(&path_expected, &path_output));

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
