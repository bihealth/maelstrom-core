/// Commonly shared data structures
pub mod bam;
pub mod bcf;
pub mod doc;
pub mod stats;

use std::str::FromStr;

use bio_types::genome::Interval as BioInterval;
use bio_types::genome::Position;
use clap::Args as ClapArgs;
use clap_verbosity_flag::Verbosity;

use crate::err::ArgError;

#[derive(ClapArgs, Debug)]
pub struct Args {
    /// Verbosity of the program
    #[clap(flatten)]
    pub verbose: Verbosity,
}

/// Helper function to prepend lines in string
pub fn prefix_lines(prefix: &str, text: &str) -> String {
    let lines: Vec<String> = text
        .lines()
        .map(|line| format!("{}{}", prefix, line))
        .collect();
    lines.join("\n")
}

/// Interval as parsed from the command line
#[derive(Debug, PartialEq, Eq, Clone, Hash)]
pub struct ArgInterval {
    internal: BioInterval,
}

impl ArgInterval {
    // pub fn new(&self, interval: &BioInterval) -> Self {
    //     return ArgInterval {
    //         internal: interval.clone(),
    //     };
    // }

    pub fn to_interval(&self) -> BioInterval {
        self.internal.clone()
    }
}

impl FromStr for ArgInterval {
    type Err = ArgError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let split = s.split(':');
        let vec = split.collect::<Vec<&str>>();
        if vec.len() != 2 {
            return Err(ArgError::IntervalInvalidFormat);
        }

        let chrom = vec[0];
        let split = vec[1].split('-');
        let vec = split.collect::<Vec<&str>>();
        if vec.len() != 2 {
            return Err(ArgError::IntervalInvalidFormat);
        }
        let start = vec[0].replace(',', "").parse::<Position>()?;
        let end = vec[1].replace(',', "").parse::<Position>()?;

        Ok(ArgInterval {
            internal: BioInterval::new(chrom.to_owned(), std::ops::Range { start, end }),
        })
    }
}
