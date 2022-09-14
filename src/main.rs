/// Main entry point for maelstrom-core
mod cmd_bam_collect_doc;
mod common;
mod err;

use clap::{Parser, Subcommand};
use console::{Emoji, Term};

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    /// Commonly used arguments
    #[clap(flatten)]
    common: common::Args,

    /// The sub command to run
    #[clap(subcommand)]
    command: Commands,
}

#[allow(clippy::large_enum_variant)]
#[derive(Subcommand, Debug)]
enum Commands {
    /// Create contigs with synthetic sequence
    BamCollectDoc(cmd_bam_collect_doc::Args),
}

fn main() -> Result<(), anyhow::Error> {
    let cli = Cli::parse();

    let term = Term::stderr();
    match &cli.command {
        Commands::BamCollectDoc(args) => {
            cmd_bam_collect_doc::run(&term, &cli.common, args)?;
        }
    }
    term.write_line(&format!("All done. Have a nice day!{}", Emoji(" 😃", "")))?;

    Ok(())
}
