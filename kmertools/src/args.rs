use clap::{Args, Parser, Subcommand};

/// Main command line tool
#[derive(Debug, Parser)]
#[command(author, version, about, long_about = None)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,

    /// Global flag to make logging silent
    #[arg(short, long)]
    pub quiet: bool,
}

/// Subcommands available
#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Processes the input file and outputs vectors
    Composition(CompositionCommand),
    /// Generates coverage histogram based on the reads
    Coverage(CoverageCommand),
}

#[derive(Debug, Args)]
pub struct CompositionCommand {
    /// Input file path
    #[arg(short, long)]
    file: String,

    /// Output vectors path
    #[arg(short, long)]
    output: String,

    /// Output type, should be one of csv, tsv, or json
    #[arg(short, long, default_value = "csv")]
    preset: String,

    /// Set k-mer size
    #[arg(short, long, default_value_t = 3)]
    k_size: usize,

    /// Set thread count
    #[arg(short, long, default_value_t = 8)]
    threads: usize,
}

#[derive(Debug, Args)]
pub struct CoverageCommand {
    /// Reads path for binning
    #[arg(short, long)]
    reads_path: String,

    /// K size for the coverage histogram
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(7..=31))]
    k_size: usize,

    /// Bin size for the coverage histogram
    #[arg(short = 's', long = "bin-size", default_value_t = 32)]
    bin_size: usize,

    /// Number of bins for the coverage histogram
    #[arg(short = 'c', long = "bin-count", default_value_t = 16)]
    bin_count: usize,

    /// Thread count for computations
    #[arg(short, long, default_value_t = 8)]
    threads: usize,

    /// Output file name
    #[arg(short, long)]
    output: String,
}
