use clap::{Args, Parser, Subcommand, ValueEnum};

const ABOUT: &str = "kmertools: DNA vectorisation

k-mer based vectorisation for DNA sequences for
metagenomics and AI/ML applications";

/// kmertools: DNA vectorisation
#[derive(Debug, Parser)]
#[command(author, version, about, long_about = ABOUT)]
pub struct Cli {
    #[command(subcommand)]
    pub command: Commands,
    // /// Turn off all loggin information
    // #[arg(short, long)]
    // pub quiet: bool,
}

// COMMON
#[derive(Debug, ValueEnum, Clone)]
pub enum Preset {
    /// Comma separated format
    Csv,
    /// Tab separated format
    Tsv,
    /// Space separated format
    Spc,
}

/// Subcommands available
#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Processes the input file and outputs vectors
    Comp {
        #[clap(subcommand)]
        command: CompositionCommands,
    },
    /// Generates coverage histogram based on the reads
    Cov(CoverageCommand),
}

// COMPOSITION
#[derive(Debug, Subcommand)]
pub enum CompositionCommands {
    /// Generate oligonucleotide frequency vectors
    Oligo(OligoCommand),
    /// Generates Chaos Game Representation
    Cgr(CGRCommand),
}

#[derive(Debug, Args)]
pub struct OligoCommand {
    /// Input file path
    #[arg(short, long)]
    pub input: String,

    /// Output vectors path
    #[arg(short, long)]
    pub output: String,

    /// Disable normalisation and output raw counts
    #[arg(short, long)]
    pub counts: bool,

    /// Set k-mer size
    #[arg(long, default_value_t = 3)]
    pub k_size: usize,

    /// Output type to write
    ///     csv: comma separated
    ///     tsv: tab separated
    ///     spc: space separated
    #[clap(value_enum, short, long, verbatim_doc_comment, default_value_t = Preset::Spc)]
    pub preset: Preset,

    /// Thread count for computations 0=auto
    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,
}

#[derive(Debug, Args)]
pub struct CGRCommand {
    /// Input file path
    #[arg(short, long)]
    pub input: String,

    /// Output vectors path
    #[arg(short, long)]
    pub output: String,

    /// Disable normalisation and output raw counts
    #[arg(short, long)]
    pub counts: bool,

    /// Set k-mer size
    #[arg(long, default_value_t = 3)]
    pub k_size: usize,

    /// Output type to write
    ///     csv: comma separated
    ///     tsv: tab separated
    ///     spc: space separated
    #[clap(value_enum, short, long, verbatim_doc_comment, default_value_t = Preset::Spc)]
    pub preset: Preset,

    /// Thread count for computations 0=auto
    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,
}

// COVERAGE
#[derive(Debug, Args)]
pub struct CoverageCommand {
    /// Input file path
    #[arg(short, long)]
    pub input: String,

    /// Input file path, for k-mer counting
    #[arg(short, long)]
    pub alt_input: Option<String>,

    /// Output vectors path
    #[arg(short, long)]
    pub output: String,

    /// K size for the coverage histogram
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(7..=31), default_value_t = 15)]
    pub k_size: u64,

    /// Output type to write
    ///     csv: comma separated
    ///     tsv: tab separated
    ///     spc: space separated
    #[clap(value_enum, short, long, verbatim_doc_comment, default_value_t = Preset::Spc)]
    pub preset: Preset,

    /// Bin size for the coverage histogram
    #[arg(short = 's', long = "bin-size", value_parser = clap::value_parser!(u64).range(5..), default_value_t = 16)]
    pub bin_size: u64,

    /// Number of bins for the coverage histogram
    #[arg(short = 'c', long = "bin-count", value_parser = clap::value_parser!(u64).range(5..), default_value_t = 16)]
    pub bin_count: u64,

    /// Disable normalisation and output raw counts
    #[arg(long)]
    pub counts: bool,

    /// Thread count for computations 0=auto
    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,
}
