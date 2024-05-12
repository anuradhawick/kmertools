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

// Presets for vector outputs
#[derive(Debug, ValueEnum, Clone)]
pub enum VecFmtPreset {
    /// Comma separated format
    Csv,
    /// Tab separated format
    Tsv,
    /// Space separated format
    Spc,
}

// Presets for minimiser outputs
#[derive(Debug, ValueEnum, Clone)]
pub enum MinFmtPreset {
    /// Conver sequences into minimiser representation
    S2m,
    /// Group sequences by minimiser
    M2s,
}

/// Subcommands available
#[derive(Debug, Subcommand)]
pub enum Commands {
    /// Generate sequence composition based features
    Comp {
        #[clap(subcommand)]
        command: CompositionCommands,
    },
    /// Generates coverage histogram based on the reads
    Cov(CoverageCommand),
    /// Bin reads using minimisers
    Min(MinimiserCommand),
}

// COMPOSITION
#[derive(Debug, Subcommand)]
pub enum CompositionCommands {
    /// Generate oligonucleotide frequency vectors
    Oligo(OligoCommand),
    // /// Generates Chaos Game Representation
    // Cgr(CGRCommand),
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
    #[arg(short, long, default_value_t = 3)]
    pub k_size: usize,

    /// Output type to write
    #[clap(value_enum, short, long, default_value_t = VecFmtPreset::Spc)]
    pub preset: VecFmtPreset,

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
    #[arg(short, long, default_value_t = 3)]
    pub k_size: usize,

    /// Output type to write
    #[clap(value_enum, short, long, default_value_t = VecFmtPreset::Spc)]
    pub preset: VecFmtPreset,

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
    #[clap(value_enum, short, long, default_value_t = VecFmtPreset::Spc)]
    pub preset: VecFmtPreset,

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

// MINIMISERS
#[derive(Debug, Args)]
pub struct MinimiserCommand {
    /// Input file path
    #[arg(short, long)]
    pub input: String,

    /// Output vectors path
    #[arg(short, long)]
    pub output: String,

    /// Minimiser size
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(7..=28), default_value_t = 10)]
    pub m_size: u64,

    /// Window size
    ///
    /// 0 - emits one minimiser per sequence (useful for sequencing reads)
    /// w_size must be longer than m_size
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(0..), verbatim_doc_comment, default_value_t = 0)]
    pub w_size: u64,

    /// Output type to write
    #[clap(value_enum, short, long, default_value_t = MinFmtPreset::S2m)]
    pub preset: MinFmtPreset,

    /// Thread count for computations 0=auto
    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,
}
