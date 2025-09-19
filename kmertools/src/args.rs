use clap::{Args, Parser, Subcommand, ValueEnum};
use composition::{cgr::CgrComputer, oligo::OligoComputer, oligocgr::OligoCgrComputer};
use coverage::CovComputer;
use ktio::fops::create_directory;
use misc::minimisers;

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
    /// Count k-mers
    Ctr(CounterCommand),
}

// COMPOSITION
#[derive(Debug, Subcommand)]
pub enum CompositionCommands {
    /// Generate oligonucleotide frequency vectors
    Oligo(OligoCommand),
    /// Generates Chaos Game Representations
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
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(3..=7), default_value_t = 3)]
    pub k_size: u64,

    /// Raw counts
    #[arg(short, long)]
    pub raw_count: bool,

    /// Output type to write
    #[clap(value_enum, short, long, default_value_t = VecFmtPreset::Spc)]
    pub preset: VecFmtPreset,

    /// Include header (with k-mer in ACGT format)
    #[clap(value_enum, short = 'H', long)]
    pub header: bool,

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

    /// Disable normalisation and output raw counts (only with k-mer mode)
    #[arg(short, long)]
    pub counts: bool,

    /// Set k-mer size or default to full sequence CGR
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(3..=7))]
    pub k_size: Option<u64>,

    /// Set vector size (output will be a square matrix with N=vecsize)
    #[arg(short, long)]
    pub vec_size: Option<u64>,

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

    /// Output directory path
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

    /// Max memory in GB
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(6..=128), default_value_t = 6)]
    pub memory: u64,

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

// COUNTER
#[derive(Debug, Args)]
pub struct CounterCommand {
    /// Input file path
    #[arg(short, long)]
    pub input: String,

    /// Output directory path
    #[arg(short, long)]
    pub output: String,

    /// k size for counting
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(10..32))]
    pub k_size: u64,

    /// Max memory in GB
    #[arg(short, long, value_parser = clap::value_parser!(u64).range(6..=128), default_value_t = 6)]
    pub memory: u64,

    /// Output ACGT instead of numeric values
    ///
    /// This requires a larger space for the final result
    /// compared to the compact numeric representation
    #[arg(short, long, verbatim_doc_comment)]
    pub acgt: bool,

    /// Thread count for computations 0=auto
    #[arg(short, long, default_value_t = 0)]
    pub threads: usize,
}

#[cfg(not(tarpaulin_include))]
pub fn cli(cli: Cli) {
    match cli.command {
        Commands::Comp { command } => match command {
            CompositionCommands::Oligo(command) => {
                let mut com = OligoComputer::new(
                    command.input,
                    command.output,
                    command.k_size as usize,
                    !command.raw_count,
                );
                if command.threads > 0 {
                    com.set_threads(command.threads);
                }
                com.set_norm(!command.counts);
                com.set_header(command.header);

                match command.preset {
                    VecFmtPreset::Csv => com.set_delim(",".to_owned()),
                    VecFmtPreset::Spc => com.set_delim(" ".to_owned()),
                    VecFmtPreset::Tsv => com.set_delim("\t".to_owned()),
                }
                if let Err(e) = com.vectorise() {
                    eprintln!("Error: {}", e);
                }
            }
            CompositionCommands::Cgr(command) => {
                if let Some(ksize) = command.k_size {
                    let vecsize = command
                        .vec_size
                        .unwrap_or((ksize as f64).powf(4.0).powf(0.5) as u64)
                        as usize;
                    let mut cgr = OligoCgrComputer::new(
                        command.input,
                        command.output,
                        ksize as usize,
                        vecsize,
                    );
                    if command.threads > 0 {
                        cgr.set_threads(command.threads);
                    }
                    cgr.set_norm(!command.counts);
                    if let Err(e) = cgr.vectorise() {
                        eprintln!("Error: {}", e);
                    }
                } else {
                    if command.counts {
                        eprintln!("Error: cannot use counts in whole sequence CGR!");
                        return;
                    }
                    let vecsize = command.vec_size.unwrap_or(1) as usize;
                    let mut cgr = CgrComputer::new(command.input, command.output, vecsize);
                    if command.threads > 0 {
                        cgr.set_threads(command.threads);
                    }
                    if let Err(e) = cgr.vectorise() {
                        eprintln!("Error: {}", e);
                    }
                }
            }
        },
        Commands::Cov(command) => {
            create_directory(&command.output).unwrap();
            let mut cov = CovComputer::new(
                command.input,
                command.output,
                command.k_size as usize,
                command.bin_size as usize,
                command.bin_count as usize,
            );
            if command.threads > 0 {
                cov.set_threads(command.threads);
            }
            if let Some(path) = command.alt_input {
                cov.set_kmer_path(path);
            }
            if command.counts {
                cov.set_norm(false);
            }
            cov.set_max_memory(command.memory as f64);
            match command.preset {
                VecFmtPreset::Csv => cov.set_delim(",".to_owned()),
                VecFmtPreset::Spc => cov.set_delim(" ".to_owned()),
                VecFmtPreset::Tsv => cov.set_delim("\t".to_owned()),
            }
            cov.build_table().unwrap();
            cov.compute_coverages();
        }
        Commands::Min(command) => {
            if command.w_size <= command.m_size && command.w_size > 0 {
                eprintln!("Window size must be longer than minimiser size!");
                return;
            }
            if command.m_size >= 31 {
                eprintln!("Minimisers longer than 30 bases not allowed!");
                return;
            }

            match command.preset {
                MinFmtPreset::M2s => minimisers::bin_sequences(
                    command.w_size as usize,
                    command.m_size as usize,
                    &command.input,
                    &command.output,
                    command.threads,
                ),
                MinFmtPreset::S2m => minimisers::seq_to_min(
                    command.w_size as usize,
                    command.m_size as usize,
                    &command.input,
                    &command.output,
                    command.threads,
                ),
            }
        }
        Commands::Ctr(command) => {
            create_directory(&command.output).unwrap();
            let mut ctr =
                counter::CountComputer::new(command.input, command.output, command.k_size as usize);
            if command.threads > 0 {
                ctr.set_threads(command.threads);
            }
            if command.acgt {
                ctr.set_acgt_output(true);
            }
            ctr.set_max_memory(command.memory as f64);
            ctr.count();
            ctr.merge(true);
        }
    }
}
