use args::CompositionCommands;
use clap::Parser;
use composition::oligo::OligoComputer;
use coverage::CovComputer;
use misc::minimisers;
mod args;

#[cfg(not(tarpaulin_include))]
fn main() {
    use ktio::fops::create_directory;

    let cli = args::Cli::parse();

    match cli.command {
        args::Commands::Comp { command } => match command {
            CompositionCommands::Oligo(command) => {
                let mut com = OligoComputer::new(command.input, command.output, command.k_size);
                if command.threads > 0 {
                    com.set_threads(command.threads);
                }
                if command.counts {
                    com.set_norm(false);
                }
                match command.preset {
                    args::VecFmtPreset::Csv => com.set_delim(",".to_owned()),
                    args::VecFmtPreset::Spc => com.set_delim(" ".to_owned()),
                    args::VecFmtPreset::Tsv => com.set_delim("\t".to_owned()),
                }
                if let Err(e) = com.vectorise() {
                    println!("Error: {}", e);
                }
            }
        },
        args::Commands::Cov(command) => {
            create_directory(&command.output).unwrap();
            let mut cov = CovComputer::new(
                command.input,
                command.output,
                command.k_size,
                command.bin_size,
                command.bin_size,
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
                args::VecFmtPreset::Csv => cov.set_delim(",".to_owned()),
                args::VecFmtPreset::Spc => cov.set_delim(" ".to_owned()),
                args::VecFmtPreset::Tsv => cov.set_delim("\t".to_owned()),
            }
            cov.build_table().unwrap();
            cov.compute_coverages();
        }
        args::Commands::Min(command) => {
            if command.w_size <= command.m_size && command.w_size > 0 {
                println!("Window size must be longer than minimiser size!");
                return;
            }
            if command.m_size >= 31 {
                println!("Minimisers longer than 30 bases not allowed!");
                return;
            }

            match command.preset {
                args::MinFmtPreset::M2s => minimisers::bin_sequences(
                    command.w_size as usize,
                    command.m_size as usize,
                    &command.input,
                    &command.output,
                    command.threads,
                ),
                args::MinFmtPreset::S2m => minimisers::seq_to_min(
                    command.w_size as usize,
                    command.m_size as usize,
                    &command.input,
                    &command.output,
                    command.threads,
                ),
            }
        }
        args::Commands::Ctr(command) => {
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
