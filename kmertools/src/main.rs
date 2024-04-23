use args::CompositionCommands;
use clap::Parser;
use composition::oligo::OligoComputer;
use coverage::CovComputer;
mod args;

fn main() {
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
                    args::Preset::Csv => com.set_delim(",".to_owned()),
                    args::Preset::Spc => com.set_delim(" ".to_owned()),
                    args::Preset::Tsv => com.set_delim("\t".to_owned()),
                }
                if let Err(e) = com.vectorise() {
                    println!("Error: {}", e);
                }
            }
            CompositionCommands::Cgr(_command) => {
                todo!("")
            }
        },
        args::Commands::Cov(command) => {
            let mut cov = CovComputer::new(
                command.input,
                command.output,
                command.k_size as usize,
                command.bin_size as usize,
                command.bin_size as usize,
            );
            if command.threads > 0 {
                cov.set_threads(command.threads);
            }
            if let Some(path) = command.alt_input {
                cov.set_kmer_path(path);
            }
            match command.preset {
                args::Preset::Csv => cov.set_delim(",".to_owned()),
                args::Preset::Spc => cov.set_delim(" ".to_owned()),
                args::Preset::Tsv => cov.set_delim("\t".to_owned()),
            }
            if let Err(e) = cov.vectorise() {
                println!("Error: {}", e);
            }
        }
    }
}
