use args::CompositionCommands;
use clap::Parser;
use composition::oligo::OligoComputer;
mod args;

fn main() {
    let cli = args::Cli::parse();

    match cli.command {
        args::Commands::Comp { command } => match command {
            CompositionCommands::Oligo(command) => {
                let mut com = OligoComputer::new(&command.input, &command.output, command.k_size);
                if command.threads > 0 {
                    com.set_threads(command.threads);
                }
                if command.counts {
                    com.set_norm(false);
                }
                if let Err(e) = com.vectorise() {
                    println!("Error: {}", e);
                }
            }
            CompositionCommands::Cgr(command) => {}
        },
        args::Commands::Cov(command) => {}
    }
}
