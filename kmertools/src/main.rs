use clap::Parser;

mod args;

fn main() {
    let cli = args::Cli::parse();

    println!("{:#?}", cli.command);
}
