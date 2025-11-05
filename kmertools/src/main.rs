use anyhow::Result;
use args::Cli;
use clap::Parser;
mod args;

#[cfg(not(tarpaulin_include))]
fn main() -> Result<()> {
    let cli = Cli::parse();
    args::cli(cli)
}
