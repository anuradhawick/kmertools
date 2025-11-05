pub mod error;
pub mod fops;
pub mod mmap;
pub mod seq;

#[cfg(test)]
mod error_tests;

// Re-export error types for convenience
pub use error::{KmertoolsError, Result};
