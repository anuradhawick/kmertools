use thiserror::Error;

/// Core error types for kmertools operations.
/// 
/// This enum defines the 5 main error categories used throughout kmertools:
/// - IoError: File I/O and stream operations
/// - ParseError: Sequence format and data parsing
/// - InvalidInputError: Invalid sequences or nucleotides
/// - ConfigError: Invalid parameters and configuration
/// - SystemError: Memory mapping and resource allocation
#[derive(Error, Debug)]
pub enum KmertoolsError {
    /// Error during file I/O operations (reading, writing, opening files)
    #[error("IO error: {0}")]
    IoError(#[from] std::io::Error),

    /// Error parsing sequence data or formats
    #[error("Parse error: {0}")]
    ParseError(String),

    /// Error due to invalid input data (e.g., bad nucleotides)
    #[error("Invalid input: {0}")]
    InvalidInputError(String),

    /// Error due to invalid configuration or parameters
    #[error("Configuration error: {0}")]
    ConfigError(String),

    /// Error related to system resources (memory mapping, allocation)
    #[error("System error: {0}")]
    SystemError(String),
}

/// Type alias for Results using KmertoolsError
pub type Result<T> = std::result::Result<T, KmertoolsError>;
