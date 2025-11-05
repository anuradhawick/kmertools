#[cfg(test)]
mod tests {
    use super::super::error::{KmertoolsError, Result};
    use std::io;

    #[test]
    fn test_io_error_conversion() {
        let io_err = io::Error::new(io::ErrorKind::NotFound, "file not found");
        let kmertools_err: KmertoolsError = io_err.into();
        
        match kmertools_err {
            KmertoolsError::IoError(_) => (),
            _ => panic!("Expected IoError variant"),
        }
    }

    #[test]
    fn test_parse_error() {
        let err = KmertoolsError::ParseError("Invalid format".to_string());
        assert_eq!(err.to_string(), "Parse error: Invalid format");
    }

    #[test]
    fn test_invalid_input_error() {
        let err = KmertoolsError::InvalidInputError("Bad nucleotide 'X'".to_string());
        assert_eq!(err.to_string(), "Invalid input: Bad nucleotide 'X'");
    }

    #[test]
    fn test_config_error() {
        let err = KmertoolsError::ConfigError("Invalid k-mer size".to_string());
        assert_eq!(err.to_string(), "Configuration error: Invalid k-mer size");
    }

    #[test]
    fn test_system_error() {
        let err = KmertoolsError::SystemError("Memory allocation failed".to_string());
        assert_eq!(err.to_string(), "System error: Memory allocation failed");
    }

    #[test]
    fn test_result_type() {
        fn returns_error() -> Result<i32> {
            Err(KmertoolsError::ParseError("test".to_string()))
        }

        let result = returns_error();
        assert!(result.is_err());
    }

    #[test]
    fn test_result_ok() {
        fn returns_ok() -> Result<i32> {
            Ok(42)
        }

        let result = returns_ok();
        assert_eq!(result.unwrap(), 42);
    }
}
