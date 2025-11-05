# Error Handling in kmertools

This document describes the error handling strategy used in kmertools, which utilizes `thiserror` and `anyhow` crates for consistent and ergonomic error management.

## Overview

kmertools uses a two-tier error handling approach:

1. **Library errors**: Defined using `thiserror` in the `ktio` crate
2. **Application errors**: Handled using `anyhow` in the CLI binary

## Core Error Types

The `ktio::error` module defines 5 core error types using `thiserror`:

### 1. IoError
**Purpose**: File I/O and stream operations  
**Usage**: Automatically converted from `std::io::Error` via `#[from]` attribute  
**Examples**:
- File not found
- Permission denied
- Write failures

```rust
use ktio::Result;
use std::fs::File;

fn open_file(path: &str) -> Result<File> {
    let file = File::open(path)?;  // Automatically converts io::Error
    Ok(file)
}
```

### 2. ParseError
**Purpose**: Sequence format and data parsing issues  
**Usage**: Created explicitly for parsing failures  
**Examples**:
- Invalid sequence format detection
- Malformed data in input files
- Numeric parsing failures

```rust
use ktio::{KmertoolsError, Result};

fn parse_format(data: &[u8]) -> Result<SeqFormat> {
    if data.is_empty() {
        return Err(KmertoolsError::ParseError(
            "Empty data stream".to_string()
        ));
    }
    // ... parsing logic
}
```

### 3. InvalidInputError
**Purpose**: Invalid input data such as bad nucleotides or malformed sequences  
**Usage**: Used when input data is structurally valid but semantically incorrect  
**Examples**:
- Invalid nucleotides (e.g., 'X' in DNA sequence)
- Sequences that don't meet requirements
- Out-of-range values

```rust
use ktio::{KmertoolsError, Result};

fn validate_nucleotide(nuc: u8) -> Result<()> {
    match nuc {
        b'A' | b'C' | b'G' | b'T' => Ok(()),
        _ => Err(KmertoolsError::InvalidInputError(
            format!("Invalid nucleotide '{}'", nuc as char)
        )),
    }
}
```

### 4. ConfigError
**Purpose**: Invalid configuration parameters and settings  
**Usage**: Validation of user-provided parameters  
**Examples**:
- Invalid k-mer sizes
- Conflicting options
- Out-of-range parameters

```rust
use ktio::{KmertoolsError, Result};

fn validate_ksize(k: usize) -> Result<()> {
    if k < 3 || k > 31 {
        return Err(KmertoolsError::ConfigError(
            format!("K-mer size {} out of valid range [3, 31]", k)
        ));
    }
    Ok(())
}
```

### 5. SystemError
**Purpose**: System resource issues like memory mapping and allocation  
**Usage**: Low-level system operation failures  
**Examples**:
- Memory mapping failures
- Resource allocation errors
- System limit exceeded

```rust
use ktio::{KmertoolsError, Result};
use memmap2::MmapMut;

fn create_mmap(path: &str, size: usize) -> Result<MmapMut> {
    // ... file creation
    unsafe {
        MmapOptions::new()
            .map_mut(&file)
            .map_err(|e| KmertoolsError::SystemError(
                format!("Memory mapping failed for {}: {}", path, e)
            ))
    }
}
```

## Usage Patterns

### In Library Code

Library functions should return `ktio::Result<T>`:

```rust
use ktio::Result;

pub fn process_sequences(path: &str) -> Result<()> {
    let reader = ktio::seq::get_reader(path)?;
    // ... processing logic
    Ok(())
}
```

### In Application Code

The CLI uses `anyhow::Result<T>` for top-level error handling:

```rust
use anyhow::{Context, Result};

fn main() -> Result<()> {
    let cli = Cli::parse();
    
    process_data()
        .context("Failed to process data")?;
    
    Ok(())
}
```

### Error Conversion

The `?` operator automatically converts between compatible error types:

```rust
use anyhow::Result;
use ktio;

fn cli_handler() -> Result<()> {
    // ktio::Result is automatically converted to anyhow::Result
    ktio::seq::get_reader("file.fa")?;
    Ok(())
}
```

## Benefits

1. **Type Safety**: Distinct error types prevent confusion
2. **Ergonomic**: Uses Rust's `?` operator naturally
3. **Informative**: Clear error messages with context
4. **Maintainable**: Centralized error definitions
5. **Composable**: Easy to add context at each level

## Migration Notes

When migrating code from `Result<T, String>` to the new error types:

1. Change function signatures to use `ktio::Result<T>`
2. Replace `.map_err(|_| String::from(...))` with appropriate error types
3. Use `?` operator for automatic error conversion
4. Add `.context()` in CLI code for additional error information

## Testing

Error types can be tested like any other type:

```rust
#[test]
fn test_invalid_input_error() {
    let err = KmertoolsError::InvalidInputError("Bad nucleotide 'X'".to_string());
    assert_eq!(err.to_string(), "Invalid input: Bad nucleotide 'X'");
}
```

## Future Improvements

Potential enhancements to consider:

1. Add error codes for programmatic error handling
2. Include source location information for debugging
3. Add structured error data for better error analysis
4. Implement custom error reporters for different output formats
