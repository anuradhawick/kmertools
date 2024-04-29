use std::{fs, path::Path};

pub fn delete_file_if_exists<P: AsRef<Path>>(path: P) -> std::io::Result<()> {
    let path = path.as_ref();
    if path.exists() {
        fs::remove_file(path)?; // Attempt to delete the file
    }
    Ok(())
}
