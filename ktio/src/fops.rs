use std::{fs, path::Path};

pub fn delete_file_if_exists<P: AsRef<Path>>(path: P) -> std::io::Result<()> {
    let path = path.as_ref();
    if path.exists() {
        fs::remove_file(path)?; // Attempt to delete the file
    }
    Ok(())
}

pub fn load_lines_sorted(path: &str) -> Vec<String> {
    let data = fs::read(path).unwrap();
    let text = String::from_utf8(data).unwrap().trim().to_string();
    let mut arr: Vec<String> = text
        .trim()
        .split('\n')
        .map(|s| s.trim().to_string())
        .collect();
    arr.sort();
    arr
}
