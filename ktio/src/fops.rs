use std::{fs, io, path::Path};

pub fn delete_file_if_exists<P: AsRef<Path>>(path: P) -> std::io::Result<()> {
    let path = path.as_ref();
    if path.exists() {
        fs::remove_file(path)?; // Attempt to delete the file
    }
    Ok(())
}

pub fn create_directory<P: AsRef<Path>>(path: P) -> io::Result<()> {
    fs::create_dir_all(path)
}

pub fn load_lines_sorted<P: AsRef<Path>>(path: P) -> Vec<String> {
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

#[test]
fn delete_file_if_exists_test() {
    assert!(delete_file_if_exists("../test_data/doesnotexist.txt").is_ok());
}

#[test]
fn create_directory_test() {
    assert!(create_directory("../test_data/madedirectory").is_ok());
}
