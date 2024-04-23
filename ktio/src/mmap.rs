use memmap2::{MmapMut, MmapOptions};
use std::{cell::UnsafeCell, fs::OpenOptions, ptr};

// https://stackoverflow.com/questions/65178245/how-do-i-write-to-a-mutable-slice-from-multiple-threads-at-arbitrary-indexes-wit/65182786#65182786
#[derive(Copy, Clone)]
pub struct MMWriter<'a, T> {
    slice: &'a [UnsafeCell<T>],
}
unsafe impl<'a, T: Send + Sync> Send for MMWriter<'a, T> {}
unsafe impl<'a, T: Send + Sync> Sync for MMWriter<'a, T> {}

impl<'a, T> MMWriter<'a, T> {
    pub fn new(slice: &'a mut [T]) -> Self {
        let ptr = slice as *mut [T] as *const [UnsafeCell<T>];
        Self {
            slice: unsafe { &*ptr },
        }
    }

    /// # Safety
    /// During multi threading, one must ensure i is mutually exclusive.
    /// It is undefined behaviour if two threads write to the same index without
    /// synchronization.
    pub unsafe fn write_at(&self, data: &[T], pos: usize)
    where
        T: Copy,
    {
        ptr::copy_nonoverlapping(data.as_ptr(), self.slice[pos].get(), data.len());
    }
}

pub fn mmap_file_for_writing(path: &str, size: usize) -> Result<MmapMut, String> {
    let file = OpenOptions::new()
        .read(true)
        .write(true)
        .truncate(true)
        .create(true)
        .open(path)
        .map_err(|_| format!("Could not mmap: {}", path))?;
    file.set_len(size as u64)
        .map_err(|_| format!("Could not mmap: {}", path))?;
    unsafe {
        MmapOptions::new()
            .map_mut(&file)
            .map_err(|_| format!("Could not mmap: {}", path))
    }
}
