pub struct CgrComputer<'a> {
    in_path: &'a str,
    out_path: &'a str,
    ksize: usize,
    kcount: usize,
    threads: usize,
    pos_map: Vec<usize>,
    norm: bool,
    delim: String,
}

impl<'a> CgrComputer<'a> {
    pub fn new() -> Self {
        // Self {}
        todo!()
    }

    pub fn set_threads(&mut self, threads: usize) {
        self.threads = threads;
    }

    pub fn set_norm(&mut self, norm: bool) {
        self.norm = norm;
    }

    pub fn set_delim(&mut self, delim: String) {
        self.delim = delim;
    }

    pub fn vectorise(&self) -> Result<(), String> {
        if self.in_path == "-" {
            return self.vectorise_batch();
        }
        self.vectorise_mmap()
    }

    fn vectorise_batch(&self) -> Result<(), String> {
        todo!()
    }

    fn vectorise_mmap(&self) -> Result<(), String> {
        todo!()
    }

    fn vectorise_one(&self, seq: &str) -> Vec<f64> {
        todo!()
    }
}
