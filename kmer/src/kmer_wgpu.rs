use crate::kmer::KmerGenerator;
use crate::Kmer;

const WGSL: &str = r#"
@group(0) @binding(0) var<storage, read> seq: array<u32>;
@group(0) @binding(1) var<storage, read_write> out_fwd: array<u32>;
@group(0) @binding(2) var<storage, read_write> out_rev: array<u32>;
@group(0) @binding(3) var<storage, read_write> out_valid: array<u32>;

struct Params {
  len: u32,
  k: u32,
}
@group(0) @binding(4) var<uniform> params: Params;

@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let i = gid.x;
  if (params.k == 0u || i + params.k > params.len) { return; }

  var fwd: u32 = 0u;
  var rev: u32 = 0u;
  let shift: u32 = (params.k - 1u) * 2u;

  for (var j: u32 = 0u; j < params.k; j = j + 1u) {
    let b = seq[i + j];
    if (b > 3u) {
      out_valid[i] = 0u;
      return;
    }
    fwd = (fwd << 2u) | b;
    rev = (rev >> 2u) | ((b ^ 3u) << shift);
  }

  out_fwd[i] = fwd;
  out_rev[i] = rev;
  out_valid[i] = 1u;
}
"#;

pub struct KmerGeneratorWGPU<'a> {
    seq: &'a [u8],
    ksize: usize,
    cursor: usize,
    buffer: Vec<(Kmer, Kmer)>,
    done: bool,
    batch_size: usize,
}

impl<'a> KmerGeneratorWGPU<'a> {
    pub fn new(seq: &'a [u8], ksize: usize) -> Self {
        Self::with_batch_size(seq, ksize, 1 << 20)
    }

    pub fn with_batch_size(seq: &'a [u8], ksize: usize, batch_size: usize) -> Self {
        Self { seq, ksize, cursor: 0, buffer: Vec::new(), done: false, batch_size: batch_size.max(ksize.max(1)) }
    }

    fn process_batch(&mut self) {
        if self.done { return; }
        if self.cursor >= self.seq.len() {
            self.done = true;
            return;
        }

        let start = self.cursor;
        let mut end = (start + self.batch_size).min(self.seq.len());
        if end < self.seq.len() && self.ksize > 1 {
            end = (end + self.ksize - 1).min(self.seq.len());
        }
        let slice = &self.seq[start..end];

        self.buffer = compute_kmers_wgpu(slice, self.ksize)
            .unwrap_or_else(|_| KmerGenerator::new(slice, self.ksize).collect());

        self.cursor = (start + self.batch_size).min(self.seq.len());
        if self.cursor >= self.seq.len() { self.done = true; }
    }
}

fn encode_base(b: u8) -> u32 {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 4,
    }
}

fn compute_kmers_wgpu(seq: &[u8], ksize: usize) -> Result<Vec<(Kmer, Kmer)>, String> {
    if ksize > 16 {
        return Err(format!(
            "WGPU backend does not support ksize {} (maximum supported ksize is 16)",
            ksize
        ));
    }
    if ksize == 0 || seq.len() < ksize {
        return Ok(Vec::new());
    }

    pollster::block_on(async move {
        let instance = wgpu::Instance::default();
        let adapter = instance
            .request_adapter(&wgpu::RequestAdapterOptions::default())
            .await
            .ok_or("No GPU adapter")?;
        let (device, queue) = adapter
            .request_device(&wgpu::DeviceDescriptor::default(), None)
            .await
            .map_err(|e| e.to_string())?;

        let in_data: Vec<u32> = seq.iter().map(|&b| encode_base(b)).collect();
        let out_len = seq.len() - ksize + 1;

        use wgpu::util::DeviceExt;
        let seq_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: Some("seq"),
            contents: bytemuck::cast_slice(&in_data),
            usage: wgpu::BufferUsages::STORAGE,
        });
        let out_fwd = device.create_buffer(&wgpu::BufferDescriptor { label: Some("out_fwd"), size: (out_len * 4) as u64, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false });
        let out_rev = device.create_buffer(&wgpu::BufferDescriptor { label: Some("out_rev"), size: (out_len * 4) as u64, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false });
        let out_valid = device.create_buffer(&wgpu::BufferDescriptor { label: Some("out_valid"), size: (out_len * 4) as u64, usage: wgpu::BufferUsages::STORAGE | wgpu::BufferUsages::COPY_SRC, mapped_at_creation: false });

        #[repr(C)]
        #[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
        struct Params { len: u32, k: u32 }
        let params = Params { len: seq.len() as u32, k: ksize as u32 };
        let params_buf = device.create_buffer_init(&wgpu::util::BufferInitDescriptor { label: Some("params"), contents: bytemuck::bytes_of(&params), usage: wgpu::BufferUsages::UNIFORM });

        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor { label: Some("kmer shader"), source: wgpu::ShaderSource::Wgsl(WGSL.into()) });
        let bgl = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("bgl"),
            entries: &[
                wgpu::BindGroupLayoutEntry { binding: 0, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: true }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                wgpu::BindGroupLayoutEntry { binding: 1, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                wgpu::BindGroupLayoutEntry { binding: 2, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                wgpu::BindGroupLayoutEntry { binding: 3, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Storage { read_only: false }, has_dynamic_offset: false, min_binding_size: None }, count: None },
                wgpu::BindGroupLayoutEntry { binding: 4, visibility: wgpu::ShaderStages::COMPUTE, ty: wgpu::BindingType::Buffer { ty: wgpu::BufferBindingType::Uniform, has_dynamic_offset: false, min_binding_size: None }, count: None },
            ],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor { label: Some("pl"), bind_group_layouts: &[&bgl], push_constant_ranges: &[] });
        let pipeline = device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor { label: Some("cp"), layout: Some(&pipeline_layout), module: &shader, entry_point: "main", compilation_options: Default::default() });
        let bg = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("bg"), layout: &bgl,
            entries: &[
                wgpu::BindGroupEntry { binding: 0, resource: seq_buf.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 1, resource: out_fwd.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 2, resource: out_rev.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 3, resource: out_valid.as_entire_binding() },
                wgpu::BindGroupEntry { binding: 4, resource: params_buf.as_entire_binding() },
            ],
        });

        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor::default());
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor::default());
            pass.set_pipeline(&pipeline);
            pass.set_bind_group(0, &bg, &[]);
            pass.dispatch_workgroups((out_len as u32).div_ceil(64), 1, 1);
        }

        let read_fwd = device.create_buffer(&wgpu::BufferDescriptor { label: Some("read_fwd"), size: (out_len * 4) as u64, usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ, mapped_at_creation: false });
        let read_rev = device.create_buffer(&wgpu::BufferDescriptor { label: Some("read_rev"), size: (out_len * 4) as u64, usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ, mapped_at_creation: false });
        let read_valid = device.create_buffer(&wgpu::BufferDescriptor { label: Some("read_valid"), size: (out_len * 4) as u64, usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ, mapped_at_creation: false });
        encoder.copy_buffer_to_buffer(&out_fwd, 0, &read_fwd, 0, (out_len * 4) as u64);
        encoder.copy_buffer_to_buffer(&out_rev, 0, &read_rev, 0, (out_len * 4) as u64);
        encoder.copy_buffer_to_buffer(&out_valid, 0, &read_valid, 0, (out_len * 4) as u64);
        queue.submit(Some(encoder.finish()));

        for b in [&read_fwd, &read_rev, &read_valid] {
            let slice = b.slice(..);
            let (tx, rx) = std::sync::mpsc::channel();
            slice.map_async(wgpu::MapMode::Read, move |r| { let _ = tx.send(r); });
            device.poll(wgpu::Maintain::Wait);
            rx.recv().map_err(|e| e.to_string())?.map_err(|e| e.to_string())?;
        }

        let fwd: Vec<u32> = bytemuck::cast_slice(&read_fwd.slice(..).get_mapped_range()).to_vec();
        let rev: Vec<u32> = bytemuck::cast_slice(&read_rev.slice(..).get_mapped_range()).to_vec();
        let val: Vec<u32> = bytemuck::cast_slice(&read_valid.slice(..).get_mapped_range()).to_vec();

        let mut out = Vec::with_capacity(out_len);
        for i in 0..out_len {
            if val[i] == 1 {
                out.push((fwd[i] as u64, rev[i] as u64));
            }
        }
        Ok(out)
    })
}

impl Iterator for KmerGeneratorWGPU<'_> {
    type Item = (Kmer, Kmer);
    fn next(&mut self) -> Option<Self::Item> {
        while self.buffer.is_empty() {
            if self.done { return None; }
            self.process_batch();
            if !self.buffer.is_empty() {
                self.buffer.reverse();
            }
        }
        self.buffer.pop()
    }
}
