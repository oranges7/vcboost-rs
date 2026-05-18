use pyo3::prelude::*;

fn char_to_value(c: char) -> u8 {
    match c {
        'A' => 25, 'C' => 100, 'T' => 75, 'G' => 50,
        '-' | 'N' | 'R' | 'Y' | 'W' | 'K' | 'M' | 'S' | 'B' => 0,
        _ => 0,
    }
}

fn char_to_code(c: char) -> u8 {
    match c {
        'A' => 1, 'C' => 2, 'T' => 3, 'G' => 4,
        '-' | 'N' | 'R' | 'Y' | 'W' | 'K' | 'M' | 'S' | 'B' => 0,
        _ => 0,
    }
}

fn sanitize_ref(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'N' | 'R' | 'Y' | 'W' | 'K' | 'M' | 'S' | 'B' => '-',
            other => other,
        })
        .collect()
}

fn normalize_mq(mq: u8) -> u8 {
    (100u16 * mq.min(60) as u16 / 60) as u8
}

#[pyfunction]
pub fn base_to_code(c: char) -> PyResult<u8> {
    Ok(char_to_code(c))
}

#[pyfunction]
pub fn base_onehot(code: u8) -> PyResult<Vec<f32>> {
    let mut onehot = vec![0.0f32; 5];
    if (code as usize) < 5 {
        onehot[code as usize] = 1.0;
    }
    Ok(onehot)
}

#[pyfunction]
pub fn encode_hap_data(hap0: Vec<String>, hap1: Vec<String>) -> PyResult<(Vec<Vec<f32>>, Vec<Vec<f32>>)> {
    let seq_len = 33usize;
    let mut hap_0_data = vec![vec![0.0f32; seq_len]; 80];
    let mut hap_1_data = vec![vec![0.0f32; seq_len]; 80];
    for (idx, seq) in hap0.iter().enumerate().take(40) {
        if seq.len() < seq_len { continue; }
        for (j, c) in seq.chars().take(seq_len).enumerate() {
            hap_0_data[idx][j] = char_to_value(c) as f32;
        }
    }
    for (idx, seq) in hap1.iter().enumerate().take(40) {
        if seq.len() < seq_len { continue; }
        for (j, c) in seq.chars().take(seq_len).enumerate() {
            hap_1_data[idx][j] = char_to_value(c) as f32;
        }
    }
    Ok((hap_0_data, hap_1_data))
}

#[pyfunction]
pub fn gen_data_one_hap(d: Vec<(String, String, u8)>) -> PyResult<(Vec<Vec<f32>>, Vec<Vec<f32>>, usize, Vec<Vec<f32>>)> {
    let seq_len = 33usize;
    let max_reads = 80usize;
    let mut hap_0_data = vec![vec![0.0f32; seq_len]; max_reads];
    let mut hap_1_data = vec![vec![0.0f32; seq_len]; max_reads];
    let mut mq_data = vec![vec![0.0f32; seq_len]; max_reads];
    let mut count = 0usize;
    for (_qname, seq, mq) in &d {
        if count >= max_reads { break; }
        if seq.len() < seq_len { continue; }
        for (j, c) in seq.chars().take(seq_len).enumerate() {
            hap_0_data[count][j] = char_to_value(c) as f32;
        }
        if *mq == 1 { hap_1_data[count] = vec![100.0f32; seq_len]; }
        else { hap_1_data[count] = vec![0.0f32; seq_len]; }
        let mq_val = normalize_mq(*mq) as f32;
        mq_data[count] = vec![mq_val; seq_len];
        count += 1;
    }
    Ok((hap_0_data, hap_1_data, count, mq_data))
}

#[pyfunction]
pub fn gen_data(ref_seq: &str, d: Vec<(String, String, u8)>) -> PyResult<Vec<Vec<Vec<f32>>>> {
    let ref1 = sanitize_ref(ref_seq);
    let seq_len = 33usize;
    let max_reads = 80usize;
    let ref_chars: Vec<char> = ref1.chars().collect();
    let mut ref_data = vec![vec![0.0f32; seq_len]; max_reads];
    if ref_chars.len() >= seq_len {
        for (j, c) in ref_chars.iter().take(seq_len).enumerate() {
            ref_data[0][j] = char_to_value(*c) as f32;
        }
    }
    let (hap_0_data, _hap_1_data, count, mq_data) = gen_data_one_hap(d)?;
    if count < max_reads {
        for i in 1..count { ref_data[i] = ref_data[0].clone(); }
    } else {
        for i in 1..max_reads { ref_data[i] = ref_data[0].clone(); }
    }
    let channels = 4usize;
    let mut result = vec![vec![vec![0.0f32; channels]; seq_len]; max_reads];
    for i in 0..max_reads {
        for j in 0..seq_len {
            result[i][j][0] = ref_data[i][j];
            result[i][j][1] = hap_0_data[i][j];
            result[i][j][2] = _hap_1_data[i][j];
            result[i][j][3] = mq_data[i][j];
        }
    }
    Ok(result)
}
