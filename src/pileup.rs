use pyo3::prelude::*;
use std::collections::HashSet;

fn char_to_value(c: char) -> u8 {
    match c {
        'A' => 25, 'C' => 100, 'T' => 75, 'G' => 50,
        '-' | 'N' | 'R' | 'Y' | 'W' | 'K' | 'M' | 'S' | 'B' => 0,
        _ => 0,
    }
}

fn sanitize_ref(seq: &str) -> String {
    seq.chars().map(|c| match c {
        'N' | 'R' | 'Y' | 'W' | 'K' | 'M' | 'S' | 'B' => '-',
        other => other,
    }).collect()
}

fn normalize_mq(mq: u8) -> u8 { (100u16 * mq.min(60) as u16 / 60) as u8 }

#[pyfunction]
pub fn generate_features(ref_seq: &str, reads: Vec<(String, String, u8)>, hap_reads_0: Vec<String>, hap_reads_1: Vec<String>) -> PyResult<Vec<Vec<Vec<f32>>>> {
    let hap_set_0: HashSet<String> = hap_reads_0.into_iter().collect();
    let hap_set_1: HashSet<String> = hap_reads_1.into_iter().collect();
    let ref_sanitized = sanitize_ref(ref_seq);
    let ref_chars: Vec<char> = ref_sanitized.chars().take(160).collect();
    let seq_len = 33usize;
    let mut ref_data = vec![vec![0.0f32; seq_len]; 80];
    for (j, c) in ref_chars.iter().take(seq_len).enumerate() { ref_data[0][j] = char_to_value(*c) as f32; }
    for i in 1..80 { ref_data[i] = ref_data[0].clone(); }
    let mut hap_0_data = vec![vec![0.0f32; seq_len]; 80];
    let mut hap_1_data = vec![vec![0.0f32; seq_len]; 80];
    let mut mq_data = vec![vec![0.0f32; seq_len]; 80];
    let mut count = 0usize;
    for (qname, seq, mq) in &reads {
        if count >= 80 { break; }
        if seq.len() < seq_len { continue; }
        let _is_hap0 = hap_set_0.contains(qname);
        let _is_hap1 = hap_set_1.contains(qname);
        for (j, c) in seq.chars().take(seq_len).enumerate() { hap_0_data[count][j] = char_to_value(c) as f32; }
        if *mq == 1 { hap_1_data[count] = vec![100.0f32; seq_len]; }
        else { hap_1_data[count] = vec![0.0f32; seq_len]; }
        let mq_val = normalize_mq(*mq) as f32;
        mq_data[count] = vec![mq_val; seq_len];
        count += 1;
    }
    let channels = 4;
    let mut result = vec![vec![vec![0.0f32; channels]; seq_len]; 80];
    for i in 0..80 { for j in 0..seq_len {
        result[i][j][0] = ref_data[i][j]; result[i][j][1] = hap_0_data[i][j];
        result[i][j][2] = hap_1_data[i][j]; result[i][j][3] = mq_data[i][j];
    }}
    Ok(result)
}

#[pyfunction]
pub fn run_variant(_chrom: &str, _v_pos: i64, info1: &str, info2: &str, ref_seq: &str, hap0_reads: Vec<(String, String)>, hap1_reads: Vec<(String, String)>) -> PyResult<Option<Vec<i16>>> {
    let min_cov = 2usize; let max_cov = 100usize;
    let info1_clean = info1.replace('N', "-").replace('R', "-").replace('Y', "-").replace('W', "-").replace('K', "-").replace('M', "-").replace('S', "-").replace('B', "-");
    let info2_clean = info2.replace('N', "-").replace('R', "-").replace('Y', "-").replace('W', "-").replace('K', "-").replace('M', "-").replace('S', "-").replace('B', "-");
    let info1_first = info1_clean.split(',').next().unwrap_or("");
    let info2_first = info2_clean.split(',').next().unwrap_or("");
    let ref_prefix = info2_first;
    let ref_suffix_start = info1_first.len();
    let ref_full = if ref_suffix_start < ref_seq.len() { format!("{}{}", ref_prefix, &ref_seq[ref_suffix_start..]) } else { ref_prefix.to_string() };
    let ref_trimmed: String = ref_full.chars().take(160).collect();
    let mut d0: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    for (name, seq) in &hap0_reads { d0.insert(name.clone(), seq.clone()); }
    let mut d1: std::collections::HashMap<String, String> = std::collections::HashMap::new();
    for (name, seq) in &hap1_reads { d1.insert(name.clone(), seq.clone()); }
    let (ref_data_0, data_0) = if d0.is_empty() { (vec![vec![0.0f32; 5]; 72], vec![vec![0.0f32; 5]; 72]) } else { match crate::msa::fast_msa(d0, &ref_trimmed, min_cov, max_cov) { Ok(r) => r, Err(_) => (vec![vec![0.0f32; 5]; 72], vec![vec![0.0f32; 5]; 72]) } };
    let (ref_data_1, data_1) = if d1.is_empty() { (vec![vec![0.0f32; 5]; 72], vec![vec![0.0f32; 5]; 72]) } else { match crate::msa::fast_msa(d1, &ref_trimmed, min_cov, max_cov) { Ok(r) => r, Err(_) => (vec![vec![0.0f32; 5]; 72], vec![vec![0.0f32; 5]; 72]) } };
    match crate::msa::encode_msa_output(ref_data_0, data_0, ref_data_1, data_1) { Ok(flat) => Ok(Some(flat)), Err(_) => Ok(None) }
}
