use pyo3::prelude::*;
use std::collections::HashMap;

fn char_to_msa_code(c: char) -> u8 {
    match c {
        'A' => 1, 'C' => 2, 'T' => 3, 'G' => 4, '-' | '_' => 0,
        _ => 0,
    }
}

fn sanitize_ref_msa(seq: &str) -> String {
    seq.chars()
        .map(|c| match c {
            'N' => '_',
            'R' | 'Y' | 'W' | 'K' | 'M' | 'S' | 'B' => '-',
            other => other,
        })
        .collect()
}

fn needleman_wunsch(seq_a: &[u8], seq_b: &[u8], match_score: i32, mismatch_penalty: i32, gap_penalty: i32) -> (Vec<u8>, Vec<u8>) {
    let len_a = seq_a.len();
    let len_b = seq_b.len();

    let mut score = vec![vec![0i32; len_b + 1]; len_a + 1];
    for i in 1..=len_a {
        score[i][0] = score[i - 1][0] + gap_penalty;
    }
    for j in 1..=len_b {
        score[0][j] = score[0][j - 1] + gap_penalty;
    }

    for i in 1..=len_a {
        for j in 1..=len_b {
            let s = if seq_a[i - 1] == seq_b[j - 1] { match_score } else { mismatch_penalty };
            score[i][j] = (score[i - 1][j - 1] + s)
                .max(score[i - 1][j] + gap_penalty)
                .max(score[i][j - 1] + gap_penalty);
        }
    }

    let mut aligned_a = Vec::new();
    let mut aligned_b = Vec::new();
    let mut i = len_a;
    let mut j = len_b;

    while i > 0 || j > 0 {
        if i > 0 && j > 0 {
            let s = if seq_a[i - 1] == seq_b[j - 1] { match_score } else { mismatch_penalty };
            if score[i][j] == score[i - 1][j - 1] + s {
                aligned_a.push(seq_a[i - 1]);
                aligned_b.push(seq_b[j - 1]);
                i -= 1;
                j -= 1;
                continue;
            }
        }
        if i > 0 && score[i][j] == score[i - 1][j] + gap_penalty {
            aligned_a.push(seq_a[i - 1]);
            aligned_b.push(b'-');
            i -= 1;
        } else {
            aligned_a.push(b'-');
            aligned_b.push(seq_b[j - 1]);
            j -= 1;
        }
    }

    aligned_a.reverse();
    aligned_b.reverse();
    (aligned_a, aligned_b)
}

fn progressive_msa(sequences: &[&str], reference: &str) -> (String, Vec<String>) {
    if sequences.is_empty() {
        return (reference.to_string(), Vec::new());
    }

    let ref_bytes: Vec<u8> = reference.bytes().collect();
    let mut profile: Vec<Vec<u8>> = vec![ref_bytes];
    let mut aligned_ref = reference.to_string();

    for seq in sequences {
        let seq_bytes: Vec<u8> = seq.bytes().collect();
        let last_seq: Vec<u8> = profile.last().unwrap().clone();
        let (new_last, new_seq) = needleman_wunsch(&last_seq, &seq_bytes, 2, -1, -1);

        let gap_insert_positions: Vec<usize> = {
            let mut positions = Vec::new();
            for (i, &c) in new_last.iter().enumerate() {
                if c == b'-' && (i == 0 || new_last[i - 1] != b'-') {
                    positions.push(i);
                }
            }
            positions
        };

        if !gap_insert_positions.is_empty() {
            let mut offset = 0isize;
            for &pos in &gap_insert_positions {
                let actual_pos = (pos as isize + offset) as usize;
                for existing in profile.iter_mut() {
                    existing.insert(actual_pos, b'-');
                }
                let mut ref_chars: Vec<u8> = aligned_ref.bytes().collect();
                ref_chars.insert(actual_pos, b'-');
                aligned_ref = String::from_utf8_lossy(&ref_chars).to_string();
                offset += 1;
            }
        }

        let mut adjusted_new_seq = new_seq.clone();
        let profile_len = profile.last().unwrap().len();
        if adjusted_new_seq.len() < profile_len {
            adjusted_new_seq.resize(profile_len, b'-');
        }

        profile.push(adjusted_new_seq);
    }

    let aligned_seqs: Vec<String> = profile[1..]
        .iter()
        .map(|s| String::from_utf8_lossy(s).to_string())
        .collect();

    (aligned_ref, aligned_seqs)
}

fn compute_pileup_matrix(
    aligned_seqs: &[String],
    reference: &str,
) -> (Vec<Vec<f32>>, Vec<Vec<f32>>) {
    let ref_sanitized = sanitize_ref_msa(reference);
    let ref_chars: Vec<char> = ref_sanitized.chars().collect();
    let max_cols = 128usize.min(ref_chars.len());

    let mut ref_mat = vec![vec![0.0f32; 5]; 128];
    for (i, c) in ref_chars.iter().take(max_cols).enumerate() {
        let code = char_to_msa_code(*c) as usize;
        if code < 5 {
            ref_mat[i][code] = 100.0;
        }
    }

    let mut count_mat = vec![vec![0u32; 5]; 128];
    for seq in aligned_seqs {
        let seq_chars: Vec<char> = seq.chars().collect();
        for i in 0..max_cols.min(seq_chars.len()) {
            let code = char_to_msa_code(seq_chars[i]) as usize;
            if code < 5 {
                count_mat[i][code] += 1;
            }
        }
    }

    let mut alt_mat = vec![vec![0.0f32; 5]; 128];
    for i in 0..max_cols {
        let total: u32 = count_mat[i].iter().sum();
        if total > 0 {
            for j in 0..5 {
                alt_mat[i][j] = (count_mat[i][j] as f32 / total as f32) * 100.0;
            }
        }
    }

    (ref_mat, alt_mat)
}

#[pyfunction]
pub fn fast_msa(
    sequences: HashMap<String, String>,
    reference: &str,
    min_cov: usize,
    max_cov: usize,
) -> PyResult<(Vec<Vec<f32>>, Vec<Vec<f32>>)> {
    let mut sample: Vec<String> = sequences.keys().cloned().collect();

    if sample.len() < min_cov {
        let zero_ref = vec![vec![0.0f32; 5]; 72];
        let zero_alt = vec![vec![0.0f32; 5]; 72];
        return Ok((zero_ref, zero_alt));
    }

    if sample.len() > max_cov {
        let seed = 812u64;
        let mut rng = seed;
        while sample.len() > max_cov {
            rng = rng.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            let idx = ((rng >> 33) as usize) % sample.len();
            sample.swap_remove(idx);
        }
    }

    sample.sort();

    let ref1 = sanitize_ref_msa(reference);
    let seq_strs: Vec<&str> = sample.iter().map(|name| sequences[name].as_str()).collect();
    let (_aligned_ref, aligned_seqs) = progressive_msa(&seq_strs, &ref1);
    let (ref_mat, alt_mat) = compute_pileup_matrix(&aligned_seqs, reference);

    let mut ref_out = vec![vec![0.0f32; 5]; 72];
    let mut alt_out = vec![vec![0.0f32; 5]; 72];
    let rows = ref_mat.len().min(72);
    for i in 0..rows {
        ref_out[i] = ref_mat[i].clone();
        alt_out[i] = alt_mat[i].clone();
    }

    Ok((ref_out, alt_out))
}

#[pyfunction]
pub fn encode_msa_output(
    ref_data_0: Vec<Vec<f32>>,
    data_0: Vec<Vec<f32>>,
    ref_data_1: Vec<Vec<f32>>,
    data_1: Vec<Vec<f32>>,
) -> PyResult<Vec<i16>> {
    let rows = 72usize;
    let cols = 5usize;
    let channels = 4usize;
    let mut result = vec![0i16; rows * cols * channels];
    for i in 0..rows {
        for j in 0..cols {
            let v0 = if i < ref_data_0.len() && j < ref_data_0[i].len() { ref_data_0[i][j] } else { 0.0 };
            let v1 = if i < data_0.len() && j < data_0[i].len() { data_0[i][j] } else { 0.0 };
            let v2 = if i < ref_data_1.len() && j < ref_data_1[i].len() { ref_data_1[i][j] } else { 0.0 };
            let v3 = if i < data_1.len() && j < data_1[i].len() { data_1[i][j] } else { 0.0 };
            let idx = (i * cols + j) * channels;
            result[idx] = v0 as i16;
            result[idx + 1] = v1 as i16;
            result[idx + 2] = v2 as i16;
            result[idx + 3] = v3 as i16;
        }
    }
    Ok(result)
}
