use pyo3::prelude::*;
use std::collections::HashSet;
use std::io::{BufRead, BufReader, Write};
use std::fs::File;

#[pyfunction]
pub fn filter_heterozygous_snps(vcf_path: &str, output_dir: &str, contig: &str) -> PyResult<usize> {
    let file = File::open(vcf_path)?;
    let reader = BufReader::new(file);
    let mut header = String::new();
    let mut variants: Vec<String> = Vec::new();
    let mut count = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') { header.push_str(&line); header.push('\n'); continue; }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 10 { continue; }
        let chr = columns[0];
        if !contig.is_empty() && chr != contig { continue; }
        let ref_base = columns[3]; let alt_base = columns[4];
        if ref_base.len() != 1 || alt_base.len() != 1 { continue; }
        let genotype = columns[9].split(':').next().unwrap_or("").replace('|', "/");
        if genotype == "0/1" || genotype == "1/0" { variants.push(line); count += 1; }
    }
    let output_path = format!("{}/{}.vcf", output_dir, contig);
    let mut out_file = File::create(&output_path)?;
    write!(out_file, "{}", header)?;
    for variant in &variants { writeln!(out_file, "{}", variant)?; }
    Ok(count)
}

#[pyfunction]
pub fn merge_predictions_to_vcf(vcf_path: &str, output_path: &str, filter_positions: HashSet<String>, exclude_xy: bool) -> PyResult<usize> {
    let file = File::open(vcf_path)?;
    let reader = BufReader::new(file);
    let mut out_file = File::create(output_path)?;
    let mut kept = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') { writeln!(out_file, "{}", line.trim_end())?; continue; }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.len() < 2 { writeln!(out_file, "{}", line.trim_end())?; kept += 1; continue; }
        if exclude_xy && (columns[0] == "chrX" || columns[0] == "chrY") { continue; }
        let key = format!("{} {}", columns[0], columns[1]);
        if !filter_positions.contains(&key) { writeln!(out_file, "{}", line.trim_end())?; kept += 1; }
    }
    Ok(kept)
}

#[pyfunction]
pub fn split_vcf_by_chr(vcf_path: &str, output_folder: &str, contig: &str, batch_size: usize) -> PyResult<usize> {
    let chromosomes: Vec<String> = (1..=22).map(|i| format!("chr{}", i)).collect();
    let chromosomes_xy: Vec<String> = chromosomes.iter().chain(&[String::from("chrX"), String::from("chrY")]).cloned().collect();
    let valid_chromosomes: &[String] = if contig == "chr1-22" { &chromosomes } else { &chromosomes_xy };
    let mut chr_counts: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    let file = File::open(vcf_path)?;
    let reader = BufReader::new(file);
    let mut batch_files: std::collections::HashSet<String> = std::collections::HashSet::new();
    let mut total_written = 0usize;
    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') { continue; }
        let columns: Vec<&str> = line.split('\t').collect();
        if columns.is_empty() { continue; }
        let chr = columns[0].to_string();
        if !valid_chromosomes.contains(&chr) { continue; }
        let count = chr_counts.entry(chr.clone()).or_insert(0);
        let batch_num = *count / batch_size;
        *count += 1;
        let out_file_name = format!("{}_{}", chr, batch_num);
        let out_path = format!("{}/{}", output_folder, out_file_name);
        let is_new = batch_files.insert(out_file_name.clone());
        if columns.len() >= 5 {
            let content = format!("{} {} {} {}", columns[0], columns[1], columns[3], columns[4]);
            if let Ok(mut f) = std::fs::OpenOptions::new().append(true).create(true).open(&out_path) {
                let _ = writeln!(f, "{}", content);
            }
        }
        total_written += 1;
    }
    Ok(total_written)
}
