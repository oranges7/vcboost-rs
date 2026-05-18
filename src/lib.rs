mod encoding;
mod msa;
mod pileup;
mod vcf;

use pyo3::prelude::*;

#[pymodule]
fn _vcboostrs(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(msa::fast_msa, m)?)?;
    m.add_function(wrap_pyfunction!(msa::encode_msa_output, m)?)?;
    m.add_function(wrap_pyfunction!(encoding::base_to_code, m)?)?;
    m.add_function(wrap_pyfunction!(encoding::base_onehot, m)?)?;
    m.add_function(wrap_pyfunction!(encoding::encode_hap_data, m)?)?;
    m.add_function(wrap_pyfunction!(encoding::gen_data_one_hap, m)?)?;
    m.add_function(wrap_pyfunction!(encoding::gen_data, m)?)?;
    m.add_function(wrap_pyfunction!(pileup::generate_features, m)?)?;
    m.add_function(wrap_pyfunction!(pileup::run_variant, m)?)?;
    m.add_function(wrap_pyfunction!(vcf::filter_heterozygous_snps, m)?)?;
    m.add_function(wrap_pyfunction!(vcf::merge_predictions_to_vcf, m)?)?;
    m.add_function(wrap_pyfunction!(vcf::split_vcf_by_chr, m)?)?;
    Ok(())
}
