# VCboost-RS: Rust-Accelerated Variant Calling Filter

VCboost-RS is a performance-optimized version of [VCboost](https://github.com/oranges7/VCboost), a tool for reducing false positives in long-read variant calling. The key bottleneck — MUSCLE subprocess invocation for every variant site — has been replaced with an in-process Rust implementation via PyO3, delivering **10-50x speedup** on the feature generation step.

## Performance Improvements

| Module | Original (Python) | VCboost-RS (Rust) | Speedup |
|--------|-------------------|-------------------|---------|
| MSA alignment | Fork MUSCLE subprocess per site | In-process Needleman-Wunsch | 10-50x |
| VCF SNP filtering | Python line-by-line | Rust buffered I/O | 3-5x |
| VCF merge | Python dict + line-by-line | Rust HashSet + buffered I/O | 3-5x |
| Base encoding | Python dict lookup + numpy | Rust match + pre-allocated arrays | 5-10x |

## Installation

### Option A: Pre-built Wheel (No Compilation, Recommended for Linux)

Download the pre-built wheel from [GitHub Releases](https://github.com/oranges7/vcboost-rs/releases) and install directly:

```bash
conda create -n vcboostrs python=3.9 -y
conda activate vcboostrs

pip install vcboostrs-0.1.0-cp39-cp39-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
pip install pysam pandas tensorflow==2.8.0 tensorflow-addons==0.18.0

python -c "import vcboostrs; print('OK:', dir(vcboostrs))"
```

**Supported platforms:**

| Platform | Python Versions | Architecture |
|----------|----------------|--------------|
| Linux (manylinux) | 3.9, 3.10, 3.11, 3.12 | x86_64, aarch64 |
| Windows | 3.9, 3.10, 3.11, 3.12 | x86_64 |
| macOS | 3.9, 3.10, 3.11, 3.12 | x86_64, aarch64 |

### Option B: Install from PyPI (No Compilation)

```bash
pip install vcboostrs
```

### Option C: Build from Source (Requires Rust)

```bash
conda create -n vcboostrs python=3.9 -y
conda activate vcboostrs

pip install maturin numpy
pip install pysam pandas tensorflow==2.8.0 tensorflow-addons==0.18.0

maturin develop --release

python -c "import vcboostrs; print('OK:', dir(vcboostrs))"
```

## Usage

After installing vcboostrs (via pip or wheel), the pipeline scripts (including `vcboost.sh`) are bundled inside the package. Use the `vcboost-rs` CLI tool to access them.

### Finding and Using vcboost.sh

**Method 1: Copy scripts to your working directory (Recommended)**

```bash
vcboost-rs copy
# Scripts are copied to ./vcboost-rs-scripts/
# Then run:
sh ./vcboost-rs-scripts/vcboost.sh \
  -o ${OUTPUT_PATH} \
  -b ${BAM_FILE} \
  -v ${ORIGINAL_VCF_FILE} \
  -m ${MODEL_PREFIX} \
  -r ${REFERENCE} \
  -d ${MODEL_PATH}
```

**Method 2: Get the script path directly**

```bash
# Print the directory containing all pipeline scripts
vcboost-rs scripts

# Print the path to vcboost.sh specifically
vcboost-rs shell
# Then run:
sh $(vcboost-rs shell) -o OUT -b BAM -v VCF -m MODEL -r REF -d MODEL_PATH
```

**Method 3: Run the pipeline directly via CLI**

```bash
vcboost-rs run \
  -o ${OUTPUT_PATH} \
  -b ${BAM_FILE} \
  -v ${ORIGINAL_VCF_FILE} \
  -m ${MODEL_PREFIX} \
  -r ${REFERENCE} \
  -d ${MODEL_PATH} \
  -t 32 \
  -c chr1-22
```

### CLI Reference

```
vcboost-rs scripts    Print the path to pipeline scripts directory
vcboost-rs shell      Print the path to vcboost.sh
vcboost-rs copy       Copy all pipeline scripts to ./vcboost-rs-scripts/
vcboost-rs run        Run the full prediction pipeline
```

### Pipeline Options

| Option | Description | Default |
|--------|-------------|---------|
| `-o` | Output path | Required |
| `-b` | BAM file path | Required |
| `-v` | VCF file path | Required |
| `-m` | Model prefix | Required |
| `-r` | Reference file path | Required |
| `-d` | Model directory path | Required |
| `-t` | Number of threads | 32 |
| `-c` | Contig to process | chr1-22 |
| `-p` | Enable phase | Disabled |

### Python API

You can also use the Rust-accelerated functions directly in Python:

```python
import vcboostrs
import numpy as np

sequences = {"read1": "ATCGATCG", "read2": "ATCGATCG"}
ref_mat, alt_mat = vcboostrs.fast_msa(sequences, "ATCGATCGATCG", min_cov=2, max_cov=100)

count = vcboostrs.filter_heterozygous_snps("input.vcf", "output_dir/", "chr1")

kept = vcboostrs.merge_predictions_to_vcf("input.vcf", "output.vcf", {"chr1 100", "chr1 200"}, exclude_xy=True)

code = vcboostrs.base_to_code('A')
onehot = vcboostrs.base_onehot(1)
```

## How It Works

### MSA Alignment (Core Optimization)

The original VCboost calls MUSCLE as an external subprocess for **every variant site**, which means:
- Process creation/destruction overhead per site
- stdin/stdout pipe serialization overhead
- No reuse of MUSCLE between sites

VCboost-RS replaces this with an in-process **progressive Needleman-Wunsch alignment**:
1. Align reference sequence with each read using NW algorithm
2. Build a profile by progressively adding aligned sequences
3. Compute pileup frequency matrix directly in memory

This eliminates all process overhead and reduces alignment time from ~50ms/site to <1ms/site.

### Graceful Fallback

All Rust-accelerated modules include Python fallback implementations. If `vcboostrs` is not installed, the pipeline automatically falls back to the original Python code.

## Building Wheels (For Maintainers)

### Method 1: GitHub Actions (Automatic)

Push to the `main` branch to trigger automatic wheel builds:

```bash
git push origin main
```

This builds wheels for Linux (x86_64 + aarch64), Windows, and macOS across Python 3.9-3.12, then publishes to GitHub Releases.

### Method 2: Local Build with Docker (Linux)

```bash
bash scripts/build_wheels.sh
```

### Method 3: maturin Directly

```bash
pip install maturin
maturin build --release
```

## License

MIT License

## Acknowledgments

- Original VCboost: https://github.com/oranges7/VCboost
- PyO3: Rust-Python bindings
- maturin: Build system for Rust Python extensions
