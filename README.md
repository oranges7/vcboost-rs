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

```bash
sh vcboost.sh \
  -o ${OUTPUT_PATH} \
  -b ${BAM_FILE} \
  -v ${ORIGINAL_VCF_FILE} \
  -m ${MODEL_PREFIX} \
  -r ${REFERENCE}
```

## Building Wheels (For Maintainers)

### Method 1: GitHub Actions (Automatic)

Push a version tag to trigger automatic wheel builds:

```bash
git tag v0.1.0
git push origin v0.1.0
```

### Method 2: Local Build with Docker (Linux)

```bash
bash scripts/build_wheels.sh
```

## License

MIT License

## Acknowledgments

- Original VCboost: https://github.com/oranges7/VCboost
- PyO3: Rust-Python bindings
- maturin: Build system for Rust Python extensions
