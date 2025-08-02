# W_bulk_RNA_seq_pipeline

A Snakemake-based pipeline for preprocessing, quality control, alignment, quantification, and transformation of bulk RNA-seq data. This pipeline is designed for flexibility and reproducibility, supporting both single-end and paired-end reads and producing normalized count matrices ready for downstream analysis.

## 🧬 Features

- Download and organize raw FASTQ files
- Quality control with FastQC and MultiQC
- Adapter trimming and deduplication using Fastp
- Alignment using HISAT2
- Gene-level quantification with FeatureCounts
- Aggregation and transformation of expression data (TPM, log2TPM, VST)
- Configurable support for single- or paired-end reads

## 🧰 Requirements

- [Snakemake](https://snakemake.readthedocs.io)
- Python ≥ 3.8
- External tools:
  - `fastqc`
  - `fastp`
  - `multiqc`
  - `hisat2`
  - `samtools`
  - `featureCounts` (from Subread)
- Python libraries: `pandas`, `numpy` (used in `aggregate_counts.py`)

## 📂 Directory Structure

```bash
.
├── config.yaml               # Set sample IDs, paths, paired-end flag, etc.
├── Snakefile                 # Main Snakemake workflow
├── aggregate_counts.py       # Aggregates counts and generates TPM, log2TPM, VST matrices
├── input/
│   └── sample_ids.txt        # List of SRR IDs or sample names
├── raw_data/                 # Input FASTQ files
├── trimmed_data/             # Output from Fastp
├── QC/
│   ├── raw/                  # FastQC + MultiQC for raw data
│   └── trimmed/              # FastQC + MultiQC for trimmed data
├── alignment/               # Aligned BAM files (sorted)
├── counts/                  # Count matrices (raw, TPM, log2TPM, VST)
└── ref_genomes/
    └── hg38/                # Contains HISAT2 index and GTF file
