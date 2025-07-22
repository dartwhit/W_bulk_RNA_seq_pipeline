# Snakefile

import os
import shutil


# Config / variables

SRR_IDs_file = config["sample_ids"]
OUTDIR = config.get("outdir")
SRA_PROJECT_ID = config.get("sra_project_id", "new_proj")
PAIRED_END = config.get("paired_end", True)

RAW_DIR = os.path.join(OUTDIR, "raw_data")
TRIMMED_DIR = os.path.join(OUTDIR, "trimmed_data")
QC_DIR = os.path.join(OUTDIR, "QC")


def get_samples():
    with open(SRR_IDs_file) as f:
        samples = [line.strip() for line in f if line.strip()]
    return samples



rule all:
    input:
        # FastQC reports on raw data
        expand(os.path.join(QC_DIR, "raw", "{sample}_R1_fastqc.html"), sample=get_samples()),
        *(
            expand(os.path.join(QC_DIR, "raw", "{sample}_R2_fastqc.html"), sample=get_samples())
            if PAIRED_END else []
        ),
        # Trimmed fastq files (to ensure trimming is done)
        expand(os.path.join(TRIMMED_DIR, "{sample}_R1_paired.fastq.gz"), sample=get_samples()),
        *(
            expand(os.path.join(TRIMMED_DIR, "{sample}_R2_paired.fastq.gz"), sample=get_samples())
            if PAIRED_END else []
        ),
        # FastQC reports on trimmed data
        expand(os.path.join(QC_DIR, "trimmed", "{sample}_R1_paired_fastqc.html"), sample=get_samples()),
        *(
            expand(os.path.join(QC_DIR, "trimmed", "{sample}_R2_paired_fastqc.html"), sample=get_samples())
            if PAIRED_END else []
        ),
        # MultiQC reports
        os.path.join(QC_DIR, "raw", "multiqc_raw_report.html"),
        os.path.join(QC_DIR, "trimmed", "multiqc_trimmed_report.html"),
        # Alignment
        expand(os.path.join(OUTDIR, "alignment", "{sample}.sorted.bam"), sample=get_samples()),
        # Counts
        expand(os.path.join(OUTDIR, "counts", "{sample}_counts.txt"), sample=get_samples()),
        os.path.join(OUTDIR, "counts", "counts_by_ensembl.txt"),
        os.path.join(OUTDIR, "counts", "counts_by_gene.txt"),
        os.path.join(OUTDIR, "counts", "tpm_by_ensembl.txt"), 
        os.path.join(OUTDIR, "counts", "tpm_by_gene.txt"),
        os.path.join(OUTDIR, "counts", "log2tpm_by_ensembl.txt"), 
        os.path.join(OUTDIR, "counts", "log2tpm_by_gene.txt"),
        os.path.join(OUTDIR, "counts", "vst_by_ensembl.txt"), 
        os.path.join(OUTDIR, "counts", "vst_by_gene.txt")




rule fastqc_raw:
    input:
        fq=lambda wildcards: next(
            os.path.join(RAW_DIR, f)
            for f in os.listdir(RAW_DIR)
            if re.match(rf"{re.escape(wildcards.sample)}_R{wildcards.read}(_.*)?\.fastq\.gz", f)
        )
    output:
        html=os.path.join(QC_DIR, "raw", "{sample}_R{read}_fastqc.html"),
        zip=os.path.join(QC_DIR, "raw", "{sample}_R{read}_fastqc.zip")
    resources:
        mem="64G",
        time="01:00:00"
    run:
        outdir = os.path.join(QC_DIR, "raw")
        os.makedirs(outdir, exist_ok=True)
        
        # Run fastqc with original input filename
        shell(f"fastqc {input.fq} -o {outdir}")
        
        # Find the generated output files (they keep the full input basename)
        input_basename = os.path.basename(input.fq)  # e.g. sample_R1_001.fastq.gz
        prefix = input_basename.rsplit(".", 2)[0]  # remove .fastq.gz, e.g. sample_R1_001

        # Expected actual output files
        actual_html = os.path.join(outdir, f"{prefix}_fastqc.html")
        actual_zip = os.path.join(outdir, f"{prefix}_fastqc.zip")
        
        # Desired output files (standardized without suffix)
        desired_html = output.html
        desired_zip = output.zip

        # Rename files to remove suffix from name
        if os.path.exists(actual_html):
            shutil.move(actual_html, desired_html)
        else:
            raise FileNotFoundError(f"Expected fastqc html output not found: {actual_html}")
        
        if os.path.exists(actual_zip):
            shutil.move(actual_zip, desired_zip)
        else:
            raise FileNotFoundError(f"Expected fastqc zip output not found: {actual_zip}")
rule multiqc_raw:
    input:
        expand(
            os.path.join(QC_DIR, "raw", "{sample}_R{read}_fastqc.html"),
            sample=get_samples(),
            read=["1", "2"] if PAIRED_END else ["1"]
        )
    output:
        os.path.join(QC_DIR, "raw", "multiqc_raw_report.html")
    resources:
        mem="64G",
        time="01:00:00"
    run:
        outdir = os.path.join(QC_DIR, "raw")
        os.makedirs(outdir, exist_ok=True)
        shell(f"multiqc -f {outdir} -o {outdir} -n multiqc_raw_report.html")


rule fastp:
    input:
        r1=lambda wildcards: next(
            os.path.join(RAW_DIR, f)
            for f in os.listdir(RAW_DIR)
            if f.startswith(f"{wildcards.sample}_R1") and f.endswith(".fastq.gz")
        ),
        r2=lambda wildcards: next(
            os.path.join(RAW_DIR, f)
            for f in os.listdir(RAW_DIR)
            if f.startswith(f"{wildcards.sample}_R2") and f.endswith(".fastq.gz")
        ) if PAIRED_END else None
    output:
        trimmed_r1=os.path.join(TRIMMED_DIR, "{sample}_R1_paired.fastq.gz"),
        trimmed_r2=os.path.join(TRIMMED_DIR, "{sample}_R2_paired.fastq.gz") if PAIRED_END else None
    threads: 8
    resources:
        mem="16G",
        time="02:00:00"
    run:
        os.makedirs(TRIMMED_DIR, exist_ok=True)
        if PAIRED_END:
            shell(
                "fastp -i {input.r1} -I {input.r2} "
                "-o {output.trimmed_r1} -O {output.trimmed_r2} "
                "-w {threads} --detect_adapter_for_pe "
                "--thread {threads} "
                "--trim_poly_g "
                "-p --dedup "
                "--html {TRIMMED_DIR}/{wildcards.sample}_fastp.html "
                "--json {TRIMMED_DIR}/{wildcards.sample}_fastp.json"
            )
        else:
            shell(
                "fastp -i {input.r1} -o {output.trimmed_r1} "
                "-w {threads} --detect_adapter_for_pe "
                "--thread {threads} "
                "--qualified_quality_phred 15 "
                "--length_required 36 "
                "--html {TRIMMED_DIR}/{wildcards.sample}_fastp.html "
                "--json {TRIMMED_DIR}/{wildcards.sample}_fastp.json"
            )




rule fastqc_trimmed:
    input:
        fq=os.path.join(TRIMMED_DIR, "{sample}_R{read}_paired.fastq.gz")
    output:
        html=os.path.join(QC_DIR, "trimmed", "{sample}_R{read}_paired_fastqc.html"),
        zip=os.path.join(QC_DIR, "trimmed", "{sample}_R{read}_paired_fastqc.zip")
    resources:
        mem="64G",
        time="01:00:00"
    run:
        outdir = os.path.join(QC_DIR, "trimmed")
        os.makedirs(outdir, exist_ok=True)
        shell(f"fastqc {input.fq} -o {outdir}")


rule multiqc_trimmed:
    input:
        expand(
            os.path.join(QC_DIR, "trimmed", "{sample}_R{read}_paired_fastqc.html"),
            sample=get_samples(),
            read=["1", "2"] if PAIRED_END else ["1"]
        )
    output:
        os.path.join(QC_DIR, "trimmed", "multiqc_trimmed_report.html")
    resources:
        mem="64G",
        time="01:00:00"
    run:
        outdir = os.path.join(QC_DIR, "trimmed")
        os.makedirs(outdir, exist_ok=True)
        shell(f"multiqc -f {outdir} -o {outdir} -n multiqc_trimmed_report.html")


rule hisat2_align:
    input:
        r1=os.path.join(TRIMMED_DIR, "{sample}_R1_paired.fastq.gz"),
        r2=os.path.join(TRIMMED_DIR, "{sample}_R2_paired.fastq.gz") if PAIRED_END else None,
    output:
        bam=os.path.join(OUTDIR, "alignment", "{sample}.sorted.bam")
    threads: 16
    resources:
        mem="128G",
        time="04:00:00"
    params:
        extra="",  # additional hisat2 options if needed
        index_prefix="ref_genomes/hg38/hg38"  # adjust to your actual HISAT2 index prefix
    run:
        import os
        outdir = os.path.join(QC_DIR, "alignment")        
        if PAIRED_END:
            shell(
                "hisat2 -p {threads} -x {params.index_prefix} "
                "-1 {input.r1} -2 {input.r2} "
                "{params.extra} | samtools sort -@ {threads} -o {output.bam}"
            )
        else:
            shell(
                "hisat2 -p {threads} -x {params.index_prefix} "
                "-U {input.r1} "
                "{params.extra} | samtools sort -@ {threads} -o {output.bam}"
            )


rule featurecounts:
    input:
        bam=os.path.join(OUTDIR, "alignment", "{sample}.sorted.bam"),
        gtf="ref_genomes/hg38/Homo_sapiens.GRCh38.110.gtf"  # adjust to your annotation path
    output:
        counts=os.path.join(OUTDIR, "counts", "{sample}_counts.txt")
    threads: 8
    resources:
        mem="128G",
        time="02:00:00"
    run:
        import os
        outdir = os.path.dirname(output.counts)
        os.makedirs(outdir, exist_ok=True)
        if PAIRED_END: 
            shell(
                "featureCounts -T {threads} -a {input.gtf} -o {output.counts} -s 0 -p {input.bam}"
            )
        else:
            shell(
                "featureCounts -T {threads} -a {input.gtf} -o {output.counts} {input.bam}"
            )

rule aggregate_counts:
    input:
        countfiles=expand(os.path.join(OUTDIR, "counts", "{sample}_counts.txt"), sample=get_samples()),
        gtf="ref_genomes/hg38/Homo_sapiens.GRCh38.110.gtf",
        outdir=os.path.join(OUTDIR, "counts/")
    output:
        ensembl_counts=os.path.join(OUTDIR, "counts", "counts_by_ensembl.txt"),
        gene_counts=os.path.join(OUTDIR, "counts", "counts_by_gene.txt"),
        ensembl_tpm=os.path.join(OUTDIR, "counts", "tpm_by_ensembl.txt"), 
        gene_tpm=os.path.join(OUTDIR, "counts", "tpm_by_gene.txt"),
        ensembl_logtpm=os.path.join(OUTDIR, "counts", "log2tpm_by_ensembl.txt"), 
        gene_logtpm=os.path.join(OUTDIR, "counts", "log2tpm_by_gene.txt"),
        ensembl_vst=os.path.join(OUTDIR, "counts", "vst_by_ensembl.txt"), 
        gene_vst=os.path.join(OUTDIR, "counts", "vst_by_gene.txt")

    threads: 1
    resources:
        mem="64G",
        time="01:00:00"
    shell:
        "python aggregate_counts.py "
        "-c {input.countfiles} "
        "-g {input.gtf} "
        "-o {input.outdir}"

