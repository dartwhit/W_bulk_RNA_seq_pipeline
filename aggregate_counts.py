#!/usr/bin/env python3

import pandas as pd
import re
import argparse
import os
import numpy as np

def parse_gtf(gtf_file):
    gene_map = {}
    gene_length = {}
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != "gene":
                continue
            info = fields[8]
            gene_id_match = re.search('gene_id "([^"]+)"', info)
            gene_name_match = re.search('gene_name "([^"]+)"', info)
            if gene_id_match:
                gene_id = gene_id_match.group(1)
                gene_name = gene_name_match.group(1) if gene_name_match else gene_id
                gene_map[gene_id] = gene_name

            try:
                start = int(fields[3])
                end = int(fields[4])
                gene_length[gene_id] = end - start + 1
                # print(f"{gene_id} ({gene_name}): start={start}, end={end}, length={gene_length[gene_id]}")
            except ValueError:
                print(f"Warning: Failed to parse start/end for line: {line.strip()}")
                gene_length[gene_id] = np.nan
    return gene_map, gene_length

def load_counts(file_path):
    df = pd.read_csv(file_path, sep='\t', comment='#')
    df = df.rename(columns={df.columns[0]: "EnsemblID"})
    counts_col = df.columns[-1]
    df = df[["EnsemblID", counts_col]]
    df = df.set_index("EnsemblID")
    df = df.rename(columns={counts_col: os.path.basename(file_path).replace('_counts.txt','')})
    return df

def calculate_tpm(counts_df, lengths):
    lengths_kb = lengths / 1000
    lengths_kb = lengths_kb.reindex(counts_df.index)

    # Warn about missing lengths
    missing = lengths_kb.isna()
    if missing.any():
        print(f"Warning: {missing.sum()} genes in counts have no length info. TPM will be zero for them.")

    lengths_kb = lengths_kb.replace(0, np.nan)
    rate = counts_df.div(lengths_kb, axis=0)
    rate_sum = rate.sum(axis=0)
    tpm = rate.div(rate_sum, axis=1) * 1e6
    tpm = tpm.fillna(0)

    return tpm


def vst_transform(counts_df):
    import scanpy as sc
    counts_numeric = counts_df.select_dtypes(include=['number'])
    adata = sc.AnnData(counts_numeric.T)
    sc.pp.normalize_total(adata, target_sum=1e6)
    sc.pp.log1p(adata)
    return pd.DataFrame(adata.X.T, index=counts_numeric.index, columns=counts_numeric.columns)

def main():
    parser = argparse.ArgumentParser(description="Aggregate featureCounts count files and produce counts, TPM, log2TPM, and VST matrices")
    parser.add_argument("-c", "--countfiles", nargs='+', required=True,
                        help="List of featureCounts count files (per sample)")
    parser.add_argument("-g", "--gtf", required=True, help="GTF annotation file")
    parser.add_argument("-o", "--outdir", required=True, help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    print("Parsing GTF for gene mapping and gene lengths...")
    gene_map, gene_length = parse_gtf(args.gtf)
    gene_map_df = pd.DataFrame.from_dict(gene_map, orient='index', columns=['GeneName']).reset_index()
    gene_map_df = gene_map_df.rename(columns={'index': 'EnsemblID'})
    gene_length_series = pd.Series(gene_length, name='Length')

    print("Loading count files and merging...")
    all_counts = [load_counts(f) for f in args.countfiles]
    combined_counts = pd.concat(all_counts, axis=1, sort=True).fillna(0).astype(int)
    combined_counts = combined_counts.reset_index().rename(columns={'index':'EnsemblID'})

    print("Merging counts with gene names...")
    combined_counts_annot = combined_counts.merge(gene_map_df, on='EnsemblID', how='left')
    combined_counts_annot['GeneName'] = combined_counts_annot['GeneName'].fillna(combined_counts_annot['EnsemblID'])

    # Raw counts by Ensembl ID
    print("Saving raw counts by Ensembl ID...")
    ensembl_counts = combined_counts.set_index('EnsemblID')
    ensembl_counts.to_csv(os.path.join(args.outdir, "counts_by_ensembl.txt"), sep="\t")

    # Raw counts aggregated by gene name
    print("Saving raw counts by Gene Name (aggregated)...")
    counts_by_gene = combined_counts_annot.groupby('GeneName').sum()
    counts_by_gene.drop('EnsemblID', axis=1).to_csv(os.path.join(args.outdir, "counts_by_gene.txt"), sep="\t")

    # Calculate TPM by Ensembl ID
    print("Calculating TPM by Ensembl ID...")
    tpm_ensembl = calculate_tpm(ensembl_counts, gene_length_series)
    tpm_ensembl.to_csv(os.path.join(args.outdir, "tpm_by_ensembl.txt"), sep="\t")

    # Calculate TPM by gene name by summing TPMs of Ensembl IDs per gene
    print("Calculating TPM by Gene Name by summing TPMs...")
    tpm_ensembl_reset = tpm_ensembl.reset_index().merge(gene_map_df, on='EnsemblID', how='left')
    tpm_ensembl_reset['GeneName'] = tpm_ensembl_reset['GeneName'].fillna(tpm_ensembl_reset['EnsemblID'])
    tpm_by_gene = tpm_ensembl_reset.groupby('GeneName').sum()
    tpm_by_gene.drop('EnsemblID', axis=1).to_csv(os.path.join(args.outdir, "tpm_by_gene.txt"), sep="\t")

    # Log2 transform TPM + 1 for Ensembl IDs
    print("Calculating log2(TPM + 1) for Ensembl ID...")
    log2_tpm_ensembl = pd.DataFrame(
        np.log2(tpm_ensembl.select_dtypes(include=['number']) + 1),
        index=tpm_ensembl.index,
        columns=tpm_ensembl.columns
    )
    log2_tpm_ensembl.to_csv(os.path.join(args.outdir, "log2tpm_by_ensembl.txt"), sep="\t")

    # Log2 transform TPM + 1 for Gene Names
    print("Calculating log2(TPM + 1) for Gene Name...")
    log2_tpm_gene = pd.DataFrame(
        np.log2(tpm_by_gene.select_dtypes(include=['number']) + 1),
        index=tpm_by_gene.index,
        columns=tpm_by_gene.columns
    )
    log2_tpm_gene.drop('EnsemblID', axis=1).to_csv(os.path.join(args.outdir, "log2tpm_by_gene.txt"), sep="\t")

    # VST (approximate) for Ensembl counts
    print("Calculating VST (approximate) for Ensembl IDs using scanpy...")
    vst_ensembl = vst_transform(ensembl_counts)
    vst_ensembl.to_csv(os.path.join(args.outdir, "vst_by_ensembl.txt"), sep="\t")

    # VST (approximate) for Gene counts
    print("Calculating VST (approximate) for Gene Names using scanpy...")
    counts_by_gene_numeric = counts_by_gene.select_dtypes(include=['number'])
    vst_gene = vst_transform(counts_by_gene_numeric)
    vst_gene.to_csv(os.path.join(args.outdir, "vst_by_gene.txt"), sep="\t")

    print("All files saved successfully!")

if __name__ == "__main__":
    main()
