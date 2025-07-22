import re
import numpy as np
import pandas as pd

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

def main():
    print("Parsing GTF for gene mapping and gene lengths...")
    gene_map, gene_length = parse_gtf('ref_genomes/hg38/Homo_sapiens.GRCh38.110.gtf')
    # print(gene_length)
    gene_map_df = pd.DataFrame.from_dict(gene_map, orient='index', columns=['GeneName']).reset_index()
    gene_map_df = gene_map_df.rename(columns={'index': 'EnsemblID'})
    # print(gene_map_df.head())
    gene_length_series = pd.Series(gene_length, name='Length')
    print(gene_length_series.head())
    
if __name__ == "__main__":
    main()
