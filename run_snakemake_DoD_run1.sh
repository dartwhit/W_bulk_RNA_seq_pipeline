#!/bin/bash
#SBATCH --job-name=snakemake_DoD
#SBATCH --output=snakemake_pipeline.%j.out
#SBATCH --error=snakemake_pipeline.%j.err
#SBATCH --time=12:00:00          # Max runtime hh:mm:ss
#SBATCH --cpus-per-task=4        # CPUs for the submission job itself
#SBATCH --mem=16G                # Memory for the submission job
#SBATCH --partition=standard     # Partition/queue name (adjust accordingly)
#SBATCH --mail-type=END,FAIL     # Email on job done or failed
#SBATCH --mail-user=f005c3n@dartmouth.edu


# Run snakemake with cluster submission for individual rules
snakemake \
    --configfile config_DoD_run1.yaml \
    --jobs 10 \
    --latency-wait 60 \
    --rerun-incomplete \
    --keep-going \
    --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem} --time={resources.time} --partition=standard" \
    --printshellcmds
