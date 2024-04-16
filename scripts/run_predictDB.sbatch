#!/bin/bash
#SBATCH --job-name=PredictDB_JOBNAME
#SBATCH --account=pi-haky
#SBATCH --partition=caslake
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --time=18:00:00
#SBATCH --mem=2G


module unload python
module unload R
module load python
module load R
source ~/.bashrc
conda activate my_predictdb

nextflow run PREDICTDB_REPO/main.nf \
--gene_annotation GENE_ANNOTATION \
--snp_annotation SNP_ANNOTATION \
--genotype GENOTYPE \
--gene_exp EXPRESSION \
--outdir OUTPUT_DIR \
--nfolds NFOLDS \
--keepIntermediate \
-resume \
-profile midway3

