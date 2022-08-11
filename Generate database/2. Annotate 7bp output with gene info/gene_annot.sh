#!/bin/bash

#SBATCH --job-name=gene_annot_1sub
#SBATCH --time=48:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --output fuzz_gene_annot.out


module load R/4.1.2-foss-2021b
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2

export OMP_NUM_THREADS=1


date
echo -- Start Job --
chmod -R 775 ../*

# User modifiable parameters
# input="../Output/Homo_sapiens/"
input="/theia/scratch/projects/brightcore/edir/Homo_sapiens"

Rscript ./geneAnnotate/geneAnnotate_optimizedTEST.R $input



echo
echo -- End Job --
date
echo
