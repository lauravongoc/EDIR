#!/bin/bash

#SBATCH --job-name=find_seq_1sub
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output fuzz_find_seq.out


module load R/4.1.2-foss-2021b
module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
module load Biopython/1.79-foss-2021b
module load rpy2/3.4.5-foss-2021b


export OMP_NUM_THREADS=1

echo -- Modules --
ml
echo


cd $SLURM_SUBMIT_DIR

echo submit directory: $PWD
echo jobid: $SLURM_JOB_ID
echo hostname: $HOSTNAME
echo
date
echo -- Start Job --
echo
chmod  -R 775 ../../genRepeats

# User modifiable parameters
input="./seqFind/Input/fuzzy_sequences_7bp.txt"
dist=1000


time python3 ./seqFind/find_seq_element_occurrences_in_seq.py $input $dist

echo
echo -- End Job --
date
echo
