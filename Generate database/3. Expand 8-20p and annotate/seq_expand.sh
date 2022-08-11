#!/bin/bash

#SBATCH --job-name=seq_expand_1sub
#SBATCH --time=10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --output fuzz_seq_expand2.out

module load R-bundle-Bioconductor/3.14-foss-2021b-R-4.1.2
module load Biopython/1.79-foss-2021b

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
input="expandSeq/Input/fuzzy.txt"
distance=1000

chmod 775 "$input"

time python3 ./expandSeq/sequence_expand.py $input $distance

echo
echo -- End Job --
date
echo

