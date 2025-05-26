#!/bin/bash
#SBATCH --chdir=/SCRATCH/AGR145/alejandrogb/
#SBATCH --nodes=1
#SBATCH --partition=albaicin
#SBATCH --mem=8G
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --job-name=sortedSAM_Zucker
#SBATCH --error=/SCRATCH/AGR145/alejandrogb/%j.err
#SBATCH --output=/SCRATCH/AGR145/alejandrogb/%j.out
#SBATCH --time=10:00:00

# LOAD MODULES
module load gnu9
module load SAMtools

for i in {1..16}
do
    SAM=Zucker_S${i}.sam
    BAM=Zucker_S${i}.bam
    SORTED_BAM=Zucker_S${i}_sorted.bam

    # Ordenar BAM
    samtools sort -o $SORTED_BAM $BAM

    # Indexar BAM ordenado
    samtools index $SORTED_BAM
done