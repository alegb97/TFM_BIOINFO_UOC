#!/bin/bash
#SBATCH --chdir=/SCRATCH/AGR145/alejandrogb/
#SBATCH --nodes=1
#SBATCH --partition=albaicin
#SBATCH --mem=8G
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --job-name=SAM_to_BAM_DIO
#SBATCH --error=/SCRATCH/AGR145/alejandrogb/%j.err
#SBATCH --output=/SCRATCH/AGR145/alejandrogb/%j.out
#SBATCH --time=10:00:00

# LOAD MODULES
module load gnu9
module load SAMtools

for i in {1..24}
do
    SAM=DIO_S${i}.sam
    BAM=DIO_S${i}.bam
    SORTED_BAM=DIO_S${i}_sorted.bam

    # Convertir SAM a BAM
    samtools view -bo $BAM $SAM

    # Ordenar BAM
    samtools sort -o $SORTED_BAM $BAM

    # Indexar BAM ordenado
    samtools index $SORTED_BAM
done