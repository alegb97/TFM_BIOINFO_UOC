#!/bin/bash
#SBATCH --chdir=/SCRATCH/AGR145/alejandrogb/
#SBATCH --nodes=1
#SBATCH --partition=NOParalela
#SBATCH --mem=8G
#SBATCH --ntasks=28
#SBATCH --ntasks-per-node=28
#SBATCH --job-name=Alinemaiento_Zucker
#SBATCH --error=/SCRATCH/AGR145/alejandrogb/%j.err
#SBATCH --output=/SCRATCH/AGR145/alejandrogb/%j.out
#SBATCH --time=03:00:00
/SCRATCH/AGR145/alejandrogb
 
source /home/AGR145/alejandrogb/.bashrc

# LOAD MODULES
module load gnu9

# COMMANDS:
####reference genome

IDX=/home/AGR145/alejandrogb/master/rn6/genome

for i in {1..16}
do
    R1=/home/AGR145/alejandrogb/master/Zucker/reads/Sample${i}/NGS014_24_*_S${i}_L005_R1_001.fastq.gz
    R2=/home/AGR145/alejandrogb/master/Zucker/reads/Sample${i}/NGS014_24_*_S${i}_L005_R2_001.fastq.gz
    SAM=Zucker_S${i}.sam
	hisat2 -p 28 --dta -x $IDX -1 $R1 -2 $R2 -S $SAM 
done