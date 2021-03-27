#!/bin/bash
#SBATCH -J submit
#SBATCH -p ckpt
#SBATCH -A pfaendtner-ckpt
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=28
#SBATCH --mem=242G
#SBATCH --time=05:30:00

## Sepcify the working directory for this job
cd $SLURM_SUBMIT_DIR

module load icc_19 icc_19-impi_2019 gcc/6.3.1
source /gscratch/pfaendtner/jpfaendt/codes/plumed2/sourceme.sh
source /gscratch/pfaendtner/jpfaendt/codes/gmx2020.5-cpu/bin/GMXRC

packmol < test.inp &> log.txt

exit 0
