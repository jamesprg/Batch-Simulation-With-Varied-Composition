#!/bin/bash
#SBATCH -J JNAME
#SBATCH --nodes=1
#SBATCH --mem=60GB
#SBATCH --constraint=broadwell
#SBATCH --time=36:00:00

## Sepcify the working directory for this job
cd $SLURM_SUBMIT_DIR

imodule load icc_17-impi_2017
module load cuda/9.1.85.3
source /gscratch/pfaendtner/sarah/codes/gromacs18.3/gromacs-2018.3/bin/bin/GMXRC
#module load icc_19 icc_19-impi_2019 gcc/6.3.1
#source /gscratch/pfaendtner/jpfaendt/codes/gmx2020.5-cpu/bin/GMXRC
source /gscratch/pfaendtner/jpfaendt/codes/plumed2/sourceme.sh

gmx_mpi mdrun -plumed plumed.dat -cpi -append -cpt 1 &>log.txt


exit 0
