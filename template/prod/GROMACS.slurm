#!/bin/bash
#SBATCH -J NAME
#SBATCH --nodes=1
#SBATCH --mem=120GB
#SBATCH --time=36:00:00

## Sepcify the working directory for this job
cd $SLURM_SUBMIT_DIR

imodule load icc_17-impi_2017
module load cuda/9.1.85.3
source /gscratch/pfaendtner/sarah/codes/gromacs18.3/gromacs-2018.3/bin/bin/GMXRC

gmx_mpi mdrun -cpi -append -cpt 1 &>log.txt


exit 0
