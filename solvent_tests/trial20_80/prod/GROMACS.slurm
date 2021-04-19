#!/bin/bash
#SBATCH -J T-FM20_80
#SBATCH --nodes=1
#SBATCH --mem=120GB
#SBATCH --time=36:00:00

## Sepcify the working directory for this job
cd $SLURM_SUBMIT_DIR

imodule load icc_17-impi_2017
module load cuda/9.1.85.3
source /gscratch/pfaendtner/sarah/codes/gromacs18.3/gromacs-2018.3/bin/bin/GMXRC
#source /gscratch/pfaendtner/jpfaendt/codes/plumed2/sourceme.sh
source /gscratch/pfaendtner/cdf6gc/codes/plumed/luiz/activate_plumed24.sh

mpiexec.hydra -np 28 gmx_mpi mdrun -ntomp 1 -s topol.tpr -o traj.trr -x traj.xtc -cpi restart -cpo restart -c confout.gro -e ener.edr -g md.log -plumed plumed.dat -multi 28 &> log.txt 



exit 0
