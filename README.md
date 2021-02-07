# Batch-Simulation-With-Varied-Composition
Python scripts which automate the creation of various molecular dynamics systems using packmol to initialize before proceeding to create gromacs ready files.

### Required Software

- [python](https://www.python.org/) (3.6)
- [packmol](http://m3g.iqm.unicamp.br/packmol/home.shtml) (>=18.169)
- [GROMACS](https://manual.gromacs.org/documentation/2020/download.html) (2020)

The script assumes you are using a slurm scheduler to launch your jobs

### Guide

- __Step 1__ - Determine the different components you wish to vary using `make_systems.py` and set inside the script. The original purpose of the code was to study how changing the amount of molecular components effected their propensity to interact. Furthermore, you may want to adjust things such as box size, solvent amount, solvents used, and ion concentrations based on your system of interest.

- __Step 2__ - Place structure files for each component desired to study into the `template/setup` directory for use when packmol is called.

- __Step 3__ - Add any necessary forcefield paths to `equil/topol.top`.

- __Step 4__ - Run `make_systems.py`in the directory you wish to be a home directory for each components directories


### To Do List

- __1__ - Add functionality for choosing to turn on/off slurm submission

- __2__ - COMPLETE: Add functionality for choosing amount of solvent



### Attribution

When using the code above to conduct original research, please reference this code repository and any of the open source software mentioned above. Thanks!
