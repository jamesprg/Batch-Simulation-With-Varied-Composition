#!/usr/bin/python
import os
import shutil
import subprocess
import itertools
import numpy as np
import MDAnalysis

#convert to mmols and get solvent totals
fms=['SOLNUM','MOHNUM']
fullW=1110.
fullM=494.
values=np.linspace(0,100,6)/100
systems=np.zeros((len(values),len(fms)))
for n,row in enumerate(systems.T):
    if n == 0:
        systems[:,n] = values
        systems[:,n] = systems[:,n]*fullW
    else:
        systems[:,n] = np.flip(values)
        systems[:,n] = systems[:,n]*fullM
control = np.asanyarray((0.25,0.75))
systems = np.vstack((systems,control))
systems[-1]= [systems[-1][0]*fullW,systems[-1][1]*fullM]
systems = np.round(systems)
systems=systems.astype('int32')
systems=np.append(systems[:5],[systems[-1]],axis=0)
outfile = open("system.txt", "w")
volume1=41.0
volume2=volume1/10.
#use seq to sweep through 
for i in systems:
	outfile.write('System number %s \n'%i)
	#Setup Phase
	sys_num=[np.round(i[0]/fullW,2),np.round(i[1]/fullM,2)]
	sys_num=[int(sys_num[0]*100),int(sys_num[1]*100)]
	system_n="_".join(str(x) for x in sys_num)
	if os.path.isdir('trial%s'%system_n)==False:
		outfile.write('did not find directory for trial%s \n'%system_n)
		shutil.copytree('template', 'trial%s'%system_n)	
	elif os.path.isfile('trial%s/equil/anneal/confout.gro'%system_n)==False:
		outfile.write('trial%s did not make it through Annealing step, we will restart this run \n'%system_n)
		#remove directory
		os.system('rm -r trial%s'%system_n)
		shutil.copytree('template', 'trial%s'%system_n)
	else:
		outfile.write('Already finished or running trial%s, moving to next system \n'%system_n)
		continue
	outfile.flush()
	os.chdir('trial%s/setup'%system_n)
	#if neither are 0
	if i[0] != 0 and i[1]!=0:
		for n,f in enumerate(i):
			sed='sed -i "s/'+fms[n]+'/'+str(f)+'/g" packmol.inp'
			subprocess.call([sed],shell=True)
		sed='sed -i "s/V/'+str(volume1)+'/g" packmol.inp'
		subprocess.call([sed],shell=True)
		subprocess.call(['packmol < packmol.inp'],shell=True)	
		subprocess.call(['gmx_mpi editconf -f mips.pdb -box %s %s %s -o ../equil/out.gro'%(volume2,volume2,volume2)],shell=True)
		#Equilibration Phase
		os.chdir('../equil/')
		for n,f in enumerate(i):
			sed='sed -i "s/'+fms[n]+'/'+str(f)+'/g" topol.top'
			subprocess.call([sed],shell=True)
	else:
		if i[0] ==0:
			os.system('cp packmol_nowater.inp packmol.inp')
			sed='sed -i "s/'+fms[1]+'/'+str(i[1])+'/g" packmol.inp'	
			subprocess.call([sed],shell=True)
			sed='sed -i "s/V/'+str(volume1)+'/g" packmol.inp'
			subprocess.call([sed],shell=True)
			subprocess.call(['packmol < packmol.inp'],shell=True)
			subprocess.call(['gmx_mpi editconf -f mips.pdb -box %s %s %s -o ../equil/out.gro'%(volume2,volume2,volume2)],shell=True)
			os.chdir('../equil/')
			os.system('cp topol_nowater.top topol.top')
			for n,f in enumerate(i):
				sed='sed -i "s/'+fms[n]+'/'+str(f)+'/g" topol.top'
				subprocess.call([sed],shell=True)
		else:
			os.system('cp packmol_nomeoh.inp packmol.inp')
			sed='sed -i "s/'+fms[0]+'/'+str(i[0])+'/g" packmol.inp'
			subprocess.call([sed],shell=True)
			sed='sed -i "s/V/'+str(volume1)+'/g" packmol.inp'
			subprocess.call([sed],shell=True)
			subprocess.call(['packmol < packmol.inp'],shell=True)
			subprocess.call(['gmx_mpi editconf -f mips.pdb -box %s %s %s -o ../equil/out.gro'%(volume2,volume2,volume2)],shell=True)
			os.chdir('../equil/')
			os.system('cp topol_nomeoh.top topol.top')
			for n,f in enumerate(i):
				sed='sed -i "s/'+fms[n]+'/'+str(f)+'/g" topol.top'
				subprocess.call([sed],shell=True)
	os.chdir('em/')
	subprocess.call(['gmx_mpi grompp -f em.mdp -c ../out.gro -p ../topol.top -maxwarn 1'],shell=True)
	p = subprocess.Popen(['gmx_mpi genion -s topol.tpr -o ../out.gro -p ../topol.top -np 1 -pname K'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'3\n')
	subprocess.call(['gmx_mpi grompp -f em.mdp -c ../out.gro -p ../topol.top -maxwarn 1'],shell=True)
	subprocess.call(['gmx_mpi mdrun'],shell=True)
	os.chdir('../anneal/')
	subprocess.call(['gmx_mpi grompp -f anneal.mdp -c ../em/confout.gro -p ../topol.top -maxwarn 1'],shell=True)
	subprocess.call(['gmx_mpi mdrun'],shell=True)
	os.chdir('../nptber/')
	subprocess.call(['gmx_mpi grompp -f npt_equil.mdp -c ../anneal/confout.gro -p ../topol.top -maxwarn 1'],shell=True)
	subprocess.call(['gmx_mpi mdrun'],shell=True)
	os.chdir('../../prod/')
	#Production run and submit job
	grofilename = '../equil/nptber/confout.gro'
	u = MDAnalysis.Universe(grofilename)
	pfile = open('plumed.dat', "w")
	full_toxin=u.select_atoms('resname ISO and not type H')
	full_toxin=np.array2string(full_toxin.ids, separator=',')[1:-1]
	pfile.write("###############################################\n")
	pfile.write('############## DEFINE RESIDUES ###############\n')
	pfile.write('###############################################\n')
	pfile.write('toxin_COM: COM ATOMS=%s\n'%full_toxin.replace(" ",""))
	r_set=u.select_atoms('resname STY and not type H')
	r_set=np.array2string(r_set.ids, separator=',')[1:-1]
	r_set=r_set.replace(" ","")
	r_set=r_set.replace("\n","")
	pfile.write('STY: COM ATOMS=%s\n'%r_set)
	pfile.write('###############################################\n')
	pfile.write('############## DEFINE COLVARS ###############\n')
	pfile.write('###############################################\n')
	pfile.write('d: DISTANCE ATOMS=STY,toxin_COM \n')
	pfile.write('METAD ...\n    ARG=d\n    SIGMA=0.0025\n    HEIGHT=1\n    PACE=500\n    TEMP=300\n    BIASFACTOR=6\n    LABEL=metad\n    GRID_MIN=0.1\n    GRID_MAX=%s\n    WALKERS_MPI\n ... METAD\n\nPRINT STRIDE=1500 ARG=* FILE=COLVAR'%volume2)
	sed='sed -i "s/JNAME/T-FM'+str(system_n)+'/g" GROMACS.slurm'
	subprocess.call([sed],shell=True)
	for j in range(0,28):
		subprocess.call(['gmx_mpi grompp -f nvt.mdp -c ../equil/nptber/confout.gro -p ../equil/topol.top -o topol%s.tpr'%j],shell=True)
	subprocess.call(['sbatch -p ckpt -A pfaendtner-ckpt --ntasks=28 ./GROMACS.slurm'],shell=True)
	outfile.write('submitting job for %s \n'%i)
	os.chdir('../../')
outfile.close()
