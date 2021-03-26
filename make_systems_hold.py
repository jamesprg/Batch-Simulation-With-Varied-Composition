#!/usr/bin/python
import os
import shutil
import subprocess
import numpy as np
import MDAnalysis 

os.system('rm slurm*out')

fms=['DEMNUM','HEMNUM','EGMNUM','SOLNUM']
max_fm=3
values=np.linspace(1,max_fm,max_fm)
temp=[]
temp.append([list((x, y, z)) for x in values for y in values for z in values])
temp = np.asarray(temp[0])
# add a column of zeros to it:
temp = np.hstack((temp,np.zeros((temp.shape[0],1))))
sol = (20/18.02)*1000
temp[:,-1]=round(sol+43+1200)
temp=temp.astype('int32')
volume=6.1
outfile = open("status.txt", "w")
#use seq to sweep through 
for n,i in enumerate(temp):
	fms=['DEMNUM','HEMNUM','EGMNUM','SOLNUM']
	outfile.write('System number %s \n'%i)
	#Setup Phase
	sys_num=[]
	for v in i[0:3]:
		sys_num.append(int(v))
	sys_num=''.join(str(i) for i in sys_num)
	#os.chdir('%s'%fms[count][:3])
	#now we check if the directory exists and if not, copy the template from one above	
	if os.path.isdir('d_h_e_%s'%sys_num)==False:
		outfile.write('did not find directory for d_h_e_%s \n'%sys_num)
		shutil.copytree('template', 'd_h_e_%s'%sys_num)	
	elif os.path.isfile('d_h_e_%s/equil/anneal/confout.gro'%sys_num)==False:
		outfile.write('d_h_e_%s did not make it through Annealing step, we will restart this run \n'%sys_num)
		#remove directory
		os.system('rm -r d_h_e_%s'%sys_num)
		shutil.copytree('template', 'd_h_e_%s'%sys_num)
	else:
		outfile.write('Already finished or running d_h_e_%s, moving to next system \n'%sys_num)
		continue
	outfile.flush()
	os.chdir('d_h_e_%s/setup'%sys_num)
	for fn,f in enumerate(i):
		sed='sed -i "s/'+fms[fn]+'/'+str(f)+'/g" packmol.inp'
		subprocess.call([sed],shell=True)
	subprocess.call(['packmol < packmol.inp'],shell=True)	
	subprocess.call(['gmx_mpi editconf -f mips.pdb -box %s %s %s -o ../equil/out.gro'%(volume,volume,volume)],shell=True)
	#Equilibration Phase
	os.chdir('../equil/')
	for fn,f in enumerate(i):
		sed='sed -i "s/'+fms[fn]+'/'+str(f)+'/g" topol.top'
		subprocess.call([sed],shell=True)
	os.chdir('em/')
	subprocess.call(['gmx_mpi grompp -f em.mdp -c ../out.gro -p ../topol.top -maxwarn 1'],shell=True)
	p = subprocess.Popen(['gmx_mpi genion -s topol.tpr -o ../out.gro -p ../topol.top -pname K -np 1'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'7\n')
	subprocess.call(['gmx_mpi grompp -f em.mdp -c ../out.gro -p ../topol.top -maxwarn 1'],shell=True)
	p = subprocess.Popen(['gmx_mpi genion -s topol.tpr -o ../out.gro -p ../topol.top -pname NA -np 21 -nname CL -nn 21'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'8\n')
	subprocess.call(['gmx_mpi grompp -f em.mdp -c ../out.gro -p ../topol.top -maxwarn 1'],shell=True)
	subprocess.call(['gmx_mpi mdrun'],shell=True)
	os.chdir('../anneal/')
	subprocess.call(['gmx_mpi grompp -f anneal.mdp -c ../em/confout.gro -p ../topol.top -maxwarn 1'],shell=True)
	#make plumed file
	grofilename = '../em/confout.gro'
	u = MDAnalysis.Universe(grofilename)
	f = open(f'plumed.dat', "w")
	fms = np.unique(u.select_atoms('not resname ISO and not resname K and not resname SOL and not resname NA and not resname CL').residues)
	full_toxin=u.select_atoms('resname ISO and not type H')
	full_toxin=np.array2string(full_toxin.ids, separator=',')[1:-1]
	#f.write('WHOLEMOLECULES ENTITY0=1-%s\n'%u.atoms[-1].id)
	f.write('toxin: COM ATOMS=%s\n'%full_toxin.replace(" ",""))
	for fm in fms:
		atoms=u.select_atoms('resname %s and resid %s and not type H'%(fm.resname,fm.resid))
		atoms=np.array2string(atoms.ids, separator=',')[1:-1]
		atoms=atoms.replace(" ","")
		atoms=atoms.replace("\n","")
		f.write('%s%s: COM ATOMS=%s\n'%(fm.resid,fm.resname,atoms))
		f.write('%s%s_dist: DISTANCE ATOMS=%s%s,toxin\n'
			%(fm.resid,fm.resname,fm.resid,fm.resname))
		#change per FM
		if fm.resname == 'DEM':	
			f.write('UPPER_WALLS ARG=%s%s_dist AT=0.52 KAPPA=500000.0 EXP=2 OFFSET=0 LABEL=uwall_%s%s\n'
			%(fm.resid,fm.resname,fm.resid,fm.resname))
		elif fm.resname == 'HEM':
			f.write('UPPER_WALLS ARG=%s%s_dist AT=0.48 KAPPA=500000.0 EXP=2 OFFSET=0 LABEL=uwall_%s%s\n'
			%(fm.resid,fm.resname,fm.resid,fm.resname))
		elif fm.resname == 'EGM':
			f.write('UPPER_WALLS ARG=%s%s_dist AT=0.54 KAPPA=500000.0 EXP=2 OFFSET=0 LABEL=uwall_%s%s\n'
			%(fm.resid,fm.resname,fm.resid,fm.resname))
	f.write('PRINT FILE=COLVAR ARG=* STRIDE=100')
	f.close()	
	subprocess.call(['gmx_mpi mdrun -plumed plumed.dat'],shell=True)
	os.chdir('../nptber/')
	os.system('cp ../anneal/plumed.dat .')
	subprocess.call(['gmx_mpi grompp -f npt_equil.mdp -c ../anneal/confout.gro -p ../topol.top -maxwarn 1'],shell=True)
	subprocess.call(['gmx_mpi mdrun -plumed plumed.dat'],shell=True)
	os.chdir('../../prod/')
	#Production run and submit job
	os.system('cp ../equil/anneal/plumed.dat .')
	#now grompp and run
	subprocess.call(['gmx_mpi grompp -f nvt.mdp -c ../equil/nptber/confout.gro -p ../equil/topol.top -maxwarn 1'],shell=True)
	sed='sed -i "s/JNAME/dhe_'+str(sys_num)+'/g" GROMACS.slurm'
	subprocess.call([sed],shell=True)
	subprocess.call(['sbatch -p ckpt -A pfaendtner-ckpt --ntasks=14 ./GROMACS.slurm'],shell=True)
	outfile.write('submitting job for d_h_e_%s \n'%sys_num)
	os.chdir('../../')
	break
outfile.close()
