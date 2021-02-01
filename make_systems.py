#!/usr/bin/python
import os
import shutil
import subprocess
import itertools
import numpy as np

os.system('rm slurm*out')

comps=['COMP1','COMP2','COMP3','SOLNUM']
max_fm=10
values=np.linspace(1,max_fm,max_fm)
systems=np.zeros((len(comps[:3])*len(values),len(comps)))
control = np.asanyarray((1,2,5))
systems[:,:3] = control
for n,row in enumerate(systems.T[:3]):
	systems[:,n][n*max_fm:(n+1)*max_fm] = values
#convert to mmols and get solvent totals
systems = systems*8
sol = (20/18.02)*1000
systems[:,-1]=round(sol+43+1200)
systems=systems.astype('int32')
outfile = open("system.txt", "w")
volume=6.1
count=0
control=[1,2,5]
#use seq to sweep through 
for n,i in enumerate(systems):
	outfile.write('System number %s \n'%i)
	if n == 0:
		for f in comps[:-1]:
			if os.path.isdir('%s'%f[:3])==False:
				outfile.write('did not find directory for component %s \n'%f[:3])
				os.system('mkdir %s'%f[:3])
			else:
				outfile.write('found directory for component %s \n'%f[:3])
		os.chdir('%s'%comps[count][:3])
	#find out which varible FM we are on
	else:
		if n%(float(len(systems))/float(len(comps)-1))==0:
			count+=1
			os.chdir('../%s'%comps[count][:3])	
	#Setup Phase
	sys_num=[]
	for v in i[0:3]:		
		sys_num.append(int(v/8))
	system_n=int(sys_num[count])
	#check if this is the control being repeated, if so copy from before
	if sys_num==control and n<(float(len(systems))/float(len(comps)-1)):
		outfile.write('Found original control \n')
		control_directory='%s/%s%s'%(comps[count][:3],comps[count][:3],system_n)
	if sys_num==control and n!=0:		
		print('%s%s is a duplicate control, skipping this step and copying previous')
		shutil.copytree('../%s'%control_directory, '%s%s'%(comps[count][:3],system_n))	
		continue
	#now we check if the directory exists and if not, copy the template from one above	
	if os.path.isdir('%s%s'%(comps[count][:3],system_n))==False:
		outfile.write('did not find directory for %s%s \n'%(comps[count][:3],system_n))
		shutil.copytree('../template_nosty', '%s%s'%(comps[count][:3],system_n))	
	elif os.path.isfile('%s%s/equil/anneal/confout.gro'%(comps[count][:3],system_n))==False:
		outfile.write('%s%s did not make it through Annealing step, we will restart this run \n'%(comps[count][:3],system_n))
		#remove directory
		os.system('rm -r %s%s'%(comps[count][:3],system_n))
		shutil.copytree('../template_nosty', '%s%s'%(comps[count][:3],system_n))
	else:
		outfile.write('Already finished or running %s%s, moving to next system \n'%(comps[count][:3],system_n))
		continue
	outfile.flush()
	os.chdir('%s%s/setup'%(comps[count][:3],system_n))
	for n,f in enumerate(i):
		sed='sed -i "s/'+comps[n]+'/'+str(f)+'/g" packmol.inp'
		subprocess.call([sed],shell=True)
	subprocess.call(['packmol < packmol.inp'],shell=True)	
	subprocess.call(['gmx_mpi editconf -f mips.pdb -box %s %s %s -o ../equil/out.gro'%(volume,volume,volume)],shell=True)
	#Equilibration Phase
	os.chdir('../equil/')
	for n,f in enumerate(i):
		sed='sed -i "s/'+comps[n]+'/'+str(f)+'/g" topol.top'
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
	subprocess.call(['gmx_mpi mdrun'],shell=True)
	os.chdir('../nptber/')
	subprocess.call(['gmx_mpi grompp -f npt_equil.mdp -c ../anneal/confout.gro -p ../topol.top -maxwarn 1'],shell=True)
	subprocess.call(['gmx_mpi mdrun'],shell=True)
	os.chdir('../../prod/')
	#Production run and submit job
	subprocess.call(['gmx_mpi grompp -f nvt.mdp -c ../equil/nptber/confout.gro -p ../equil/topol.top -maxwarn 1'],shell=True)
	sed='sed -i "s/JNAME/%s'%comps[count][:3]+str(system_n)+'/g" GROMACS.slurm'
	subprocess.call([sed],shell=True)
	subprocess.call(['sbatch -p ckpt -A pfaendtner-ckpt --ntasks=10 ./GROMACS.slurm'],shell=True)
	outfile.write('submitting job for %s%s \n'%(comps[count][:3],system_n))
	os.chdir('../../')
#	break
outfile.close()
