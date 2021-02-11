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
systems=systems.astype('int32')
outfile = open("rdf_system.txt", "w")
count=0
control=[1,2,5]
#use seq to sweep through 
for n,i in enumerate(systems):
	outfile.write('System number %s \n'%i)	
	if n == 0:
		for f in comps[:-1]:
			if os.path.isdir('%s'%f[:3])==False:
				outfile.write('did not find directory for component %s \n'%f[:3])
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
		print('%s%s is a duplicate control, skipping this step and copying previous'%(comps[count][:3],system_n))
		os.system('cp ../%s/prod/*xvg %s%s/prod/'%(control_directory,comps[count][:3],system_n))
		outfile.flush()
		continue	
	outfile.flush()
		
	if os.path.isdir('%s%s/'%(comps[count][:3],system_n))==True:
		outfile.write('%s%s has been found, checking for rdf\n'%(comps[count][:3],system_n))
		if os.path.isfile('%s%s/prod/sol.xvg'%(comps[count][:3],system_n))==False:
			outfile.write('did not find final sol.xvg, starting analysis\n')
		else:
			outfile.write('%s%s has finished rdfs, moving to next system \n'%(comps[count][:3],system_n))
			continue
	else:
		outfile.write('%s%s has not been found, moving to next system\n'%(comps[count][:3],system_n))
		continue
	os.chdir('%s%s/prod/'%(comps[count][:3],system_n))
	outfile.flush()
	#make index file that has the interactions of interest
	p = subprocess.Popen(['gmx_mpi make_ndx -f topol.tpr'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'2&aH6\n2&tO\n4&aH8\n3&aO1\n5&aO2\n5&aO4\n25|26\n10&tHW\nq\n')
	outfile.write('finished index file for %s%s \n'%(comps[count][:3],system_n))
	p = subprocess.Popen(['gmx_mpi rdf -f traj_comp.xtc -s confout.gro -n index.ndx -o hem.xvg'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'22\n23\n')
	os.chdir('../../')
outfile.close()	
