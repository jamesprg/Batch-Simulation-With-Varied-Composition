#!/usr/bin/python
import os
import shutil
import subprocess
import itertools
import numpy as np
import MDAnalysis

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
systems[:,-1]=round(sol+1)
systems=systems.astype('int32')
outfile = open("system.txt", "w")
volume=6.1
count=0
control=[1,2,5]
#use seq to sweep through 
for n,i in enumerate(systems):
	comps=['DEMNUM','HEMNUM','EGMNUM','SOLNUM']
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
		#os.system("rm -r %s%s"%(comps[count][:3],system_n))
		os.system("rm -r %s%s/prod/analysis"%(comps[count][:3],system_n))
		shutil.copytree('../%s/prod/analysis'%control_directory, '%s%s/prod/analysis'%(comps[count][:3],system_n))
		outfile.flush()
		continue	
	outfile.flush()
		
	if os.path.isdir('%s%s/'%(comps[count][:3],system_n))==True:
		outfile.write('%s%s has been found, checking for completed analysis \n'%(comps[count][:3],system_n))
		if os.path.isfile('%s%s/prod/analysis/driver/cn_group_COLVAR'%(comps[count][:3],system_n))==False:
			outfile.write('did not find colvar file')
			if os.path.isdir('%s%s/prod/analysis/'%(comps[count][:3],system_n))==True:
				os.system("rm -r %s%s/prod/analysis/"%(comps[count][:3],system_n))
				outfile.write('removed analysis \n')
			shutil.copytree('../template_nosty/prod/analysis', '%s%s/prod/analysis'%(comps[count][:3],system_n))	
		else:
			outfile.write('%s%s has not been started, moving to next system \n'%(comps[count][:3],system_n))
			continue
	os.chdir('%s%s/prod/'%(comps[count][:3],system_n))
	outfile.flush()
	#make centered.gro and use for the rest of the script
	p = subprocess.Popen(['gmx_mpi trjconv -f traj_comp.xtc -center -pbc mol -o centered.xtc -s topol.tpr -dump 0'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'0\n0\n')
	p = subprocess.Popen(['gmx_mpi trjconv -f centered.xtc -center -pbc mol -o centered.gro -s topol.tpr -dump 0'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'0\n0\n')
	p = subprocess.Popen(['gmx_mpi trjconv -f traj_comp.xtc -center -pbc mol -o centered.xtc -dump 0'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'0\n0\n')
	p = subprocess.Popen(['gmx_mpi trjconv -f centered.xtc -center -pbc mol -o centered.gro -s -dump 0'],stdin=subprocess.PIPE,shell=True)
	p.communicate(b'0\n0\n')				
	#move into dumpass
	os.chdir('analysis/driver/dumpmass/')
	subprocess.call(['gmx_mpi grompp -f npt.mdp -c ../../../confout.gro -p ../../../../equil/topol.top -maxwarn 1'],shell=True)
	subprocess.call(['gmx_mpi mdrun -plumed dumpmass.dat'],shell=True)
	os.chdir('../')
	grofilename = '../../centered.gro'
	u = MDAnalysis.Universe(grofilename)
	#switching based
        #CN and grouped FMs
	pfile = open('plumed_cn_group.dat', "w")
	compsp = np.unique(u.select_atoms('not resname ISO and not resname K and not resname CL and not resname NA').residues.resnames)
	full_toxin=u.select_atoms('resname ISO and not type H')
	full_toxin=np.array2string(full_toxin.ids, separator=',')[1:-1]
	pfile.write("###############################################\n")
	pfile.write('############## DEFINE RESIDUES ###############\n')
	pfile.write('###############################################\n')
	pfile.write('toxin: GROUP ATOMS=%s\n'%full_toxin.replace(" ",""))
	for fm in compsp:
		fm_set=u.select_atoms('resname %s'%fm)
		temp,temp_COM=[],[]
		for r in np.unique(fm_set.resids):
			r_set=u.select_atoms('resname %s and resid %s and not type H'%(fm,r))
			r_set=np.array2string(r_set.ids, separator=',')[1:-1]
			r_set=r_set.replace(" ","")
			r_set=r_set.replace("\n","")
			pfile.write('%s%s: GROUP ATOMS=%s\n'%(fm,r,r_set))
			if r != np.unique(fm_set.resids)[-1]:
				temp.append('%s%s,'%(fm,r))
			else:
				temp.append('%s%s\n'%(fm,r))
		pfile.write('%s_f: GROUP ATOMS=%s'%(fm,''.join(temp)))
	pfile.write('###############################################\n')
	pfile.write('############## DEFINE COLVARS ###############\n')
	pfile.write('###############################################\n')
	for n,fm in enumerate(compsp):
#		pfile.write('%s_coord_COM: COORDINATION GROUPA=%s_COM GROUPB=toxin_COM R_0=0.4 NN=60 MM=120 \n'%(fm,fm))
		pfile.write('%s_coord_strong: COORDINATION GROUPA=%s_f GROUPB=toxin R_0=0.4 NN=60 MM=120 \n'%(fm,fm))
		pfile.write('%s_coord: COORDINATION GROUPA=%s_f GROUPB=toxin R_0=0.4  \n'%(fm,fm))
		for x in range(len(compsp)):
			if x != n:
				pfile.write('%s_%s_coord: COORDINATION GROUPA=%s_f GROUPB=%s_f R_0=0.4\n'%(fm,compsp[x],fm,compsp[x]))
	pfile.write('PRINT FILE=cn_group_COLVAR ARG=* STRIDE=5')
	pfile.close()

	#COM based min distance
	pfile = open('plumed_com_min_dist.dat', "w")
	compsp = np.unique(u.select_atoms('not resname ISO and not resname K and not resname CL and not resname NA').residues.resnames)
	full_toxin=u.select_atoms('resname ISO and not type H')
	full_toxin=np.array2string(full_toxin.ids, separator=',')[1:-1]
	pfile.write("###############################################\n")
	pfile.write('############## DEFINE RESIDUES ###############\n')
	pfile.write('###############################################\n')
	pfile.write('toxin: COM ATOMS=%s\n'%full_toxin.replace(" ",""))
	for fm in compsp:
		fm_set=u.select_atoms('resname %s'%fm)
		temp=[]
		for r in np.unique(fm_set.resids):
			r_set=u.select_atoms('resname %s and resid %s and not type H'%(fm,r))
			r_set=np.array2string(r_set.ids, separator=',')[1:-1]
			r_set=r_set.replace(" ","")
			r_set=r_set.replace("\n","")
			pfile.write('%s%s: COM ATOMS=%s\n'%(fm,r,r_set))
			if r != np.unique(fm_set.resids)[-1]:
			    temp.append('%s%s,'%(fm,r))
			else:
			    temp.append('%s%s\n'%(fm,r))
		pfile.write('%s_f: GROUP ATOMS=%s'%(fm,''.join(temp)))
	pfile.write('###############################################\n')
	pfile.write('############## DEFINE COLVARS ###############\n')
	pfile.write('###############################################\n')
	for n,fm in enumerate(compsp):
	 	pfile.write('%s_full_d: DISTANCES GROUPA=%s_f GROUPB=toxin MIN={BETA=100.}\n'%(fm,fm))
 		for x in range(len(compsp)):
 			if x != n:
 				pfile.write('%s_%s_d: DISTANCES GROUPA=%s_f GROUPB=%s_f MIN={BETA=100.}\n'%(fm,compsp[x],fm,compsp[x]))
	pfile.write('PRINT FILE=com_min_dist ARG=*.min STRIDE=5')
	pfile.close()
	subprocess.call(['plumed driver --mc dumpmass/mcfile --mf_xtc ../../traj_comp.xtc --timestep 2 --plumed plumed_com_min_dist.dat'],shell=True)
	subprocess.call(['plumed driver --mc dumpmass/mcfile --mf_xtc ../../traj_comp.xtc --timestep 2 --plumed plumed_cn_group.dat'],shell=True)
	outfile.write('finished plumed for %s%s \n'%(comps[count][:3],system_n))
	os.chdir('../../../../')
#	break
outfile.close()	
