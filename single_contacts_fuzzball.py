#Usage: time ipython ~/Desktop/prj/bsff/bsff_scripts/single_contacts_fuzzball.py
#intended to be run on the pdb files generated in the fragment_x directories
#within Transformed_Aligned_pdbs dir from the 'align'command in BSFF protocol
#returns a single pdb file containing all residues with at least one atom
#within 4.0 angstroms of the ligand fragment
#cd desktop/prj/bsff/compounds/pac/Transformed_Aligned_PDBs/Fragment_1

import sys
import math
import os
import pandas as pd
import collections

#get the names of the pdb files in working directory
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']


#get the ligand coordinates and hetatm lines from one of pdb files
fuzzball_lines=[]
fragment_coords=[]
f=open(pdbfiles[0],'r')
for line in f.readlines():
    if line[0:6]=='HETATM':
        fuzzball_lines.append(line)
        x=float(line[30:38].strip()) ;y=float(line[38:46].strip()) ;z=float(line[46:54].strip())
        fragment_coords.append((x,y,z))
f.close()

#function to return only digits from string
def only_numerics(seq):
    seq_type= type(seq)
    return seq_type().join(filter(seq_type.isdigit, seq))

#set of functions to get distance between two points
def displace(p1,p2):
	x = p1[0] - p2[0]
	y = p1[1] - p2[1]
	z = p1[2] - p2[2]
	return (x,y,z)
def norm(x):
    return math.sqrt(sum(i**2 for i in x))
def dist(p1,p2):
	v = displace(p1,p2)
	return norm(v)

#okay now do the things
for file in pdbfiles:
    f=open(file,'r')
    atomlines=[line for line in f.readlines() if line[0:4]=='ATOM'];f.close()
    f.close();seen_resnums=[];residue_indices=[]#start and end indices of each residue in atomlines list
    for index, line in enumerate(atomlines): #go through here and store the indices
        resnum=only_numerics(line[22:29])
        if resnum not in seen_resnums:
            seen_resnums.append(resnum)
            start=index
            try:
                next_resnum=only_numerics(atomlines[index+1][22:29])
                if resnum!=next_resnum:
                    last=index+1
                    residue_indices.append((start,last,resnum))
            except:
                last=len(atomlines)
                residue_indices.append((start,last,resnum))
        else:
            try:
                next_resnum=only_numerics(atomlines[index+1][22:29])
                if resnum!=next_resnum:
                    last=index+1
                    residue_indices.append((start,last,resnum))
            except:
                last=len(atomlines)
                residue_indices.append((start,last,resnum))
    for a,b,c in residue_indices: #go through each residue and check for atoms w/in 4 angstroms
        current_residue_lines=[i for i in atomlines[a:b]]
        try:
            current_res_coords=[(float(line[31:38].strip()),float(line[39:46].strip()),float(line[47:55].strip())) for line in current_residue_lines]
            distances=[]
            for fragment_atom_coords in fragment_coords:
                for current_res_atom_coords in current_res_coords:
                    d=dist(fragment_atom_coords,current_res_atom_coords)
                    distances.append(d)
            if min(distances)<4.0:
                for i in current_residue_lines:
                    fuzzball_lines.append(i)
        except:
            print('problem parsing file: '+ str(file) + ' residue: ' + str(c))


n_files=str(len(os.listdir(os.getcwd())))
print(n_files+' files parsed for contacts')

##############    re order fuzzball residues from 1 to whatever
############# this is a pain in the ass because contact residues from different
############# files can have the same residue number
fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
starts=[];lasts=[]
for index, line in enumerate(fuzzball_lines): #collecting start/end indices for unique res
        resnum=only_numerics(line[22:29])
        resname=line[17:20]
        try:
            lastresnum=only_numerics(fuzzball_lines[index-1][22:29])
            lastresname=fuzzball_lines[index-1][17:20]
            if resnum!=lastresnum or resname!=lastresname:
                start=index
                starts.append(start)
        except:
            start=index
            starts.append(start)
        try:
            nextresname=fuzzball_lines[index+1][17:20]
            next_resnum=only_numerics(fuzzball_lines[index+1][22:29])
            if resnum!=next_resnum or resname!=nextresname:
                last=index+1
                lasts.append(last)
        except:
            last=len(fuzzball_lines)
            lasts.append(last)

for index,start in enumerate(starts): #put the indices together for each res
    fuzzball_residue_indices.append((start,lasts[index]))
fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices

clean_fuzzball_lines=[] #the renumbered lines
for line in fuzzball_lines: #add ligand lines
    if line[0:6]=='HETATM':
        clean_fuzzball_lines.append(line)
new_residue_number=1
for (a,b) in fuzzball_residue_indices: #go through each residue and edit line w new resnum
    current_residue_lines=[i for i in fuzzball_lines[a:b]]
    for line in current_residue_lines:
        if new_residue_number<10:
            newline=line[0:23]+str(new_residue_number)+'       '+line[31:]
            clean_fuzzball_lines.append(newline)
        elif 9<new_residue_number<100:
            newline=line[0:23]+str(new_residue_number)+'      '+line[31:]
            clean_fuzzball_lines.append(newline)
        elif 99<new_residue_number<1000:
            newline=line[0:23]+str(new_residue_number)+'    '+line[31:]
            clean_fuzzball_lines.append(newline)
        elif 999<new_residue_number<10000:
            newline=line[0:23]+str(new_residue_number)+'   '+line[31:]
            clean_fuzzball_lines.append(newline)
        elif 9999<new_residue_number<100000:
            newline=line[0:23]+str(new_residue_number)+'  '+line[31:]
            clean_fuzzball_lines.append(newline)
        else:
            print(str(new_residue_number))
    new_residue_number+=1



#write fuzzball pdb file
fuzzball_pdb=open('contact_fuzzball.pdb','w')
for line in clean_fuzzball_lines:
    fuzzball_pdb.write(line)
fuzzball_pdb.close()
