#Usage: time ipython ~/Desktop/prj/bsff/bsff_scripts/scf2.py fragment2_contact_fuzzball.pdb
#intended to be run on the pdb files generated in the fragment_x directories
#within Transformed_Aligned_pdbs dir from the 'align'command in BSFF protocol
#returns a single pdb file containing all residues with at least one atom
#within 4.0 angstroms of the ligand fragment
#cd desktop/prj/bsff/compounds/pac/Transformed_Aligned_PDBs/Fragment_2

import sys
import math
import os
import pandas as pd
import collections
from scipy.spatial import distance
#get the names of the pdb files in working directory
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']

#storing ligand coordinates from one of the pdb files in list fragment_coords
#storing full ligand lines in list liglines
#storing contact residue lines in list fuzzball_lines
liglines=[]
fragment_coords=[]
f=open(pdbfiles[0],'r') #get lig lines
for line in f.readlines():
    if line[0:6]=='HETATM':
        liglines.append(line)
        x=float(line[31:38]) ;y=float(line[38:46]) ;z=float(line[46:54])
        fragment_coords.append((x,y,z))
f.close()
fuzzball_lines=[]
for file in pdbfiles: #
    f=open(file,'r')
    atomlines=[line for line in f.readlines() if line[0:4]=='ATOM'];f.close()
    for i in atomlines:
        fuzzball_lines.append(i)

#first gotta get indices of unique residues in fuzzball_lines
fuzzball_residue_indices=[];fuzzball_seen_resnums=[]
starts=[];lasts=[]
for index, line in enumerate(fuzzball_lines): #collecting start/end indices for unique res
        resnum=int(fuzzball_lines[index][22:26])
        resname=line[17:20]
        try:
            lastresnum=int(fuzzball_lines[index-1][22:26])
            lastresname=fuzzball_lines[index-1][17:20]
            if resnum!=lastresnum or resname!=lastresname:
                start=index
                starts.append(start)
        except:
            start=index
            starts.append(start)
        try:
            nextresname=fuzzball_lines[index+1][17:20]
            next_resnum=int(fuzzball_lines[index+1][22:26])
            if resnum!=next_resnum or resname!=nextresname:
                last=index+1
                lasts.append(last)
        except:
            last=len(fuzzball_lines)
            lasts.append(last)
for index,start in enumerate(starts): #put the indices together for each res
    fuzzball_residue_indices.append((start,lasts[index]))
fuzzball_residue_indices=list(set(fuzzball_residue_indices)) #ensure no redundant indices
fuzzball_residue_indices=sorted(fuzzball_residue_indices, key=lambda first: first[0])
og_n_res=str(len(fuzzball_residue_indices))


#filter out residues that suck
def filters(fuzzball_lines,fuzzball_residue_indices):
    startnres=len(fuzzball_residue_indices)
    for a,b in fuzzball_residue_indices:
        l=b-a
        if l<4:
            fuzzball_residue_indices.remove((a,b))
    #now extend to requirements for individual residues to have complete sidechains
    for a,b in fuzzball_residue_indices:
        resname=fuzzball_lines[a][17:20]
        rlines=fuzzball_lines[a:b]
        for line in rlines:
            if line[11:16].split()[0][0]=='H' or line[11:16].split()[0]=='OXT':
                rlines.remove(line)
        l=len(rlines) #this is to make sure that im only looking at non hydrogens
        if resname=='GLY':
            pass
        elif resname=='ALA':
            if l<5:
                fuzzball_residue_indices.remove((a,b))
        elif resname=='SER' or resname=='CYS':
            if l<6:
                fuzzball_residue_indices.remove((a,b))
        elif resname=='VAL' or resname=='THR' or resname=='PRO':
            if l<7:
                fuzzball_residue_indices.remove((a,b))
        elif resname=='ILE' or resname=='LEU' or resname=='MET' or resname=='ASN' or resname=='ASP':
            if l<8:
                fuzzball_residue_indices.remove((a,b))
        elif resname=='GLN' or resname=='LYS' or resname=='GLU':
            if l<9:
                fuzzball_residue_indices.remove((a,b))
        elif resname=='HIS':
            if l<10:
                fuzzball_residue_indices.remove((a,b))
        elif resname=='PHE' or resname=='ARG':
            if l<11:
                fuzzball_residue_indices.remove((a,b))
        elif resname=='TYR':
            if l<12:
                fuzzball_residue_indices.remove((a,b))
        elif resname=='TRP':
            if l<14:
                fuzzball_residue_indices.remove((a,b))
    #remove residues with zero occupancy or high b factor (>60 A^2)
    for a,b in fuzzball_residue_indices:
        lines=fuzzball_lines[a:b]
        for line in lines:
            occupancy=line[54:60]
            temp_factor=line[60:66]
            if float(occupancy)==0. or float(temp_factor)>60.0:
                fuzzball_residue_indices.remove((a,b))
                break
    end_nres=len(fuzzball_residue_indices)
    for a,b in fuzzball_residue_indices: #remove res without atom win 4 angs
        lines=fuzzball_lines[a:b]
        res_dists=[]
        for line in lines:
            atom_coords=(float(line[31:39].strip()),float(line[39:47].strip()),float(line[47:54].strip()))
            for frag_atom in fragment_coords:
                d=distance.euclidean(atom_coords,frag_atom)
                res_dists.append(d)
        if min(res_dists)<=4.1:
            pass
        else:
            fuzzball_residue_indices.remove((a,b))
    return startnres, end_nres

#for some reason the filters remove more and more residues
#upon additional application, so i am applying them over and over until the
#number of remaining res converges
#idk why this is happening and it is extremekly upsetting
#because it could mean something is horribly wrong here
#HOWEVER the fact that the number of res does converge
#and it converges to the same number every time
#seems to imply that for whatever reason some residues just make it past
#the filters on their first go even though they shouldnt
#so hopefully at the very least I can trust that these remaining residues do
#indeed satisfy the desired criteria
for i in range(100):
    s,n = filters(fuzzball_lines,fuzzball_residue_indices)
    print(s,n)
    if s==n:
        break


#reorder the residue numbers of the filtered contact residues
#it is useful to have this unique identifier for residues
clean_fuzzball_lines=[] #the renumbered lines
for line in liglines:
    newline=line[0:23]+'1'+'     '+line[29:]
    clean_fuzzball_lines.append(newline)
new_residue_number=2
for (a,b) in fuzzball_residue_indices: #go through each residue and edit line w new resnum
    current_residue_lines=[i for i in fuzzball_lines[a:b]]
    for line in current_residue_lines:
        if new_residue_number<10:
            newline=line[0:22]+' '+str(new_residue_number)+'     '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9<new_residue_number<100:
            newline=line[0:22]+' '+str(new_residue_number)+'    '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99<new_residue_number<1000:
            newline=line[0:22]+' '+str(new_residue_number)+'   '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 999<new_residue_number<10000:
            newline=line[0:22]+' '+str(new_residue_number)+'  '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9999<new_residue_number:
            newline=line[0:22]+' '+str(new_residue_number)+' '+line[29:]
            clean_fuzzball_lines.append(newline)
        else:
            print(str(new_residue_number))
    new_residue_number+=1

#print this stuff cus why not
n_files=str(len(os.listdir(os.getcwd())))
print(n_files+' files parsed for contacts')
print(str(len(fuzzball_residue_indices))+' quality contact residues found')
print('out of '+og_n_res+' total residues')

#write the output pdb fuzzball
ofile=open(sys.argv[1],'w')
for line in clean_fuzzball_lines:
    ofile.write(line)
ofile.close()
