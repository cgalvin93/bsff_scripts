#Usage: time ipython ~/Desktop/prj/bsff/bsff_scripts/scf2.py fragment2_contact_fuzzball.pdb
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

##############    re order fuzzball residues from 1 to whatever
############# this is a pain in the ass because contact residues from different
############# files can have the same residue number
liglines=[]
fragment_coords=[]
f=open(pdbfiles[0],'r') #get lig lines
for line in f.readlines():
    if line[0:6]=='HETATM':
        liglines.append(line)
        x=float(line[30:38]) ;y=float(line[38:46]) ;z=float(line[46:54])
        fragment_coords.append((x,y,z))
f.close()
fuzzball_lines=[]
for file in pdbfiles: #
    f=open(file,'r')
    atomlines=[line for line in f.readlines() if line[0:4]=='ATOM'];f.close()
    for i in atomlines:
        fuzzball_lines.append(i)
#first gotta get indices of unique residues
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

clean_fuzzball_lines=[] #the renumbered lines
new_residue_number=1
for (a,b) in fuzzball_residue_indices: #go through each residue and edit line w new resnum
    current_residue_lines=[i for i in fuzzball_lines[a:b]]
    for line in current_residue_lines:
        if new_residue_number<10:
            newline=line[0:23]+str(new_residue_number)+'     '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9<new_residue_number<100:
            newline=line[0:23]+str(new_residue_number)+'    '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 99<new_residue_number<1000:
            newline=line[0:23]+str(new_residue_number)+'   '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 999<new_residue_number<10000:
            newline=line[0:23]+str(new_residue_number)+'  '+line[29:]
            clean_fuzzball_lines.append(newline)
        elif 9999<new_residue_number:
            newline=line[0:23]+str(new_residue_number)+' '+line[29:]
            clean_fuzzball_lines.append(newline)
        else:
            print(str(new_residue_number))
    new_residue_number+=1


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


def pdb_to_df(pdblist):
    l=[]
    for line in pdblist:
        if line[0:4]=='ATOM':
            d={}
            d['recordname']=line[0:6]
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=line[20:22]
            d['resnum']=line[22:29]
            d['achar']=line[29:30]
            d['x']=line[30:38]
            d['y']=line[38:46]
            d['z']=line[46:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']=line[78:80]
            l.append(d)
    df=pd.DataFrame(l)
    f.close()
    return df

def df_to_pdb(dataframe,ofile):
    of=open(ofile,'w')
    for i in range(df.shape[0]):
        line=''.join(df.iloc[i])
        of.write(line)
    of.close()

# mdf=pd.DataFrame(columns=['recordname',
#                           'atomnumber',
#                           'atomname',
#                           'altloc',
#                           'resname',
#                           'chain',
#                           'resnum',
#                           'achar',
#                           'x',
#                           'y',
#                           'z',
#                           'occupancy',
#                           'temp_factor',
#                           'seg',
#                           'element',
#                           'q']
#                     )

# for file in pdbfiles:
#     odf=pdb_to_df(file)
#     mdf=mdf.append(odf,ignore_index=True)

#okay load pdbfile data frame with reordered numbers
df=pdb_to_df(clean_fuzzball_lines)
#find residues with at least one contact atom
contact_dict=collections.defaultdict(list)
for i in range(df.shape[0]):
    ptn_atom_coords=( float(df.iloc[i][8]), float(df.iloc[i][9]), float(df.iloc[i][10]) )
    for ligand_atom_coords in fragment_coords:
        d=dist(ptn_atom_coords, ligand_atom_coords)
        if d<=4.0:
            resnum=only_numerics(df.iloc[i][6])
            resname=df.iloc[i][4]
            atom_name=df.iloc[i][2]
            contact_atom_info=(resname,atom_name,ptn_atom_coords[0],ptn_atom_coords[1],ptn_atom_coords[2])
            if contact_atom_info not in contact_dict[resnum]:
                contact_dict[resnum].append(contact_atom_info)

#write contact lines to df
#why is this so damn slow???
contact_df=pd.DataFrame(columns=['recordname',
                                 'atomnumber',
                                 'atomname',
                                 'altloc',
                                 'resname',
                                 'chain',
                                 'resnum',
                                 'achar',
                                 'x',
                                 'y',
                                 'z',
                                 'occupancy',
                                 'temp_factor',
                                 'seg',
                                 'element',
                                 'q']
                         )

rns=[]
for i in df['resnum']:
    rns.append(int(i))
n_cont=0
for key in contact_dict.keys():
    resnum=int(key)
    rows=[i for i, e in enumerate(rns) if e == resnum]
    if len(rows)>=4:
        contact_df=contact_df.append(df.iloc[min(rows):max(rows)],ignore_index=True)
        n_cont+=1

#write to import pdb
df_to_pdb(contact_df,sys.argv[1])

n_files=str(len(os.listdir(os.getcwd())))
print(n_files+' files parsed for contacts')
print(str(n_cont)+' contact residues found')
