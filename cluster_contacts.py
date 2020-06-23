
# ~Desktop/prj/BSFF/bsff_scripts/cluster_contacts.py
# ~Desktop/prj/BSFF/output/contactfuzzballs/contactsubset.pdb
# USAGE: time python cluster_contacts.py path/to/contact_fuzzball.pdb

import sys
import pandas as pd
import math
import collections

#function to return only digits from string
def only_numerics(seq):
    seq_type= type(seq)
    return seq_type().join(filter(seq_type.isdigit, seq))

#store ligand atom coordinates in list,
#store protein residue info in pandas df resnum, resname, atom name, coords(xyz)
ligand_atoms=[]
protein_atoms=[]
contact_fuzzball=open('Fragment2_contacts.pdb','r') #################testing########
# contact_fuzzball=open(sys.argv[1],'r')
for line in contact_fuzzball.readlines():
    if line[0:6]=='HETATM':
        x=float(line[31:38].strip()) ;y=float(line[39:46].strip()) ;z=float(line[47:55].strip())
        ligand_atoms.append((x,y,z))
    elif line[0:4]=='ATOM':
        d={}
        x=float(line[31:38].strip()) ;y=float(line[39:46].strip()) ;z=float(line[47:55].strip())
        d['res_num']=only_numerics(line[23:29])
        d['res_name']=line[17:20]
        d['atom_name']=line[13:17].strip()
        d['coordinates']=(x,y,z)
        protein_atoms.append(d)
contact_fuzzball.close()
df=pd.DataFrame(protein_atoms) #this is the dataframe with residue info


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

#create a dictionary for contact atoms for each residue = resnum:[(atom1,x1,y1,z1),(...)]
contact_dict=collections.defaultdict(list)
row_count=0
for ptn_atom_coords in list(df['coordinates']):
    for ligand_atom_coords in ligand_atoms:
        d=dist(ptn_atom_coords, ligand_atom_coords)
        if d<4.0:
            resname=df.iloc[row_count][1]
            atom_name=df.iloc[row_count][2]
            contact_atom_info=(resname,atom_name,ptn_atom_coords[0],ptn_atom_coords[1],ptn_atom_coords[2])
            if contact_atom_info not in contact_dict[df.iloc[row_count][0]]: #single ptn atom can be win cutoff of mult lig atoms
                contact_dict[df.iloc[row_count][0]].append(contact_atom_info)
    row_count+=1


######looking at residues in contact fuzzball that now appear
######to actually not have a contact atom...why then are they added  by scf.py? 
for i in range(9298):
    if not contact_dict[str(i)]:
        rowc=df[df['res_num'] == str(i)].index
        coords=list(df['coordinates'][rowc])
        for coord in coords:
            for ligand_atom_coords in ligand_atoms:
                d=dist(coord, ligand_atom_coords)
                print(d)
