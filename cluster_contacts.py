
# ~Desktop/prj/BSFF/bsff_scripts/cluster_contacts.py
# ~Desktop/prj/BSFF/output/contactfuzzballs/contactsubset.pdb
# USAGE: time python cluster_contacts.py path/to/contact_fuzzball.pdb

import sys
import pandas as pd
import math
import collections
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

#function to return only digits from string
def only_numerics(seq):
    seq_type= type(seq)
    return seq_type().join(filter(seq_type.isdigit, seq))

'''
# contact_fuzzball=open(sys.argv[1],'r')
def pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:4]=='ATOM':
            d={}
            d['recordname']=line[0:6]
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=line[17:20]
            d['chain']=line[20:22]
            d['resnum']=line[22:26]
            d['achar']=line[26:31]
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
    return df
df=pdb_to_df('Fragment2_contacts.pdb')

def df_to_pdb(dataframe,ofile):
    of=open(ofile,'w')
    for i in range(df.shape[0]):
        line=''.join(df.iloc[i])
        of.write(line)
    of.close()


'''
#store ligand atom coordinates in list,
#store protein residue info in pandas df resnum, resname, atom name, coords(xyz)
ligand_atoms=[]
protein_atoms=[]
contact_fuzzball=open('Fragment2_contacts.pdb','r') #################testing########

contact_fuzzball.close()
df=pd.DataFrame(protein_atoms) #this is the dataframe with residue info

for line in contact_fuzzball.readlines():
    if line[0:6]=='HETATM':
        x=float(line[31:38].strip()) ;y=float(line[39:46].strip()) ;z=float(line[47:54].strip())
        ligand_atoms.append((x,y,z))
    elif line[0:4]=='ATOM':
        d={}
        x=float(line[31:38].strip()) ;y=float(line[39:46].strip()) ;z=float(line[47:54].strip())
        d['res_num']=only_numerics(line[23:29])
        d['res_name']=line[17:20]
        d['atom_name']=line[12:17].strip()
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

'''
######looking at residues in contact fuzzball that now appear
######to actually not have a contact atom...why then are they added  by scf.py?
######looking at the distances printed from this block
######it seems they are mostly much much bigger than 4... need to revisit scf.py
# In [14]: len(contact_dict)
# Out[14]: 4237

for i in range(9298):
    if not contact_dict[str(i)]:
        rowc=df[df['res_num'] == str(i)].index
        coords=list(df['coordinates'][rowc])
        for coord in coords:
            for ligand_atom_coords in ligand_atoms:
                d=dist(coord, ligand_atom_coords)
                print(d)
'''

#lets break down to bb only vs sc only vs both
bb=[]
sc=[]
scbb=[]
for key in contact_dict.keys():
    nbb=0
    nsc=0
    for contact_atom_info in contact_dict[key]:
        atom_name=contact_atom_info[1]
        if atom_name=='C' or atom_name=='CA'or atom_name=='O' or atom_name=='N':
            nbb+=1
        elif atom_name!='C' and atom_name!='CA'and atom_name!='O' and atom_name!='N':
            nsc+=1
    if nbb>0 and nsc==0:
        bb.append(key)
    if nsc>0 and nbb==0:
        sc.append(key)
    if nbb>0 and nsc>0:
        scbb.append(key)
        sc.append(key) #counting scbb amongst sc contacts

print('There are '+str(len(contact_dict))+' total contact residues')
print(str(len(bb))+' are backbone contacts')
print(str(len(sc))+' are sidechain contacts')


'''
# gonna look at residue frequency dists for bb, sc...
# which one does both resemble more?

amino_acids=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY',
             'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
             'THR','TRP','TYR','VAL']
def return_res_frequencies(keys, data_dict):
    frequency_dict={}
    contact_resnames=[]
    total_res=float(len(keys))
    for key in keys:
        contact_resname=data_dict[key][0][0]
        contact_resnames.append(contact_resname)
    for resname in amino_acids:
        count=0
        for contact_resname in contact_resnames:
            if contact_resname==resname:
                count+=1
        frequency_dict[resname]=count/total_res
    return frequency_dict.keys(),frequency_dict.values()

from matplotlib.backends.backend_pdf import PdfPages
pdf = PdfPages('resfreqs.pdf')

bbx,bby=return_res_frequencies(bb, contact_dict)
plt.bar(bbx, bby)
plt.xticks(rotation='vertical')
plt.xlabel('Residue')
plt.ylabel('Frequency')
plt.title('bb')
pdf.savefig()
plt.clf()

scx,scy=return_res_frequencies(sc, contact_dict)
plt.bar(scx, scy)
plt.xticks(rotation='vertical')
plt.xlabel('Residue')
plt.ylabel('Frequency')
plt.title('sc')
pdf.savefig()
plt.clf()

scbbx,scbby=return_res_frequencies(scbb, contact_dict)
plt.bar(scbbx, scbby)
plt.xticks(rotation='vertical')
plt.xlabel('Residue')
plt.ylabel('Frequency')
plt.title('sc+bb')
pdf.savefig()
plt.clf()

pdf.close()

# GOING TO TREAT BBSC CONTACTS AS SC CONTACTS
# I REASON THAT IF THERE IS A SC ATOM THAT CLOSE TO THE LIGAND THEN
# ITS IDENTITY AND ORIENTATION MUST CERTAINLY BE RELEVANT
# AND THERE ARE SUFFICIENT BB EXCLUSIVE CONTACTS THAT I CAN USE THIS
# EXCLUSIVE GROUP TO EXTRACT ENOUGH INFO ABOUT BB PREFERRED GEOMETRIES AND STUFF
# IE I DONT 'NEED' THE INFO ABOUT BB CONTACTS IN THE BBSC CONTACTS

# next group sc by chemistry
# 0 | 1	Hydrophobic, aliphatic (AILV) contact (1) or not (0)
# 0 | 1	Hydrophobic, aromatic (FWY) contact (1) or not (0)
# 0 | 1	Polar (NCQMST) contact (1) or not (0)
# 0 | 1	Charged, Acidic (DE) contact (1) or not (0)
# 0 | 1	Charged, Basic (HKR) contact (1) or not (0)
# 0 | 1	Glycine contact (1) or not (0)
# 0 | 1	Proline contact (1) or not (0)
#
# and bb by contact atom
# 0 | 1	Backbone carbonyl contact (1) or not (0)
# 0 | 1	Backbone amino contact (1) or not (0)
# 0 | 1	Backbone C/CA contact (1) or not (0)



'''

amino_acids=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY',
             'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
             'THR','TRP','TYR','VAL']

#lists to hold keys for residues of given contact chemistry
aliphatic=[]
aromatic=[]
polar=[]
charged_acidic=[]
charged_basic=[]
glycines=[]
prolines=[]
methionines=[]
'''
not sure what to do with methionine, its more nonpolar right?
but feels weird to call it aliphatic hphobe
also i am inclined to double count tyrosine as both an aromatic and a polar
interaxn, for now keeping it aromatic exclusive
also what to do with his? for now keeping it charge basic as did james
'''
#these are for the res frequency stats for all sc contacts
sc_frequency_dict={}
sc_contact_resnames=[]
for key in sc:
        contact_resname=contact_dict[key][0][0]
        sc_contact_resnames.append((key,contact_resname))
        if contact_resname=='ALA' or contact_resname=='ILE' or contact_resname=='LEU' or contact_resname=='VAL':
           aliphatic.append(key)
        elif contact_resname=='PHE' or contact_resname=='TYR' or contact_resname=='TRP':
            aromatic.append(key)
        elif contact_resname=='SER' or contact_resname=='THR' or contact_resname=='ASN' or contact_resname=='GLN' or contact_resname=='CYS':
             polar.append(key)
        elif contact_resname=='GLU' or contact_resname=='ASP':
            charged_acidic.append(key)
        elif contact_resname=='ARG' or contact_resname=='LYS' or contact_resname=='HIS':
            charged_basic.append(key)
        elif contact_resname=='GLY':
            glycines.append(key)
        elif contact_resname=='PRO':
            prolines.append(key)
        elif contact_resname=='MET':
            methionines.append(key)
#bar plot contact chemistry frequencies
total_res=float(len(contact_dict))
contact_chem_categories=['bb','aliph','arom','polar','q-',
                         'q+','gly','pro','met']
contact_chem_freqs=[]
nbb=len(bb)/total_res;contact_chem_freqs.append(nbb)
naliphatic=len(aliphatic)/total_res;contact_chem_freqs.append(naliphatic)
naromatic=len(aromatic)/total_res;contact_chem_freqs.append(naromatic)
npolar=len(polar)/total_res;contact_chem_freqs.append(npolar)
ncharged_acidic=len(charged_acidic)/total_res;contact_chem_freqs.append(ncharged_acidic)
ncharged_basic=len(charged_basic)/total_res;contact_chem_freqs.append(ncharged_basic)
nglycines=len(glycines)/total_res;contact_chem_freqs.append(nglycines)
nprolines=len(prolines)/total_res;contact_chem_freqs.append(nprolines)
nmethionines=len(methionines)/total_res;contact_chem_freqs.append(nmethionines)
#
pdf = PdfPages('contact_statistics.pdf')
#
plt.bar(contact_chem_categories, contact_chem_freqs)
plt.xticks(rotation='vertical')
plt.xlabel('Contact Chemistry')
plt.ylabel('Frequency')
plt.title('Contact Chemistry Frequencies')
pdf.savefig()
plt.clf()


#amino acid distribution for sidechain contacts
sc_total_res=float(len(sc))
res_key_dict=collections.defaultdict(list)
for resname in amino_acids:
    count=0
    for key,contact_resname in sc_contact_resnames:
        if contact_resname==resname:
            count+=1
            res_key_dict[resname].append(key)
    sc_frequency_dict[resname]=count/sc_total_res
scx,scy=sc_frequency_dict.keys(),sc_frequency_dict.values()
plt.bar(scx, scy)
plt.xticks(rotation='vertical')
plt.xlabel('Residue')
plt.ylabel('Frequency')
plt.title('Amino Acid Frequencies for Sidechain Contacts')
pdf.savefig()
plt.clf()

pdf.close()


'''
okay next I want to cluster based on geometry for sc res
there is some ambiguity as to whether i should do this on the level
of individual residues or on the level of contact chemistry
that is, res within same contact chem category are probably playing the same role
ie hb donor/acceptor w a specific lig atom
OH BUT REMEMBER, FIRST I JUST WANNA LOOK AT DISTRIBUTION OF GEOMETRIES
STILL, FIRST QUESTION IS TO DEFINE THESE DISTRIBUTIONS
aliphatic=[]
aromatic=[]
polar=[]
charged_acidic=[]
charged_basic=[]
glycines=[]
prolines=[]
methionines=[]
A FRAG LIKE PHENYL RING ACTUALLY HAS NO HB DONOR OR ACCEPTOR
SO CANT DEFINE HBOND PROPS



'''
#well one way or another, I need access to all atom coords for individual residues
#in a convenient way
#since not all have hydrogens, I'm gonna exlude hydrogens
residue_coords_dict=collections.defaultdict(list)
for i in range(df.shape[0]):
        resnum=df.iloc[i][0]
        atom_name=df.iloc[i][2]
        ptn_atom_coords=df.iloc[i][3]
        if atom_name[0]!='H' and atom_name!='OXT' and atom_name!='C' and atom_name!='CA'and atom_name!='O' and atom_name!='N':
            contact_atom_info=(float(ptn_atom_coords[0]),float(ptn_atom_coords[1]),float(ptn_atom_coords[2]))
            if contact_atom_info not in residue_coords_dict[resnum]:
                residue_coords_dict[resnum].append(contact_atom_info)
'''
#now I also need all the keys corresponding to individual res
#going back and putting this into res freq block

checking that res keys sum to same as sc keys:
In [5]: ns=[]
   ...: for key in res_key_dict.keys():
   ...:     x=len(res_key_dict[key])
   ...:     ns.append(x)
   ...: print(sum(ns))
3298

In [6]: len(sc)
Out[6]: 3298

LOOKS GOOD

now I want to do things with coordinates of individual residues
for aliphatic and aromatics at least,
i am gonna start by just iding clusters of res with low rmsd
'''
#return the rmsd between two sets of cartesian coordinates
def calc_rmsd(v1,v2):
    diff = np.array(v1) - np.array(v2)
    N = len(v1)
    return np.sqrt((diff * diff).sum() / N)


#I think the following code will create clusters for aliphatic and
#aromatic contacts based off of the rmsd of the sidechains
cluster_dict=collections.defaultdict(list)
for key in res_key_dict.keys():
    if key=='ALA' or key=='ILE' or key=='LEU' or key=='VAL' or key=='PHE' or key=='TRP' or key=='TYR':
        already_clustered=[]
        for index,strc1 in enumerate(res_key_dict[key]):
            if strc1 not in already_clustered:
                already_clustered.append(strc1)
                cluster_dict[strc1]
                current_cluster=[]
                for strc2 in res_key_dict[key]:
                    if strc2!=strc1 and strc2 not in already_clustered:
                        x=len(residue_coords_dict[strc1])
                        y=len(residue_coords_dict[strc2])
                        if x!=y:
                            print('different number of atoms between residues '+str(strc1)+' and '+str(strc2))
                            print('cannot calculate rmsd :(')
                        else:
                            rmsd=calc_rmsd(residue_coords_dict[strc1],residue_coords_dict[strc2])
                            if rmsd<2.5:
                                already_clustered.append(strc2)
                                cluster_dict[strc1].append(strc2)
                            else:
                                continue
    else:
        pass
print(cluster_dict)
print(str(len(cluster_dict)))
'''
cutoff 1.5
ns=[]
for key in res_key_dict.keys():
    if key=='ALA' or key=='ILE' or key=='LEU' or key=='VAL' or key=='PHE' or key=='TRP' or key=='TYR':
        x=len(res_key_dict[key])
        ns.append(x)
print(str(sum(ns)))
1848
ns=[]
ns.append(len(cluster_dict))
for key in cluster_dict.keys():
    x=len(cluster_dict[key])
    ns.append(x)
print(str(sum(ns)))
1848

cutoff 2.5
In [70]: np.median(ns)
Out[70]: 1.0
In [71]: np.mean(ns)
Out[71]: 3.85
'''


'''
need something to write the cluster fuzzballs


'''
