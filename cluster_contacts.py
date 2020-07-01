#this script takes a contact fuzzball pdb file as input
#such as the one generated from contact_fuzzball.py script
#returns a pdf showing some stats on the contact chemistry, amino acid
#distribution, and cluster population distribution
#also returns pdbs of clusters with at lest X members, where you can enter
#X around line 340 (if b>X)
#clustering is done based off of sidechain rmsd being <1.5 angs between 2 res
#
# USAGE: time python cluster_contacts.py path/to/contact_fuzzball.pdb frag_3letter_ID
# for me: time ipython ~/Desktop/prj/bsff/bsff_scripts/cluster_contacts.py fragment1_contact_fuzzball.pdb cx8

#imports
import sys
import os
import pandas as pd
import math
import collections
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

###
###defining a bunch of functions
###

#function to return only digits from string
def only_numerics(seq):
    seq_type= type(seq)
    return seq_type().join(filter(seq_type.isdigit, seq))

#convert pdb file to pandas dataframe
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
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']=line[78:80]
            l.append(d)
    df=pd.DataFrame(l)
    return df


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

#convert dataframe to pdb file
def df_to_pdb(dataframe,ofile):
    of=open(ofile,'w')
    for i in range(dataframe.shape[0]):
        line=''.join(dataframe.iloc[i])
        of.write(line)
    of.close()

#convert ligand pdb to dataframe
def lig_pdb_to_df(pdbfile):
    l=[]
    f=open(pdbfile,'r')
    for line in f.readlines():
        if line[0:6]=='HETATM':
            d={}
            d['recordname']='ATOM'+'  '
            d['atomnumber']=line[6:11]
            d['atomname']=line[11:16]
            d['altloc']=line[16:17]
            d['resname']=sys.argv[2]
            d['chain']=' X'
            d['resnum']=line[22:29]
            d['achar']=line[29:31]
            d['x']=line[31:39]
            d['y']=line[39:47]
            d['z']=line[47:54]
            d['occupancy']=line[54:60]
            d['temp_factor']=line[60:66]
            d['seg']=line[66:76]
            d['element']=line[76:78]
            d['q']=line[78:80]
            l.append(d)
        else:
            break
    df=pd.DataFrame(l)
    return df

#return the rmsd between two sets of cartesian coordinates
def calc_rmsd(v1,v2):
    diff = np.array(v1) - np.array(v2)
    N = len(v1)
    return np.sqrt((diff * diff).sum() / N)


#############
#now doing stuff
##############


#load contact fuzzball as pandas dataframe
df=pdb_to_df(sys.argv[1])

#load ligand coords and df
ligand_atoms=[]
contact_fuzzball=open(sys.argv[1],'r')
ligdf=lig_pdb_to_df(sys.argv[1])
for line in contact_fuzzball.readlines():
    if line[0:6]=='HETATM':
        x=float(line[31:39].strip()) ;y=float(line[39:47].strip()) ;z=float(line[47:54].strip())
        ligand_atoms.append((x,y,z))
contact_fuzzball.close()

#create a dictionary where key is resnumber
#value is contact atom information (w/in 4.1 angs of lig fragment)
contact_dict=collections.defaultdict(list)
for i in range(df.shape[0]):
    ptn_atom_coords=( float(df.iloc[i][8]), float(df.iloc[i][9]), float(df.iloc[i][10]) )
    for ligand_atom_coords in ligand_atoms:
        d=dist(ptn_atom_coords, ligand_atom_coords)
        if d<4.1:
            resname=df.iloc[i][4]
            atom_name=str(df.iloc[i][2].split()[0])
            contact_atom_info=(resname,atom_name,ptn_atom_coords[0],ptn_atom_coords[1],ptn_atom_coords[2])
            if contact_atom_info not in contact_dict[df.iloc[i][0]]:
                contact_dict[df.iloc[i][6]].append(contact_atom_info)

#lets break down contacts to bb only vs sc only vs both
bb=[]
sc=[]
scbb=[]
for key in contact_dict.keys():
    nbb=0
    nsc=0
    for contact_atom_info in contact_dict[key]:
        atom_name=str(contact_atom_info[1])
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

#sc_contact_resnames holds res identities of sc contacts
#sc_frequency_dict will eventually have key=resname,value=n contacts for that res
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
pdfname=sys.argv[1][:-4]+'_statistics.pdf'
pdf = PdfPages(pdfname)
#
sortlist=[]
for i,e in enumerate(contact_chem_freqs):
    sortlist.append((contact_chem_categories[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
contact_chem_categories=[];contact_chem_freqs=[]
for a,b in sortlist:
    contact_chem_categories.append(a);contact_chem_freqs.append(b)
#
plt.bar(contact_chem_categories, contact_chem_freqs)
plt.xticks(rotation='vertical')
plt.xlabel('Contact Chemistry')
plt.ylabel('Frequency')
plt.title('Contact Chemistry Frequency')
pdf.savefig()
plt.clf()


#plot amino acid distribution for sidechain contacts
sc_total_res=float(len(sc))
res_key_dict=collections.defaultdict(list)
for resname in amino_acids:
    count=0
    for key,contact_resname in sc_contact_resnames:
        if contact_resname==resname:
            count+=1
            res_key_dict[resname].append(key)
    sc_frequency_dict[resname]=count/sc_total_res
scx,scy=list(sc_frequency_dict.keys()),sc_frequency_dict.values()
sortlist=[]
for i,e in enumerate(scy):
    sortlist.append((scx[i],e))
sortlist=sorted(sortlist, reverse=True, key=lambda nmem: nmem[1])
scx=[];scy=[]
for a,b in sortlist:
    scx.append(a);scy.append(b)
plt.bar(scx, scy)
plt.xticks(rotation='vertical')
plt.xlabel('Residue')
plt.ylabel('Frequency')
plt.title('Amino Acid Frequencies for Sidechain Contacts')
pdf.savefig()
plt.clf()




#I need access to all atom coords for individual residues
#in a convenient way to find rmsd for clustering
#since not all have hydrogens, I'm gonna exlude hydrogens (also oxt atoms)
residue_coords_dict=collections.defaultdict(list)
for i in range(df.shape[0]):
        resnum=df.iloc[i][6]
        atom_name=str(df.iloc[i][2].split()[0])
        ptn_atom_coords=( float(df.iloc[i][8]), float(df.iloc[i][9]), float(df.iloc[i][10]) )
        if atom_name[0]!='H' and atom_name!='OXT' and atom_name!='C' and atom_name!='CA'and atom_name!='O' and atom_name!='N':
            contact_atom_info=(atom_name,float(ptn_atom_coords[0]),float(ptn_atom_coords[1]),float(ptn_atom_coords[2]))
            if contact_atom_info not in residue_coords_dict[resnum]:
                residue_coords_dict[resnum].append(contact_atom_info)


#creater clusters based off of sidechain rmsd
cluster_dict=collections.defaultdict(list)
for key in res_key_dict.keys():
    if key!='GLY' and key!='PRO' and key!='ALA' and key!='CYS': #excluding these res from clustering, as did james
        already_clustered=[]
        for index,strc1 in enumerate(res_key_dict[key]):
            if strc1 not in already_clustered:
                already_clustered.append(strc1)
                current_cluster=[]
                for strc2 in res_key_dict[key]:
                    if strc2!=strc1 and strc2 not in already_clustered:
                        vs1=[];vs2=[]
                        for a,b,c,d in residue_coords_dict[strc1]:
                            for e,f,g,h in residue_coords_dict[strc2]:
                                if a==e:
                                    vs1.append((b,c,d))
                                    vs2.append((f,g,h))
                                else:
                                    pass
                        x=len(vs1)
                        y=len(vs2)
                        if x!=y:
                            print('different number of atoms between residues '+str(strc1)+' and '+str(strc2))
                            print('cannot calculate rmsd :(')
                        else:
                            rmsd=calc_rmsd(vs1,vs2)
                            if rmsd<1.5:
                                already_clustered.append(strc2)
                                cluster_dict[strc1].append(strc2)
                            else:
                                continue
    else:
        pass


#okay so next thing is to create fuzzball pdbs of clusters
##############################
#compile and organize cluster populations
cluster_populations=[]
for key in cluster_dict.keys():
    n_members=len(cluster_dict[key])+1
    cluster_populations.append((key,int(n_members)))
cluster_populations=sorted(cluster_populations, reverse=True, key=lambda nmem: nmem[1])
#plot cluster populations
ncm=[y for x,y in cluster_populations]
clabs=[i+1 for i in range(len(cluster_populations))]
plt.bar(clabs, ncm)
plt.xlabel('Cluster ID')
plt.ylabel('Number of Members')
plt.title('Cluster Population Distribution')
pdf.savefig()
plt.clf()
pdf.close()

#make a nice home for stats and cluster pdbs
newdir=pdfname[:-4]+'cluster_results'
os.mkdir(newdir)
os.rename(pdfname,newdir+'/'+pdfname)

#export cluster fuzzballs
rns=[]
for i in df['resnum']:
    rns.append(int(i))
for a,b in cluster_populations:
    cluster_df=pd.DataFrame()
    if b>=0:
        resnumbers=[int(i) for i in cluster_dict[a]]
        resnumbers.append(int(a))
        for resnum in resnumbers:
            rows=[i for i,e in enumerate(rns) if e == resnum]
            res_df=df.iloc[min(rows):max(rows)+1]
            cluster_df=cluster_df.append(res_df,ignore_index=True)
        cluster_df=cluster_df.append(ligdf,ignore_index=True)
        resnamee=str(res_df.iloc[0][4])
        ofilename=str(a.strip())+'_'+resnamee+'_cluster_'+str(b)+'.pdb'
        df_to_pdb(cluster_df,ofilename)
        os.rename(ofilename,newdir+'/'+ofilename)
    else:
        break
