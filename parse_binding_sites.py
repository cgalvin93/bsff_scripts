#takes the csv with binding site solns, + the fuzzball pdb as input
#returns an analagous csv with binding sites containing at least one polar residue
#also returns statistics on binding site composition, which at this point
#only includes a bar graph showing the amino acid distribution for the binding sites
#USAGE: python parse_binding_sites.py path/to/binding_site_solns.csv path/to/fuzzball.pdb outfilename

import pandas as pd
import sys
import matplotlib.pyplot as plt

#load binding site solns csv as pandas dataframe
bsdf = pd.read_csv(sys.argv[1])


#load fuzzball and store atom lines
fuzzball_file=open(sys.argv[2])
fuzzball_residues=[] #to store relevant info on residues
for line in fuzzball_file.readlines():
    if line[0:4]=='ATOM':
        resname=line[17:20]
        resnum=line[22:29]
        fuzzball_residues.append((resnum.strip(),resname))
fuzzball_file.close()



#function to filter binding sites on given criteria. at the moment only criteria
#is includes polar residue
polar_list=['SER','THR','CYS','ASN','GLN','TYR','ASP','GLU','HIS'] #including charged res and his
def filter_binding_sites(bs_indices, fuzzball_pdb_list):
    n_polar=0
    for bs_residue in bs_indices:
        for resnum,resname in fuzzball_pdb_list:
            if bs_residue==resnum:
                if resname in polar_list:
                    n_polar+=1
    if n_polar>0:
        return True

#iterate through binding sites and check if they pass filters
rows_list=[] #stores info for rows corresponding to binding sites that pass filter
for x in bsdf['Residue_indicies']:
    row_count=0
    if filter_binding_sites(x,fuzzball_residues)==True:
        d={}
        d['Residue_indicies']=bsdf.iloc[row_count][0]
        d['Obj_score']=bsdf.iloc[row_count][1]
        d['Conformer']=bsdf.iloc[row_count][2]
        rows_list.append(d)
    row_count+=1

#output binding sites containing polar residues to csv
output_df=pd.DataFrame(rows_list) #construct dataframe from accepted rows
output_df.to_csv(sys.argv[3]+'.csv')


#get distribution of amino acids from fuzzball_residues list
all_resnames=[b for a,b in set(fuzzball_residues)]
amino_acids=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY',
             'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
             'THR','TRP','TYR','VAL']

total_res=float(len(all_resnames))
aa_dist_dict={}
for i in amino_acids:
    count=0
    for x in all_resnames:
        if i==x:
            count+=1
    aa_dist_dict[i]=count/total_res

keys=aa_dist_dict.keys()
values=aa_dist_dict.values()
plt.bar(keys, values)
plt.xticks(rotation='vertical')
plt.xlabel('Residue')
plt.ylabel('Frequency')
plt.savefig(sys.argv[3]+'.pdf')
