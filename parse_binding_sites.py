#takes the csv with binding site solns, + the fuzzball pdb as input
#returns an analagous csv with binding sites containing at least one polar residue
#also returns statistics on binding site composition, which at this point
#only includes a bar graph showing the amino acid distribution for the binding sites
#USAGE: python ~/Desktop/prj/bsff/bsff_scripts/parse_binding_sites.py ~/Desktop/Files_IB2/IB2_0052-iter_0-fuzz_0-3_residue_solutions.csv ~/Desktop/Files_IB2/IB2_0052-iter_0-fuzz_0.pdb testIB2

import pandas as pd
import sys
import matplotlib.pyplot as plt

#load binding site solns csv as pandas dataframe
bsdf = pd.read_csv(sys.argv[1])

##function to return only digits in str
def only_numerics(seq):
    seq_type= type(seq)
    return seq_type().join(filter(seq_type.isdigit, seq))
#load fuzzball and store atom lines
fuzzball_file=open(sys.argv[2])
fuzzball_residues=[] #to store relevant info on residues
for line in fuzzball_file.readlines():
    if line[0:4]=='ATOM':
        resname=line[17:20]
        resnum=line[22:29]
        fuzzball_residues.append((resnum.strip(),resname))
fuzzball_file.close()
fuzzball_residues=list(set(fuzzball_residues))



#function to filter binding sites on given criteria. at the moment only criteria
#is includes polar residue
polar_list=['SER','THR','CYS','ASN','GLN','TYR','ASP','GLU','HIS'] #including charged res and his
def filter_binding_sites(bs_indices, fuzzball_pdb_list):
    n_polar=0
    for bs_residue in bs_indices: #bs indices should be the list of n res numbers
        for resnum,resname in fuzzball_pdb_list:
            if int(only_numerics(bs_residue))==int(only_numerics(resnum)):
                if resname in polar_list:
                    n_polar+=1
    if n_polar>0:
        return True



#iterate through binding sites and check if they pass filters
rows_list=[] #stores info for rows corresponding to binding sites that pass filter
row_count=0
for x in list(bsdf['Residue_indicies']):
    bs_res_indices=[only_numerics(i) for i in x.split(',')]
    print('checking residues '+str(bs_res_indices)+' for polar contacts')
    if filter_binding_sites(bs_res_indices,fuzzball_residues)==True:
        d={}
        d['Residue_indicies']=bsdf.iloc[row_count][0]
        d['Obj_score']=bsdf.iloc[row_count][1]
        d['Conformer']=bsdf.iloc[row_count][2]
        rows_list.append(d)
    row_count+=1


#output binding sites containing polar residues to csv
output_df=pd.DataFrame(rows_list) #construct dataframe from accepted rows
print(output_df)
output_df.to_csv(sys.argv[3]+'.csv')
n_filtered=output_df.shape[0]
print(str(n_filtered)+' binding sites found containing polar contact')

#get distribution of amino acids from fuzzball_residues list
all_resnames=[b for a,b in set(fuzzball_residues)]
amino_acids=['ALA','ARG','ASN','ASP','CYS','GLU','GLN','GLY',
             'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER',
             'THR','TRP','TYR','VAL']

print('there are '+str(len(all_resnames))+ ' total residues in the fuzzball')
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
