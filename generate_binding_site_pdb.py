#takes binding site csv, fuzzball pdb, and row number of desired output pdb (from first column of csv)
#returns a pdb of that binding site
#USAGE: python generate_binding_site_pdb.py path/to/binding_site_solns.csv path/to/fuzzball.pdb row_num

#
import sys
import pandas as pd

#
bsdf = pd.read_csv(sys.argv[1])
desired_pdb=int(sys.argv[3])

#function to return only digits in str
def only_numerics(seq):
    seq_type= type(seq)
    return seq_type().join(filter(seq_type.isdigit, seq))

#acquire and clean indices of residues in desired binding site
indices=bsdf['Residue_indicies'][desired_pdb].split(',')
for i in indices:
    indices[indices.index(i)]=only_numerics(i)
for i in indices:
    indices[indices.index(i)]=int(i)


output_pdb=[] #store desired lines for pdb file
#get the ligand lines from fuzzball pdb
fuzzball_file=open(sys.argv[2])
for line in fuzzball_file.readlines():
    if line[0:6]=='HETATM':
        output_pdb.append(line)
fuzzball_file.close()

#find the lines for desired bs res in fuzzball pdb
fuzzball_file=open(sys.argv[2])
for line in fuzzball_file.readlines():
    if line[0:4]=='ATOM':
        resnum=only_numerics(line[22:29])
        for k in indices:
            if int(resnum)==k:
                output_pdb.append(line)
fuzzball_file.close()

#output pdb
pdb_filename=sys.argv[1]+'_row_'+sys.argv[3]+'.pdb'
f=open(pdb_filename,'w')
for i in output_pdb:
    f.write(i)
f.close()
