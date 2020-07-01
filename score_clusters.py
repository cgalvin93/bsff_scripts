#USAGE: time ipython score_clusters.py

import collections
import os
import pandas as pd
from pyrosetta import *
init()
from pyrosetta.toolbox import atom_pair_energy

sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec
sf.set_weight(fa_atr, 0.423)
sf.set_weight(fa_rep, 0.100)
sf.set_weight(fa_sol, 0.372)
sf.set_weight(hbond_sc, 0.245)
sf.set_weight(fa_elec, 0.026)

#rosetta ligand score function
#sf = create_score_function('ligand')

# rawterms=list(sf.get_nonzero_weighted_scoretypes())
# scoreterms=[str(i).split('.')[1] for i in rawterms]
###########FOR SOME REASON THIS ONLY WORKS WITH THESE FOUR SCORETERMS ....
scoreterms=['fa_atr','fa_rep','fa_sol','fa_elec']
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']
cluster_data=[]
for file in pdbfiles:
    cluster_dict={}
    l = ['arl.params']
    p = Pose()
    generate_nonstandard_residue_set(p,l)
    pose_from_file(p, file)
    ligand_resnum=p.total_residue()
    term_scores=collections.defaultdict(list)
    cluster_dict['cluster_population']=p.total_residue()-1
    for term in scoreterms:
        try:
            all_res_energies=[]
            for i in range(p.total_residue()-1):
                resnum=i+1
                pwise_energies=list(atom_pair_energy._atom_pair_energy_table(sf, str(term), p.residue(ligand_resnum), p.residue(1), threshold=0.0))
                atom_energies=[]
                for atom_list in pwise_energies:
                    for entries in atom_list[1:]:
                        atom_energies.append(entries[1])
                res_energy=sum(atom_energies)
                all_res_energies.append(res_energy)
            term_scores[term].append(all_res_energies)
        except:
            print('Could not evaluate score term: '+term)
    avg_E=0 #this is the average residue-ligand interaction energy for the cluster
    for term in term_scores.keys():
        l=[]
        for res_scores in list(term_scores[term]):
            for res_score in res_scores:
                l.append(res_score)
        terme=sum(l)/float(len(l))
        cluster_dict[term]=terme
        avg_E+=terme
    cluster_dict['total_score']=avg_E
    cluster_data.append(cluster_dict)
df=pd.DataFrame(cluster_data)

####################
#alright now I wanna figure out how to look at hydrogen bond energies if they exist 
hbs=pyrosetta.rosetta.core.scoring.hbonds.HBondSet(p,True)
list(hbs.hbonds())
hbs.nhbonds()
hbs.residue_hbonds(1)
