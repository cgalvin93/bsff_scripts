#USAGE: time ipython score_clusters.py

import collections
import os
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


scoreterms=['fa_atr','fa_rep','fa_sol','fa_elec'] #wont let me do hbond_sc for some reason!!!!!!!!!!!!!
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']
cluster_data=[]
for file in pdbfiles:
    l = ['arl.params']
    p = Pose()
    generate_nonstandard_residue_set(p,l)
    pose_from_file(p, file)
    ligand_resnum=p.total_residue()
    term_scores=collections.defaultdict(list)
    for term in scoreterms:
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
    avg_E=0 #this is the average residue-ligand interaction energy for the cluster
    for term in term_scores.keys():
        l=[]
        for res_scores in list(term_scores[term]):
            for res_score in res_scores:
                l.append(res_score)
        terme=sum(l)/float(len(l))
        print(terme)
        avg_E+=terme
    cluster_data.append((p.total_residue()-1,avg_E))
