#USAGE: time ipython ~/Desktop/prj/bsff/bsff_scripts/score_clusters.py XXX.params
#USAGE: time ipython ~/Desktop/prj/bsff/bsff_scripts/score_clusters.py cx8.params

import sys
import collections
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyrosetta import *
init()
from pyrosetta.toolbox import atom_pair_energy

#using the terms and weights that james used to filter individual contacts
sf = ScoreFunction()
from pyrosetta.rosetta.core.scoring import fa_atr, fa_rep, fa_sol,hbond_sc, fa_elec, hbond_bb_sc
sf.set_weight(fa_atr, 1)
sf.set_weight(fa_rep, .55)
sf.set_weight(fa_sol, 1)
sf.set_weight(hbond_sc, 1)
sf.set_weight(fa_elec, 1)
sf.set_weight(hbond_bb_sc,1)

#rosetta ligand score function
#sf = create_score_function('ligand')


pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']
cluster_data=[]
for file in pdbfiles:
    frag = [sys.argv[1]]
    p = Pose()
    generate_nonstandard_residue_set(p,frag)
    pose_from_file(p, file)
    ligand_resnum=p.total_residue()
    cluster_dict={}
    cluster_dict['cluster_population']=p.total_residue()-1
    term_scores=collections.defaultdict(list) #all scores for res in cluster, for each term
    for i in range(p.total_residue()-1): #HERE ITERATING THROUGH RES IN CLUSTER
        resnum=i+1
        contact_pose=Pose()
        ligand_pose=p.residue(p.total_residue()).clone()
        res_pose=p.residue(resnum).clone()
        contact_pose.append_residue_by_jump(res_pose, 1)
        contact_pose.append_residue_by_jump(ligand_pose, 1)
        # rosetta.core.pack.optimizeH(contact_pose, sf)
        sf(contact_pose)
        w=contact_pose.energies().energy_graph().find_energy_edge(1,2)
        w.fill_energy_map()
        hbsc=w[rosetta.core.scoring.hbond_sc];term_scores['hbond_sc'].append(hbsc)
        hbbbsc=w[rosetta.core.scoring.hbond_bb_sc];term_scores['hbond_bb_sc'].append(hbbbsc)
        faatr=w[rosetta.core.scoring.fa_atr];term_scores['fa_atr'].append(faatr)
        farep=w[rosetta.core.scoring.fa_rep];term_scores['fa_rep'].append(farep)
        fasol=w[rosetta.core.scoring.fa_sol];term_scores['fa_sol'].append(fasol)
        faelec=w[rosetta.core.scoring.fa_elec];term_scores['fa_elec'].append(faelec)
    avg_E=0 #this is the average residue-ligand interaction energy for the cluster
    for term in term_scores.keys():
        term_mean=np.mean(list(term_scores[term]))
        term_var=np.var(list(term_scores[term]))
        cluster_dict[term]=term_mean
        avg_E+=term_mean
    cluster_dict['total_score']=avg_E
    cluster_data.append(cluster_dict)

#export correlation matrix
df=pd.DataFrame(cluster_data)
cdf=df.corr()
cdf.to_csv('cluster_correlations.csv')

#export plot of cluster population vs total score
xs=[];ys=[]
for i in range(df.shape[0]):
    xs.append(df.iloc[i][7])
    ys.append(df.iloc[i][0])
plt.scatter(xs,ys)
plt.xlabel('Average Score of Cluster Contact')
plt.ylabel('Cluster Size')
plt.title('Cluster Score vs. Frequency')
plt.savefig('score_v_pop.pdf')
plt.close()

'''

# rawterms=list(sf.get_nonzero_weighted_scoretypes())
# scoreterms=[str(i).split('.')[1] for i in rawterms]
###########FOR SOME REASON THIS ONLY WORKS WITH THESE FOUR SCORETERMS ....
scoreterms=['fa_atr','fa_rep','fa_sol','fa_elec']
pdbfiles=[file for file in os.listdir(os.getcwd()) if file[-3:]=='pdb']
cluster_data=[]
for file in pdbfiles:
    cluster_dict={}
    l = [sys.argv[1]]
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
# l=['cx8.params']
# p = Pose()
# generate_nonstandard_residue_set(p,l)
# pose_from_file(p, '102_ARG_cluster_62.pdb')
# hbs=pyrosetta.rosetta.core.scoring.hbonds.HBondSet(p,True)
# list(hbs.hbonds())
# hbs.nhbonds()
# hbs.residue_hbonds(63)
#######################
# Okay well none of this is working and I think if I want to look at the hbonds
# Im gonna have to do it myself by extracting coordinates of relevant atoms
# the good thing is pyrosetta automatically adds hydrogens and I think they have standard
# names so it shouldnt be hard to get the coordinates and operate on them
# just gonna have to go through a lot to define donor and acceptor atoms for each case
# and like working out which is which depending on which is closest for like,
# for example carboxylate has two oxygens, arg several donor H etc.....

'''
