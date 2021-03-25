# -*- coding: utf-8 -*-
"""
Created on Sun Feb 21 14:05:04 2021

@author: lucia
"""

# Import scaffoldgraph
import scaffoldgraph as sg

# Import networkx
import networkx as nx

# Import plotting tools
import matplotlib.pyplot as plt

import numpy as np

# Import rdkit
from rdkit.Chem import Draw
from rdkit import Chem
import rdkit
import random
import os
from collections import Counter
os.chdir('C:/Users/lucia/OneDrive/Desktop/ScaffoldGraph-master')

#%%
sdf_file = ('limonin.sdf') # Example SDF file (200 PubChem compounds)
#smi_file = ('smiles_plateIA.smi')
supplier = Chem.SDMolSupplier(sdf_file)
network = sg.ScaffoldNetwork.from_sdf(sdf_file, progress=True)

# We can access the number of molecule nodes and scaffold nodes in the graph
n_scaffolds_net = network.num_scaffold_nodes
n_molecules_net = network.num_molecule_nodes

print('\nGenerated scaffold network from {} molecules with {} scaffolds\n'.format(n_molecules_net, n_scaffolds_net))
#%%
# Molecules are stored in the network with their _Name property as a key
# When using SDF format rdkit assigns the _Name property from the TITLE section
# In this case this refers to a PubChem ID

molecules = list(network.get_molecule_nodes())

# Calculate number of scaffolds and singletons
frag_to_save = []
for i in range(len(molecules)):
    mol_id = i
    pubchem_id = molecules[mol_id]
    smiles_iter = network.nodes[pubchem_id]['smiles']
    frags_iter = sg.get_all_murcko_fragments(Chem.MolFromSmiles(smiles_iter))
    num_rings = []
    for mol in frags_iter:
        num_rings.append(rdkit.Chem.rdMolDescriptors.CalcNumRings(mol))
    m = max(num_rings)    
    position_max = [i for i, j in enumerate(num_rings) if j == m]
    frag_to_save.append(frags_iter[position_max[0]]) 
list_smiles = []
for frag_saved in frag_to_save:
    list_smiles.append(Chem.MolToSmiles(frag_saved))

smiles_no_duplicate = set(list_smiles)

count_unique_smiles = Counter(list_smiles)
number_of_unique_smiles= []

for key_count in count_unique_smiles.keys():
    if count_unique_smiles[key_count]==1:
        number_of_unique_smiles.append(count_unique_smiles[key_count])
# Number of scaffold
N = len(smiles_no_duplicate)
Nsing = sum(number_of_unique_smiles)

print("Number of scaffolds is",N)
print("Number of singletons is",Nsing)
#%%
# We can generate a scaffold tree from the SDF file just like before

#tree = sg.ScaffoldTree.from_smiles_file(smi_file, progress=True)
tree = sg.ScaffoldTree.from_sdf(sdf_file, progress=True)

# We can access the number of molecule nodes and scaffold nodes in the graph
n_scaffolds_tree = tree.num_scaffold_nodes
n_molecules_tree = tree.num_molecule_nodes

#print('\nGenerated scaffold tree from {} molecules with {} scaffolds\n'.format(n_molecules_tree, n_scaffolds_tree))

# The output is a forest structure (multiple trees)

#print('Graph is a Forest:', nx.is_forest(tree))
#%%
# We can get the number of scaffolds in each hierarchy easily (The numbers are different to the network)

counts_tree = tree.get_hierarchy_sizes()
lists = sorted(counts_tree.items())
x, y = zip(*lists)

# Plot sizes as bar chart

plt.figure(figsize=(8, 6))
plt.bar(x, y)
plt.xlabel('Hierarchy')
plt.ylabel('Scaffold Count')
plt.title('Number of Scaffolds per Hierarchy (Tree)')
plt.show()
#%%
# frag_all_to_save = []
# num_rings_all = []
# num_frag_per_mol = []
# for i in range(len(molecules)):
#     mol_id = i
#     pubchem_id = molecules[mol_id]
#     smiles_iter = network.nodes[pubchem_id]['smiles']
#     frags_iter = sg.get_all_murcko_fragments(Chem.MolFromSmiles(smiles_iter))
#     num_frag_per_mol.append(len(frags_iter))
#     frag_all_to_save.append(frags_iter)
# all_frags = [item for sublist in frag_all_to_save for item in sublist]
# for mol in all_frags:
#     num_rings_all.append(rdkit.Chem.rdMolDescriptors.CalcNumRings(mol))
# #all_frags = [item for sublist in frag_all_to_save for item in sublist]

# list_all_smiles = []
# for frag_iteration in all_frags:
#     list_all_smiles.append(Chem.MolToSmiles(frag_iteration))

# count_all_smiles = Counter(list_all_smiles)

# perc_mol = []
# num_mol  = []
# for key_smile in count_all_smiles.keys():
#     perc_mol.append(count_all_smiles[key_smile]/len(molecules))
#     num_mol.append(count_all_smiles[key_smile])

# perc_mol_sorted = 1-np.sort(np.array(perc_mol))[::-1]

# list_smiles = []
# for frag_saved in frag_to_save:
#     list_smiles.append(Chem.MolToSmiles(frag_saved))

# smiles_no_duplicate = set(list_smiles)

# count_unique_smiles = Counter(list_smiles)
# number_of_unique_smiles= []

# for key_count in count_unique_smiles.keys():
#     if count_unique_smiles[key_count]==1:
#         number_of_unique_smiles.append(count_unique_smiles[key_count])
        
# def ecdf(data):
#     """ Compute ECDF """
#     x = np.sort(data)
#     n = x.size
#     y = np.arange(1, n+1) / n
#     return(x,y)
