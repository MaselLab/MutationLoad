#!/usr/bin/python

#import necessary packages
import tskit
import math
import sys
import io
import numpy as np
import pandas as pd
import os
import argparse

# Define command-line arguments
parser = argparse.ArgumentParser(description="Output of coalescent Ne of simulation")
parser.add_argument("popsize", type=int, help="Input census population size:")
parser.add_argument("burnin", type=float, help="Input generation of burn-in end:")

# Parse command-line arguments
args = parser.parse_args()

#read tables for tskit
with open('sitetable.txt') as f:
    sites = f.read()

with open('nodetable.txt') as f:
    nodes = f.read()

with open('mutationtable.txt') as f:
    mutations = f.read()

with open('edgetable.txt') as f:
    edges = f.read()

#load in the tree sequence data
ts = tskit.load_text(
    nodes = io.StringIO(nodes),
    edges = io.StringIO(edges),
    sites = io.StringIO(sites),
    mutations = io.StringIO(mutations),
    strict = False)

#ts_2 = TableCollection.tree_sequence("tables.trees")

num_samples = ts.get_sample_size()
print(f"the size of the sample of text-based ts is {num_samples}.")
#num_samples = ts_2.get_sample_size()
#print(f"the size of the sample of direct load ts is {num_samples}.")

N = args.popsize
burnin = args.burnin
G = float(N*100) #change this if you are running generations as N*10 or N*100

threshold_gen = (G-burnin)*N #this is the burnin generation threhold to discard anything that might've fixed during the burnin

#output the fixation time table for fixed mutations in final sample
with open('fixed_mut.txt', 'w', encoding='utf-8') as f:
    print('mutation', sep='\t', file=f)
    for tree in ts.trees(root_threshold=2):
        for root in tree.roots:
            # Create a list to store mutations that meet the time threshold
            filtered_mutations = []
            for mutation in tree.mutations():
                # Get the associated node ID for the mutation
                node_id = mutation.node

                # Access the time attribute of the node
                node_time = tree.time(node_id)
                
                if node_time <= threshold_gen:
                    filtered_mutations.append(mutation)
            
            # Iterate through the filtered mutations to find those where root == mutation.node
            for mutation in filtered_mutations:
                if root == mutation.node:
                    print(mutation.derived_state, sep='\t', file=f)
#fixed mutation table
mutfixtable = pd.read_csv('fixed_mut.txt', sep="\t", header=0)
#seperate mutations into beneficial and deleterious subsets
mut_ben = mutfixtable[mutfixtable["mutation"] > 0]
mut_del = mutfixtable[mutfixtable["mutation"] < 0]

#printing fluxes
print("Beneficial flux is: ", f'{((2 * (mut_ben["mutation"] + 1).apply(math.log)).sum())/(G-burnin):.10f}')
print("Deleterious flux is: ", f'{((2 * (mut_del["mutation"] + 1).apply(math.log)).sum())/(G-burnin):.10f}')

#Calculating Average branch length between pair of sample nodes
print("Calculating coalescent Ne: ", f'{(ts.diversity(mode="branch"))/(2*N)}')