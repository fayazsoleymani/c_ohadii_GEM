from cobra.io import read_sbml_model, load_matlab_model
import csv
from itertools import combinations

model= read_sbml_model('../data/models/without_complexes/iCO1515_auto_100.xml')

interactions= dict()
with open('../../data/genes/string_interactions.tsv', 'r') as file:
    reader= csv.reader(file, delimiter= '\t')
    next(reader, None)
    for row in reader:
        interactions[(row[2], row[3])]= row[-1]


mapping= dict()
with open('../../data/genes/string_mapping.tsv', 'r') as file:
    reader= csv.reader(file, delimiter= '\t')
    next(reader, None)
    for row in reader:
        if row[2] not in mapping:
            mapping[row[2]] = [row[1]]
        else:
            mapping[row[2]].append(row[1])



interactions_ohadii= dict()
for pair, score in interactions.items():
    for gene_1 in mapping[pair[0]]:
        for gene_2 in mapping[pair[1]]:
            interactions_ohadii[(gene_1, gene_2)]= float(score)
len(interactions_ohadii)


import networkx as nx

confidence_threshold = 0.95
complex_counts= 0

G = nx.Graph()
for (g1, g2), score in interactions_ohadii.items():
    if score >= confidence_threshold:
        G.add_edge(g1, g2)

for rxn in model.reactions:

    genes = [g.strip() for g in rxn.gene_reaction_rule.split(' or ')]
    subgraph = G.subgraph(genes)

    used_genes = set()
    gpr_parts = []


    complex_found= False
    for component in nx.connected_components(subgraph):
        
        component = sorted(component)
        if len(component) == 1:
            gpr_parts.append(component[0])
        else:
            gpr_parts.append(' and '.join(component))
            complex_counts +=1
            complex_found = True
        used_genes.update(component)

    # Handle genes not present in the interaction graph at all
    leftovers = [g for g in genes if g not in used_genes]
    gpr_parts.extend(leftovers)

    # Assign updated GPR rule
    rxn.gene_reaction_rule = ' or '.join(gpr_parts)

    if complex_found:
        print(rxn.id, rxn.gene_reaction_rule)
        print('-' * 120)
print(complex_counts)


from cobra.io import save_matlab_model, write_sbml_model

save_matlab_model(model, '../../data/models/iCO1515_hetero_gpr.mat')
write_sbml_model(model, '../../data/models/iCO1515_hetero_gpr.xml')
