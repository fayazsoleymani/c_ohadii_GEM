from cobra.io import read_sbml_model
from cobra.sampling.optgp import OptGPSampler
import os
from cobra.sampling import sample

iCre1355= read_sbml_model('iCre1355_auto.xml');

samples_iCre1355= sample(iCre1355, n= 5*10e3, method= 'optgp', thinning= 1000, processes= os.cpu_count()-4, seed= 42)
samples_iCre1355.to_csv('../../data/rxns/iCre1355_n5000_t1000.csv', sep=',',
                        encoding='utf-8', index=False, header=True)
active_rxns_iCre1355= samples_iCre1355.columns[(samples_iCre1355 != 0).any(axis=0)].tolist()


model_ours= read_sbml_model('../../data/models/iCO1515_auto_100.xml')
samples_ours= sample(model_ours, n= 5*10e3, method= 'optgp', thinning= 1000, processes= os.cpu_count()-4, seed= 42)
samples_ours.to_csv('../../data/rxns/iCO1515_n5000_t1000.csv',
                    sep=',', encoding='utf-8', index=False, header=True)
active_rxns_ours= samples_ours.columns[(samples_ours != 0).any(axis=0)].tolist()


import csv
rxn2growth_correlation= dict()
with open('../../data/rxns/iCO1515_irrev_correlations.csv', 'r') as file:
    reader= csv.reader(file, delimiter= ',')
    for row in reader:
        rxn= row[0]
        if rxn.endswith('_r'):
            rxn= rxn.replace('_r', '')
        rxn2growth_correlation[rxn]= float(row[1])


from scipy.stats import ttest_ind
import numpy as np
from statsmodels.stats.multitest import multipletests

significants_greater= []
significants_greater_short= []
significants_less= []
significants_less_short= []
count_significants_greater=0
for gene1, gene2 in alignments_ours_iCre1355:
    
    ours_gene= model_ours.genes.get_by_id(gene1)
    ours_rxns= [rxn.id for rxn in ours_gene.reactions]
    
    iCre1355_gene= iCre1355.genes.get_by_id(iCre1355_fasta2gene[gene2])
    iCre1355_rxns= [rxn.id for rxn in iCre1355_gene.reactions]
    
    pairs= [(x, y) for x in ours_rxns for y in iCre1355_rxns]
    
    significance_greater= True
    p_values_greater= []
    significance_less= True
    p_values_less= []
    
    rxn_pairs_greater= []
    rxn_mean_pairs_greater= []
    
    rxn_pairs_less= []
    rxn_mean_pairs_less= []

    correlations_with_growth= []
    
    for pair in pairs:
        rxn_ours, rxn_iCre1355= pair[0], pair[1]
        if rxn_ours in active_rxns_ours and rxn_iCre1355 in active_rxns_iCre1355:
        
            flux_sample_ours= np.array(samples_ours[rxn_ours])
            flux_sample_iCre1355= np.array(samples_iCre1355[rxn_iCre1355])
            
            _, p_value_greater= ttest_ind(flux_sample_ours, flux_sample_iCre1355,
                                     alternative="greater", equal_var=False)
            _, p_value_less= ttest_ind(flux_sample_ours, flux_sample_iCre1355,
                                       alternative="less", equal_var=False)
            
            
            if p_value_greater > 0.01:
                significance_greater = False
            else:
                p_values_greater.append(p_value_greater)
                rxn_pairs_greater.append((rxn_ours, rxn_iCre1355))
                rxn_mean_pairs_greater.append((np.mean(flux_sample_ours).item(), np.mean(flux_sample_iCre1355).item()))
                correlations_with_growth.append(rxn2growth_correlation[rxn_ours])
                
                
            if p_value_less > 0.01:
                significance_less = False
            else:
                p_values_less.append(p_value_greater)
                rxn_pairs_less.append((rxn_ours, rxn_iCre1355))
                rxn_mean_pairs_less.append((np.mean(flux_sample_ours), np.mean(flux_sample_iCre1355)))
    
    
    if len(p_values_greater) >= 1 and significance_greater:
        adj_p_values_greater = multipletests(p_values_greater, method='fdr_bh')[1][0]
        if adj_p_values_greater < 0.01:
            count_significants_greater +=1
            significants_greater.append([ours_gene.id, iCre1355_gene.id,
                                         rxn_pairs_greater, rxn_mean_pairs_greater,
                                         correlations_with_growth, adj_p_values_greater])
            significants_greater_short.append([ours_gene.id, iCre1355_gene.id,
                                               adj_p_values_greater, rxn2growth_correlation[rxn_ours]])
                 
                
    if len(p_values_less) >= 1 and significance_less:
        adj_p_values_less = multipletests(p_values_less, method='fdr_bh')[1][0]
        if adj_p_values_greater < 0.01:
            for i in range(len(p_values_less)):
                if i == 0:
                    significants_less.append([ours_gene.id, iCre1355_gene.id,
                                              rxn_pairs_less[i], rxn_mean_pairs_less[i],
                                              rxn2growth_correlation[rxn_ours], p_values_less[i]])
                else:
                    significants_less.append(["", "",
                                              rxn_pairs_less[i], rxn_mean_pairs_less[i],
                                              rxn2growth_correlation[rxn_ours], p_values_less[i]])
            significants_less_short.append([ours_gene.id, iCre1355_gene.id,
                                            adj_p_values_less, rxn2growth_correlation[rxn_ours]])

            
        
print(len(significants_greater_short), '/', len(alignments_ours_iCre1355))
print(len(significants_less), '/', len(alignments_ours_iCre1355))


with open('../../data/genes/significant_genes_iCO1515_iCre1355_shortened_new.csv', 'w') as file:
    writer= csv.writer(file, delimiter= ',')
    writer.writerow(["iCO1515Gene", "iCre1355Gene", "adjusted-p-value", "CorrelationWithGrowth"])
    for element in significants_greater_short:
        writer.writerow(element)
