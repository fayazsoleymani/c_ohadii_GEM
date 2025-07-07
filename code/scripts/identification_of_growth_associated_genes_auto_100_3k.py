from cobra.io import read_sbml_model
from cobra.sampling import sample
import os

model_100= read_sbml_model('../../data/models/iCO1515_auto_100.xml')
samples_100= sample(model_100, n= 5*10e3, method= 'optgp', thinning= 1000, processes= os.cpu_count()-4, seed= 42)
samples_100.to_csv('../../data/rxns/iCO1515_auto_100_n5000_t1000.csv',
                   sep=',', encoding='utf-8', index=False, header=True)
active_rxns_100= samples_100.columns[(samples_100 != 0).any(axis=0)].tolist()


model_3k= read_sbml_model('../../data/models/iCO1515_auto_3k.xml')
samples_3k= sample(model_3k, n= 5*10e3, method= 'optgp', thinning= 1000, processes= os.cpu_count()-4, seed= 42)
samples_3k.to_csv('../../data/rxns/iCO1515_auto_3k_n5000_t1000.csv',
                  sep=',', encoding='utf-8', index=False, header=True)
active_rxns_3k= samples_3k.columns[(samples_3k != 0).any(axis=0)].tolist()


from scipy.stats import ttest_ind
import numpy as np

significants= []
significant_genes= []
for gene in shared:
    
    gene_100= model_100.genes.get_by_id(gene)
    rxns_100= set([rxn.id for rxn in gene_100.reactions])
    
    gene_3k= model_3k.genes.get_by_id(gene)
    rxns_3k= set([rxn.id for rxn in gene_3k.reactions])
    
    rxns_intersect= rxns_100.intersection(rxns_3k)
    
    significance= True
    sample_means= []
    p_values= []
    for rxn in rxns_intersect:
        flux_sample_100= np.array(samples_100[rxn])
        flux_sample_3k= np.array(samples_3k[rxn])
            
        stat, p_value= ttest_ind(flux_sample_3k, flux_sample_100,
                                 alternative="greater")
        if p_value < 0.01:
            p_values.append(p_value.item())
            sample_means.append((np.mean(flux_sample_3k).item(), np.mean(flux_sample_100).item()))
        else:
            significance = False
            
    if significance:
        significants.append([gene_3k.id, rxns_intersect,sample_means,p_values])
        significant_genes.append(gene_3k.id)
print(len(significants), "/", len(shared))



with open('../../data/genes/significant_genes_auto_100_3k_new.csv', 'w') as file:
    writer= csv.writer(file, delimiter= ',')
    writer.writerow(['gene', "rxns", "meanFlux(auto_3k, auto_100)", "p-values"])
    for element in significants:
        writer.writerow(element)
