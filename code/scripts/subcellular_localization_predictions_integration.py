from scipy.io import loadmat
import csv


def dict_values_probs(input_dict):
    total= sum(input_dict.values())
    for key, value in input_dict.items():
        input_dict[key]= value/total
    return input_dict

def clean_deeploc_data(input_file):
    deeploc_data= dict()
    with open(input_file, 'r') as file:
        reader= csv.reader(file, delimiter= ',')
        header= next(reader, None)
        for row in reader:
            loc_dict= {
                'cytoplasm': float(row[3]),
                'extracellular': float(row[5]),
                'mitochondrion': float(row[7]),
                'thylakoid/plastid': float(row[8]),
                'endoplasmic_reticulum': float(row[9]),
                'lysosome/vacuole': float(row[10]),
                'peroxisome': float(row[12]),
                'chloroplast': 0
            }
            deeploc_data[row[0]]= dict_values_probs(loc_dict)
    return deeploc_data


def clean_targetp_data(input_file):
    targetp_data= dict()
    with open(input_file, 'r') as file:
        reader= csv.reader(file, delimiter= '\t')
        next(reader, None)
        next(reader, None)
        for row in reader:
            other_pred= float(row[2])
            
            loc_dict= {
                'cytoplasm': other_pred / 5,
                'mitochondrion': float(row[4]),
                'extracellular': other_pred / 5,
                'chloroplast':float(row[5]),
                'thylakoid/plastid':float(row[6]),
                'endoplasmic_reticulum': other_pred / 5,
                'lysosome/vacuole': other_pred / 5,
                'peroxisome': other_pred / 5
            }
            targetp_data[row[0]]= dict_values_probs(loc_dict)
    return targetp_data


def clean_localizer_data(input_file):
    localizer_data=dict()
    with open(input_file, 'r') as file:
        reader= csv.reader(file, delimiter= '\t')
        for _ in range(4):
            next(reader, None)
        for row in reader:
            if len(row) == 4:
                loc_dict= {
                    'chloroplast': float(row[1].split(' ')[1].replace('(', '')) if row[1].startswith('Y') else 0,
                    'mitochondrion': float(row[2].split(' ')[1].replace('(', '')) if row[2].startswith('Y') else 0
                }
                loc_dict = {k: v for k, v in loc_dict.items() if v != 0}
                other_pred= 1 - sum(loc_dict.values())
                if other_pred < 0:
                    other_pred= 0
                remained = set(subcellular_locations)- set(loc_dict.keys())
                for remain in remained:
                    loc_dict[remain]= other_pred / len(remained)    
                
                localizer_data[row[0].split(' ')[0]] = dict_values_probs(loc_dict)
    return localizer_data


def clean_mulocdeep_data(input_file):
    mulocdeep_data= dict()
    with open('../data/models/localization/MULocDeep/MULocDeep_all.txt', 'r') as file:
        reader= csv.reader(file, delimiter='\t')
        for row in reader:
            loc_dict= {
                'cytoplasm': float(row[2].split(':')[1]),
                'mitochondrion': float(row[4].split(':')[1]),
                'thylakoid/plastid': float(row[7].split(':')[1]),
                'endoplasmic_reticulum': float(row[6].split(":")[1]),
                'lysosome/vacuole': float(row[9].split(':')[1]),
                'peroxisome': float(row[10].split(':')[1])
            }
            loc_dict['chloroplast']= 0
            loc_dict['extracellular']= 0 
            mulocdeep_data[row[0].split(' ')[0].replace('>', '')] = dict_values_probs(loc_dict)

    return mulocdeep_data


def clean_busca_data(input_file):

    busca_data= dict()
    subcellular_loc_set= set()
    with open(input_file, 'r') as file:
        reader= csv.reader(file, delimiter= ',')
        header= next(reader, None)
        for row in reader:
            loc= row[2].split(':')[1]
            subcellular_loc_set.add(loc)
            loc_dict= dict()
            if loc in subcellular_locations:
                loc_dict[loc]= float(row[3])
            elif loc == 'extracellular space':
                loc_dict['extracellular']= float(row[3])
            elif loc == 'chloroplast thylakoid lumen':
                loc_dict['thylakoid/plastid']= float(row[3])
            elif loc in ['nucleus', 'plasma membrane', 'chloroplast outer membrane', 'endomembrane system',
                        'mitochondrial membrane', 'organelle membrane', 'anchored component of plasma membrane',
                        'chloroplast thylakoid membrane']:
                loc_dict = dict()
            other_pred= 1- sum(loc_dict.values())
            remained= set(subcellular_locations)- set(loc_dict.keys())
            for remain in remained:
                loc_dict[remain]= other_pred/len(remained)
            busca_data[row[0]]= loc_dict
    return busca_data



def integrate_all(genes, subcellular_locations, deeploc_data, targetp_data, localizer_data, mulocdeep_data, busca_data):

    gene2loc_pred= dict()
    for gene in genes:
        temp_dict= dict()
        for loc in subcellular_locations:
            temp_dict[loc]= deeploc_data[gene][loc] + targetp_data[gene][loc] + localizer_data[gene][loc] + mulocdeep_data[gene][loc] + busca_data[gene][loc]
        gene2loc_pred[gene]= dict_values_probs(temp_dict)

    gene2loc_pred_max_one= dict()
    for gene, loc2pred in gene2loc_pred.items():
        gene2loc_pred_max_one[gene]= dict()
        largest_pred = max(loc2pred.values())
        for loc in loc2pred:
            gene2loc_pred_max_one[gene][loc]= loc2pred[loc]/largest_pred


    return gene2loc_pred_max_one


def main():

    model= loadmat('../../data/models/c_ohadii_draft.mat')['draft']
    genes= [entry[0][0] for entry in model[0][0][13]]

    gene2loc_count= dict()
    subcellular_locations= ['cytoplasm', 'mitochondrion', 'extracellular', 'chloroplast', 'thylakoid/plastid', 'endoplasmic_reticulum', 'lysosome/vacuole', 'peroxisome']

    deeploc_data= clean_deeploc_data('../../data/genes/deeploc_all.csv')
    targetp_data= clean_targetp_data('../../data/genes/targetP.txt')
    localizer_data= clean_localizer_data('../../data/genes/Results_localizer.txt')
    mulocdeep_data= clean_mulocdeep_data('../../data/genes/MULocDeep_all.txt')
    busca_data= clean_busca_data('../../data/genes/BUSCA_all.csv')

    gene2loc_pred_max_one= integrate_all(genes, subcellular_locations, deeploc_data, targetp_data, localizer_data, mulocdeep_data, busca_data)
    

    with open('../../data/genes//GSS_file_new.csv', 'w') as file:
        writer= csv.writer(file, delimiter= ',')
        writer.writerow(['protein', 'cytoplasm', 'extracellular', 'mitochondrion', 'thylakoid/plastid', 
                         'endoplasmic_reticulum', 'lysosome/vacuole', 'peroxisome', 'chloroplast'])
        for protein, loc2pred in gene2loc_pred_max_one.items():
            writer.writerow([protein, loc2pred['cytoplasm'], loc2pred['extracellular'], loc2pred['mitochondrion'],
                            loc2pred['thylakoid/plastid'], loc2pred['endoplasmic_reticulum'],
                            loc2pred['lysosome/vacuole'], loc2pred['peroxisome'], loc2pred['chloroplast']])
    

if __name__ == "__main__":
    main()






