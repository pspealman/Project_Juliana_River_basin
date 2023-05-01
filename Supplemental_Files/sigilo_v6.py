# -*- coding: utf-8 -*-
"""
Created on Sun Feb  9 23:03:48 2020

sigilo - for the generation and analysis of PICRUSt2 output

python sigilo --generate_heatmap -i pred_metagenome_unstrat.tsv -sig significant_objects_file -o metagenome

python sigilo.py --ko_enrichment -c pred_metagenome_contrib.tsv -t taxonomy.tsv -k ko00001.keg -pct 0.1 -pval 0.05 -o metagenome_contrib

python sigilo.py --versatility -c pred_metagenome_contrib.tsv -t taxonomy.tsv -k ko00001.keg -pct 0.1 -pval 0.05 -o metagenome_contrib

python sigilo.py --pathway_enrichment -c pred_metagenome_contrib.tsv -t taxonomy.tsv -k ko00001.keg -pct 0.05 -pval 0.05 -o metagenome_contrib_CPP

python sigilo.py --otu2taxa -f feature-table.biom.txt -t taxonomy.tsv -o feature_taxonomy.tab

Version: Public 0.1 (I refuse to give up because I haven't tried all possible ways.)
Version: Beta 0.4 (Error Construct): added pathway_enrichment
Version: Beta 0.5 (Restaurant Dismiss): added versatility
Version: Beta 0.6 (Obscure Threat): added summative taxonomic levels

@author: pspea
"""
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-hm',"--generate_heatmap", action='store_true')
parser.add_argument('-i',"--input_abundance_file")
parser.add_argument('-sig',"--significant_objects_file")

parser.add_argument('-ko',"--ko_enrichment", action='store_true')
parser.add_argument('-c',"--contrib_file")
parser.add_argument('-t',"--taxonomy")
parser.add_argument('-select', '--select_taxa_list')
parser.add_argument('-pct', '--pct_threshold')
parser.add_argument('-pval', '--pval_threshold')
parser.add_argument('-k',"--kegg_file")

parser.add_argument('-plot_ko',"--plot_ko_enrichment", action='store_true')
parser.add_argument('-s',"--taxa_source")

parser.add_argument('-v',"--versatility", action="store_true")


parser.add_argument('-pathway',"--pathway_enrichment", action='store_true')
parser.add_argument('-plot_pathway',"--plot_pathway_enrichment", action='store_true')

parser.add_argument('-o2t',"--otu2taxa", action='store_true')
parser.add_argument('-f',"--feature_file")

parser.add_argument('-o',"--output_file")

args = parser.parse_args()

#import packages:
import numpy as np
import scipy.stats as stats
import pickle

if args.taxa_source:
    taxa_source = args.taxa_source
else:
    taxa_source = 'silva'

if taxa_source.lower() == 'silva' or taxa_source.lower() == 's':
    prefixe = ['D_0__','D_1__','D_2__','D_3__','D_4__','D_5__','D_6__','D_7__','D_8__','D_9__','D_10__','D_11__','D_12__','D_13__','D_14__']
   
if 'green' in taxa_source.lower() or taxa_source.lower() == 'gg':
    prefixe = ['k__','p__','c__','o__','f__','g__','s__']
    
rank_order = ['species','genus','family','order','class','phylum','kingdom']
convert_taxa_to_rank = {'kingdom':0, 'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5, 'species': 6}

rank_order = ['species']
#convert_taxa_to_rank = {'kingdom':0, 'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5, 'species': 6}


def load_sig_obj(so_file_name = args.significant_objects_file):
    so_file = open(so_file_name)
    sig_object_set = set()
    
    for line in so_file:
        ko = line.split(',')[0]
        if len(ko) > 0:
            ko = ko.replace('"','')
            sig_object_set.add(ko)
    so_file.close()
    
    return(sig_object_set)
    
def load_ko(kegg_file_name):
    universal_ko_dict = {}
    universal_ko_lookup = {}
    
    keg_file = open(kegg_file_name)
    
    for line in keg_file:
        if line[0] =='A':
            #A09100 Metabolism
            a_ko = line[1:6]
            universal_ko_lookup[a_ko]=line[7:].strip()
    
        if line[0] == 'B':
            #B  09102 Energy metabolism
            if len(line.strip()) > 1:
                b_ko = line[3:8]
                universal_ko_lookup[b_ko]=line[8:].strip()
    
        if line[0] == 'C':
            #C    00010 Glycolysis / Gluconeogenesis [PATH:ko00010]
            c_ko = line[5:10]
            universal_ko_lookup[c_ko]=line[10:].strip()
    
        if line[0] == 'D':
            #D      K00844  HK; hexokinase [EC:2.7.1.1]
            d_ko = line[7:13]
            universal_ko_lookup[d_ko]=line[13:].strip()
            
            if d_ko not in universal_ko_dict:
                universal_ko_dict[d_ko] = {'A':set(), 'B':set(), 'C':set()}
                
            universal_ko_dict[d_ko]['A'].add(a_ko)
            universal_ko_dict[d_ko]['B'].add(b_ko)
            universal_ko_dict[d_ko]['C'].add(c_ko)
    
    keg_file.close()  
    
    return(universal_ko_dict, universal_ko_lookup)
    
def ko_to_pathway(universal_ko_dict, universal_ko_lookup):
    
    ko_super_set_type = {'Photosynthesis [PATH:ko00195]':set(), 
                         'Sulfur metabolism [PATH:ko00920]':set(), 
                         'Nitrogen metabolism [PATH:ko00910]': set(), 
                         'Methane metabolism [PATH:ko00680]': set(), 
                         'Pentose phosphate pathway [PATH:ko00030]': set(),
                         'Carbon fixation in photosynthetic organisms [PATH:ko00710]':set(),
                         'Carbon fixation pathways in prokaryotes [PATH:ko00720]':set(),
                         'CPP':set(),
                         'CNM':set(),
                         'CN':set()}
       
    ko_to_path_lookup = {}
    
    for ko in universal_ko_dict:                
        for each in universal_ko_dict[ko]['C']:
            if universal_ko_lookup[each] == 'Carbon fixation pathways in prokaryotes [PATH:ko00720]':
                ko_super_set_type['Carbon fixation pathways in prokaryotes [PATH:ko00720]'].add(ko)
                ko_to_path_lookup[ko]='Carbon fixation pathways in prokaryotes'
                
            if universal_ko_lookup[each] == 'Carbon fixation in photosynthetic organisms [PATH:ko00710]':
                ko_super_set_type['CPP'].add(ko)
                ko_to_path_lookup[ko]='CPP'
                
            if universal_ko_lookup[each] == 'Photosynthesis [PATH:ko00195]':
                ko_super_set_type['CPP'].add(ko)
                ko_to_path_lookup[ko]='CPP'
            
            if universal_ko_lookup[each] == 'Methane metabolism [PATH:ko00680]':
                ko_super_set_type['Methane metabolism [PATH:ko00680]'].add(ko)
                ko_to_path_lookup[ko]='Methane'

            if universal_ko_lookup[each] == 'Sulfur metabolism [PATH:ko00920]':
                ko_super_set_type['Sulfur metabolism [PATH:ko00920]'].add(ko)
                ko_to_path_lookup[ko]='Sulfur metabolism'
                
            if universal_ko_lookup[each] == 'Nitrogen metabolism [PATH:ko00910]':
                ko_super_set_type['Nitrogen metabolism [PATH:ko00910]'].add(ko)
                ko_to_path_lookup[ko]='Nitrogen metabolism [PATH:ko00910]'
                
            if universal_ko_lookup[each] == 'Pentose phosphate pathway [PATH:ko00030]':
                ko_super_set_type['CPP'].add(ko)
                ko_to_path_lookup[ko]='CPP'
                
    return(ko_super_set_type, ko_to_path_lookup)

def parse_line(line, runmode='log'):
    global universal_ko_lookup
    
    print(line)
    if runmode == 'log':
        KO = line.split('\t')[0].strip()
        #
        P1 = np.log(max(float(line.split('\t')[1]),1))
        P2 = np.log(max(float(line.split('\t')[2]),1))
        P3 = np.log(max(float(line.split('\t')[3]),1))
        #
        S1 = np.log(max(float(line.split('\t')[4]),1))
        S2 = np.log(max(float(line.split('\t')[5]),1))
        S3 = np.log(max(float(line.split('\t')[6]),1))        
        #
        V2 = np.log(max(float(line.split('\t')[7]),1))
        V3 = np.log(max(float(line.split('\t')[8]),1))
        #
        if KO in universal_ko_lookup:
            description = universal_ko_lookup[KO]
        else:
            description = KO
    else:
        KO = line.split('\t')[0].strip()
        #
        P1 = float(line.split('\t')[1])
        P2 = float(line.split('\t')[2])
        P3 = float(line.split('\t')[3])
        #
        S1 = float(line.split('\t')[4])
        S2 = float(line.split('\t')[5])
        S3 = float(line.split('\t')[6])        
        #
        V2 = float(line.split('\t')[7])
        V3 = float(line.split('\t')[8])
        #
        if KO in universal_ko_lookup:
            description = universal_ko_lookup[KO]
        else:
            description = KO
    
    return(KO, P1, P2, P3, S1, S2, S3, V2, V3, description)
   
def go_rep_heatmap(ko_master_dict, output_filename, max_value):
    global ko_dict
        
    title_is = ('Heatmap of Significant KOs').format()

    site_list = ["Source_1", "Source_2", "Source_3", "Valley_1", "Valley_2", "Mangrove_1", "Mangrove_2", "Mangrove_3"]
    ko_list = []
    value_array = []
    
    for ko, data in ko_master_dict.items():
        print(ko, data)
        value_array.append(data)
        ko_line = ('{}: {}').format(ko, ko_dict[ko])
        ko_list.append(ko_line)

    fig = go.Figure(data=go.Heatmap(
                    z=value_array,
                    x=site_list,
                    y=ko_list,
                    hoverongaps = False,
                    zmin=0, zmax=max_value))
    
    fig.update_layout(
            title=title_is,
            autosize=False,
            width=3600,
            height=3600,
            margin=go.layout.Margin(
                    l=50,
                    r=50,
                    b=100,
                    t=100,
                    pad=4
                    )
            )

    fig.show()
    fig.write_image(output_filename)
        
def go_median_heatmap(ko_master_dict, output_filename, max_value):
    global ko_dict
          
    title_is = ('Heatmap of Significant KOs').format()

    site_list = ["Source", "Valley", "Mangrove"]
    ko_list = []
    value_array = []
    
    for ko, data in ko_master_dict.items():
        data=ko_master_dict[ko]
        print(ko, data)
        value_array.append(data)
        ko_line = ('{}: {}').format(ko, ko_dict[ko])
        ko_list.append(ko_line)

    fig = go.Figure(data=go.Heatmap(
                    z=value_array,
                    x=site_list,
                    y=ko_list,
                    hoverongaps = False,
                    zmin=0, zmax=max_value))
    fig.update_layout(
            title=title_is,
            autosize=False,
            width=3600,
            height=3600,
            margin=go.layout.Margin(
                    l=50,
                    r=50,
                    b=100,
                    t=100,
                    pad=4
                    )
            )
    
    fig.show()
    fig.write_image(output_filename)
    
def return_log10(each_set):
    new_set = []
    for each_obs in each_set:
        if each_obs == 0:
            each_obs = 1
        else:
            each_obs = np.log10(each_obs)

        new_set.append(each_obs)
    
    return(new_set)
        
def with_sig_object():
    import plotly.graph_objects as go
    #universal_ko_dict, universal_ko_lookup = load_ko('C:/Gresham/Project_Gravimondo/Project_Osun/sigilo/ko00001.keg')
    universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
    ko_super_set_type, ko_to_path_lookup = ko_to_pathway(universal_ko_dict, universal_ko_lookup)
    so_file = (args.significant_objects_file)
    #so_file = ('C:/Gresham/Project_Gravimondo/Project_Osun/Aldex_results/aldex_significant_CNPS_KO.csv')
    ko_set = load_sig_obj(so_file)
    
    is_val_file = open('isval.tab','w')
        
    for super_condition, ko_select in ko_super_set_type.items():
        print(super_condition)
        if len(ko_select)>0:
            condition = ko_to_path_lookup[list(ko_select)[0]]
            ko_dict = {}
            rep_sig_ko_dict = {}
            median_sig_ko_dict = {}
            rep_all_ko_dict = {}
            median_all_ko_dict = {}
            sig_max_value = 0
            all_max_value = 0
            
            all_rep_heatmap_name = ('{}_all_replicates_heatmap.pdf').format(condition)
            sig_rep_heatmap_name = ('{}_0.05_replicates_heatmap.pdf').format(condition)
            all_median_heatmap_name = ('{}_all_median_heatmap.pdf').format(condition)
            sig_median_heatmap_name = ('{}_0.05_median_heatmap.pdf').format(condition)
            
            ko_list = []
            for ko in ko_select:
                ko_list.append(ko)
                
            #ko_import_file.close()
            abun_val_file = open('C:/Gresham/Project_Gravimondo/Project_Osun/Picrust2_Results/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv')
        
            for line in abun_val_file:
                #KO	P1.1	P1.2	P1.3	P2.1	P2.2	P2.3	P3.1	P3.2	P3.3	description
                if line.split('\t')[0] != 'function':
                    KO, P1, P2, P3, S1, S2, S3, V2, V3, description = parse_line(line)
                    
                    if KO in ko_set and KO in ko_list:
                        sig_max_value  = max([P1, P2, P3, S1, S2, S3, V2, V3, sig_max_value])
                                
                        if KO not in rep_sig_ko_dict:
                            ko_dict[KO]=description.strip()
                            rep_sig_ko_dict[KO]=[S1, S2, S3, V2, V3, P1, P2, P3]
                            median_sig_ko_dict[KO]=[np.median([S1, S2, S3]),np.median([V2, V3]), np.median([P1, P2, P3])]
                        else:
                            print('Error: Duplicate KO Identified')
                            quit()
                            
                    if KO in ko_list:
                        all_max_value  = max([P1, P2, P3, S1, S2, S3, V2, V3, all_max_value])
                                
                        if KO not in rep_all_ko_dict:
                            ko_dict[KO]=description.strip()
                            rep_all_ko_dict[KO]=[S1, S2, S3, V2, V3, P1, P2, P3]
                            median_all_ko_dict[KO]=[np.median([S1, S2, S3]),np.median([V2, V3]), np.median([P1, P2, P3])]
                        else:
                            print('Error: Duplicate KO Identified')
                            quit()
            
            abun_val_file.close()
            
            go_rep_heatmap(rep_sig_ko_dict, sig_rep_heatmap_name, sig_max_value)
        
            go_median_heatmap(median_sig_ko_dict, sig_median_heatmap_name, sig_max_value)
            
            go_rep_heatmap(rep_all_ko_dict, all_rep_heatmap_name, all_max_value)
        
            go_median_heatmap(median_all_ko_dict, all_median_heatmap_name, all_max_value)
            
            outline = ('For {}, {} significant of {} detected, out of {} total.').format(condition, len(rep_sig_ko_dict), len(rep_all_ko_dict), len(ko_list))
            is_val_file.write(outline)
        
    is_val_file.close()


def without_sig_object():
    #abun_val_file = open(args.input_abundance_file)
    abun_val_file = open('C:/Gresham/Project_Gravimondo/Project_Osun/Picrust2_results/KO_metagenome_out/pred_metagenome_unstrat.tsv/pred_metagenome_unstrat.tsv')
    ko_dict = {}
    master_ko_dict = {}
    median_ko_dict = {}

    max_value = 0
   
    for line in abun_val_file:
        print(line)
        #KO	P1.1	P1.2	P1.3	P2.1	P2.2	P2.3	P3.1	P3.2	P3.3	description
        if line.split('\t')[0] != 'function':
            KO, M1, M2, M3, P1, P2, P3, S1, S2, S3, V2, V3, description = parse_line(line, 'nolog')
            hV1 = np.median([V2,V3])
            #TODO make ko not list based                
            try:
                obs = np.array([[M1, M2, M3,], [P1, P2, P3], [S1, S2, S3], [hV1, V2, V3]])
                chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)

                if pval <= 10:
                    KO, M1, M2, M3, P1, P2, P3, S1, S2, S3, V2, V3, description = parse_line(line)
                    max_value  = max([M1, M2, M3, P1, P2, P3, S1, S2, S3, V2, V3, max_value])
                    
                    if KO not in ko_dict:
                        ko_dict[KO]=description.strip()
                        
                        master_ko_dict[KO]=[M1, M2, M3, P1, P2, P3, S1, S2, S3, V2, V3]
                        median_ko_dict[KO]=[np.mean([M1, M2, M3]), np.mean([P1, P2, P3]),np.mean([S1, S2, S3]),np.mean([V2, V3])]
                        print(KO, np.mean([M1, M2, M3]), np.mean([P1, P2, P3]),np.mean([S1, S2, S3]),np.mean([V2, V3]))
                    else:
                        print('Error: Duplicate KO Identified')
                        quit()
            except:
                continue
    
    abun_val_file.close()
    
    go_master_heatmap(master_ko_dict)
    
    go_median_heatmap(median_ko_dict)
    

## TODO parse fuck
def parse_taxonomy(taxa_file_name, taxa_cutoff_num):
    taxa_file = open(taxa_file_name)
    
    taxa_set = set()
    taxa_to_otu_dict = {}
    otu_to_taxa_dict = {}
        
    for line in taxa_file:
        if line[0]!='#':
            line = line.replace('"','')
            line = line.strip()
            otu = line.split('\t')[0]
            taxa = line.split('\t')[1]
            for each in prefixe:
                taxa = taxa.replace(each,'')

            taxa = taxa.replace(';','_')
            while taxa[-1] == '_':
                taxa = taxa[:-1]
                                
            if taxa.count('_') >= taxa_cutoff_num:
                taxa_set.add(taxa)
                
                if taxa.count('_') > taxa_cutoff_num:
                    taxa_list = taxa.split('_')[:taxa_cutoff_num+1]
                    taxa = ''
                    for each in taxa_list:
                        taxa+=str(each)+'_'
                    
                    if taxa[-1] == '_':
                        taxa = taxa[:-1]
                               
                if taxa not in taxa_to_otu_dict:
                    taxa_to_otu_dict[taxa] = [otu]
                else:
                    taxa_to_otu_dict[taxa].append(otu)
                    
                if otu not in otu_to_taxa_dict:
                    otu_to_taxa_dict[otu] = taxa
    
                else:
                    print('err')
                
    taxa_file.close()
        
    return(taxa_to_otu_dict, otu_to_taxa_dict, taxa_set)

    
def convert_sample_to_site(sample):
    if sample in ['P1', 'P2', 'P3']:
        return('P')
    if sample in ['S1', 'S2', 'S3']:
        return('S')
    if sample in ['V1', 'V2', 'V3']:
        return('V')
        
def parse_pathway_contrib(contrib_name,  taxa_dict, level, convert_to_site=True):
    global ordered_list, ko_to_path_lookup
    contrib = open(contrib_name)
    ko_dict = {'total':{}}
#    taxonomic_level_total = {}
    ko_tfa_dict = {}
    pathway_totals_dict = {}
    
    for line in contrib:
        if 'function' not in line:
            if convert_to_site:
                site = convert_sample_to_site(line.split('\t')[0])
            else:
                site = line.split('\t')[0]
                
            ko = line.split('\t')[1]

            if ko in ordered_list:
                otu = line.split('\t')[2]
                
                taxon_function_abun = float(line.split('\t')[6])
                
                if otu in taxa_dict:
                    taxon = taxa_dict[otu]
                                                                                                                
                    if site not in ko_dict:
                        ko_dict[site] = {}
                    if level not in ko_dict[site]:
                        ko_dict[site][level] = {}    
                    if ko not in ko_dict[site][level]:
                        ko_dict[site][level][ko] = {}                        
                    if taxon not in ko_dict[site][level][ko]:
                        ko_dict[site][level][ko][taxon] = 0
                        
                    ko_dict[site][level][ko][taxon] += taxon_function_abun
                        
                    pathway = ko_to_path_lookup[ko]
                    
                    #pathway_totals_dict[level][site][pathway]
                    if level not in pathway_totals_dict:
                        pathway_totals_dict[level] = {} 
                    if site not in pathway_totals_dict[level]:
                        pathway_totals_dict[level][site] = {} 
                    if pathway not in pathway_totals_dict[level][site]:
                        pathway_totals_dict[level][site][pathway] = 0
                                                                                          
                    pathway_totals_dict[level][site][pathway] += taxon_function_abun
                            
                    if ko not in ko_tfa_dict:
                        ko_tfa_dict[ko] = set()
                    else:
                        ko_tfa_dict[ko].add(site)
                    
                    taxon = taxon.rsplit(';',1)[0]
                         
    contrib.close()
    
    return(ko_dict, ko_tfa_dict, pathway_totals_dict)

#def parse_contrib(contrib_name,  taxa_dict, level, convert_to_site=True):
#    global ordered_list
#    contrib = open(contrib_name)
#    ko_dict = {}
##    taxonomic_level_total = {}
#    ko_tfa_dict = {}
#    
#    for line in contrib:
#        if 'function' not in line:
#            if convert_to_site:
#                site = convert_sample_to_site(line.split('\t')[0])
#            else:
#                site = line.split('\t')[0]
#                
#            ko = line.split('\t')[1]
#
#            if ko in ordered_list:
#                otu = line.split('\t')[2]
#                
#                taxon_function_abun = float(line.split('\t')[6])
#                
#                if otu in taxa_dict:
#                    taxon = taxa_dict[otu]
#                                                                                                                
#                    if site not in ko_dict:
#                        ko_dict[site] = {}
#                    if level not in ko_dict[site]:
#                        ko_dict[site][level] = {}    
#                    if ko not in ko_dict[site][level]:
#                        ko_dict[site][level][ko] = {}
#                    if taxon not in ko_dict[site][level][ko]:
#                        ko_dict[site][level][ko][taxon] = 0
#                    if 'total' not in ko_dict[site][level][ko]:
#                        ko_dict[site][level][ko]['total'] = 0
#                    
#                    ko_dict[site][level][ko][taxon] += taxon_function_abun
#                    ko_dict[site][level][ko]['total'] += taxon_function_abun
#                            
#                    if ko not in ko_tfa_dict:
#                        ko_tfa_dict[ko] = set()
#                    else:
#                        ko_tfa_dict[ko].add(site)
#                    
#                    taxon = taxon.rsplit(';',1)[0]
#                         
#    contrib.close()
#    
#    return(ko_dict, ko_tfa_dict)

def check_metabolism(ko):
    global universal_ko_dict
    
    if ko in universal_ko_dict:
        #print(universal_ko_dict[ko])
        meta_string = ('{},{},{}').format(universal_ko_dict[ko]['A'], universal_ko_dict[ko]['B'], universal_ko_dict[ko]['C'])
    else:
        meta_string = ('{}').format(ko)
        
    return(meta_string)
        
#    
    
def run_mwu(x_set, y_set):
    if sum(x_set) > 30 or sum(y_set) > 30:
        _w, p_xvy = stats.mannwhitneyu(x_set, y_set)
        return(p_xvy)
    else:
        return(1)
        
def fill_out(inlist, length):
    
    while len(inlist) < length:
        inlist.append(0)
        
    return(inlist)
    
def try_kruskal(m_is, p_is, s_is, v_is):
    if (sum(m_is) >= 5) and (sum(p_is) >= 5) and (sum(s_is) >= 5) and (sum(v_is) >= 5):
        if (m_is != p_is) and (p_is != s_is) and (s_is != v_is):
            is_list = [m_is, p_is, s_is, v_is]
            for x in range(len(is_list)):
                for y in range(len(is_list)):
                    if x != y:
                        if set(is_list[x]) == set(is_list[y]):
                            return(False)
            return(True)
            
    return(False)
    
def run_chi2(x_num, x_den, y_num, y_den):
    if (x_num) > 5 or (y_num) > 5:
        obs = np.array([[max(x_num,1), x_den], [max(y_num,1), y_den]])
        chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
        return(pval)
    else:
        return(1)
        
def run_chi2x3(x_num, x_den, y_num, y_den, z_num, z_den):
    if (x_num) > 5 or (y_num) > 5 or (z_num) > 5:
        obs = np.array([[max(x_num,1), x_den], [max(y_num,1), y_den], [max(z_num,1), z_den]])
        chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
        return(pval)
    else:
        return(1)
        
def pathway_build(taxa_cutoff_name, pct_threshold=0.05, pval_threshold=0.05):
    
    outfile_name = ('{}_{}.log').format(args.output_file, taxa_cutoff_name)
    outfile = open(outfile_name,'w')
    
    correction = 0
    pval_ct = 0
    global ko_to_path_lookup, ko_dict, pathway_totals_dict
    
    correction = len(ko_to_path_lookup)

    pathway_round = {}
    plot_round = {}
    pct_round = {}
    
    for site in ko_dict:
        if site != 'totals':
            for level in ko_dict[site]:
                for ko in ko_dict[site][level]:                      
                    for taxon, taxa_fa in ko_dict[site][level][ko].items():
                        if ko in ko_to_path_lookup:
                            pathway = ko_to_path_lookup[ko]
                           
                            if level not in pathway_round:
                                pathway_round[level] = {}
                            if taxon not in pathway_round[level]:
                                pathway_round[level][taxon] = {}
                            if pathway not in pathway_round[level][taxon]:
                                pathway_round[level][taxon][pathway] = {}
                            if site not in pathway_round[level][taxon][pathway]:
                                pathway_round[level][taxon][pathway][site] = {}
                            if ko not in pathway_round[level][taxon][pathway][site]:
                                pathway_round[level][taxon][pathway][site][ko] = []
                            
                            pathway_round[level][taxon][pathway][site][ko].append(taxa_fa) 
                                                    
    for level in pathway_round:
        for taxon in pathway_round[level]:
            for pathway in pathway_round[level][taxon]:
                temp_path = {}
                temp_total = {}
                for site in pathway_round[level][taxon][pathway]:
                    temp_path[site]=[]
                    for ko in pathway_round[level][taxon][pathway][site]:
                        for taxa_fa in pathway_round[level][taxon][pathway][site][ko]:
                            temp_path[site].append(taxa_fa)   
                    
                    #print(pathway_totals_dict)
                    temp_total[site] = pathway_totals_dict[level][site][pathway]
                
                uname = ('{}_{}_{}').format(level, taxon, pathway)
                p_is = [0]
                s_is = [0]
                v_is = [0]
                
                p_den = 1
                s_den = 1
                v_den = 1
                
                if 'P' in temp_path:
                    p_is = temp_path['P']
                    p_den = temp_total['P']
                if 'S' in temp_path:
                    s_is = temp_path['S']
                    s_den = temp_total['S']
                if 'V' in temp_path:
                    v_is = temp_path['V']
                    v_den = temp_total['V']
                    
                p_is = fill_out(p_is, max([len(p_is), len(s_is), len(v_is)]))
                s_is = fill_out(s_is, max([len(p_is), len(s_is), len(v_is)]))
                v_is = fill_out(v_is, max([len(p_is), len(s_is), len(v_is)]))
                                
                runmode = 'chi2'
                if runmode == 'kruskal':
                    if try_kruskal(p_is, s_is, v_is):
                        _stat, pval = stats.kruskal(p_is, s_is, v_is)                    
                    else:
                        pval = 1
                else:
                    pval = run_chi2x3(sum(p_is), p_den, sum(s_is), s_den, sum(v_is), v_den)

                    
                if pval*correction <= pval_threshold:
                    valcheck = set()
                    site_set = set()
                    for site in temp_path:
                        if len(temp_path[site]) >= 3:
                            site_set.add(site)
                            
                    for x_site in site_set:
                        for y_site in site_set:
                            if x_site != y_site:
                                x_set = temp_path[x_site]
                                y_set = temp_path[y_site]
                                
                                x_den = temp_total[x_site]
                                y_den = temp_total[y_site]
                                
                                if (sum(x_set) >= pct_threshold*x_den) or (sum(y_set) >= pct_threshold*y_den):
                                
                                    if runmode == 'kruskal':
                                        pval_2 = run_mwu(x_set, y_set)
                                    else:
    
                                        pval_2 = run_chi2(sum(x_set), x_den, sum(y_set), y_den)
                                    
                                    if (np.median(y_set) == 0):
                                        x_ratio = np.median(x_set)
                                    else:
                                        x_ratio = np.median(x_set)/np.median(y_set)
                                    
                                    if (np.median(x_set) == 0):
                                        y_ratio = np.median(y_set)
                                    else:
                                        y_ratio = np.median(y_set)/np.median(x_set)
                                        
                                    if ((pval_2*correction) <= pval_threshold) and ((x_ratio >= 1+pct_threshold) or (y_ratio >= 1+pct_threshold)):
                                        print(pval, pval_2, x_ratio, y_ratio)
                                        checkname = ('{}v{}').format(min(x_site,y_site),max(x_site,y_site))
                                        valcheck.add(checkname)
                        
                    if len(valcheck) >= 2:
                        p_pct = sum(p_is)/p_den
                        s_pct = sum(s_is)/s_den
                        v_pct = sum(v_is)/v_den                                    
                        outline = ('{level}\t{taxon}\t{pathway}\t{pval}\t{p_is}\t{p_den}\t{s_is}\t{s_den}\t{v_is}\t{v_den}\t{ko}\n').format(level=level, taxon=taxon, pathway=pathway, pval=(pval*correction), p_is=sum(p_is), p_den=p_den, s_is=sum(s_is), s_den=s_den, v_is=sum(v_is), v_den=v_den, ko=list(pathway_round[level][taxon][pathway][site].keys()))
 
                        outfile.write(outline)
                        pval_ct+=1
                        
                        if uname in plot_round:
                            print('uname duplicate', plot_round[uname])
                        
                        plot_round[uname]={'P':p_is, 'S': s_is, 'V': v_is, 'ko':list(pathway_round[level][taxon][pathway][site].keys())}
                        pct_round[uname]={'P':p_pct, 'S': s_pct, 'V': v_pct, 'ko':list(pathway_round[level][taxon][pathway][site].keys())}
                        
    
    pickle_name = ('pct_{}_{}.p').format(args.output_file, taxa_cutoff_name)
    pickle.dump(pct_round, open(pickle_name, 'wb'))
    
    pickle_name = ('plot_{}_{}.p').format(args.output_file, taxa_cutoff_name)
    pickle.dump(plot_round, open(pickle_name, 'wb'))
        
    print(pval_ct)
    outfile.close()
    
def clean_up_name(name):
    name = name.strip()
    
    while name[-1] == ',':
        name=name[:-1]
        name = name.strip()
    
    return(name)
#    
def make_plot_figures():
    for rank in ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']:
        file_name = ('plot_metagenome_contrib_{}.p').format(rank)
    
        file = open(file_name,'rb')
        plot_round = pickle.load(file)
        file.close()
        
        for uname in plot_round:
            p_tag = ('Pristine')
            s_tag = ('Source')
            v_tag = ('Valley')
            x_data = p_tag, s_tag, v_tag
            
            p_list = return_log10(plot_round[uname]['P'])
            s_list = return_log10(plot_round[uname]['S'])
            v_list = return_log10(plot_round[uname]['V'])
            y_data = p_list, s_list, v_list

            print(uname)
            uname = clean_up_name(uname)
          
            fig = go.Figure()
            colors = 'rgba(44, 160, 101, 0.5)', 'rgba(93, 164, 214, 0.5)', 'rgba(128, 0, 128, 0.5)'
    
            outfile_name = ('funmaps/{}_enrichment.pdf').format(uname)
            print(outfile_name)
            
            for xd, yd, cls in zip(x_data, y_data, colors):
                    fig.add_trace(go.Box(
                        #,
                        y=yd,
                        name=xd,
                        
                        boxpoints='all',
                        notched=True,
                        jitter=0.5,
                        whiskerwidth=0.2,
                        fillcolor=cls,
                        line_color=cls,
                        marker_size=5,
                        line_width=1,
                        showlegend=False)
                    )
                    
            fig.update_layout(
                title=uname.replace(',', '\n'),
                xaxis_title="Sample Site",
                yaxis_title="Log10(Relative Functional Abundance)",
                        font=dict(
                            family="Courier New, monospace",
                            size=10,
                            color="#7f7f7f"
                        )
            )
            
                    
            fig.show()
            fig.write_image(outfile_name)
    
def apply_first_round(pct_threshold=0.2, pval_threshold=0.05):
    total_ct = 0
    pct_ct = 0
    pval_ct = 0
    global ko_dict, ko_tfa_dict, taxa_gfc

    second_round = {}
    
    for site in ko_dict:
        for level in ko_dict[site]:
            if level == 4:
                for ko in ko_dict[site][level]:
                    correction = len(ko_dict[site][level])*len(ko_dict[site])*len(ko_dict)
                    for taxon, taxa_fa in ko_dict[site][level][ko].items():
                        if taxon != 'total':
                            total_ct += 1                               
                            taxon_function_abun = ko_dict[site][level][ko]['total']
        
                            if (taxa_fa/taxon_function_abun) >= pct_threshold:
                                pval = stats.binom_test(taxa_fa, n=(taxon_function_abun))*correction
                                pct_ct += 1
                                if pval <= pval_threshold:
                                    print(ko, site, level, taxon)
                                    print(taxa_fa, taxon_function_abun, pval, correction)
                                    
                                    pval_ct+=1
                                    for ko_is in universal_ko_dict[ko]['C']:
                                        meta_string = universal_ko_lookup[ko_is]
                                    
                                        if taxon not in second_round:
                                            second_round[taxon] = {}
                                            
                                        if meta_string not in second_round[taxon]:
                                            second_round[taxon][meta_string] = set()
                
                                        second_round[taxon][meta_string].add(ko)                     
    print(total_ct, pct_ct, pval_ct)
    return(second_round)

def calc_sites(level, ko, taxon):
    # m_fa, p_fa, s_fa, v_fa = calc_sites(level, ko, taxon)
    global ko_dict
    
    m_fa, p_fa, s_fa, v_fa = 0, 0, 0, 0
    for site in ko_dict:
        try:
            taxa_fa = ko_dict[site][level][ko][taxon]
            #taxon_function_abun = ko_dict[site][level][ko]['total']
            
            if site == 'M':
                m_fa = taxa_fa
                
            if site == 'P':
                p_fa = taxa_fa
                
            if site == 'S':
                s_fa = taxa_fa
                
            if site == 'V':
                v_fa = taxa_fa
        
        except:
            taxa_fa = 0
                    
    return(m_fa, p_fa, s_fa, v_fa)
        
def apply_second_round(measures=3, pct_threshold=0.2, pval_threshold=0.05, select_taxa=False):
    global ko_dict, second_round, universal_ko_lookup
    
    if select_taxa:
        outname = ('{}.tsv').format(select_taxa.rsplit('__',1)[1])
    else:
        outname = ('global_criteria_{}_{}.tsv').format(pct_threshold, pval_threshold)
    out_temp = open(outname,'w')
       
    for site in ko_dict:
        for level in ko_dict[site]:
            if level == 4:
                for ko in ko_dict[site][level]:
                    for taxon, taxa_fa in ko_dict[site][level][ko].items():
                        if taxon != 'total':
                            correction = len(ko_dict[site][level])*len(ko_dict[site])*len(ko_dict)       
                            process = True
                            
                            if select_taxa:
                                process = False
                                if select_taxa in taxon:
                                    process = True
                            
                            if process:
                                for ko_is in universal_ko_dict[ko]['C']:
                                    meta_string = universal_ko_lookup[ko_is]
                                
                                    if taxon in second_round:
                                        if meta_string in second_round[taxon]:
                                            if len(second_round[taxon][meta_string]) >= measures:
                                                taxon_function_abun = ko_dict[site][level][ko]['total']

                                                if (taxa_fa/taxon_function_abun) >= pct_threshold:
                                                    print(ko, taxon, site, taxa_fa, taxon_function_abun) 
                                                    pval = stats.binom_test(taxa_fa, n=(taxon_function_abun))*correction
                                                    
                                                    if pval <= 0.05:
                                                        if (taxa_fa/taxon_function_abun) >= pct_threshold:                                    
                                                            function = universal_ko_lookup[ko]
                                                            m_fa, p_fa, s_fa, v_fa = calc_sites(level, ko, taxon)
                                                            outline = ('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n').format(site, meta_string, level, taxon, ko, function, (taxa_fa/taxon_function_abun), m_fa, p_fa, s_fa, v_fa, pval, taxa_fa, taxon_function_abun)
                                                            print(outline)
                                                            out_temp.write(outline)                                           
    out_temp.close()
    
def parse_select_taxa():
    select_taxa = set()
    select_file = args.select_taxa_list
    
    for line in select_file:
        select_taxa.add(line.strip())
    
    select_file.close()
    return(select_taxa)  
    
def check_for_null(taxa):
    taxa_list = taxa.split(';')
    
    taxa_string = ''
    
    for each in ['D_0__', 'D_1__', 'D_2__', 'D_3__', 'D_4__', 'D_5__', 'D_6__']:
        for every in taxa_list:
            if each in every:
                new_taxa, test = every.split(each)
                print(test)
                
                if len(test) > 0:
                    taxa_string += every + ';'
                else:
                    return(taxa_string[:-1])
                    
    return(taxa_string[:-1])
    
def clean_up_taxa(taxa):
    if 'Unassigned' in taxa:
        return(False)
    
    if 'D_7__' in taxa:
        taxa = taxa.split('D_7__')[0]
        
    taxa = check_for_null(taxa)
    
    taxa = taxa.replace('D_0__','k__')
    taxa = taxa.replace('D_1__','p__')
    taxa = taxa.replace('D_2__','c__')    
    taxa = taxa.replace('D_3__','o__')
    taxa = taxa.replace('D_4__','f__')
    taxa = taxa.replace('D_5__','g__')
    taxa = taxa.replace('D_6__','s__')
    
    return(taxa)
    
def check_sum(abundance_list):
    total = 0

    for each in abundance_list:
        total += float(each)
        
    return(total)
    
if args.generate_heatmap:
    if args.significant_objects_file:
        with_sig_object()
    else:
        without_sig_object()
        
if args.pathway_enrichment:
    if args.pct_threshold:
        pct_threshold = float(args.pct_threshold)
    else:
        pct_threshold = 1
        
    if args.pval_threshold:
        pval_threshold = float(args.pval_threshold)
    else:
        pval_threshold = 0.05
    
    ordered_list = []
    
    universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
    ko_super_set_type, ko_to_path_lookup = ko_to_pathway(universal_ko_dict, universal_ko_lookup)
    
    taxonomy_name = args.taxonomy
    contrib_name = args.contrib_file
    
    for superset, ko_list in ko_super_set_type.items():
        for ko in ko_list:
            ordered_list.append(ko)
    

    
    for taxa_cutoff_name in rank_order:
        taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]
        
        print('For ', taxa_cutoff_name)
        print('... running parse_taxonomy')
        taxa_to_otu_dict, taxa_dict, taxa_set = parse_taxonomy(taxonomy_name, taxa_cutoff_num)
        print('... running parse_contrib')
        ko_dict, ko_tfa_dict, pathway_totals_dict = parse_pathway_contrib(contrib_name, taxa_dict, taxa_cutoff_name)
        print('... running pathway_build')
        pathway_build(taxa_cutoff_name, pct_threshold, pval_threshold)
        
#if args.ko_enrichment:
#    #TODO fix 'select_taxa'
#    #select_taxa_set = parse_select_taxa()
#        
#    universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
#
#    taxonomy_name = args.taxonomy
#    path = args.path_to_ko_files
#    contrib_name = args.contrib_file
#    
#    pct_threshold = float(args.pct_threshold)
#    pval_threshold = float(args.pval_threshold)
#                            
#    ko_super_set_type = {'Sulfur metabolism [PATH:ko00920]':set(), 'Nitrogen metabolism [PATH:ko00910]': set(), 'Carbohydrate metabolism': set(), 'Pentose phosphate pathway [PATH:ko00030]': set()}
#    ko_to_path_lookup = {}
#    
#    for ko in universal_ko_dict:
##        for each in universal_ko_dict[ko]['B']:
##            if universal_ko_lookup[each] == 'Carbohydrate metabolism':
##                ko_super_set_type['Carbohydrate metabolism'].add(ko)
##                ko_to_path_lookup[ko]='Carbohydrate metabolism'
#        for each in universal_ko_dict[ko]['C']:
#            if universal_ko_lookup[each] == 'Sulfur metabolism [PATH:ko00920]':
#                ko_super_set_type['Sulfur metabolism [PATH:ko00920]'].add(ko)
#                ko_to_path_lookup[ko]='Sulfur metabolism'
#                
#            if universal_ko_lookup[each] == 'Nitrogen metabolism [PATH:ko00910]':
#                ko_super_set_type['Nitrogen metabolism [PATH:ko00910]'].add(ko)
#                ko_to_path_lookup[ko]='Nitrogen metabolism'
#                
#            if universal_ko_lookup[each] == 'Pentose phosphate pathway [PATH:ko00030]':
#                ko_super_set_type['Pentose phosphate pathway [PATH:ko00030]'].add(ko)
#                ko_to_path_lookup[ko]='Pentose phosphate pathway'
#    
#    ordered_list = []
#    
#    for superset, ko_list in ko_super_set_type.items():
#        for ko in ko_list:
#            ordered_list.append(ko)
#        
#    #TODO replace ko_to_object_dict with universal_ko_lookup
#    #metabolism_dict, ordered_list, ko_to_object_dict = map_kos(path)
#    #universal_ko_dict, universal_ko_lookup = load_ko()
#    print('running  parse_taxonomy')
#    taxa_to_otu_dict, taxa_dict, taxa_set = parse_taxonomy(taxonomy_name)
#    print('running parse_contrib')
#    ko_dict, ko_tfa_dict = parse_contrib(contrib_name)
#
#    print('running first_build')
#    second_round = apply_first_round(pct_threshold, pval_threshold)
#       
#    print('running second_round')
#    apply_second_round(3, pct_threshold, pval_threshold)
#
#    # Todo fix select_taxa
#    for select_taxa in select_taxa_set:
#        apply_second_round(3, pct_threshold, pval_threshold, select_taxa)
#        
if args.otu2taxa:                
    taxonomy = open(args.taxonomy)
    
    otu_dict = {}
    
    for line in taxonomy:
        if 'Feature' not in line:
            otu, taxa, _confidence = line.split('\t')
            
            if otu not in otu_dict:
                otu_dict[otu] = clean_up_taxa(taxa)
            else:
                print('something wrong')
    taxonomy.close()
    
    feature = open(args.feature_file)
    
    taxa_count_check = {}
    
    for line in feature:
        if line[0]!='#':
            line = line.strip()
            # OTU_ID	P1-1	P1-2	P1-3	P2-1	P2-2	P2-3	P3-1	P3-2	P3-3
            otu, P1, P2, P3, S1, S2, S3, V1, V2 = line.split('\t')
            taxa = otu_dict[otu.strip()]
            
            if taxa not in taxa_count_check:
                taxa_count_check[taxa] = check_sum([P1, P2, P3, S1, S2, S3, V1, V2])
            else:
                taxa_count_check[taxa] += check_sum([P1, P2, P3, S1, S2, S3, V1, V2])
    
    feature.close()
    
    feature = open(args.feature_file)
    outfile = open(args.output_file, 'w')
    
    header = ('otu\ttaxa\tP1\tP2\tP3\tS1\tS2\tS3\tV1\tV2\n')
    
    outfile.write(header)
    for line in feature:
        if line[0]!='#':
            line = line.strip()
            # OTU_ID	P1-1	P1-2	P1-3	P2-1	P2-2	P2-3	P3-1	P3-2	P3-3
            otu, P1, P2, P3, S1, S2, S3, V1, V2 = line.split('\t')
            taxa = otu_dict[otu.strip()]

            if taxa_count_check[taxa]>=5:
                if taxa:
                    outline = ('{otu}\t{taxa}\t{P1}\t{P2}\t{P3}\t{S1}\t{S2}\t{S3}\t{V1}\t{V2}\n').format(
                            otu=otu, taxa=taxa, 
                            P1=P1, P2=P2, P3=P3, 
                            S1=S1, S2=S2, S3=S3,
                            V1=V1, V2=V2)
                    
                    outfile.write(outline)
                    
    outfile.close()
    feature.close()
    
    
##
if args.plot_ko_enrichment:
    import plotly.graph_objects as go
    
    def mod_null_set(isset):
        if len(isset) <1:
            isset.append(0)
        return(isset)
            
    convert_taxa_to_rank = {'kingdom':0, 'phylum':1, 'class':2, 'order':3, 'family':4, 'genus':5, 'species': 6}
    convert_rank_to_taxa = {0:'kingdom', 1:'phylum', 2:'class', 3:'order', 4:'family', 5:'genus', 6:'species'}
    enrich_dict = {}
    infile = open('C:/Gresham/Project_Gravimondo/Project_Osun/sigilo/global_criteria_new_0.2_0.05.tsv')
    
    for line in infile:
        if line[0] != '#':
        #site, meta_string, level, taxon, ko, function, (taxa_fa/taxon_function_abun), m_fa, p_fa, s_fa, v_fa, pval, taxa_fa, taxon_function_abun
            site = line.split('\t')[0]
            path = line.split('\t')[1]
            level = line.split('\t')[2]
            taxa = line.split('\t')[3]
            ko = line.split('\t')[4]
            m_fa = float(line.split('\t')[7])
            p_fa = float(line.split('\t')[8])
            s_fa = float(line.split('\t')[9])
            v_fa = float(line.split('\t')[10])
                           
            if taxa not in enrich_dict:
                enrich_dict[taxa] = {}
                
            if path not in enrich_dict[taxa]:
                enrich_dict[taxa][path] = {}
                
            if ko not in enrich_dict[taxa][path]:
                enrich_dict[taxa][path][ko]=[m_fa, p_fa, s_fa, v_fa]
                
            else:
                m_o, p_o, s_o, v_o = enrich_dict[taxa][path][ko]
                m_o += m_fa
                p_o += p_fa
                s_o += s_fa
                v_o += v_fa
                enrich_dict[taxa][path][ko] = m_o, p_o, s_o, v_o
    
    outfile = open('C:/Gresham/Project_Gravimondo/Project_Osun/sigilo/binned_criteria_new_0.2_0.05.tsv', 'w')

    ct = 0        
    for taxa in enrich_dict:
        for path in enrich_dict[taxa]:
            ko_set = enrich_dict[taxa][path]
                        
            if len(ko_set) >= 3:    
                if '/' in path:
                    path = path.replace('/','_or_')

                if ':' in path:
                    path = path.replace(':','_')

                if '/' in taxa:
                    taxa = taxa.replace('/','_or_')
                
                taxa_cutoff_name = convert_rank_to_taxa[taxa.count(';')]
    
                i_tag = ('Impacted').format(path)
                p_tag = ('Pristine').format(path)
                s_tag = ('Source').format(path)
                v_tag = ('Valley').format(path)
                x_data = i_tag, p_tag, s_tag, v_tag
                
                i_list = []
                p_list = []
                s_list = []
                v_list = []
                
                for ko, set_list in ko_set.items():
                    log_list = (return_log10(set_list))
                                               
                    i_list.append(log_list[0])
                    p_list.append(log_list[1])
                    s_list.append(log_list[2])
                    v_list.append(log_list[3])
                    
                i_list = mod_null_set(i_list)
                p_list = mod_null_set(p_list)
                s_list = mod_null_set(s_list)
                v_list = mod_null_set(v_list)
                
                y_data = i_list, p_list, s_list, v_list
            
                outline = ('{}\t{}\t{}\t{}\t{}\n').format(level, site, taxa, path, ko_set)
                outfile.write(outline)
                print(outline)
                ct+=1
                
                
                fig = go.Figure()
                colors = 'rgba(255, 144, 14, 0.5)', 'rgba(44, 160, 101, 0.5)', 'rgba(93, 164, 214, 0.5)', 'rgba(128, 0, 128, 0.5)'


                outfile_name = ('taxa_enrichment_{}_{}_enrichment.pdf').format(taxa_cutoff_name, path)
                print(outfile_name)
                
                for xd, yd, cls in zip(x_data, y_data, colors):
                        fig.add_trace(go.Box(
                            #,
                            y=yd,
                            name=xd,
                            
                            boxpoints='all',
                            notched=False,
                            jitter=0.5,
                            whiskerwidth=0.2,
                            fillcolor=cls,
                            line_color=cls,
                            marker_size=5,
                            line_width=1,
                            showlegend=False)
                        )
                        
                fig.update_layout(
                    title=taxa.replace(';','\n'),
                    xaxis_title="Sample Site",
                    yaxis_title="Log10(Relative Functional Abundance)",
#                    font=dict(
#                        family="Courier New, monospace",
#                        size=18,
#                        color="#7f7f7f"
#                    )
                )
                        
                fig.show()
                #1/0
                fig.write_image(outfile_name)
                #1/0
                        #fig.close()
    
    print(ct)
    outfile.close()
    
def define_type(s_taxa, p_taxa, v_taxa, raw_s_fa, raw_p_fa, raw_v_fa):
    coord_dict = {'S':[0,0], 'V':[0,0], 'P':[0,0]} 
    s_fa = np.median(raw_s_fa)
    p_fa = np.median(raw_p_fa)
    v_fa = np.median(raw_v_fa)
    
    coord_dict = {'S':[s_taxa, s_fa], 'V':[v_taxa, v_fa], 'P':[p_taxa, p_fa]}
    
    return(coord_dict)      

def plot_funbubbles():
    import plotly.graph_objects as go
    import pickle 

    file = open('C:/Gresham/Project_Gravimondo/Project_Osun/CPP/pct_metagenome_contrib_CPP_family.p','rb')
    fun_dict = pickle.load(file)
    file.close()
    
    condition_dict = {}
    site_abundance_dict = {}
    
    taxon_set = set()
        
    for each in fun_dict:
        taxa = each.split('_')
        condition = taxa[-1]
        if condition not in condition_dict:
            condition_dict[condition] = set()
        
        taxa_list = taxa[1:-1]
        taxon = ''
        for taxa in taxa_list:
            taxon += taxa +'_'
        taxon = taxon[:-1]
        
        condition_dict[condition].add(taxon)
        taxon_set.add(taxon)
        
        for site in fun_dict[each]:
            if site != 'ko':
                if condition not in site_abundance_dict:
                    site_abundance_dict[condition] = {}
                if taxon not in site_abundance_dict[condition]:
                    site_abundance_dict[condition][taxon] = {'S':0, 'V':0, 'P':0}
                  
#                val = np.median(fun_dict[each][site])
#                
#                if val == 0:
#                    val = 1 
#                
#                val = np.log10(val)
#                
#                if val < 0:
#                    val = 0
                
                site_abundance_dict[condition][taxon][site] += fun_dict[each][site]
                        
    for condition, taxa_set in condition_dict.items():
        x_compound_dict = {'S':[], 'V':[], 'P':[]}
        y_compound_dict = {'S':[], 'V':[], 'P':[]}
    
        for taxon in taxa_set:
            for site in site_abundance_dict[condition][taxon]:
                val = site_abundance_dict[condition][taxon][site]

                x_compound_dict[site].append(val)
                y_compound_dict[site].append(taxon)
            
        outfile_name = ('family_{}.pdf').format(condition)
        
        import plotly.graph_objects as go
        
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x=x_compound_dict['S'],
            y=y_compound_dict['S'],
            marker=dict(color='rgba(93, 164, 214, 0.5)', size=6),
            mode="markers",
            name="Source",
        ))
        
        fig.add_trace(go.Scatter(
            x=x_compound_dict['V'],
            y=y_compound_dict['V'],
            marker=dict(color='rgba(128, 0, 128, 0.5)', size=6),
            mode="markers",
            name="Valley",
        ))
        
        fig.add_trace(go.Scatter(
            x=x_compound_dict['P'],
            y=y_compound_dict['P'],
            marker=dict(color='rgba(44, 160, 101, 0.5)', size=6),
            mode="markers",
            name="Mangrove",
        ))

        fig.update_layout(title=condition,
                          xaxis_title="Percent Relative Abundance",
                          yaxis_title="Taxa",
                          font_size=10,
                          width=1500,
                          height=600)
        
        fig.show()

        fig.write_image(outfile_name)
        
        print(len(taxon_set))
        
def apply_versatility():
    global ko_dict
    
    versa_dict = {}

    for site in ko_dict:
        for level in ko_dict[site]:
            if level not in versa_dict:
                versa_dict[level] = {}
            if site not in versa_dict[level]:
                # a set of all taxa in a site at a level
                #   this is for normalization ot the KO specific numbers of taxa
                versa_dict[level][site] = {'total_taxa': set()}
                
            for ko in ko_dict[site][level]:
                if ko not in versa_dict[level][site]:
                    # the KO of every site should have:
                    #   how many taxa have the ko:
                    #       versa_dict[level][site][ko]['taxa_num'].add(taxon)
                    #   a taxa specific functional abundance:
                    #       versa_dict[level][site][ko][taxon] += (ko_dict[site][level][ko][taxon])
                    #   a distribution of all functional abundances
                    #       versa_dict[level][site][ko]['fa_dist'].append(versa_dict[level][site][ko][taxon])
                    versa_dict[level][site][ko] = {'taxa_num':set(), 'fa_dist':[]}
                
                for taxon in ko_dict[site][level][ko]:
                    versa_dict[level][site][ko]['taxa_num'].add(taxon)
                    versa_dict[level][site]['total_taxa'].add(taxon)
                    
                    if taxon not in versa_dict[level][site][ko]:
                        versa_dict[level][site][ko][taxon] = 0
                    
                    #The KO functional abundance of the a specific taxa
                    #   not normalised
                    #   ko_dict[site][level][ko]['total']
                    versa_dict[level][site][ko][taxon] += (ko_dict[site][level][ko][taxon])
                    versa_dict[level][site][ko]['fa_dist'].append(versa_dict[level][site][ko][taxon])

    fet_ct = 0
    versatile_dict = {}
    
    for level in versa_dict:
        temp_site_dict = {}
        site_list = {}
        
        for site in versa_dict[level]:
            # a set of all taxa in a site at a level
            total_taxa = len(versa_dict[level][site]['total_taxa'])

            for ko in versa_dict[level][site]:
                #filters 'total_taxa'
                if 'taxa_num' in versa_dict[level][site][ko]:
                    #only consider those sites with ko in each replicate
                    if ko not in site_list:
                        site_list[ko] = set()
                        
                    site_list[ko].add(site)
                    
                    #   how many taxa have the ko:
                    taxa_num = len(versa_dict[level][site][ko]['taxa_num'])
                    fa_dist = versa_dict[level][site][ko]['fa_dist']
                    
                    if ko not in temp_site_dict:
                        temp_site_dict[ko]={}
                        
                    if site not in temp_site_dict[ko]:
                        temp_site_dict[ko][site]={'taxa_num': 0, 'total_taxa': total_taxa, 'fa_dist': []}
                        
                    temp_site_dict[ko][site]['taxa_num'] += taxa_num
                    if len(fa_dist) > 1:
                        for each in fa_dist:
                            temp_site_dict[ko][site]['fa_dist'].append(each)
                    else:
                        temp_site_dict[ko][site]['fa_dist'].append(fa_dist[0])                                       
                    #temp_site_dict[ko][site]['total_taxa'] = total_taxa
                        
    
        for ko in temp_site_dict:
            process =  True
            #ko_sites = site_list[ko]
#            print(ko_sites)
            #for each in ['S', 'V', 'P']:
#            for each in ['S1', 'S2', 'S3', 'V1', 'V2', 'V3', 'P1', 'P2', 'P3']:
#                if each not in ko_sites:
#                    process = False
                    
            if process:
                uname = ('{}_{}').format(level, ko)
                print(uname)
#                for x_sample_set in [['M1', 'M2', 'M3'], ['P1', 'P2', 'P3']]:
#                    for y_sample_set in [['M1', 'M2', 'M3'], ['P1', 'P2', 'P3']]:
#                        if x_sample_set != y_sample_set:
                x_set = []
                x_t_set = []
                x_fa = []
                x_ratio = []
                
                y_set = []
                y_t_set = []
                y_fa = []
                y_ratio = []
                
                z_set = []
                z_t_set = []
                z_fa = []
                z_ratio = []
                
                for x in ['S1', 'S2', 'S3']:
                    if x in temp_site_dict[ko]:
                        x_num = temp_site_dict[ko][x]['taxa_num']
                        x_total = temp_site_dict[ko][x]['total_taxa']
                        x_fa_dist = temp_site_dict[ko][x]['fa_dist']
                        
                    else:
                        x_num = 0
                        x_total = 1
                        x_fa_dist = [1]
                    
                    x_set.append(x_num)
                    x_t_set.append(x_total)
                    
                    if len(x_fa_dist) > 1:
                        for each in x_fa_dist:
                            x_fa.append(each)
                            x_ratio.append((x_num/x_total)/each)
                    else:
                        x_ratio.append((x_num/x_total)/fa_dist[0])
                        x_fa.append(fa_dist[0])
                        
                    
                for y in ['V1', 'V2']:
                    if y in temp_site_dict[ko]:
                        y_num = temp_site_dict[ko][y]['taxa_num']
                        y_total = temp_site_dict[ko][y]['total_taxa']
                        y_fa_dist = temp_site_dict[ko][y]['fa_dist']
                    else:
                        y_num = 0
                        y_total = 1
                        y_fa_dist = [1]
                        
                    y_set.append(y_num)
                    y_t_set.append(y_total)
                    
                    if len(y_fa_dist) > 1:
                        for each in y_fa_dist:
                            y_fa.append(each)
                            y_ratio.append((y_num/y_total)/each)
                    else:
                        y_ratio.append((y_num/y_total)/fa_dist[0])
                        y_fa.append(fa_dist[0])
                        
                for z in ['P1', 'P2', 'P3']:
                    if z in temp_site_dict[ko]:
                        z_num = temp_site_dict[ko][z]['taxa_num']
                        z_total = temp_site_dict[ko][z]['total_taxa']
                        z_fa_dist = temp_site_dict[ko][z]['fa_dist']
                    else:
                        z_num = 0
                        z_total = 1
                        z_fa_dist = [1]
                    
                    z_set.append(z_num)
                    z_t_set.append(z_total)
                    
                    if len(z_fa_dist) > 1:
                        for each in z_fa_dist:
                            z_fa.append(each)
                            z_ratio.append((z_num/z_total)/each)
                    else:
                        z_ratio.append((z_num/z_total)/fa_dist[0])
                        z_fa.append(fa_dist[0])
                        
                print(x_ratio)
                print(y_ratio)
                print(z_ratio)
                
                _u, mwu_pval_1 = run_mwu(x_ratio, y_ratio) 
                _u, mwu_pval_2 = run_mwu(y_ratio, z_ratio)
                _u, mwu_pval_3 = run_mwu(x_ratio, z_ratio)

                if min(mwu_pval_1, mwu_pval_2, mwu_pval_3) <= 0.05:
#                obs = np.array([x_set, y_set, z_set])
#                chi2, pval, dof, expected = stats.chi2_contingency(obs, correction=True)
                
#                if pval <= 0.05:
#                    if uname not in versatile_dict:
#                        versatile_dict[uname] = {}
#                        
#                    versatile_dict[uname]['M_fa'] = x_fa
#                    versatile_dict[uname]['P_fa'] = y_fa
#                    versatile_dict[uname]['M_num'] = sum(x_set)
#                    versatile_dict[uname]['P_num'] = sum(y_set)
#                    versatile_dict[uname]['M_norm'] = sum(x_t_set)
#                    versatile_dict[uname]['P_norm'] = sum(y_t_set)
                                                    
                #_od, fet_pval= stats.fisher_exact([[sum(x_set), sum(x_t_set)-sum(x_set)],[sum(y_set), sum(y_t_set)-sum(y_set)]])
                #if pval <= 0.05:
                    print(ko)
                    if uname not in versatile_dict:
                        versatile_dict[uname] = {}
                        
                    #coord_dict = define_type((sum(x_set)/sum(x_t_set)), (sum(y_set)/sum(y_t_set)), x_fa, y_fa)
                        
                    versatile_dict[uname]['S_fa'] = x_fa
                    versatile_dict[uname]['V_fa'] = y_fa
                    versatile_dict[uname]['P_fa'] = z_fa
                    versatile_dict[uname]['S_num'] = sum(x_set)
                    versatile_dict[uname]['V_num'] = sum(y_set)
                    versatile_dict[uname]['P_num'] = sum(z_set)
                    versatile_dict[uname]['S_norm'] = sum(x_t_set)
                    versatile_dict[uname]['V_norm'] = sum(y_t_set)
                    versatile_dict[uname]['P_norm'] = sum(z_t_set)
                    #versatile_dict[uname]['istype'] = istype
                    #versatile_dict[uname]['ratio_taxa'] = ratio_taxa
                    #versatile_dict[uname]['ratio_fa'] = ratio_fa
                    versatile_dict[uname]['S_ratio'] = x_ratio
                    versatile_dict[uname]['V_ratio'] = y_ratio
                    versatile_dict[uname]['P_ratio'] = z_ratio                    
#                print(x_fa)
#                print(y_fa)

                    fet_ct+=1
#                    if level not in versatile_dict:
#                        versatile_dict[level] = set()
#                    versatile_dict[level].add(ko)
#                    
                #print(obs, pval, chi_ct, fet_pval, fet_ct)
                    
    pickle_name = ('plot_versatility.p').format()
    pickle.dump(versatile_dict, open(pickle_name, 'wb'))
                     
    return(versa_dict)

def make_versatility_figures():
    import plotly.express as px
    convert_rank_to_taxa = {0:'kingdom', 1:'phylum', 2:'class', 3:'order', 4:'family', 5:'genus', 6:'species'}
    
    file = open('plot_versatility.p','rb')
    versatile_dict = pickle.load(file)
    file.close()
    
    for uname in versatile_dict:
        i_tag = ('Impacted, n = {} of {}: {}. Median = {}').format(versatile_dict[uname]['M_num'], versatile_dict[uname]['M_norm'], round(versatile_dict[uname]['M_num']/versatile_dict[uname]['M_norm'],2), np.median(versatile_dict[uname]['M_fa']))
        p_tag = ('Pristine, n = {} of {}: {}. Median = {}').format(versatile_dict[uname]['P_num'], versatile_dict[uname]['P_norm'], round(versatile_dict[uname]['P_num']/versatile_dict[uname]['P_norm'],2), np.median(versatile_dict[uname]['P_fa']))
        #s_tag = ('Source').format(path)
        #v_tag = ('Valley').format(path)
        x_data = i_tag, p_tag
        
        i_list = return_log10(versatile_dict[uname]['M_fa'])
        p_list = return_log10(versatile_dict[uname]['P_fa'])
        #s_list = return_log10(plot_round[uname]['S'])
        #v_list = return_log10(plot_round[uname]['V'])
        y_data = i_list, p_list
    
        ko = uname.split('_')[1]
        
        taxa_level = convert_rank_to_taxa[int(uname.split('_')[0])]
                        
        fig = go.Figure()
        colors = 'rgba(255, 144, 14, 0.5)', 'rgba(44, 160, 101, 0.5)'

        outfile_name = ('{}_{}_versatility.pdf').format(taxa_level, ko)
        print(outfile_name)
        
        for xd, yd, cls in zip(x_data, y_data, colors):
                fig.add_trace(go.Box(
                    #,
                    y=yd,
                    name=xd,
                    boxpoints='all',
                    notched=True,
                    jitter=0.5,
                    whiskerwidth=0.2,
                    fillcolor=cls,
                    line_color=cls,
                    marker_size=5,
                    line_width=1,
                    showlegend=False)
                )
                
        fig.update_layout(
            title=ko,
            xaxis_title="Sample Site",
            yaxis_title="Log10(Relative Functional Abundance)",
                    font=dict(
                        family="Courier New, monospace",
                        size=10,
                        color="#7f7f7f"
                    )
        )
        
                
        fig.show()
        fig.write_image(outfile_name)

    
    x_list = []
    y_list = []
    size_list = []
    istype_list = []
    uname_list = []
    type_dict = {}    
    for uname in versatile_dict:
        if 'K02591' in uname:
            size_list.append(int(uname.split('_')[0])**3)
            x_list.append(np.log10(versatile_dict[uname]['ratio_taxa']))
            y_list.append(np.log10(versatile_dict[uname]['ratio_fa']))
            
            istype = versatile_dict[uname]['istype']
            istype_list.append(istype)
            uname_list.append(uname)
            
            if istype not in type_dict:
                type_dict[istype] = []
            type_dict[istype].append(uname)
        
    fig = px.scatter(x=x_list, y=y_list, color=istype_list, size=size_list)
    fig.show()
    fig.write_image('Nitrogen_example.pdf')
    
    print(type_dict)

if args.versatility:
    #TODO fix 'select_taxa'
    #select_taxa_set = parse_select_taxa()
    taxonomy_name = args.taxonomy        
    universal_ko_dict, universal_ko_lookup = load_ko(args.kegg_file)
    contrib_name = args.contrib_file
    
    pct_threshold = float(args.pct_threshold)
    pval_threshold = float(args.pval_threshold)
    
    ko_super_set_type, ko_to_path_lookup = ko_to_pathway(universal_ko_dict, universal_ko_lookup)
    
    ordered_list = []
    for superset, ko_list in ko_super_set_type.items():
        for ko in ko_list:
            ordered_list.append(ko)
                                    
    for taxa_cutoff_name in rank_order:
            taxa_cutoff_num = convert_taxa_to_rank[taxa_cutoff_name]
            
            print('For ', taxa_cutoff_name)
            print('... running parse_taxonomy')
            taxa_to_otu_dict, taxa_dict, taxa_set = parse_taxonomy(taxonomy_name, taxa_cutoff_num)
            print('... running parse_contrib')
            ko_dict, ko_tfa_dict, pathway_totals_dict = parse_pathway_contrib(contrib_name, taxa_dict, taxa_cutoff_name, False)
    
#    print('running parse_contrib')
#    ko_dict, ko_tfa_dict = parse_pathway_contrib(contrib_name, False)

            print('running apply_versatility')
            versa_dict = apply_versatility()
       

        