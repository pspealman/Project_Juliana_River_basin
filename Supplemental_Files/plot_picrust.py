# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 21:36:17 2023

@author: pspea
"""
import pandas as pd
import numpy as np

import plotly.io as pio
pio.renderers.default = "browser"
import plotly.graph_objects as go

def make_figure(ko_set, threshold, ko_infile_name, ko_pred_file_name, output_figure_name):
    
    ko_lookup_dict = {}
    
    ko_infile = open(ko_infile_name)
    for line in ko_infile:
        ko = line.split('\t')[0].split('ko:')[1]
        if ';' in line.split('\t')[1]:
            name = line.split('\t')[1].split(';')[0]
            desc = line.split('\t')[1].split(';')[1].strip()
            
        else:
            desc = line.split('\t')[1].strip()
            name = 'desc_' + desc[0:4]
        
        ko_lookup_dict[ko] = {'name':name,
                              'desc':desc}
        
    ko_infile.close()
    
        
    counts_df = pd.read_csv(ko_pred_file_name, sep='\t', index_col=0)
    counts_dict = counts_df.to_dict(orient='index')
    
    for ko in ko_lookup_dict:
        if ko in counts_dict:
            for sample in counts_dict[ko]:
                simple_sample = sample[0]
                
                if simple_sample not in ko_lookup_dict[ko]:
                    ko_lookup_dict[ko][simple_sample] = []
                    
                ko_lookup_dict[ko][simple_sample].append(counts_dict[ko][sample])
                
    ko_list = []
    for ko in ko_lookup_dict:
        ismin = np.inf
        ismax = -1*np.inf
        for simple_sample in ['S', 'V', 'P']:
            val = np.log2(np.median(ko_lookup_dict[ko][simple_sample]))
            if val < ismin:
                ismin = val
            if val > ismax:
                ismax = val
        if abs(ismax - ismin) > threshold:
            ko_list.append(ko)
            
    #ko_list = list(ko_lookup_dict.keys())            
    ko_list.sort()
    
    x_list = []
    z_list = []
    
    for ko in ko_list:
        x = ko_lookup_dict[ko]['name']
        z = []
        for simple_sample in ['S', 'V', 'P']:
            val = np.log2(np.median(ko_lookup_dict[ko][simple_sample]))
            z.append(val)
            
        x_list.append(x)
        z_list.append(z)
    
    fig = go.Figure(data=go.Heatmap(
                       z=z_list,
                       y=x_list,
                       x=['Spring', 'Valley', 'Mangrove'],
                       hoverongaps = False,
                       zmin = -100,
                       zmax = 20))
    
    
    fig.update_layout(
        title=ko_set+'_threshold_'+str(threshold),
        xaxis_title="Site",
        yaxis_title="KO name",
        legend_title="log2 Median",
    )
    
    fig.update_layout(
    autosize=False,
    width=250,
    height=(len(ko_list)*20)+200,
    )
    
    fig.show()
    
    fig.write_image(output_figure_name)

base_dir = ('C:/Gresham/tiny_projects/Project_Osun/picrust/')
ko_pred_file_name = ('{}/ko_picrust2_feature-table.tsv').format(base_dir)

for threshold in set([15]):


    ko_set = ('ko01504_Antimicrobial_resistance_genes')
    ko_infile_name = ('{}/{}.txt').format(base_dir, ko_set)
    output_figure_name = ('{}/{}_threshold_{}.pdf').format(base_dir, ko_set, threshold)
    
    make_figure(ko_set, threshold, ko_infile_name, ko_pred_file_name, output_figure_name)
    
    ko_set = ('ko02042_Bacterial_toxins')
    ko_infile_name = ('{}/{}.txt').format(base_dir, ko_set)
    output_figure_name = ('{}/{}_threshold_{}.pdf').format(base_dir, ko_set, threshold)
    
    make_figure(ko_set, threshold, ko_infile_name, ko_pred_file_name, output_figure_name)
    
    ko_set = ('ko_metabolism')
    ko_infile_name = ('{}/{}.txt').format(base_dir, ko_set)
    output_figure_name = ('{}/{}_threshold_{}.pdf').format(base_dir, ko_set, threshold)
    
    make_figure(ko_set, threshold, ko_infile_name, ko_pred_file_name, output_figure_name)
    
    ko_set = ('ko_ko01110_Biosynthesis_of_secondary_metabolites')
    ko_infile_name = ('{}/{}.txt').format(base_dir, ko_set)
    output_figure_name = ('{}/{}_threshold_{}.pdf').format(base_dir, ko_set, threshold)
    
    make_figure(ko_set, threshold, ko_infile_name, ko_pred_file_name, output_figure_name)
    
    ko_set = ('ko01120_Microbial_metabolism_in_diverse_environments')
    ko_infile_name = ('{}/{}.txt').format(base_dir, ko_set)
    output_figure_name = ('{}/{}_threshold_{}.pdf').format(base_dir, ko_set, threshold)
    
    make_figure(ko_set, threshold, ko_infile_name, ko_pred_file_name, output_figure_name)



