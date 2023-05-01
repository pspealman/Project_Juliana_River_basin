# -*- coding: utf-8 -*-
"""
Created on Tue Apr 18 16:01:02 2023

@author: pspea
"""
import pandas as pd
import numpy as np

import plotly.io as pio
pio.renderers.default = "browser"
import plotly.graph_objects as go
import plotly.express as px

from scipy.stats import sem

base_dir = ('C:/Gresham/tiny_projects/Project_Osun/qiime_results/')

sample_set = ['ancom_bc_Spring_Mangrove',
              'ancom_bc_Spring_Valley',
              'ancom_bc_Valley_Mangrove']

sig_taxa = set()

for sample in sample_set:
    filename = ('{}/{}/{}_q_val_slice.csv').format(base_dir, sample, sample)
    
    infile = open(filename)
    
    for line in infile:
        if '(Intercept)' not in line:
            taxa, _i, qval = line.split(',')
            
            if float(qval) <= 0.05:
                sig_taxa.add(taxa)
        

infile_name = ('C:/Gresham/tiny_projects/Project_Osun/qiime_results/osun_genus_counts.tsv')
counts_df = pd.read_csv(infile_name, sep='\t', index_col=0)
counts_dict = counts_df.to_dict(orient='index')

total_dict = {}

for taxa in counts_dict:
    for sample in counts_dict[taxa]:
        if sample not in total_dict:
            total_dict[sample] = 0
            
        total_dict[sample] += counts_dict[taxa][sample]
        
        
rel_counts_dict = {}

for taxa in sig_taxa:
    
    if taxa not in rel_counts_dict:
        rel_counts_dict[taxa] = {}
    
    for sample in counts_dict[taxa]:
        rel_count = counts_dict[taxa][sample]/total_dict[sample]
        simple_sample = sample[0]
        
        if simple_sample not in rel_counts_dict[taxa]:
            rel_counts_dict[taxa][simple_sample] = []

        rel_counts_dict[taxa][simple_sample].append(rel_count)
        
distance_dict = {}
    
for taxa in rel_counts_dict:
    distance_dict[taxa] = 0
    for simple1 in rel_counts_dict[taxa]:
        for simple2 in rel_counts_dict[taxa]: 
            if simple1 != simple2:
                val1 = np.median(rel_counts_dict[taxa][simple1])
                val2 = np.median(rel_counts_dict[taxa][simple2])
                
                distance_dict[taxa] += abs(val1-val2) 
    
distance_dict_sorted = dict(sorted(distance_dict.items(), key=lambda item: item[1],reverse=True))

high_scoring_taxa = list(distance_dict_sorted.keys())
high_scoring_taxa = high_scoring_taxa[0:10]

distance_dict = {}
    
for taxa in rel_counts_dict:
    distance_dict[taxa] = 0
    for simple in rel_counts_dict[taxa]:
        distance_dict[taxa] += np.median(rel_counts_dict[taxa][simple])
    
distance_dict_sorted = dict(sorted(distance_dict.items(), key=lambda item: item[1],reverse=True))

high_scoring_taxa = list(distance_dict_sorted.keys())
high_scoring_taxa = high_scoring_taxa[0:10]

cor_ct = 0
cor_limit = len(px.colors.qualitative.Vivid)
taxa_to_color = {}
for taxa in high_scoring_taxa:
    #if taxa in major_set:
        cor = px.colors.qualitative.Vivid[cor_ct]
        cora = cor.replace(')',',0.1)').replace('rgb','rgba')
        taxa_to_color[taxa] = {'cor':cor, 'cora':cora}
        
        cor_ct += 1
        #loop color if need be
        if cor_ct >= cor_limit:
            cor_ct = 0

output_figure_name = ('C:/Gresham/tiny_projects/Project_Osun/qiime_results/top10_pop_ancom_results.pdf')

fig = go.Figure()

simple_to_int = {'S':0, 'V':1, 'P':2}

def return_pretty_name(taxa):
    taxa_count = taxa.count(';')
    
    for index in range(taxa_count, 0, -1):
        taxa_level = taxa.split(';')[index]            
        pretty_name = taxa_level.rsplit('__',index)[1]
        
        if pretty_name not in ['', 'uncultured', 'group']:
            return(pretty_name)

for taxa in high_scoring_taxa:
    if taxa in rel_counts_dict:
        pretty_name = return_pretty_name(taxa)
        
        spring_list = rel_counts_dict[taxa]['S']
        valley_list = rel_counts_dict[taxa]['V']
        mangrove_list = rel_counts_dict[taxa]['P']
                   
        #get color from dict
        cor = taxa_to_color[taxa]['cor']
        cora = taxa_to_color[taxa]['cora']
        
        #x-axis
        x_forward = []
        x_reverse = []
        # calculate y-axis values
        y = []
        y_upper = []
        y_lower = []
        
        for simple in ['S','V','P']:
            islist = rel_counts_dict[taxa][simple]
            
            x_forward.append(simple_to_int[simple])
            
            y_med = np.median(islist)
            y_sem = sem(islist)
            
            y.append(y_med)
            y_upper.append(y_med+y_sem)
            y_lower.append(y_med-y_sem)
            
        y_lower = y_lower[::-1]
        x_reverse = x_forward[::-1]
                    
        #add std shades, except background
        fig.add_trace(go.Scatter(
            x=x_forward+x_reverse,
            y=y_upper+y_lower,
            fill='toself',
            fillcolor=cora,
            line_color='rgba(255,255,255,0)',
            showlegend=False,
            name=pretty_name))
        
        #add solid median lines
        fig.add_trace(go.Scatter(
            x=x_forward+x_reverse, 
            y=y,
            line_color=cor,
            name=pretty_name))
        
    # fig.update_traces(line_shape='spline',
    #                   line_smoothing=1.3,
    #                   selector=dict(type='scatter'))
    

fig.show()
fig.write_image(output_figure_name)
