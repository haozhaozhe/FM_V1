import scanpy as sc
import sys
sys.path.append('./scctools/scctools')
from scctools import *

# # Load Data

ad = sc.read('./FM27_cell_133454_wk.h5')

import matplotlib.patheffects as pe
matplotlib.rcParams.update({'font.size': 12})
fig = plt.figure()

gs = GridSpec(1,3, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])

ax1 = sc.pl.umap(ad, color = 'cell_label',ax=ax1,frameon = False, title=' ',show=False,legend_loc = None, size = 0.4)
ax2 = sc.pl.umap(ad,color = 'OSTN',ax=ax2,frameon = False,size = 0.4, title='OSTN',show=False,legend_loc = None, cmap = 'RdGy_r')
ax3 = sc.pl.umap(ad,color = 'LGI2',ax=ax3,frameon = False,size = 0.3,  title='LGI2',show=False,legend_loc = None, cmap = 'RdGy_r')

ax1.spines['bottom'].set_color('#000000')
ax1.spines['left'].set_color('#000000')
ax1.spines['top'].set_color('#FFFFFF')
ax1.spines['right'].set_color('#FFFFFF') 

ax1.xaxis.label
ax1.xaxis.label.set_size(12)
ax1.yaxis.label.set_size(12)

plt.savefig("./Fig4_OSTN_marker_genes.pdf")


# # Violin Plot for OSTN

admk_OSTN = ad[ad.obs['subclass']=='IT OSTN',:]

admm = sc.read('./mm_SS.h5')
admm_rorb = admm[admm.obs['subclass_label']=='L4/5 IT CTX',:]

import matplotlib.patheffects as pe
matplotlib.rcParams.update({'font.size': 12})
fig = plt.figure()
gs = GridSpec(1,2, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])

ax1 = sc.pl.violin(admk_OSTN, keys = ['OSTN'], ax = ax1,use_raw=False, layer= 'raw', show = False, jitter = False)
ax2 = sc.pl.violin(admm_rorb, keys = ['Ostn'], ax = ax2, use_raw = False, layer = 'raw',show = False, jitter = False)

ax1.xaxis.label.set_size(12)
ax1.yaxis.label.set_size(12)

ax2.xaxis.label.set_size(12)
ax2.yaxis.label.set_size(12)

ax1.set_xlabel('Monkey\nIT OSTN')
ax2.set_xlabel('Mouse\nL4/5 IT')

ax1.set_ylabel('Nomalized Expression')
ax2.set_ylabel('')

ax2.yaxis.set_ticklabels([])

ax1.set_ylim([-0.1, 3])
ax2.set_ylim([-0.1, 3])

ax1.spines['bottom'].set_color('#000000')
ax1.spines['left'].set_color('#000000')
ax1.spines['top'].set_color('#FFFFFF')
ax1.spines['right'].set_color('#FFFFFF') 

ax2.spines['bottom'].set_color('#FFFFFF')
ax2.spines['left'].set_color('#FFFFFF')
ax2.spines['top'].set_color('#FFFFFF')
ax2.spines['right'].set_color('#FFFFFF') 
ax2.set_yticks( [])
plt.savefig("./Fig_OSTN_violin_mkmm.pdf")


# # plot violin

glut = ad[ad.obs['class']=='Exc',:]
basic_gene = ['SNAP25','SLC17A7','RORB']
OSTN_special = ['OSTN','LGI2','SMYD2','SPON1','NSG1','HAPLN4','CERK']
intersect_gene = basic_gene+OSTN_special
glut_color = dict(zip(list(glut.obs['cell_label'].cat.categories)),list((glut.uns['cell_label_colors']))))
sc.pl.stacked_violin(glut, groupby='cell_label',figsize = (12,6),linewidth = 0.2, var_names= intersect_gene, jitter=False, swap_axes = True,dendrogram = False,layer = 'CPM',palette=glut_color,save = '/OSTN_violin_cluster.pdf')

# # Plot GO

def plot_enrich_square(data, n_terms=20, cmap = 'Reds_r', vmax_scale = 2, vmin_scale = 2, save=False,):
# Plotting GO enrichment terms, only keep p value, remove intersection
# data,  Pandas Dataframe output by gprofiler
# n_terms, number of terms to plot
# save, the path to save the file, if empty, not save anything
# vmax_scale, vmin_scale, factors to adjust display color

    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import seaborn as sb
    from matplotlib import colors
    from matplotlib import rcParams
   
    def scale_data_5_75(data):
        mind = np.min(data)
        maxd = np.max(data)
    
        if maxd == mind:
            maxd=maxd+1
            mind=mind-1
        
        drange = maxd - mind
        return ((((data - mind)/drange*0.70)+0.05)*100)
    # Test data input
    if not isinstance(data, pd.DataFrame):
        raise ValueError('Please input a Pandas Dataframe output by gprofiler.')
        
    if not np.all([term in data.columns for term in ['p_value', 'name', 'intersection_size']]):
        raise TypeError('The data frame {} does not contain enrichment results from gprofiler.'.format(data))
    
    data_to_plot = data.iloc[:n_terms,:].copy()
    data_to_plot['go.id'] = data_to_plot.index
    data_to_plot = data_to_plot.sort_values(by=['p_value'])
    
    min_pval = data_to_plot['p_value'].min()
    max_pval = data_to_plot['p_value'].max()
    
    # Scale intersection_size to be between 5 and 75 for plotting
    #Note: this is done as calibration was done for values between 5 and 75
    #data_to_plot['scaled.overlap'] = scale_data_5_75(data_to_plot['intersection_size'])
    
    norm = colors.LogNorm(min_pval, max_pval)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])

    rcParams.update({'font.size': 14, 'font.weight': 'bold'})

    #sb.set(style="whitegrid")

    path = plt.scatter(x=[1] *n_terms, y="name", c='p_value', cmap=cmap, 
                       norm=colors.LogNorm(min_pval, max_pval),
                       vmax = max_pval*vmax_scale, vmin = min_pval/vmin_scale,
                       data=data_to_plot, linewidth=1, edgecolor="k",marker="s", s = 300)
    ax = plt.gca()
    ax.invert_yaxis()

    if save: ax.set_title(save.split('/')[-1])
    ax.set_ylabel('')
    #ax.set_xlabel('Gene ratio', fontsize=14, fontweight='bold')
    ax.xaxis.grid(False)
    ax.yaxis.grid(False)
    
    ax.spines['bottom'].set_color('#FFFFFF')
    ax.spines['left'].set_color('#FFFFFF')
    ax.spines['top'].set_color('#FFFFFF')
    ax.spines['right'].set_color('#FFFFFF') 

    ax.set_xticks([])
    ax.tick_params(axis='y', length=0)

    
    # Shrink current axis by 20%
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    # Get tick marks for this plot
    #Note: 6 ticks maximum
    min_tick = np.floor(np.log10(min_pval)).astype(int)
    max_tick = np.ceil(np.log10(max_pval)).astype(int)
    tick_step = np.ceil((max_tick - min_tick)/6).astype(int)
    
    # Ensure no 0 values
    if tick_step == 0:
        tick_step = 1
        min_tick = max_tick-1
    
    ticks_vals = [10**i for i in range(max_tick, min_tick-1, -tick_step)]
    ticks_labs = ['$10^{'+str(i)+'}$' for i in range(max_tick, min_tick-1, -tick_step)]

    #Colorbar
    fig = plt.gcf()
    cbaxes = fig.add_axes([1.0, 0.15, 0.1, 0.2])
    cbar = ax.figure.colorbar(sm, ticks=ticks_vals, shrink=0.5, anchor=(0,0.1), cax=cbaxes)
    cbar.ax.set_yticklabels(ticks_labs,fontsize=12,fontweight='normal')
    cbar.set_label("Adjusted p-value", fontsize=12,fontweight='normal')
    cbar.ax.invert_yaxis() 
    if save:
        plt.savefig(save+'.pdf', dpi=600, format='pdf')
        plt.savefig(save+'.png', dpi=600, format='png')
    plt.show()

OSTN_UP_KEGG_results = pd.read_csv('./Supp_6_RORB_OSTN_DOWN_KEGG.csv')OSTN_UP_BP_results = pd.read_csv('./Supp_6_RORB_OSTN_DOWN_KEGG.csv')
OSTN_DOWN_KEGG_results = pd.read_csv('./Supp_6_RORB_OSTN_DOWN_KEGG.csv')
OSTN_DOWN_BP_results = pd.read_csv('./Supp_6_RORB_OSTN_DOWN_BP.csv')

OSTN_UP_results = OSTN_UP_KEGG_results.append(OSTN_UP_BP_results)
OSTN_DOWN_results = OSTN_DOWN_KEGG_results.append(OSTN_DOWN_BP_results)

plot_list = [1,7,9,12,13, 17,18,20,21,19,22,24, 45 ,63,38]
OSTN_UP_results_plot = OSTN_UP_results.iloc[plot_list,:]
OSTN_UP_results_plot = OSTN_UP_results_plot.sort_values(by=['p_value'])
plot_enrich_square(OSTN_UP_results_plot,n_terms=15, save = './Fig4_OSTN_UP_GO_square.pdf')

plot_list = [0,1,2,5,11,13,14,16,17,22,28,29,31,32,36]
OSTN_DOWN_results_plot = OSTN_DOWN_results.iloc[plot_list,:]
OSTN_DOWN_results_plot = OSTN_DOWN_results_plot.sort_values(by=['p_value'])
plot_enrich_square(OSTN_DOWN_results_plot,n_terms=15,cmap = 'Blues_r', save = './Fig4_OSTN_DOWN_GO_square.pdf')


# # Heatmap

heatmap_genes = ['OSTN','CRYM', 'ELAVL2','KCNIP2','MET','EFHD2','ESRRG','SPON1','LMO4','NECAB1','NGEF','CCK','KCNH3','RPH3A']
glut_sub= glut[glut.obs['cell_label'].isin(['Exc L4-6 RORB OSTN','Exc L1-6 RORB FAM198B']),:]
matplotlib.rcParams['figure.dpi'] = 600
matplotlib.rcParams.update({'font.size': 5})
sc.pl.heatmap(glut_sub, heatmap_genes, groupby='subclass', figsize=(2,3),log = True,standard_scale = 'var',
              cmap='seismic',var_group_rotation=90, save = "/Fig4_OSTN_DE_heatmap.pdf")

