
import scanpy as sc
import sys
sys.path.append('./scctools_0.4')
from scctools import *

# # load data

ad = sc.read('./FM27_cell_133454_wk.h5')
glut = FM27_cell[FM27_cell.obs['class']=='Exc',:]


# # violin plot

key = glut
matplotlib.rcParams.update({'font.size': 14})
basic_gene = ['SLC17A7','HPCAL1','NPY']
L23NPY_special = ['HTR6','KCNN4','KCNK12','NOV','ALCAM','DRD3','GNG4']
intersect_gene = basic_gene+L23NPY_special
g  = sc.pl.stacked_violin(key, groupby='cell_label',figsize = (7,3),linewidth = 0.1,
                     var_names= intersect_gene, jitter=False,palette=glut.uns['cell_label_colors'],
                     swap_axes = True,dendrogram = False,layer = 'CPM',save = '/Fig3_Exc_NPY_violin.pdf')


# # GO plot

def plot_enrich_square(data, n_terms=20, cm = 'Reds_r', save=False):
# Plotting GO enrichment terms, only keep p value, remove intersection
# data,  Pandas Dataframe output by gprofiler
# n_terms, number of terms to plot
# save, the path to save the file, if empty, not save anything
    
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
    data_to_plot = data.iloc[:n_terms,:].copy()
    data_to_plot['go.id'] = data_to_plot.index

    min_pval = data_to_plot['p_value'].min()
    max_pval = data_to_plot['p_value'].max()
    
    # Scale intersection_size to be between 5 and 75 for plotting
    #Note: this is done as calibration was done for values between 5 and 75
    #data_to_plot['scaled.overlap'] = scale_data_5_75(data_to_plot['intersection_size'])
    
    norm = colors.LogNorm(min_pval, max_pval)
    sm = plt.cm.ScalarMappable(cmap=cm, norm=norm)
    sm.set_array([])

    rcParams.update({'font.size': 14, 'font.weight': 'normal'})

    #sb.set(style="whitegrid")

    path = plt.scatter(x=[1] *n_terms, y="name", c='p_value', cmap=cm, 
                       norm=colors.LogNorm(min_pval, max_pval),
                       vmax = max_pval*2, vmin = min_pval/2,
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
    cbar.set_label("Adjusted p-value", fontsize=12)
    cbar.ax.invert_yaxis() 
    if save:
        plt.savefig(save+'.pdf', dpi=600, format='pdf')
        plt.savefig(save+'.png', dpi=600, format='png')
    #plt.show()

L23L23npy_UP_enrich_results_KEGG=pd.read_csv('./Sup_2_L23ExcNPY_UP_KEGG.csv')

figsize(0.8,6)
plot_enrich_square(L23L23npy_UP_enrich_results_KEGG, n_terms=15, cm = 'Reds_r',save = './ExcL23npy_UP_square' )

# # plot UMAP

matplotlib.rcParams.update({'font.size': 14})
#fig = plt.figure(constrained_layout=True)
fig = plt.figure()

gs = GridSpec(1,5, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[0, 3])
ax5 = fig.add_subplot(gs[0, 4])

ax1 = sc.pl.umap(FM27_cell, ax = ax1, color = ['cell_label'], s = 0.5, show = False, frameon = False, cmap = 'RdGy_r', legend_loc = False) 

ax2 = sc.pl.umap(FM27_cell, ax = ax2, color = ['NPY'], s = 0.5, show = False, frameon = False, cmap = 'RdGy_r')

ax3 = sc.pl.umap(FM27_cell, ax = ax3, color = ['DRD3'], s =7, show = False, frameon = False, cmap = 'RdGy_r')

ax4 = sc.pl.umap(FM27_cell, ax = ax4, color = ['SLC17A7'],s = 1, show = False, frameon = False, cmap = 'RdGy_r')

ax5 = sc.pl.umap(FM27_cell, ax = ax5, color = ['GAD1'], s = 1, show = False, frameon = False, cmap = 'RdGy_r') 

fig.tight_layout(w_pad=1)
plt.savefig("./Fig3_Exc_NPY_UMAP_all.pdf")

