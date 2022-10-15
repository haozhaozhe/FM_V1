import scanpy as sc
import pandas as pd

# # Fig1

# # prepare data
#load data
ad =sc.read('data/Figure/FM27_cell_133454_wk.h5')

# # UMAP

sc.pl.umap(ad, color = 'cell_label')
plt.savefig('./Fig1_umap_cell_label.pdf')

sc.pl.umap(ad, color = 'Layer')
plt.savefig('./Fig1_umap_layer.pdf')


# # Cell counts bar graph

cell_cat = ad.obs['cell_label'].cat.categories
cell_counts = pd.DataFrame(ad.obs['cell_label'].value_counts())
cell_counts = cell_counts.reindex(cell_cat)

barWidth = 1
cell_type_all = list(cell_counts.index)

bars_original = np.array(cell_counts['cell_label'])
bars1 = np.log10(bars_original)

r1 = np.arange(len(bars1))
    
fig = plt.figure(figsize = (12,0.8))
plt.rcParams['font.size'] = '2'
    
gs = GridSpec(1, 1, figure=fig)
sns.set(font_scale=1)
    
    
ax = fig.add_subplot(gs[0, 0])

ax.bar(r1, bars1, width = barWidth,color = ad.uns['cell_label_colors'])

ax.set_yticks([0,2,4]) 
ax.set_xticklabels(cell_type_all,rotation=90)
ax.set_xticks(np.arange(len(cell_type_all)))

ax.set_xlim([-1,67])

ax.set_facecolor('white')

ax.spines['bottom'].set_color('#000000')
ax.spines['left'].set_color('#000000') 
    
plt.savefig("./Fig1_cell_counts.png", dpi = 600)
   


# # distribution map

layer_dist_plot = pd.crosstab(ad.obs['Layer'], ad.obs['cell_label'])

layer_dist_plot = layer_dist_plot/sum(layer_dist_plot, axis =0) 

g = sns.clustermap(layer_dist_plot, cmap = 'Reds',xticklabels = 1, yticklabels = 1,                figsize = (12,0.8),               cbar_pos = (0, 0, 0.2, 0.2),annot_kws={"size": 5},row_cluster=False,col_cluster=False)
g.gs.update(left=-0.5)

plt.savefig("./Fig1_Layer_distribution.pdf", dpi = 600)


# # dotplot and dendrogram

marker_genes = ['RORB','FEZF2','HPCAL1','THEMIS','CD63','OSTN','KAT2A','SNCB','SYP',                'CNTNAP5', 'UNC80','CORO6','SMARCA2','JAM2','CPLX2','HTR2C','CRIM1','BTBD11','PDE1A',                'SYT6','NXPH4','NPY','CBLN4','TAGLN2','WFDC2','VAMP1','CCK','RGS12']

plot_genes = ['SLC17A7','RORB','FEZF2','HPCAL1','GAD1','LHX6','NR2F2','SST','PVALB','VIP','LAMP5','MOG','PDGFRA','AQP4','FLT1','C1QB',]
sc.pl.dotplot(ad, var_names = plot_genes, groupby = 'cell_label', swap_axes = True, 
                 layer = 'raw',dendrogram = True, save = '/Fig_1_dotplot_.pdf')

sc.pl.dendrogram(ad, groupby = 'cell_label', save = '/Fig1_dendro_NC.pdf')

