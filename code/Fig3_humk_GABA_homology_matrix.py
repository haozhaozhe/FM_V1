import scanpy as sc
import harmonypy as hm

import sys
sys.path.append('./scctools)
from scctools import *


# # Load and check data

#monkey
ad=sc.read('data/Figure/FM27_cell_133454_wk.h5')
gaba = ad[ad.obs['class']=='Inh',:]
gaba

#human
adhu=sc.read("./gaba_huMK_500Marker1.h5")

#match .obs

adhu.obs['class'] = adhu.obs['class_label']
adhu.obs['subclass'] = adhu.obs['subclass_label']
adhu.obs['cell_label'] = adhu.obs['cluster_label']

admk = data_subset(gaba,key = 'subclass', n=500, frac=1)


# # Matching dataset

adhu_m,admk_m =  match_gene(adhu, admk,ref_csv_path = "data/python/human2monkey.csv",
                Ln_col_name = 'Crab-eating macaque gene name',
                Ln_col_ID = 'Crab-eating macaque gene stable ID',
                Tr_species = 'Human',Ln_species = 'Monkey')


# # select HVG and merge data

adTr = adhu_m.copy()
adLn = admk_m.copy()
n_HVG = 2000
Ln_species = 'Monkey'
Tr_species = 'Human'

adLn.var[Ln_species] = adLn.var_names
adLn.var_names = adLn.var[Tr_species]

if scipy.sparse.issparse(adTr.layers['raw']):
        adTr.X = adTr.layers['raw'].todense().copy()
else: adTr.X = adTr.layers['raw'].copy()     
        
if scipy.sparse.issparse(adLn.layers['raw']):
        adLn.X = adLn.layers['raw'].todense().copy()
else: adLn.X = adLn.layers['raw'].copy()
    
adTr_selected_gene = adTr.var_names[np.array(np.std(adTr.X, axis=0).argsort())[0][-n_HVG:][::-1]]
        

adLn_selected_gene = adLn.var_names[np.array(np.std(adLn.X, axis=0).argsort())[-n_HVG:][::-1]]

select_gene  = unique(list(adTr_selected_gene)+list(adLn_selected_gene))

adTr_sele = adTr[:,select_gene]
adTr_sele.obs['source']= Tr_species

adLn_sele = adLn[:,select_gene]
adLn_sele.obs['source']= Ln_species

TrLn = adTr_sele.concatenate(adLn_sele)

def sele_dataset(adTr, adLn, n_HVG = 2000, Ln_species = 'Monkey', Tr_species = 'Human'):
#generate new dataset with HVG only
#input adTr: trainging set, full dataset; adLn:testing set, full dataset; n_HVG, number of HVG
#output Tr_sub, Ln_sub: splited Training and testing dataset, with combine calculated PCA
    adLn.var[Ln_species] = adLn.var_names
    adLn.var_names = adLn.var[Tr_species]
    
    if scipy.sparse.issparse(adTr.layers['raw']):
        adTr.X = adTr.layers['raw'].todense()
    else: adTr.X = adTr.layers['raw'].copy()     
        
    if scipy.sparse.issparse(adLn.layers['raw']):
        adLn.X = adLn.layers['raw'].todense()
    else: adLn.X = adLn.layers['raw'].copy()
    

    adTr_selected_gene = adTr.var_names[np.array(np.std(adTr.X, axis=0).argsort())[0][-n_HVG:][::-1]]
        

    adLn_selected_gene = adLn.var_names[np.array(np.std(adLn.X, axis=0).argsort())[-n_HVG:][::-1]]
    
    select_gene  = unique(list(adTr_selected_gene[0])+list(adLn_selected_gene))
    
    adTr_sele = adTr[:,select_gene]
    adTr_sele.obs['source']= Tr_species

    adLn_sele = adLn[:,select_gene]
    adLn_sele.obs['source']= Ln_species

    TrLn = adTr_sele.concatenate(adLn_sele)
     
    return TrLn

humk = TrLn

sc.pp.filter_cells(humk, min_genes=10)
sc.pp.filter_genes(humk, min_cells=2) 

humk.obs['n_counts'] = np.matrix(humk.layers['raw'].sum(axis=1)).A1
humk = run_PCA_Harmony(humk, run_Harmony = True, theta = 5,rep = 5)
humk.obs['species'] = humk.obs['source']

#  # plot matrix
plot_matrix = two_species_heatmap(humk, species_1 = 'Monkey', species_2 = 'Human',species_1_key = 'cell_label', species_2_key = 'cell_label',louvain = 2.5,figure_path = 'Figure/Fig3_gaba_heatmap.png')

plot_matrix_reorder_gaba = plot_matrix_gaba.copy()
plot_matrix_reorder_gaba = plot_matrix_reorder_gaba.reindex([
       'Inh L1-3 PAX6 DBI','Inh L1-3 PAX6 RAMP1','Inh L1-3 PAX6 KLHL4', 'Inh L1-3 PAX6 SYT10',\
       'Inh L1-3 ADARB2 MSMO1',\
        'Inh L1-3 VIP PRKCG','Inh L1-3 VIP HTR2C', 
       'Inh L1-3 VIP USP36', 'Inh L1-3 VIP SOCS2', 'Inh L1-3 VIP SCGN',\
       'Inh L1-3 VIP ABCA8','Inh L1-3 VIP RPS15A', 'Inh L1-3 VIP VIP', 
        'Inh L1-3 LAMP5 RELN','Inh L1-6 LAMP5 NKX2-1','Inh L1-6 SST CHODL','Inh L1-3 SST PENK',                                 
        'Inh L1-3 SST FUOM','Inh L1-3 SST PKIA', 'Inh L1-3 SST CPLX2',
        'Inh L4-6 SST TRDN', 
        'Inh L1-3 PVALB FILIP1','Inh L1-3 PVALB UNC5B', 
        'Inh L4-6 PVALB OSTN', ])


plot_matrix_reorder_gaba = plot_matrix_reorder_gaba[['Inh L1 PAX6 CA4', 'Inh L1 PAX6 GRIP2',
                                                     'Inh L1-6 VIP RCN1','Inh L1-3 PAX6 NABP1', \
        'Inh L1-5 VIP KCNJ2', 'Inh L1-6 VIP PENK',  'Inh L3-6 VIP KCTD13',  
        'Inh L2-4 PVALB C8orf4','Inh L5-6 PVALB STON2', 
	'Inh L3-4 PVALB HOMER3', 'Inh L3 VIP CBLN1',                    
       'Inh L1 LAMP5 NDNF',
       'Inh L5-6 LAMP5 SFTA3', 'Inh L1-6 LAMP5 CA13', 
       'Inh L6 SST NPY','Inh L3-5 SST MAFB','Inh L4-5 PVALB TRIM67',  'Inh L4-6 SST MTHFD2P6',
        'Inh L5-6 SST TH',
       'Inh L1-6 PVALB SCUBE3',  
       'Inh L1-3 PVALB WFDC2', ]]

graph = sns.heatmap(plot_matrix_reorder_gaba, cmap=cm, cbar=True, xticklabels=1,yticklabels=1, linewidth = 0.2, linecolor = 'gray')

fig = graph.get_figure()
fig.savefig('./Fig_gaba_heatmap.pdf')

# # plot matrix

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(constrained_layout=True,figsize=(5.5,4.5))
# # plot umap
gs = GridSpec(1, 1, figure=fig)
#ax1 = fig.add_subplot(gs[0, 0:])
ax1 = fig.add_subplot(gs[0, 0])
ax1 = sc.pl.umap(humk , color = 'species',show = False,
                 frameon = False, s = 20, ax = ax2)

plt.savefig("./Fig_heatmap_umap_gaba_harmony.pdf")

matplotlib.rcParams.update({'font.size': 14})
fig = plt.figure(constrained_layout=True,figsize=(4.5,4.5))
#sns.despine()

gs = GridSpec(1, 1, figure=fig)
ax1 = fig.add_subplot(gs[0, 0:])
ax2= sc.pl.umap(gaba, color=['subclass'],ax=ax1, title='Inh. Neurons', show=False, size=40, legend_loc = 'on data',legend_fontoutline = 2)

ax2.spines['bottom'].set_color('#000000')
ax2.spines['left'].set_color('#000000') 

ax2.spines['top'].set_color('#FFFFFF')
ax2.spines['right'].set_color('#FFFFFF') 

ax2.set_facecolor('white')

ax2.xaxis.label.set_size(12)
ax2.yaxis.label.set_size(12)

plt.savefig(".Fig_heatmap_umap_gaba_subclass.pdf")

