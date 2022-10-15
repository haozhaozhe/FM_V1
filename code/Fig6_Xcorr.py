
import pickle
import scanpy as sc
from matplotlib import ticker

import sys
sys.path.append('./scctools')
from scctools import *
from adjustText import adjust_text


# # load data

ad = sc.read('./FM27_cell_133454_wk.h5')
all_gene = list(ad.var_names)

es_mk=pd.read_csv('./res_mk.csv', index_col = 0)
res_mm=pd.read_csv('./res_mm.csv', index_col = 0)
res_hu=pd.read_csv('./res_hu.csv', index_col = 0)


# # crosscorr

with open('/media/zhe/zhe3/MonkeyPFC/data/python/GO_CrossCorr/go_all_dict_20210417.pckl', 'rb') as f:
    go_all_dict = pickle.load(f)
go_all_dict['long-term synaptic potentiation and its regulation'] = unique(list(go_all_dict['long-term synaptic potentiation']) + list(go_all_dict['regulation of long-term synaptic potentiation']))

# # finalize caluclation

# ### prepare loop

key_list_hu = (res_hu.index)
key_list_mk = (res_mk.index)
key_list_mm = (res_mm.index)

cluster_list = list(set(key_list_hu).intersection(set(key_list_mk)).intersection(set(key_list_mm)))
geneset_list = list(go_all_dict.keys()) 

# ### Loop
humk_corr_rho = pd.DataFrame()
humk_corr_pval = pd.DataFrame()

mmmk_corr_rho = pd.DataFrame()
mmmk_corr_pval = pd.DataFrame()

mmhu_corr_rho = pd.DataFrame()
mmhu_corr_pval = pd.DataFrame()

for clust in cluster_list:
    print(' ')
    print(clust)
    for geneset_name in geneset_list:
        gene_list_all = list(go_all_dict[geneset_name])
        gene_list_overlap = (list(set(gene_list_all) & set(all_gene)))
        print(geneset_name)
        print(len(gene_list_overlap))
        
        if len(gene_list_overlap) <40: continue # skip the geneset with low gene numbers
            
        res_mk_sele = res_mk.loc[clust,gene_list_overlap]
        res_hu_sele = res_hu.loc[clust,gene_list_overlap]
        res_mm_sele = res_mm.loc[clust,gene_list_overlap]
        
        rho_humk, pval_humk = stats.spearmanr(res_mk_sele, res_hu_sele)
        rho_mmmk, pval_mmmk = stats.spearmanr(res_mk_sele, res_mm_sele)
        rho_mmhu, pval_mmhu = stats.spearmanr(res_mm_sele, res_hu_sele)
        
        humk_corr_rho.loc[clust,geneset_name] = rho_humk
        humk_corr_pval.loc[clust,geneset_name] = pval_humk
        
        mmmk_corr_rho.loc[clust,geneset_name] = rho_mmmk
        mmmk_corr_pval.loc[clust,geneset_name] = pval_mmmk 
        
        mmhu_corr_rho.loc[clust,geneset_name] = rho_mmhu
        mmhu_corr_pval.loc[clust,geneset_name] = pval_mmhu 

# # PAUL custom
# using customer picked genes by Paul
# from paper: Paul Huang 2017 Cell

# # HGNC

# # prepare HGNC

# In[27]:


with open('./HGNC_zh.pckl','rb') as f:
    HGNC_genesets = pickle.load(f)
HGNC_df = pd.DataFrame(list(HGNC_genesets.items()),columns = ['Category','gene_list']) 
HGNC_df.set_index('Category', drop=True, append=False, inplace=True, verify_integrity=False)

HGNC_df['n_genes'] = ''

for cat in HGNC_df.index:
    n_gene = len(HGNC_df.loc[cat].gene_list)
    HGNC_df.loc[cat,'n_genes'] = n_gene
# # HGNC ion channels
all_cat = list(HGNC_df.index)
channel_cat = [s for s in all_cat if "channel" in s]
HT5_cat = [s for s in all_cat if "5-hydroxytryptamine" in s]

add_list = ['Anoctamins','Bestrophins']
for add_one in add_list:
    channel_cat.append(add_one)
len(channel_cat)

#channel_cat.remove('Acid sensing ion channel subunits') 
notwant_list = ['Acid sensing ion channel subunits','Cation channels sperm associated ',                'Volume regulated anion channel subunits','Transmembrane channel like family',               'Two pore segment channels','Zinc activated channels',]
for notwant in notwant_list:
    channel_cat.remove(notwant) 

channel_gene = list()
for cat in channel_cat:
    cat_gene = HGNC_df.loc[cat,'gene_list']
    #print(HGNC_df.loc[cat,'n_genes'])
    #print(cat)
    if type(cat_gene) == str:
        channel_gene.append(cat_gene)
    else:
        channel_gene = channel_gene + cat_gene
len(channel_gene)

serotonin_cat =['5-hydroxytryptamine receptors, G protein-coupled', '5-hydroxytryptamine receptors, ionotropic ']

serotonin_gene = list()
for cat in serotonin_cat:
    cat_gene = HGNC_df.loc[cat,'gene_list']
    #print(HGNC_df.loc[cat,'n_genes'])
    print(cat)
    if type(cat_gene) == str:
        serotonin_gene.append(cat_gene)
    else:
        serotonin_gene = serotonin_gene + cat_gene
len(serotonin_gene)

receptor_cat = [ 'Vasoactive intestinal peptide receptor family','Angiotensin receptors',                  'Peptide receptors','Arginine vasopressin and oxytocin receptors',                'Bradykinin receptors','Calcitonin receptors',                'Cholinergic receptors muscarinic', 'Cholinergic receptors nicotinic subunits',                'Chemerin receptor','Cannabinoid receptors','Dopamine receptors',                'Gamma-aminobutyric acid type B receptor subunits','Gamma-aminobutyric acid type A receptor subunits',                'Glycine receptors','Glutamate ionotropic receptor AMPA type subunits',                'Glutamate ionotropic receptor delta type subunits', 'Glutamate ionotropic receptor kainate type subunits',                'Glutamate ionotropic receptor NMDA type subunits','Glutamate metabotropic receptors',                '5-hydroxytryptamine receptors, G protein-coupled','5-hydroxytryptamine receptors, ionotropic ',                'Melanocortin receptors','Melatonin receptors','Neuropeptide receptors','Neurotensin receptors',               'Opsin receptors', 'Opioid receptors','Purinergic receptors P2X',                 'P2Y receptors','Somatostatin receptors',]

receptor_gene = list()
for cat in receptor_cat:
    cat_gene = HGNC_df.loc[cat,'gene_list']
    #print(HGNC_df.loc[cat,'n_genes'])
    print(cat)
    if type(cat_gene) == str:
        receptor_gene.append(cat_gene)
    else:
        receptor_gene = receptor_gene + cat_gene
len(receptor_gene)

GO_HGNC_dict = go_all_dict

GO_HGNC_dict['ion channels'] = channel_gene
GO_HGNC_dict['synapes receptor'] = receptor_gene
GO_HGNC_dict['serotonin receptors'] = serotonin_gene

GO_HGNC_dict.keys()

with open('data/python/GO_CrossCorr/Paul_custom_gene.pckl', 'rb') as f:   
    geneset_Paul = pickle.load(f)
geneset_merge = {**GO_HGNC_dict, **geneset_Paul}
# # create plots for figure

# ## scatter plots

cluster_key = 'GABAergic'
gene_set_key = 'nervous system development'

gene_list_all = list(go_all_dict[gene_set_key])

gene_list_overlap = list(set(gene_list_all) & set(list(res_mk.columns)))
    #len(gene_list_overlap) 

res_mk_sele = res_mk.loc[cluster_key,gene_list_overlap]
res_hu_sele = res_hu.loc[cluster_key,gene_list_overlap]
res_mm_sele = res_mm.loc[cluster_key,gene_list_overlap]

rho_humk, pval_humk = stats.spearmanr(res_mk_sele, res_hu_sele)
rho_mmmk, pval_mmmk = stats.spearmanr(res_mk_sele, res_mm_sele)
rho_mmhu, pval_mmhu = stats.spearmanr(res_hu_sele, res_mm_sele)

fig = plt.figure(figsize = (8,4))
gs = GridSpec(1, 2, figure=fig)
plt.tight_layout()

sns.set(font_scale=1)
ax1 = fig.add_subplot(gs[0, 0])
ax3 = fig.add_subplot(gs[0, 1])

ax1.scatter(res_mk_sele, res_hu_sele, color = '#3333ff', s = 5, alpha = 0.5)
ax3.scatter(res_mm_sele, res_hu_sele, color = '#000000', s = 5, alpha = 0.5)
    
ax1.text (-1.8,1.8,'rho = '+str(format(rho_humk,'.3f')), color = 'black')
ax1.text (-1.8,1.5,'p = '+str(format(pval_humk,'.1e')), color = 'black')

ax3.text (-1.8,1.8,'rho = '+str(format(rho_mmhu,'.3f')), color = 'black')
ax3.text (-1.8,1.5,'p = '+str(format(pval_mmhu,'.1e')), color = 'black')

ax1.title.set_text(f"{cluster_key}\n{gene_set_key}")
ax3.title.set_text(f"{cluster_key}\n{gene_set_key}")


x1 = np.array(res_mk_sele,dtype=float)
y1 = np.array(res_hu_sele,dtype=float)
z1 = numpy.polyfit(x1, y1, 1,)
p1 = numpy.poly1d(z1)
ax1.plot(x1,p1(x1),'-', color = "#3333ff")


dif = y1-p1(x1)
max10_ind = np.argpartition(dif, -3)[-3:]
max_dot_mk = res_mk_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax1.scatter(max_dot_mk,max_dot_hu,facecolors='none', edgecolors='#3333ff')
h1_txt = []
for gene in list(max_dot_hu.index):
    h1_txt.append(ax1.annotate(gene,(max_dot_mk[gene],max_dot_hu[gene]), size = 8))   

dif = y1-p1(x1)
max10_ind = np.argpartition(-dif, -3)[-3:]
max_dot_mk = res_mk_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]

ax1.scatter(max_dot_mk,max_dot_hu,facecolors='none', edgecolors='#3333ff')

h2_txt = []
for gene in list(max_dot_hu.index):
    h2_txt.append(ax1.annotate(gene,(max_dot_mk[gene],max_dot_hu[gene]), size = 8))   

x3 = np.array(res_mm_sele,dtype=float)
y3 = np.array(res_hu_sele,dtype=float)
z3 = numpy.polyfit(x3, y3, 1,)
p3 = numpy.poly1d(z3)
ax3.plot(x3,p3(x3),"k-")


dif = y3-p3(x3)
max10_ind = np.argpartition(dif, -3)[-3:]
max_dot_mm = res_mm_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax3.scatter(max_dot_mm,max_dot_hu,facecolors='none', edgecolors='k')
h3_txt = []
for gene in list(max_dot_hu.index):
    h3_txt.append(ax3.annotate(gene,(max_dot_mm[gene],max_dot_hu[gene]), size = 8))
    
dif = y3-p3(x3)
max10_ind = np.argpartition(-dif, -3)[-3:]
max_dot_mm = res_mm_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax3.scatter(max_dot_mm,max_dot_hu,facecolors='none', edgecolors='k')
h4_txt = []
for gene in list(max_dot_hu.index):
    h4_txt.append(ax3.annotate(gene,(max_dot_mm[gene],max_dot_hu[gene]), size = 8))
#adjust_text(h_txt)

ax1.spines['bottom'].set_color('#000000')
ax1.spines['left'].set_color('#000000') 
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')

ax3.spines['bottom'].set_color('#000000')
ax3.spines['left'].set_color('#000000') 
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')

ax1.set_xlabel(f'Monkey')
ax1.set_ylabel(f'Human')
#ax2.set_xlabel(f'Monkey {cluster_key}')
#ax2.set_ylabel(f'Mouse {cluster_key}')
ax3.set_xlabel(f'Mouse')
ax3.set_ylabel(f'Human')
ax1.set_xlim([-2,2])
ax1.set_ylim([-2,2])
ax1.set_facecolor('white')
ax3.set_facecolor('white')
plt.savefig(f'./neuron_{cluster_key}_{gene_set_key}_scatter.pdf',dpi=600)


cluster_key = 'Glutamatergic'
gene_set_key = 'nervous system development'

gene_list_all = list(go_all_dict[gene_set_key])

gene_list_overlap = list(set(gene_list_all) & set(list(res_mk.columns)))
    #len(gene_list_overlap) 

res_mk_sele = res_mk.loc[cluster_key,gene_list_overlap]
res_hu_sele = res_hu.loc[cluster_key,gene_list_overlap]
res_mm_sele = res_mm.loc[cluster_key,gene_list_overlap]

rho_humk, pval_humk = stats.spearmanr(res_mk_sele, res_hu_sele)
rho_mmmk, pval_mmmk = stats.spearmanr(res_mk_sele, res_mm_sele)
rho_mmhu, pval_mmhu = stats.spearmanr(res_hu_sele, res_mm_sele)

fig = plt.figure(figsize = (8,4))
gs = GridSpec(1, 2, figure=fig)
plt.tight_layout()

sns.set(font_scale=1)
ax1 = fig.add_subplot(gs[0, 0])
ax3 = fig.add_subplot(gs[0, 1])
ax1.scatter(res_mk_sele, res_hu_sele, color = '#3333ff', s = 5, alpha = 0.5)
ax3.scatter(res_mm_sele, res_hu_sele, color = '#000000', s = 5, alpha = 0.5)
ax1.text (-0.45,0.4,'rho = '+str(format(rho_humk,'.3f')), color = 'black')
ax1.text (-0.45,0.3,'p = '+str(format(pval_humk,'.1e')), color = 'black')
ax3.text (-0.45,0.4,'rho = '+str(format(rho_mmhu,'.3f')), color = 'black')
ax3.text (-0.45,0.3,'p = '+str(format(pval_mmhu,'.1e')), color = 'black')
ax1.title.set_text(f"{cluster_key}\n{gene_set_key}")
ax3.title.set_text(f"{cluster_key}\n{gene_set_key}")

x1 = np.array(res_mk_sele,dtype=float)
y1 = np.array(res_hu_sele,dtype=float)
z1 = numpy.polyfit(x1, y1, 1,)
p1 = numpy.poly1d(z1)
ax1.plot(x1,p1(x1),"-",color = "#3333ff")

x3 = np.array(res_mm_sele,dtype=float)
y3 = np.array(res_hu_sele,dtype=float)
z3 = numpy.polyfit(x3, y3, 1,)
p3 = numpy.poly1d(z3)
ax3.plot(x3,p3(x3),"k-")

ax1.spines['bottom'].set_color('#000000')
ax1.spines['left'].set_color('#000000') 
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')

ax3.spines['bottom'].set_color('#000000')
ax3.spines['left'].set_color('#000000') 
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')

ax1.set_xlabel(f'Monkey')
ax1.set_ylabel(f'Human')
ax3.set_xlabel(f'Mouse')

ax1.set_xlim([-1,1])
ax1.set_ylim([-1,1])

ax3.set_xlim([-1,1])
ax3.set_ylim([-1,1])
  
ax1.set_facecolor('white')
ax3.set_facecolor('white')

dif = y1-p1(x1)
max10_ind = np.argpartition(dif, -2)[-2:]
max_dot_mk = res_mk_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax1.scatter(max_dot_mk,max_dot_hu,facecolors='none', edgecolors='#3333ff')
for gene in list(max_dot_hu.index):
    ax1.annotate(gene,(max_dot_mk[gene],max_dot_hu[gene]), size = 8,)


dif = y1-p1(x1)
max10_ind = np.argpartition(-dif, -2)[-2:]
max_dot_mk = res_mk_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]

ax1.scatter(max_dot_mk,max_dot_hu,facecolors='none', edgecolors='#3333ff')
for gene in list(max_dot_hu.index):
    ax1.annotate(gene,(max_dot_mk[gene],max_dot_hu[gene]), size = 8)

dif = y3-p3(x3)
max10_ind = np.argpartition(dif, -3)[-3:]
max_dot_mm = res_mm_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax3.scatter(max_dot_mm,max_dot_hu,facecolors='none', edgecolors='k')
for gene in list(max_dot_hu.index):
    ax3.annotate(gene,(max_dot_mm[gene],max_dot_hu[gene]), size = 8)

dif = y3-p3(x3)
max10_ind = np.argpartition(-dif, -3)[-3:]
max_dot_mm = res_mm_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax3.scatter(max_dot_mm,max_dot_hu,facecolors='none', edgecolors='k')
for gene in list(max_dot_hu.index):
    ax3.annotate(gene,(max_dot_mm[gene],max_dot_hu[gene]), size = 8)    

potentiation_sets = ['long-term synaptic potentiation', 'regulation of long-term synaptic potentiation']

cluster_key = 'Glutamatergic'
gene_set_key = 'long_term_positive_pententiation_related'

go_all = []
for gene_set in potentiation_sets:
    go_all = go_all + list(go_all_dict[gene_set])
len(go_all)

gene_list_all = unique(go_all)
gene_list_overlap = list(set(gene_list_all ) & set(list(res_mk.columns)))

res_mk_sele = res_mk.loc[cluster_key,gene_list_overlap]
res_hu_sele = res_hu.loc[cluster_key,gene_list_overlap]
res_mm_sele = res_mm.loc[cluster_key,gene_list_overlap]

rho_humk, pval_humk = stats.spearmanr(res_mk_sele, res_hu_sele)
rho_mmmk, pval_mmmk = stats.spearmanr(res_mk_sele, res_mm_sele)
rho_mmhu, pval_mmhu = stats.spearmanr(res_hu_sele, res_mm_sele)

fig = plt.figure(figsize = (8,4))
gs = GridSpec(1, 2, figure=fig)

sns.set(font_scale=1)
ax1 = fig.add_subplot(gs[0, 0])
ax3 = fig.add_subplot(gs[0, 1])

ax1.scatter(res_mk_sele, res_hu_sele, color = '#3333ff', s = 5, alpha = 0.5)
ax3.scatter(res_mm_sele, res_hu_sele, color = '#000000', s = 5, alpha = 0.5)


ax1.text (-0.25,0.3,'rho = '+str(format(rho_humk,'.3f')), color = 'black')
ax1.text (-0.25,0.26,'p = '+str(format(pval_humk,'.1e')), color = 'black')

ax3.text (-0.45,0.3,'rho = '+str(format(rho_mmhu,'.3f')), color = 'black')
ax3.text (-0.45,0.26,'p = '+str(format(pval_mmhu,'.1e')), color = 'black')

ax1.title.set_text(f"{cluster_key}\n{gene_set_key}")
ax3.title.set_text(f"{cluster_key}\n{gene_set_key}")


x1 = np.array(res_mk_sele,dtype=float)
y1 = np.array(res_hu_sele,dtype=float)
z1 = numpy.polyfit(x1, y1, 1,)
p1 = numpy.poly1d(z1)
ax1.plot(x1,p1(x1),"-",color = "#3333ff")

x3 = np.array(res_mm_sele,dtype=float)
y3 = np.array(res_hu_sele,dtype=float)
z3 = numpy.polyfit(x3, y3, 1,)
p3 = numpy.poly1d(z3)
ax3.plot(x3,p3(x3),"k-")

ax1.spines['bottom'].set_color('#000000')
ax1.spines['left'].set_color('#000000') 
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')

ax3.spines['bottom'].set_color('#000000')
ax3.spines['left'].set_color('#000000') 
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')

ax1.set_xlabel(f'Monkey')
ax1.set_ylabel(f'Human')

ax3.set_xlabel(f'Mouse')
ax3.set_ylabel(f'Human')

ax1.set_facecolor('white')
ax3.set_facecolor('white')

dif = y1-p1(x1)
max10_ind = np.argpartition(dif, -3)[-3:]
max_dot_mk = res_mk_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax1.scatter(max_dot_mk,max_dot_hu,facecolors='none', edgecolors='#3333ff')
for gene in list(max_dot_hu.index):
    ax1.annotate(gene,(max_dot_mk[gene],max_dot_hu[gene]), size = 8,)


dif = y1-p1(x1)
max10_ind = np.argpartition(-dif, -3)[-3:]
max_dot_mk = res_mk_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]

ax1.scatter(max_dot_mk,max_dot_hu,facecolors='none', edgecolors='#3333ff')

for gene in list(max_dot_hu.index):
    ax1.annotate(gene,(max_dot_mk[gene],max_dot_hu[gene]), size = 8)   
    
dif = y3-p3(x3)
max10_ind = np.argpartition(dif, -3)[-3:]
max_dot_mm = res_mm_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax3.scatter(max_dot_mm,max_dot_hu,facecolors='none', edgecolors='k')
for gene in list(max_dot_hu.index):
    ax3.annotate(gene,(max_dot_mm[gene],max_dot_hu[gene]), size = 8)

dif = y3-p3(x3)
max10_ind = np.argpartition(-dif, -3)[-3:]
max_dot_mm = res_mm_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax3.scatter(max_dot_mm,max_dot_hu,facecolors='none', edgecolors='k')
for gene in list(max_dot_hu.index):
    ax3.annotate(gene,(max_dot_mm[gene],max_dot_hu[gene]), size = 8)    

plt.savefig(f'./neuron_{cluster_key}_{gene_set_key}_scatter.png',dpi=600)
plt.savefig(f'./neuron_{cluster_key}_{gene_set_key}_scatter.pdf',dpi=600)


potentiation_sets = ['long-term synaptic potentiation', 'regulation of long-term synaptic potentiation',]

cluster_key = 'GABAergic'
gene_set_key = 'long_term_positive_pententiation_related'

go_all = []
for gene_set in potentiation_sets:
    go_all = go_all + list(go_all_dict[gene_set])
len(go_all)

gene_list_all = unique(go_all)
gene_list_overlap = list(set(gene_list_all ) & set(list(res_mk.columns)))

res_mk_sele = res_mk.loc[cluster_key,gene_list_overlap]
res_hu_sele = res_hu.loc[cluster_key,gene_list_overlap]
res_mm_sele = res_mm.loc[cluster_key,gene_list_overlap]

rho_humk, pval_humk = stats.spearmanr(res_mk_sele, res_hu_sele)
rho_mmmk, pval_mmmk = stats.spearmanr(res_mk_sele, res_mm_sele)
rho_mmhu, pval_mmhu = stats.spearmanr(res_hu_sele, res_mm_sele)

fig = plt.figure(figsize = (8,4))
gs = GridSpec(1, 2, figure=fig)
plt.tight_layout()

sns.set(font_scale=1)
ax1 = fig.add_subplot(gs[0, 0])
ax3 = fig.add_subplot(gs[0, 1])

ax1.scatter(res_mk_sele, res_hu_sele, color = '#3333ff', s = 5, alpha = 0.5)
ax3.scatter(res_mm_sele, res_hu_sele, color = '#000000', s = 5, alpha = 0.5)


ax1.text (-0.7,2,'rho = '+str(format(rho_humk,'.3f')), color = 'black')
ax1.text (-0.7,1.5,'p = '+str(format(pval_humk,'.1e')), color = 'black')

ax3.text (-1,2,'rho = '+str(format(rho_mmhu,'.3f')), color = 'black')
ax3.text (-1,1.5,'p = '+str(format(pval_mmhu,'.1e')), color = 'black')

ax1.title.set_text(f"{cluster_key}\n{gene_set_key}")
#ax2.title.set_text(f"{cluster_key}\n{gene_set_key}")
ax3.title.set_text(f"{cluster_key}\n{gene_set_key}")

x1 = np.array(res_mk_sele,dtype=float)
y1 = np.array(res_hu_sele,dtype=float)
z1 = numpy.polyfit(x1, y1, 1,)
p1 = numpy.poly1d(z1)
ax1.plot(x1,p1(x1),"-",color = "#3333ff")

x3 = np.array(res_mm_sele,dtype=float)
y3 = np.array(res_hu_sele,dtype=float)
z3 = numpy.polyfit(x3, y3, 1,)
p3 = numpy.poly1d(z3)
ax3.plot(x3,p3(x3),"k-")

ax1.spines['bottom'].set_color('#000000')
ax1.spines['left'].set_color('#000000') 
ax1.spines['top'].set_color('none')
ax1.spines['right'].set_color('none')

ax3.spines['bottom'].set_color('#000000')
ax3.spines['left'].set_color('#000000') 
ax3.spines['top'].set_color('none')
ax3.spines['right'].set_color('none')

ax1.set_xlabel(f'Monkey')
ax1.set_ylabel(f'Human')

ax3.set_xlabel(f'Mouse')
ax3.set_ylabel(f'Human')

ax1.set_facecolor('white')
ax3.set_facecolor('white')


dif = y1-p1(x1)
max10_ind = np.argpartition(dif, -3)[-3:]
max_dot_mk = res_mk_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax1.scatter(max_dot_mk,max_dot_hu,facecolors='none', edgecolors='#3333ff')
for gene in list(max_dot_hu.index):
    ax1.annotate(gene,(max_dot_mk[gene],max_dot_hu[gene]), size = 8,)

dif = y1-p1(x1)
max10_ind = np.argpartition(-dif, -3)[-3:]
max_dot_mk = res_mk_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]

ax1.scatter(max_dot_mk,max_dot_hu,facecolors='none', edgecolors='#3333ff')

for gene in list(max_dot_hu.index):
    ax1.annotate(gene,(max_dot_mk[gene],max_dot_hu[gene]), size = 8)

dif = y3-p3(x3)
max10_ind = np.argpartition(dif, -3)[-3:]
max_dot_mm = res_mm_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax3.scatter(max_dot_mm,max_dot_hu,facecolors='none', edgecolors='k')
for gene in list(max_dot_hu.index):
    ax3.annotate(gene,(max_dot_mm[gene],max_dot_hu[gene]), size = 8)

dif = y3-p3(x3)
max10_ind = np.argpartition(-dif, -3)[-3:]
max_dot_mm = res_mm_sele[max10_ind]
max_dot_hu = res_hu_sele[max10_ind]
ax3.scatter(max_dot_mm,max_dot_hu,facecolors='none', edgecolors='k')
for gene in list(max_dot_hu.index):
    ax3.annotate(gene,(max_dot_mm[gene],max_dot_hu[gene]), size = 8)    

plt.savefig(f'./neuron_{cluster_key}_{gene_set_key}_scatter.png',dpi=600)
plt.savefig(f'./neuron_{cluster_key}_{gene_set_key}_scatter.pdf',dpi=600)

# # barplot

def compare_two_rho(rho1, rho2, n):

# from https://www.ibm.com/support/pages/differences-between-correlations
# * H0: R1 = R2; r1 & r2 are sample corr of x,y for groups 1 & 2 .
# * n1 and n2 are sample sizes for groups 1 and 2.
# compute z1 = .5*ln((1+r1)/(1-r1)).
# compute z2 = .5*ln((1+r2)/(1-r2)).
# compute sezdiff = sqrt(1.06/(n1 - 3) + 1.06/(n2-3)).
# compute ztest = (z1 - z2)/sezdiff.
# COMPUTE alpha = 2*(1 - cdf.normal(abs(ztest),0,1)).
   
#compare the two rhos from the spearsman's test
#rho1, rho2, cross-corr for rho1 and rho2
#n sample size
#compare the cross correlation rho
    import scipy
    
    z1 = .5*np.log((1+rho1)/(1-rho1))
    z2 = .5*log((1+rho2)/(1-rho2))
    sezdiff = sqrt(1.06/(n - 3) + 1.06/(n-3))
    ztest = (z1 - z2)/sezdiff
    alpha = 2*(1 - scipy.stats.norm.cdf(abs(ztest)))
    return (alpha)

key_list_hu = (res_hu.index)
key_list_mk = (res_mk.index)
key_list_mm = (res_mm.index)

cluster_list = list(set(key_list_hu).intersection(set(key_list_mk)).intersection(set(key_list_mm)))
geneset_list = list(geneset_merge.keys())
cluster_list

humk_HGNC_corr_rho = pd.DataFrame()
humk_HGNC_corr_pval = pd.DataFrame()
humk_HGNC_corr_p_dif = pd.DataFrame()

mmmk_HGNC_corr_rho = pd.DataFrame()
mmmk_HGNC_corr_pval = pd.DataFrame()
mmmk_HGNC_corr_p_dif = pd.DataFrame()

mmhu_HGNC_corr_rho = pd.DataFrame()
mmhu_HGNC_corr_pval = pd.DataFrame()
mmhu_HGNC_corr_p_dif = pd.DataFrame()

for clust in cluster_list:
    print(' ')
    print(clust)
    for geneset_name in geneset_list:
        gene_list_all = list(geneset_merge[geneset_name])
        gene_list_overlap = (list(set(gene_list_all) & set(all_gene)))
        print(geneset_name)
        print(len(gene_list_overlap))
        
        n_gene = len(gene_list_overlap)
        if n_gene <10: continue # skip the geneset with low gene numbers
        
        res_mk_sele = res_mk.loc[clust,gene_list_overlap]
        res_hu_sele = res_hu.loc[clust,gene_list_overlap]
        res_mm_sele = res_mm.loc[clust,gene_list_overlap]
        
        rho_humk, pval_humk = stats.spearmanr(res_mk_sele, res_hu_sele)
        rho_mmmk, pval_mmmk = stats.spearmanr(res_mk_sele, res_mm_sele)
        rho_mmhu, pval_mmhu = stats.spearmanr(res_mm_sele, res_hu_sele)
        
        humk_HGNC_corr_rho.loc[clust,geneset_name] = rho_humk
        humk_HGNC_corr_pval.loc[clust,geneset_name] = pval_humk        
        
        mmmk_HGNC_corr_rho.loc[clust,geneset_name] = rho_mmmk
        mmmk_HGNC_corr_pval.loc[clust,geneset_name] = pval_mmmk 
        
        mmhu_HGNC_corr_rho.loc[clust,geneset_name] = rho_mmhu
        mmhu_HGNC_corr_pval.loc[clust,geneset_name] = pval_mmhu 
        
        #compare the two correlations
        humk_HGNC_corr_p_dif.loc[clust,geneset_name] = compare_two_rho(rho_humk, rho_mmhu, n_gene)


# In[64]:


def Xcorr_perm(X,Y, n_perm):
#try to find the significance of Xcorr using permutation
#X, array of x
#Y, array of y
#n_perm, number of permutation
#returns all the rhos, ps and the p-value of the perm
    
    #calculate real rho
    rho_real, p_real  = stats.spearmanr(a,b)
    
    #permutation test
    import random
    shf_y = Y
    rho_perm = []
    p_perm = []
    for i in range(0,n_perm):
        random.shuffle(shf_y)
        rho_i, pval_i = stats.spearmanr(X, shf_y)
        rho_perm.append(rho_i)
        p_perm.append(pval_i)
    
    #calculate P    
    p_all = 1-sum(rho_perm<rho_real)/n_perm
    
    return(rho_perm, p_all)


# # plot only humk > humm terms

# In[34]:


gene_sets_picked = [
             'axon','dendrite','synapes receptor','ion channels', \
             'ionotropic glutamate receptor complex',\
             'AMPA glutamate receptor complex',
             'serotonin receptors','Neuropeptide receptors', \
             'calcium ion binding','All Trimeric G_proteins_alpha+beta+gamma',\
             'cadherin binding','Nrxn+Nlgn+Nphx+Dystroglycan+Cbln','Semaphorin+Plexin',         
             'nervous system development', 'long-term synaptic potentiation and its regulation']

cluster_key = 'GABAergic'

# width of the bars
#scale_factor = 0.5
barWidth = 0.3
gene_sets_all = list(humk_HGNC_corr_rho.keys())
    # Choose the height of the blue bars
gene_sets = gene_sets_picked
bars1 = humk_HGNC_corr_rho.loc[cluster_key,gene_sets] #*scale_factor
#bars2 = mmmk_corr_rho.loc[key,:]
bars3 = mmhu_HGNC_corr_rho.loc[cluster_key,gene_sets]#*scale_factor
p = humk_HGNC_corr_p_dif.loc[cluster_key,gene_sets] 

# The x position of bars
r1 = np.arange(len(bars1))
#r2 = [x + barWidth for x in r1]
r3 = [x + barWidth + 0 for x in r1]
    
fig = plt.figure(figsize = (8,3))
plt.rcParams['font.size'] = '14'
    
gs = GridSpec(1, 1, figure=fig)
sns.set(font_scale=1)
        
ax = fig.add_subplot(gs[0, 0])
    
ax.bar(r1, bars1, width = barWidth, color = '#3333ff', edgecolor = 'white',linewidth=0.5,label='HS vs FM')
ax.bar(r3, bars3, width = barWidth, color = '#000000', edgecolor = 'white',linewidth=0.5, label='HS vs MM')
 
ax.set_xticklabels(gene_sets,rotation=90)
ax.set_ylabel('rho')
ax.set_xticks(np.arange(0,len(r1))+barWidth)    
ax.legend(frameon=False)
 
plt.ylim([-0.2,1.0])
ax.set_xlim([-0.3,15])
# Show graphic
#plt.show()
ax.set_facecolor('white')
ax.set_title(f'{cluster_key} Cross Correlation of selected genesets')

ax.spines['bottom'].set_color('#000000')
ax.spines['left'].set_color('#000000') 
plt.savefig(f'./neuron_{cluster_key}_barplot_picked.pdf',dpi=600)


cluster_key = 'Glutamatergic'

# width of the bars
#scale_factor = 0.5
barWidth = 0.3
gene_sets_all = list(humk_HGNC_corr_rho.keys())
    # Choose the height of the blue bars
gene_sets = gene_sets_picked
bars1 = humk_HGNC_corr_rho.loc[cluster_key,gene_sets] #*scale_factor
bars3 = mmhu_HGNC_corr_rho.loc[cluster_key,gene_sets]#*scale_factor
p = humk_HGNC_corr_p_dif.loc[cluster_key,gene_sets] 

r1 = np.arange(len(bars1))
r3 = [x + barWidth + 0 for x in r1]

fig = plt.figure(figsize = (8,3))
plt.rcParams['font.size'] = '14'
    
gs = GridSpec(1, 1, figure=fig)
sns.set(font_scale=1)
     
ax = fig.add_subplot(gs[0, 0])
    
# Create blue bars
ax.bar(r1, bars1, width = barWidth, color = '#3333ff', edgecolor = 'white',linewidth=0.5,label='HS vs FM')
ax.bar(r3, bars3, width = barWidth, color = '#000000', edgecolor = 'white',linewidth=0.5, label='HS vs MM')

ax.set_xticklabels(gene_sets,rotation=90)
ax.set_ylabel('rho')
ax.set_xticks(np.arange(0,len(r1))+barWidth)    
ax.legend(frameon=False)
 
plt.ylim([-0.2,1.0])
ax.set_xlim([-0.3,15])
# Show graphic
#plt.show()
ax.set_facecolor('white')
ax.set_title(f'{cluster_key} Cross Correlation of selected genesets')

ax.spines['bottom'].set_color('#000000')
ax.spines['left'].set_color('#000000') 
plt.savefig(f'./neuron_{cluster_key}_barplot_picked.pdf',dpi=600)


# # Distribution

#plot 442 HGNC genesets from Paul paper
with open('data/python/GO_CrossCorr/Rho_distribution_Paul_442.pckl', 'rb') as f:
    geneset = pickle.load(f)
geneset.keys()

distr_humk_corr_rho = pd.DataFrame()
distr_humk_corr_pval = pd.DataFrame()

dristr_mmmk_corr_rho = pd.DataFrame()
dristr_mmmk_corr_pval = pd.DataFrame()

distr_mmhu_corr_rho = pd.DataFrame()
distr_mmhu_corr_pval = pd.DataFrame()

#prepare for loop
key_list_hu = (res_hu.index)
key_list_mk = (res_mk.index)
key_list_mm = (res_mm.index)

cluster_list = list(set(key_list_hu).intersection(set(key_list_mk)).intersection(set(key_list_mm)))
geneset_list = list(geneset.keys()) 
cluster_list

for clust in cluster_list:
    print(' ')
    print(clust)
    for geneset_name in geneset_list:
        gene_list_all = list(geneset[geneset_name])
        gene_list_overlap = (list(set(gene_list_all) & set(all_gene)))
        #print(geneset_name)
        #print(len(gene_list_overlap))
        
        if len(gene_list_overlap) <3: continue # skip the geneset with low gene numbers
            
        res_mk_sele = res_mk.loc[clust,gene_list_overlap]
        res_hu_sele = res_hu.loc[clust,gene_list_overlap]
        res_mm_sele = res_mm.loc[clust,gene_list_overlap]
        
        rho_distr_humk, pval_distr_humk = stats.spearmanr(res_mk_sele, res_hu_sele)
        rho_dristr_mmmk, pval_dristr_mmmk = stats.spearmanr(res_mk_sele, res_mm_sele)
        rho_distr_mmhu, pval_distr_mmhu = stats.spearmanr(res_mm_sele, res_hu_sele)
        
        distr_humk_corr_rho.loc[clust,geneset_name] = rho_distr_humk
        distr_humk_corr_pval.loc[clust,geneset_name] = pval_distr_humk
        
        dristr_mmmk_corr_rho.loc[clust,geneset_name] = rho_dristr_mmmk
        dristr_mmmk_corr_pval.loc[clust,geneset_name] = pval_dristr_mmmk 
        
        distr_mmhu_corr_rho.loc[clust,geneset_name] = rho_distr_mmhu
        distr_mmhu_corr_pval.loc[clust,geneset_name] = pval_distr_mmhu 

sns.distplot(distr_humk_corr_rho.T["GABAergic"], label="HS Vs FM Inh.", hist = False, #kde=True, 
             kde_kws={"color": "#3333ff", "lw": 2, 'linestyle':'-'})
sns.distplot(distr_humk_corr_rho.T["Glutamatergic"], label="HS Vs FM Exc.",hist = False,# kde=True,
            kde_kws={"color": "#3333ff", "lw": 2, 'linestyle':'--'})
sns.distplot(distr_mmhu_corr_rho.T["GABAergic"],  label="HS Vs MM Inh.", hist = False, #kde=True,
            kde_kws={"color": "#000000", "lw": 2, 'linestyle':'-'})
sns.distplot(distr_mmhu_corr_rho.T["Glutamatergic"], label="HS Vs MM Exc.",  hist = False, #kde=True,
            kde_kws={"color": "#000000", "lw": 2, 'linestyle':'--'})

plt.legend()
plt.xlabel('rho')
plt.ylabel('Kernel density estimation')
plt.savefig('./Density.pdf')

