import scanpy as sc
import sys
sys.path.append('./scctools')
from scctools import *


# # Load and check data
#monkey
ad=sc.read('data/Figure/FM27_cell_133454_wk.h5')
glut = ad[ad.obs['class']=='Exc',:]

#human
adhu=sc.read("./glut_huMK_500Marker1.h5")

#match .obs

adhu.obs['class'] = adhu.obs['class_label']
adhu.obs['subclass'] = adhu.obs['subclass_label']
adhu.obs['cell_label'] = adhu.obs['cluster_label']


# # load admk
admk = data_subset(glut,key = 'subclass', n=500, frac=1)

# # Matching dataset

adhu_m,admk_m =  match_gene(adhu, admk,ref_csv_path = "data/python/human2monkey.csv",
                Ln_col_name = 'Crab-eating macaque gene name',
                Ln_col_ID = 'Crab-eating macaque gene stable ID',
                Tr_species = 'Human',Ln_species = 'Monkey')

adhu_m.obs['subclass_label'].value_counts()


# # select HVG and merge data

adTr = adhu_m.copy()
adLn = admk_m.copy()
n_HVG = 2000
Ln_species = 'Monkey'
Tr_species = 'Human'

marker1 = ['FEZF2','THEMIS','RORB','HPCAL1']
marker2_mk = list(unique(admk.obs['marker2']))
marker2_hu = list(unique(adhu_m.obs['cell_label'].str.split(' ').str[3]))

markerGenes = unique(marker1 + marker2_mk + marker2_hu)
markerGenes = list(set(list(admk.var_names)).intersection(set(list(markerGenes))))

adLn.var[Ln_species] = adLn.var_names
adLn.var_names = adLn.var[Tr_species]

if scipy.sparse.issparse(adTr.layers['raw']):
        adTr.X = adTr.layers['raw'].todense().copy()
else: adTr.X = adTr.layers['raw'].copy()     
        
if scipy.sparse.issparse(adLn.layers['raw']):
        adLn.X = adLn.layers['raw'].todense().copy()
else: adLn.X = adLn.layers['raw'].copy()
    
adTr_selected_gene = adTr.var_names[np.array(np.std(adTr.X, axis=0).argsort())[0][-n_HVG:][::-1]]
adLn_selected_gene = adLn.var_names[np.array(np.std(adLn.X, axis=0).argsort())[0][-n_HVG:][::-1]]

select_gene  = unique(list(adTr_selected_gene)+list(adLn_selected_gene)+list(markerGenes))

select_gene = list(set(list(select_gene)).intersection(set(list(adTr.var_names))))
select_gene = list(set(list(select_gene)).intersection(set(list(adLn.var_names))))

adTr_sele = adTr[:,select_gene]
adTr_sele.obs['source']= Tr_species

adLn_sele = adLn[:,select_gene]
adLn_sele.obs['source']= Ln_species

TrLn = adTr_sele.concatenate(adLn_sele)
TrLn

humk = TrLn

sc.pp.filter_cells(humk, min_genes=10)
sc.pp.filter_genes(humk, min_cells=2) 


humk.obs['n_counts'] = np.matrix(humk.layers['raw'].sum(axis=1)).A1

def  find_heatmap_para(humk,theta,rep, louvain, savepath):
    #humk = run_PCA_Harmony(humk, run_Harmony = True, theta = 10,rep = 8)
    
    plot_matrix = two_species_heatmap(humk, species_1 = 'Monkey', species_2 = 'Human',species_1_key = 'cell_label', species_2_key = 'cell_label',                        louvain =louvain,figure_path = 'test_heatmap.png')
   
    print(f"theta = {theta}, rep = {rep}, louvain = {louvain}")
    
    import matplotlib.colors as mc

    def NonLinCdict(steps, hexcol_array):
        cdict = {'red': (), 'green': (), 'blue': ()}
        for s, hexcol in zip(steps, hexcol_array):
            rgb =matplotlib.colors.hex2color(hexcol)
            cdict['red'] = cdict['red'] + ((s, rgb[0], rgb[0]),)
            cdict['green'] = cdict['green'] + ((s, rgb[1], rgb[1]),)
            cdict['blue'] = cdict['blue'] + ((s, rgb[2], rgb[2]),)
        return cdict

        #hc = ['#e5e5ff', '#acacdf', '#7272bf', '#39399f', '#000080']
    hc = ['#ffffff', '#cccccc', '#aaaaaa', '#666666', '#000000']
    print('here')
    th = [0, 0.4, 0.45,0.5, 1]
    #th = [0, 0.2, 0.3,0.5, 1]

    cdict = NonLinCdict(th, hc)
    cm = mc.LinearSegmentedColormap('test', cdict)

    #figsize(12,6)
    graph = sns.heatmap(plot_matrix, cmap=cm, cbar=True, xticklabels=1,yticklabels=1, linewidth = 0.004, linecolor = 'gray')
    fig = graph.get_figure()
    save_name = f'{savepath}_theta{theta}_rep{rep}_louvain{louvain}'
    print(save_name)
    fig.savefig(save_name)

humk.obs['species'] = humk.obs['source']

humk = run_PCA_Harmony(humk, run_Harmony = True, theta = 8,rep = 4)


def two_species_heatmap(ad, species_1 = 'FM', species_2 = 'HS',species_1_key = 'subclass', species_2_key = 'marker2',                        louvain = 0,figure_path = 'test_heatmap.png'):
#generate hodge fig 5d heatmap
#input: ad: pre_merged, harmony corrected two species data, with .obs['species'] mark the species
#input : species_key, which .observe to generate the heatmap
#input: louvain, the louvain resolution for cluster the merged data
#input: figure_path, the path to save the heatplot, if figure_path = 0, do not save figure

    #prepare data
    if not louvain ==0:
        sc.pp.neighbors(ad,  metric='euclidean',use_rep = 'X_harmony' )
        sc.tl.louvain(ad, resolution = louvain,key_added = 'louvain')

    sc.pl.umap(ad, color = ['species','louvain'])
    
    ad_1 = ad[ad.obs['species'].isin([species_1]),:] 
    sc.pl.umap(ad_1, color = [species_1_key,'louvain'], legend_loc = 'on data')
    
    ad_2 = ad[ad.obs['species'].isin([species_2]),:]
    sc.pl.umap(ad_2, color = [species_2_key,'louvain'], legend_loc = 'on data')
    
    df_1 = pd.crosstab(ad_1.obs[species_1_key], ad_1.obs['louvain'], normalize ='index')
    df_2 = pd.crosstab(ad_2.obs[species_2_key], ad_2.obs['louvain'], normalize ='index')


    #generate heatmap matrix
    mk_cluster_all = df_1.index
    hu_cluster_all = df_2.index 
    low_sum_matrix = pd.DataFrame(index=mk_cluster_all, columns = hu_cluster_all)

    low_sum_matrix.index.name = species_1 +'_cluster'
    low_sum_matrix.columns.name = species_2 +'_cluster'
    
    
    low_sum_matrix = pd.DataFrame(index=mk_cluster_all, columns = hu_cluster_all)
    for mk_cluster in mk_cluster_all:
        for hu_cluster in hu_cluster_all:
        
            two_row = np.column_stack([df_1.loc[mk_cluster],df_2.loc[hu_cluster]])
            low_sum = sum(two_row.min(axis=1))
            low_sum_matrix.loc[mk_cluster, hu_cluster] = low_sum
            
    low_sum_matrix.to_csv('data/python/Homology_map/temp.csv')
    del(low_sum_matrix)
    low_sum_matrix =pd.read_csv('data/python/Homology_map/temp.csv')
    low_sum_matrix.set_index((species_1+'_cluster'), inplace = True)
    # for some reason, direct use matrix does not work.... save and read did the trick...
    
    low_sum_matrix_sort = low_sum_matrix.copy()

    ss = sns.clustermap(low_sum_matrix_sort, cmap = 'Greys',  cbar = True, col_cluster = False,  xticklabels=1, yticklabels=1,  method =  'single')
    
    hu_cluster_all = low_sum_matrix_sort.columns
    low_sum_matrix_sort1 = low_sum_matrix_sort
    low_sum_matrix_sort2 = low_sum_matrix_sort1.iloc[ss.dendrogram_row.reordered_ind]  
    
    low_sum_matrix1 =  low_sum_matrix_sort2
    for hu_cluster in hu_cluster_all:
        line = list(low_sum_matrix1.loc[:,hu_cluster])
        ss = np.array(line)
        tmp = []
        for kk in  range(len(ss)):
            if kk < len(ss) - 5:
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3+ss[kk+3]/4+ss[kk+4]/5)
            elif kk < len(ss) - 4: 
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3+ss[kk+3]/4)
            elif kk < len(ss) - 3: 
                tmp.append( ss[kk]+ss[kk+1]/2+ss[kk+2]/3 )
            elif kk < len(ss) - 2: 
                tmp.append( ss[kk]+ss[kk+1]/2 )
            elif kk < len(ss) - 1: 
                tmp.append( ss[kk] )
            elif kk < len(ss): 
                tmp.append(ss[kk]+ss[kk]/2+ss[kk]/3+ss[kk]/4+ss[kk]/5+ss[kk]/6)
        low_sum_matrix_sort2.loc['max_order',hu_cluster]= float(tmp.index(max(tmp)))
        del tmp
    low_sum_matrix_sort2 = low_sum_matrix_sort2.sort_values(by = ['max_order'],axis=1)
    plt.figure(figsize = (15,15))
    graph = sns.heatmap(low_sum_matrix_sort2[0:-1], cmap="Greys", cbar=True, xticklabels=1,yticklabels=1, linewidth = 0.01, linecolor = 'gray')
    
    plot_matrix = low_sum_matrix_sort2[0:-1]
    
    if not figure_path == 0:
        fig = graph.get_figure()
        fig.savefig(figure_path)
        
    return(plot_matrix)    

plot_matrix = two_species_heatmap(humk, species_1 = 'Monkey', species_2 = 'Human',species_1_key = 'cell_label', species_2_key = 'cell_label',                        louvain =4,figure_path = 'test_heatmap.png')


# # plot heatmap

figsize(6,5)
plot_matrix_reorder = plot_matrix.copy()
plot_matrix_reorder = plot_matrix_reorder.reindex(['Exc L1-3 HPCAL1 GPR83', 'Exc L1-3 HPCAL1 CBLN4','Exc L1-3 HPCAL1 WFDC2',                                          'Exc L1-3 HPCAL1 NPY',                                          'Exc L1-3 RORB SMARCA2','Exc L1-3 RORB CACNA1G','Exc L1-3 RORB VAMP1',                                          'Exc L1-6 RORB SMYD2','Exc L1-6 RORB FAM198B','Exc L4-6 RORB OSTN',                                          'Exc L4-6 RORB CD63','Exc L4-6 RORB CNTNAP5',                                          'Exc L4-6 RORB TNK2','Exc L4-6 RORB ARC','Exc L4-6 RORB PHLDB2','Exc L4-6 RORB SNCB',                                           'Exc L4-6 FEZF2 SLC24A2', 'Exc L4-6 FEZF2 SEC11C',                                          'Exc L4-6 FEZF2 BTBD11','Exc L4-6 FEZF2 MKX',                                            'Exc L4-6 FEZF2 SYT6','Exc L4-6 FEZF2 HTR2C','Exc L4-6 HPCAL1 NXPH4',
                                          'Exc L4-6 THEMIS CCK',
                                          'Exc L4-6 THEMIS KRT17',
                                          
                                          ])

plot_matrix_reorder = plot_matrix_reorder[['Exc L3 LINC00507 PSRC1','Exc L2-3 LINC00507 RPL9P17',                                           'Exc L3-4 RORB SEMA6D', 'Exc L3 RORB CARTPT','Exc L2-4 RORB GRIK1','Exc L4 RORB BHLHE22',                                           'Exc L4-5 RORB ASCL1','Exc L4-5 RORB AIM2', 'Exc L4 RORB CACNG5',                                           'Exc L6 FEZF2 KRT17', 'Exc L6 FEZF2 TBC1D26','Exc L6 FEZF2 VWA2',                                           'Exc L6 FEZF2 FAM95C', 'Exc L5-6 FEZF2 ANKRD20A1', 'Exc L6 FEZF2 CPZ',
                                            'Exc L4-5 RORB RPL31P31',\
                                           'Exc L4-5 RORB LCN15','Exc L4-5 RORB HNRNPA1P46', 'Exc L5 RORB SNHG7', \
                                          'Exc L6 FEZF2 TBCC', 'Exc L4-6 RORB HPCA','Exc L5-6 FEZF2 MYBPHL','Exc L6 FEZF2 P4HA3',\
                                           'Exc L6 FEZF2 SLITRK6','Exc L5-6 THEMIS TMEM233','Exc L5-6 THEMIS GPR21'\
                                          ]]

graph = sns.heatmap(plot_matrix_reorder, cmap='greys', cbar=True, xticklabels=1,yticklabels=1, linewidth = 0.2, linecolor = 'gray')

fig = graph.get_figure()
fig.savefig('./Fig_glut_heatmap.pdf')


# # plot umaps

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(constrained_layout=True,figsize=(5,4.5))

gs = GridSpec(1, 1, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])

ax1 = sc.pl.umap(humk, color=['species'],ax=ax1,           title='Exc. Neurons', show=False, size=20)

ax2.spines['bottom'].set_color('#000000')
ax2.spines['left'].set_color('#000000') 

ax2.spines['top'].set_color('#FFFFFF')
ax2.spines['right'].set_color('#FFFFFF') 

ax2.set_facecolor('white')

ax2.xaxis.label.set_size(12)
ax2.yaxis.label.set_size(12)

plt.savefig("./Fig_heatmap_umap_glut_harmony.pdf")


glut_mk = humk[humk.obs['species']=='Monkey',:]
glut_mk.obs['subclass_new'] = glut_mk.obs['cell_label'].str.split(' ').str[2]    
glut_hu = humk[humk.obs['species']=='Human',:]
glut_hu.obs['subclass_new'] = glut_hu.obs['cell_type_alias_label'].str.split(' ').str[2]    

import matplotlib.patheffects as pe

matplotlib.rcParams.update({'font.size': 10})
fig = plt.figure(constrained_layout=True,figsize=(6,3))

gs = GridSpec(1, 2, figure=fig)

ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])

ax1 = sc.pl.umap(glut_mk, color=['subclass_new'],ax=ax1,legend_fontoutline = 3,           title='Macaque', show=False, size=20, legend_loc = 'on data', frameon = False)

ax2= sc.pl.umap(glut_hu, color=['subclass_new'],ax=ax2,legend_fontoutline = 3,           title='Human', show=False, size=20, legend_loc = 'on data', frameon = False)

plt.savefig("./Fig_heatmap_umap_glut_subclass_split.pdf")

