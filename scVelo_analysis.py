import anndata as ann
import os 
import scvelo as scv
import pandas as pd
import scipy as scp
import numpy as np
import time as time

os.chdir("D:\\rstudio\\Opus\\redo_pilot")

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.settings.set_figure_params('scvelo')  # for beautified visualization

time.sleep(600)

d1 = scv.read("replicate1.loom", cache=True)
d2 = scv.read("replicate2.loom", cache=True)
d3 = scv.read("replicate3.loom", cache=True)

d1.obs = pd.DataFrame(d1.obs.index.values, columns=["cell_id"])
d2.obs = pd.DataFrame(d2.obs.index.values, columns=["cell_id"])
d3.obs = pd.DataFrame(d3.obs.index.values, columns=["cell_id"])


d1.var = pd.DataFrame(d1.var.index.values, columns=["gene_id"])
d2.var = pd.DataFrame(d2.var.index.values, columns=["gene_id"])
d3.var = pd.DataFrame(d3.var.index.values, columns=["gene_id"])


df1 = pd.DataFrame(d1.layers['spliced'].toarray())
df2 = pd.DataFrame(d2.layers['spliced'].toarray())
df3 = pd.DataFrame(d3.layers['spliced'].toarray())

spliced = pd.concat([df1,df2,df3])

df1 = pd.DataFrame(d1.layers['unspliced'].toarray())
df2 = pd.DataFrame(d2.layers['unspliced'].toarray())
df3 = pd.DataFrame(d3.layers['unspliced'].toarray())

unspliced = pd.concat([df1,df2,df3])

#df1 = pd.DataFrame(d1.layers['matrix'].toarray())
#df2 = pd.DataFrame(d2.layers['matrix'].toarray())
#df3 = pd.DataFrame(d3.layers['matrix'].toarray())

#my_matrix = pd.concat([df1,df2,df3])



obsnum = np.array(pd.read_csv("obsnum.csv")).flatten() - 1 #because python iteraters from 0
varnum = np.array(pd.read_csv("varnum.csv")).flatten() - 1 
myumap = pd.read_csv("umap.csv")
vargene = pd.read_csv("vargene.csv")
obsID = pd.read_csv("obsID.csv")
clusters = pd.read_csv("labels.csv").iloc[:,1]
topvargenes= pd.read_csv("topvargenes.csv")


my_matrix = pd.read_csv("countsarray.csv", index_col=0).transpose()
my_matrix.shape
spliced = spliced.iloc[obsnum,varnum]
unspliced = unspliced.iloc[obsnum,varnum]
#my_matrix = my_matrix.iloc[obsnum,varnum]

#my_matrix = my_matrix.iloc[:,np.array(topvargenes).flatten()]
#spliced = spliced.iloc[:,np.array(topvargenes).flatten()]
#unspliced = unspliced.iloc[:,np.array(topvargenes).flatten()]
#vargene = vargene.iloc[np.array(topvargenes).flatten()]

my_matrix = scp.sparse.csr_matrix(my_matrix.values)
spliced = scp.sparse.csr_matrix(spliced.values)
unspliced = scp.sparse.csr_matrix(unspliced.values)


alldata = ann.AnnData(my_matrix)
alldata.layers['spliced'] = spliced
alldata.layers['unspliced'] = unspliced
alldata.obs = pd.concat([obsID,clusters], axis=1).rename(columns={'obsID':'CellID','myelo.labels':'clusters'})
alldata.var = vargene
scv.pp.filter_and_normalize(alldata, min_shared_counts=20)
scv.pp.filter_genes(alldata, min_shared_cells=20)
#scv.pp.normalize_per_cell(alldata)
#scv.pp.filter_genes_dispersion(alldata)
#scv.pp.log1p(alldata)
####
scv.pp.moments(alldata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(alldata)
scv.tl.velocity_graph(alldata)
alldata.obsm['X_wnnUMAP'] = np.array(myumap)
scv.pl.velocity_embedding_stream(alldata, title='', n_neighbors=40, arrow_size=0.7, smooth=0.8, min_mass=3, max_length=50, basis='wnnUMAP', color='clusters', size=120, legend_fontsize=6)
scv.pl.velocity_embedding(alldata, basis='wnnUMAP', color='clusters', arrow_length=4, arrow_size=3, dpi=150)
scv.tl.recover_dynamics(alldata)
scv.tl.velocity(alldata, mode='dynamical')
scv.tl.velocity_graph(alldata)
scv.tl.latent_time(alldata, vkey='velocity')
scv.pl.scatter(alldata, basis='wnnUMAP', color='latent_time', color_map='gnuplot', size=80, colorbar=True)

scv.tl.recover_dynamics(alldata)
top_genes = alldata.var['fit_likelihood'].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(alldata, var_names=top_genes, groupby='labels')
scv.tl.velocity(alldata, diff_kinetics=True)
scv.tl.velocity_graph(alldata)
scv.pl.velocity_embedding_stream(alldata, title='', n_neighbors=30, arrow_size=0.7, smooth=0.8, min_mass=3, max_length=50, basis='wnnUMAP', color='labels', size=120, legend_fontsize=6)
scv.tl.latent_time(alldata)
scv.pl.scatter(alldata, basis='wnnUMAP', color='latent_time', color_map='gnuplot', size=80, colorbar=True)
scv.pl.velocity_embedding(alldata, basis='wnnUMAP', color='labels', arrow_length=4, arrow_size=3, dpi=150)


scv.pl.scatter(alldata, basis=pd.DataFrame(alldata.var.Gene[np.array(top_genes[:15])]), ncols=5, add_outline='fit_diff_kinetics')

