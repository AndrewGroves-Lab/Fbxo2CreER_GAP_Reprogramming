```python
# Run this in a R kernel

# Step 1: Convert data from Seurat to Python / anndata

library(Seurat)
library(ggplot2)
library(Matrix)



# save metadata table:
seurat_obj  <- readRDS(file = "fbxo2_filtered_analysis_Ishwar.RDS")

seurat_obj$barcode <- colnames(seurat_obj)
colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]

write.csv(seurat_obj@meta.data, file='metadata.csv', quote=F, row.names=F)
colnames(seurat_obj@meta.data)

Idents(seurat_obj)  <- seurat_obj$`cluster_names`


# write expression counts matrix
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0('counts.mtx'))
dim(counts_matrix)

# write dimesnionality reduction matrix, in this example case pca matrix
write.csv(seurat_obj@reductions$pca@cell.embeddings, file = "pca.csv", quote = F, row.names = F)# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file='gene_names.csv',
  quote=F,row.names=F,col.names=F
)




```


```python
# Continue with the next steps in a python kernel
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd

# load sparse matrix:
X = io.mmread("scVelo/counts.mtx")
X

# load cell metadata:
cell_meta = pd.read_csv("scVelo/metadata.csv")

# load gene names:
with open("scVelo/gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()
    
gene_names[1:10] #sanity check

# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names

# load dimensional reduction:
pca = pd.read_csv("scVelo/pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T


# plot a UMAP colored by sampleID to test:
sc.pl.umap(adata, color=["seurat_clusters"], frameon=False, save=True, alpha = 0.5)
adata

# save dataset as anndata format
adata.write('/media/ishwar/DATA1/GrovesLab/Melissa_collab/Fbxo2CreER/new_fbxo2_analysis/scVelo/fbxo2_pc25_res15.h5ad')

```


```python
# Constructing and setting up the spliced/unspliced matrices
# Load libraries

import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad
import loompy as lp

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
#cr.settings.verbosity = 2

adata.X #check

gene_name = adata.var
# gene_name = gene_name.to_numpy()
type(gene_name)
# print(gene_name)
# gene_list = gene_name.tolist()

# print(gene_list)
gene_list = gene_name.index
gene_list = gene_list.to_list()
type(gene_list)
print(gene_list[:10])


#load the spliced and unspliced matrices
exon = pd.read_table("/Exon_UMI_matrix_032222_withSox2.txt")
exon = exon.drop_duplicates(subset=["gene_symbol"])
exon = exon.loc[exon["gene_symbol"].isin(gene_list)]
spliced = exon.iloc[:, 2:343]
spliced = spliced.transpose()


intron = pd.read_table("/Intron_UMI_matrix_051722.txt")
intron = intron.drop_duplicates(subset=["gene_symbol"])
intron = intron.loc[intron["gene_symbol"].isin(gene_list)]
unspliced = intron.iloc[:,2:343]
unspliced = unspliced.transpose()

#Add the splcied and unpspliced layers to the anndata object
adata.layers["spliced"] = spliced
adata.layers["unspliced"] = unspliced


sc.pl.umap(adata, color= "cell_type",alpha = 0.7)

# Computing RNA velocity usign scVelo
scv.pl.proportions(adata, groupby="cell_type", save = "/scVelo/fbx_fil_grp_clust_norm.png", dpi = 300 )

#preprocessing steps
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

#computing velocity
#using the dynamic model
scv.tl.recover_dynamics(adata, n_jobs = 8 )
scv.tl.velocity(adata, mode = "dynamical")
scv.tl.velocity_graph(adata)
scv.tl.recover_latent_time(adata)


#save the dynamical modeling results for quick loading later
adata.write('data/fbxo2_dyn_mdl.h5ad', compression='gzip')
#adata = scv.read('data/fbxo2_dyn_mdl.h5ad')


# Plot the velocity stream
scv.pl.velocity_embedding_stream(adata, basis='umap', color=["cluster_names"], 
                                 palette={'Ctrl Cells 1': 'chartreuse3', 'Ctrl Cells 2': 'red', 'Ctrl Cells 3': 'skyblue', 'Ctrl Supp Cells': 'royalblue', 'Rprg HC': '#DA70D6', 'Rprg SC New': 'darkgreen', 'Rprg SC Orig': 'gold', 'Unknwn 1':'gray', 'Unknwn 2':'gray'}, 
                                 min_mass=4, frameon = False, legend_loc = False, alpha = 0.7, smooth = 1, arrow_style= "-|>",
                                  figsize=(20,15), linewidth=4, arrow_size = 7, integration_direction = "both", size=4000, save="fbxo2_stream_vel_dyn_mdl_no_labels_both_direction.tiff", dpi=300, )


scv.pl.velocity_embedding(adata, arrow_length=2, arrow_size=2, dpi=300, color="cluster_names")              


scv.pl.velocity_embedding_stream(adata, basis="umap", color="latent_time",frameon=False ,alpha=0.8,arrow_size=5, color_map="plasma", save='scVelo-umap-latent_time.png', figsize=(20,15), linewidth=3)


sc.pl.violin(adata, keys='latent_time',groupby="cluster_names", color = "cluster_names" ,save='scVelo-violin-latent_time.png', dpi = 300)

#Top-likelihood genes
#Driver genes display pronounced dynamic behavior and are systematically detected via their characterization by high likelihoods in the dynamic model.
kwargs = dict(cbar_pos=(0, 0.58, .03, .4))
top_genes= adata.var['fit_likelihood'].sort_values(ascending=False).index[:200]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time',col_color=['cell_type'], n_convolve=100, row_cluster=True, colorbar=True, save = "likelihood_genes_heatmap.png", col_cluster=False, figsize=(11,17), **kwargs)
print(top_genes)


# Downstream analysis
# Identify cluster specific highly dynamic genes
scv.tl.rank_velocity_genes(adata, groupby='cluster_names', min_corr=.3)
df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()

scv.tl.rank_dynamical_genes(adata, groupby='cluster_names')
df1 = scv.get_df(adata, 'rank_dynamical_genes/names')
df1.head(5)


```
