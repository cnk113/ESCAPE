def loom2h5ad():
    import scanpy as sc
    import loompy
    adata = sc.read_loom('seuratObj.loom')
    adata.write('seuratObj.h5ad')
