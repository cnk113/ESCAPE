# Run CellBender
CellBender = function(h5, cells, total, calls, dim=200, layers=1000, epochs=300){
  command = paste("cellbender remove-background --input", eval(h5),
                 "--output output.h5",
                 "--cuda",
                 "--expected_cells", eval(cells),
                 "--total-droplets-included", eval(total),
                 "--z-dim", eval(dim)
                 "--z-layers", eval(layers)
                 "--epochs", eval(epochs))
  count <- Read10X_h5('output.h5')
  count
}

# Remove ambient RNA, empty droplets, and runs doublet removals tools
preprocess = function(h5, exp.cells, total.cells, cell.calls, z.dim, z.layers, epochs, fdr=0.01){
  count <- CellBender(h5, cells=exp.cells, total=total.cells, calls=cell.calls, dim=z.dim, layers=z.layers, epochs=epochs)
  seuratObj <- QuickCB2(h5file=h5, FDR_threshold=fdr, AsSeurat=T)
  consensus.barcodes <- intersect(colnames(count), colnames(seuratObj))
  count <- count[,consensus.barcodes]
  source('detectDoublets.py') 
  pred = detectDoublets(count, exp.cells, cell.calls) # Runs scrublet and doubletdetection 
  pred[[3]] <- DoubletDecon(count)
  pred[[4]] <- DoubletFinder(count, prop)
  consensus
  count[,consensus]
}


# Run DoubletFinder
DoubletFinder = function(count, prop=0.1){
  seuratObject = CreateSeuratObject(count)
  seuratObject = SCTransform(seuratObject)
  seuratObject = RunPCA(seuratObject)
  seuratObject = RunUMAP(seuratObject, dims=1:30)
  sweep = paramSweep_v3(seuratObject, PCs=1:30, sct=F)
  sweep.stats = summarizeSweep(sweep, GT=F)
  sweep.pk = find.pK(sweep.stats)
  homotypic = modelHomotypic(seuratObject@meta.data$seurat.clusters)
  nExp_poi <- round(prop*length(colnames(seuratObject)))
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
  seuratObject <- doubletFinder_v3(seuratObject, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  seuratObject <- doubletFinder_v3(seuratObject, PCs = 1:30, pN = 0.25, pK = 0.09, nExp = nExp_poi.adj, reuse.pANN = "pANN_0.25_0.09_913", sct = FALSE)
}



# Run DoubletDecon
DoubletDecon = function(count){
  seuratObject <- CreateSeuratObject(count)
  seuratObject <- SCTransform(seuratObject)
  seuratObject <- RunPCA(seuratObject)
  seuratObject <- FindNeighbors(seuratObject, dims=1:30)
  seuratObject <- FindClusters(seuratObject, algorithm=4)
  newFiles = Improved_Seurat_Pre_Process(seuratObject, num_genes=50, write_files=FALSE)
  results = Main_Doublet_Decon(rawDataFile=newFiles$newExpressionFile,
                             groupsFile=newFiles$newGroupsFile,
                             filename=filename,
                             location=location,
                             fullDataFile=NULL,
                             removeCC=FALSE,
                             species="hsa",
                             rhop=1.1,
                             write=TRUE,
                             PMF=TRUE,
                             useFull=FALSE,
                             heatmap=FALSE,
                             centroids=TRUE,
                             num_doubs=100,
                             only50=FALSE,
                             min_uniq=4)

