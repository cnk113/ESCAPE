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
  cellBenderCount <- CellBender(h5, cells=exp.cells, total=total.cells, calls=cell.calls, dim=z.dim, layers=z.layers, epochs=epochs)
  cb2Count <- QuickCB2(h5, FDR_threshold=fdr, AsSeurat=T) # Run CB2
  seuratObj <- Read10X_h5(h5)
  emptyDropsCount <- emptyDrops(seuratObj@assays$RNA@counts) # Run EmptyDrops
  calls <- list(colnames(cellBenderCount), colnames(cb2Count), colnames(emptyDropsCount))
  consensus.barcodes <- Reduce(intersect, calls)
  seuratObj <- subset(seuratObj, cells=consensus.barcodes)
  as.loom(seuratObj, filename="seuratObj.loom")
  source('loom2h5ad.py')
  loom2h5ad()


  source('detectDoublets.py') 
  pred = detectDoublets(seuratObj@assays$RNA@counts, exp.cells, cell.calls) # Runs scrublet and doubletdetection 
  pred[[3]] <- DoubletDecon(seuratObj)
  pred[[4]] <- DoubletFinder(seuratObj, prop)
  pred[[5]] <- Solo()
  consensus
  count[,consensus]
}


# Run Solo
Solo = function(h5){
  command = paste("cellbender remove-background --input", eval(h5),
		  "--output output.h5",)
  
}


# Run DoubletFinder
DoubletFinder = function(seuratObject, prop=0.1){
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
DoubletDecon = function(seuratObject){
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
