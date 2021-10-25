################# EPO_analysis helper methods ############
# This script contains helper methods for the analysis
# of our 10X experiment on HSCs treated with EPO
# @author: Jason Cosgrove (jason.cosgrove@curie.fr)
# @date:   02/10/19
##########################################################


#returns a list of booleans stating TRUE if a given variable is NOT IN the list
"%ni%" <- Negate("%in%")


# do PCA diagnostic plots to see how well our normalisation is working
# i.e. compare normalised data to non-normalised,
# in terms of how well it deals with library sizes, and feature size effects on the
# variance of the data. You need to pass in the precise assay slot you wish to plot, 
# and also the seurat object
checkNormalisation <- function(assay,seuratobject){
  
  #convert the seurat object to a single cell experiment object, as their plots
  #are better for this kind of an analysis
  if(assay == "SCT"){
    lsks.sce <- as.SingleCellExperiment(seuratobject,assay = "SCT")
  }
  else if(assay == "RNA"){
    lsks.sce <- as.SingleCellExperiment(seuratobject,assay = "RNA")
  }
  
  p1 <- plotPCASCE(
    lsks.sce,
    colour_by = "nCount_RNA",
    #by_exprs_values= "logcounts",
    size_by="nFeature_RNA",
    rerun=TRUE,#add_ticks=F
    run_args = list(exprs_values= "counts")
  ) + ggtitle("raw data")
  
  p2 <- plotPCASCE(
    lsks.sce,
    colour_by = paste("nCount_",assay, sep = ""),
    #by_exprs_values= "logcounts",
    size_by=paste("nFeature_",assay, sep = ""),
    rerun=TRUE,#add_ticks = F,
    run_args = list(exprs_values= "logcounts")
  ) + ggtitle(paste("normalised_",assay, sep = ""))
  
  return(list(p1 = p1, p2 = p2))
  
}


# for two vectors of strings, calculate the jaccard index (overlap)
# between them. Inputs two lists of genes, and outputs a number between 0 and 1
# quantifying overlap between the two gene-sets. 1 suggests perfect overlap and 0
# is no overlap
calculateJaccard <- function(a,b){
  length(intersect(a, b))/length(union(a, b))
}



#install and load required packages
usePackage <- function(p) 
{
  if (!is.element(p, installed.packages()[,1]))
    install.packages(p, dep = TRUE)
  require(p, character.only = TRUE)
}


#get rid of psudogenes, and ribosomal genes from downstream analyses. 
removeNonInformativeGenes <- function(dataset){
  genes.to.remove <- c(rownames(dataset)[grepl("Gm",rownames(dataset))],
                       rownames(dataset)[grepl("Rik",rownames(dataset))],
                       rownames(dataset)[grepl("Rp",rownames(dataset))])
  
  dataset <- as.matrix(dataset)
  dataset <- dataset[rownames(dataset) %ni% genes.to.remove,]
  return(dataset)
  
}


generateSeuratObject <- function(){
  dataset.epo <- Read10X_h5("count_bc_T1/outs/filtered_gene_bc_matrices_h5.h5")
  dataset.ctrl <- Read10X_h5("count_bc_T2/outs/filtered_gene_bc_matrices_h5.h5") 
  
  colnames(dataset.epo) <-  paste(colnames(dataset.epo), "Epo", sep = "")
  colnames(dataset.ctrl) <- paste(colnames(dataset.ctrl), "Ctrl", sep = "")
  dataset <- cbind(dataset.epo,dataset.ctrl)
  
  #remove genes which are not informative for the analysis
  dataset <- removeNonInformativeGenes(dataset)
  
  #add the annotations for EPO treated and control epostatus
  condition <- c(rep("EPO_treated",ncol(dataset.epo)),rep("Control",ncol(dataset.ctrl)))
  cell_anns <- data.frame(condition = condition)
  rownames(cell_anns) <- colnames(dataset)
  
  # now we have the dataset and the metadata we can make the seurat object
  # each gene must occur in at least 20 cells, and each cell must have at least 750 unique genes expressed.
  # We will revisit these filters further downstream in the analysis
  SlamHSCs <- CreateSeuratObject(counts= dataset, min.cells = 20, min.features = 750,
                                 meta.data = cell_anns, project = "EPO")
  
  return(SlamHSCs)
}



clusterRobustnessAnalysis <- function(SlamHSCs,min_resolution, max_resolution){
  
  SlamHSCs <- FindNeighbors(object = SlamHSCs, dims = 1:10, verbose = FALSE,reduction = "pca")
  for(i in seq(from=min_resolution, to=max_resolution, by=0.1)){
    SlamHSCs <- FindClusters(object = SlamHSCs, resolution = i)
  }
  return(SlamHSCs)
  
}



customGeneSetScore <- function(SlamHSCs, gene.set.names, gene.set.fp = "genesets/published_genesets.csv" ){
  
  gene.sets <- read.csv(fp)
  
  for(i in 1:length(gene.set.names)){
    genes.of.interest <- paste(gene.sets[,gene.set.names[i]])
    genes.of.interest <- genes.of.interest[genes.of.interest != ""]
    genes.of.interest <- intersect(genes.of.interest,rownames(SlamHSCs@assays$SCT))
    SlamHSCs <- AddModuleScore(SlamHSCs, features = list(genes.of.interest),name = gene.set.names[i])
  }
  
  return(SlamHSCs)
  
}

volcanoPlot <- function(res){
  res$gene <- rownames(res)
  
  res$sign <- 0
  res$sign[which(res$p_val < 0.05 & res$avg_logFC > 0.1)] <- 2
  res$sign[which(res$p_val < 0.05 & res$avg_logFC < -0.1)] <- 1
  
  p <- ggplot(data=res, aes(x=avg_logFC, y=-log10(p_val), colour=as.factor(sign))) + geom_point( size=2) +
    scale_color_manual(name="", values=c("4" = "orange","3" = "blue" ,"2" = "red","1"=
                                           "grey30", "0"=rgb(220/255,220/255, 220/255,0.2))) +  
    theme(legend.position = "none") + xlim(-0.78,0.78) + 
    xlab("log2 fold change") + ylab("-log10 adj pvalue") + 
    geom_vline(xintercept=c(-0.1, 0.1), linetype=2) + 
    geom_hline(yintercept=-log10(0.05), linetype=2)   
  #geom_text(aes(label=ifelse(res$gene %in% almut.genes, as.character(res$gene),'')),hjust=1,vjust=1, colour = "black")
  
  p  + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "grey")
  )
  
  return(p)
  
}


#this method will create a random binary variable within the metadataslot. The group sizes will be the same as EPO vs untreated
generateRandomGroups <- function(sobj,seed){
  
  set.seed(seed * 4342)
  random.index <- sample(1:ncol(sobj),size =  table(sobj@meta.data$condition)[2][[1]]) 
  
  #assign a binary variable to cells in the random.index group or not
  ident <- rep(0,1,ncol(sobj))
  names(ident) <- seq(1,ncol(sobj))
  ident[random.index] <- 1
  
  names(ident) <- colnames(sobj)
  sobj<- AddMetaData(object = sobj, metadata = ident, col.name = "random")
  return(sobj)
}

performDEA <- function(seed){
  
  SlamHSCs <- generateRandomGroups(SlamHSCs,seed)
  #create some random group assignments
  SlamHSCs@meta.data$ident <- as.factor(SlamHSCs@meta.data$random)
  names(SlamHSCs@meta.data$ident) <- colnames(SlamHSCs@assays$RNA@data)
  Idents(SlamHSCs) <- SlamHSCs@meta.data$ident
  
  
  return(nrow(FindAllMarkers(SlamHSCs, verbose = F,
                             test.use="LR",logfc.threshold = 0.05,
                             only.pos=F,return.thresh = 0.05)))
  
}



setEPOStatus <- function(SlamHSCs, threshold = 0.9){
  #cassign a threhsold for the expressin of our composite score and call all cells above this threshold EPOresponsive
  SlamHSCs@meta.data$eporesponder <- SlamHSCs@meta.data$EPOnet > quantile(SlamHSCs@meta.data$EPOnet ,c(threshold))
  SlamHSCs@meta.data$epostatus <- SlamHSCs@meta.data$condition
  levels(SlamHSCs@meta.data$epostatus) <- c(levels(SlamHSCs@meta.data$condition), "EPO_responsive", "EPO_nonresponsive")
  SlamHSCs@meta.data$epostatus[SlamHSCs@meta.data$eporesponder == TRUE] <- "EPO_responsive"
  SlamHSCs@meta.data$epostatus[SlamHSCs@meta.data$epostatus == "EPO_treated"] <- "EPO_nonresponsive"
  return(SlamHSCs)
}




generate_EPOResponseSignature <- function(SlamHSCs,markers){
  
  markers <- FindAllMarkers(SlamHSCs, verbose = F,
                            test.use="LR",logfc.threshold = 0.05,
                            only.pos=T,return.thresh = 0.05)
  
  
  #filter the signature to take only genes that are significantly enriched in the EPO treated group
  EPO.sig.filtered.up <- markers[markers$cluster == "EPO_treated" & markers$p_val_adj < 0.05,]
  EPO.sig.filtered.down <- markers[markers$cluster =="Control" & markers$p_val_adj < 0.05,]
  
  
  #create a composite score for all of these genes and overlay expression onto our UMAP visualisation
  SlamHSCs <- AddModuleScore(SlamHSCs, features = list(rownames(EPO.sig.filtered.up)),name = "EPOup")
  SlamHSCs <- AddModuleScore(SlamHSCs, features = list(rownames(EPO.sig.filtered.down)),name = "EPOdown")
  
  SlamHSCs@meta.data$EPOnet <-SlamHSCs@meta.data$EPOup1 - SlamHSCs@meta.data$EPOdown1 
  return(SlamHSCs)
  
}



nearestNeighbourMapping <- function(referenceDataset,queryDataset,plotColor = "purple"){
  # run the PCA
  #pca <- prcomp(t(dahlin@assays$RNA@data[genes,]), tol = 0.1)
  #save(pca,file = "dahlin_pcrompPCA.Rda")
  #this step is computationally costly so we just load in one we have performed before
  load("datasets/Dahlin/dahlin_pcrompPCA.Rda")
  
  #lets do external positive controls 
  
  #Use the loadings from the ref dataset PCA to get PCA coords for the query dataset
  pca.scaled <- scale(t(queryDataset@assays$RNA@data[genes,]), pca$center, pca$scale) %*% pca$rotation 
  
  #run the nn algo on the PCA coords
  nn <- nn2(pca$x[,1:10],pca.scaled[,1:10], k = 20)
  
  #get the mean umap coords for all nearest neighbours
  umap.coords <- matrix(0, nrow = ncol(queryDataset),ncol = 2)                    
  rownames(umap.coords) <- colnames(queryDataset)
  for(i in 1:ncol(queryDataset)){
    umap.coords[i,] <- colMeans(dahlin@reductions$umap@cell.embeddings[nn$nn.idx[i,],][1:10,])
  }
  
  
  #plot the results
  x <- DimPlot(object = dahlin,reduction = "umap") + NoLegend()
  g <- ggplot_build(x)
  umap.coords <- umap.coords[colnames(queryDataset),]
  
  plot(x$data[,1],x$data[,2],col = "grey",pch = 20, cex = 0.5) + 
    points(x = umap.coords[,1], y = umap.coords[,2] ,cex = 1, pch =20,col = plotColor)
  
}



prepWilsonDataset <- function(){
  load("datasets/Wilson/wilson_processed.Rda") # Wilson dataset
  wilson <- CreateSeuratObject(counts = counts(wilson), project = "hscs", min.cells = 3, min.features = 200)
  wilson <- FindVariableFeatures(wilson, selection.method = "vst", nfeatures = 2000)
  wilson <- NormalizeData(wilson)
  save(wilson,file = "datasets/Wilson/wilson_seurat.rda")
  
  
}