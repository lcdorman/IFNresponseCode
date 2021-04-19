library(Seurat)
library(ggplot2)

#Load MG only (LD) file called "mgAVM02"
load(file = "~/Desktop/LD_AVM02/Data/mgonlyAVM02.RData")

#Load all cells (LD) file called "MG"
load("/Users/whippoorwill/Desktop/LD_AVM02/Data/MG_umap.RData")

#Load LD+ Stevens dataset, pre-clustered called "both"
load("/Users/whippoorwill/Desktop/LD_AVM02/Data/LD_Scombined.RData")


#Write a function that will print out a graph given a metadata or feature name, seurat object, and split.by/group.by parameters
#sobject = seurat object name
#graphtype = "violin","feature","dim"
#feature = gene name or metadata column to feature
#group = only for violin; metadata to group by
#split = metadata column to split on
#cellnames = vector: metadata column, subset to include (any length); required for violin plot
#namecard = token to include in the plot name (usually the name of the sobject)
setwd("~/Desktop/plots")

PrintSeuratGraph = function(namecard = "a",sobject,graphtype = "feature",feature = NULL,group = NULL,split=NULL,cellnames=NULL){
  if (!is.null(cellnames)){
    Idents(both) = cellnames[1]
    cells = colnames(both)[Idents(both) %in% cellnames[2:length(cellnames)]]} 
  else {cells = cellnames}
  if (graphtype == "feature"){
    graph = FeaturePlot(sobject,features = feature,split.by = split, cells = cells)
  }
  if (graphtype == "violin"){
    graph = VlnPlot(sobject,features = feature, pt.size = 0.1, idents = cellnames[2:length(cellnames)],group.by = group, split.by = split)
  }
  if (graphtype == "dim"){
    graph = DimPlot(sobject,cells = cells, group.by = group, split.by = split)
    
  }
  name = paste0(feature,"_",graphtype,namecard,".eps")
  graph
  setEPS()
  postscript(name)
  print(graph)
  dev.off()
}

#should print a graph called "Ptprc_featureboth.RData" showing expression of ptprc split by id (SvL) using every cell

features = c("nCount_RNA","percent.mito","Mertk","Axl","Tyrobp","P2ry12","P2ry13","Tmem119","Trem2","Cx3cr1","Spp1","Hexb","Fcrls","C3","Ptprc","Itgam","Csf1r","Cd68","Ifitm3","Ctsb","Ctsd")
all(features %in% rownames(both))
for(feature in features){
  PrintSeuratGraph(namecard = "MGonly",sobject=mgAVM02,graphtype = "feature",feature = feature)
}

for(feature in features){
  PrintSeuratGraph(namecard = "both",sobject=both,graphtype = "feature",feature = feature,split = "id")
}

for(feature in features){
  PrintSeuratGraph(namecard = "LD",sobject=MG,graphtype = "feature",feature = feature)
}

groups = c("age","sex","sample_description","celltypecluster","seurat_clusters")
for(group in groups){
  PrintSeuratGraph(namecard = "LD",sobject=MG,graphtype = "dim",group = group,feature = group)
}

for(group in groups){
  PrintSeuratGraph(namecard = "MGonly",sobject=mgAVM02,graphtype = "dim",group = group,feature = group,split = "sample_description")
}

for(group in groups){
  PrintSeuratGraph(namecard = "both",sobject=both,graphtype = "dim",group = group,feature = group,split = "id")
}

for(feature in features){
  PrintSeuratGraph(namecard = "MGonly",sobject=mgAVM02,graphtype = "vln",group = "celltypecluster",feature = feature,split = "sample_description",cellnames = c("age","P5"))
}

i=1
feature = features[i]
setEPS()
postscript(paste0(feature,"P5vln_MGonly.eps"))
VlnPlot(mgAVM02,feature,group.by = "celltypecluster",split.by = "sample_description",idents = "P5",pt.size = 0.1)
dev.off()


r = table(mgAVM02$celltypecluster,mgAVM02$sample_description)
Control_P5 = r[,1]*100/colSums(r)[1]
Control_P7 = r[,2]*100/colSums(r)[2]
Deprived_P5 = r[,3]*100/colSums(r)[3]
Deprived_P7 = r[,4]*100/colSums(r)[4]
relconpct = cbind(Control_P5,Deprived_P5,Control_P7,Deprived_P7)
relconpct

setEPS()
postscript("clustercomppctall.eps")
barplot(relconpct, main="Cluster composition by percent of sample",
        xlab="Cluster", ylab = "Number of Cells", ylim = c(0,100), col=c("darkblue","lightblue","red","orange","purple","lightgrey","darkgrey","magenta","yellow","black","red3","cornflowerblue"),axisnames = T,
        width = .2,xlim = c(0,2),legend = rownames(relconpct), space = 0.6,cex.names = 0.7,axis.lty = 1)
dev.off()


Idents(mgAVM02) = "seurat_clusters"
mgAVM02= BuildClusterTree(mgAVM02)
tree = mgAVM02@tools$BuildClusterTree
setEPS()
postscript("newtreeMGonly.eps")
plot.phylo(tree, use.edge.length = T, direction = "rightwards")
dev.off()

FeaturePlot(mgAVM02,features = "percent.mito")
FeaturePlot(mgAVM02,features = "nCount_RNA")
