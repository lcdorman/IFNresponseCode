library(Seurat)

library(ggplot2); theme_set(theme_classic())

#Load MG only (LD) file called "mgAVM02"
#load("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Sequencing/LD_AVM02/Data/mgAVM02_filtered_Feb2020.RData")


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
    Idents(sobject) = cellnames[1]
    cells = colnames(sobject)[Idents(sobject) %in% cellnames[2:length(cellnames)]]} 
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
all(features %in% rownames(mgAVM02))
for(feature in features){
  PrintSeuratGraph(namecard = "MGonly",sobject=mgAVM02,graphtype = "feature",feature = feature)
}

for(feature in features){
  PrintSeuratGraph(namecard = "newthresholds_allcells",sobject=MG,graphtype = "feature",feature = feature)
}
for(feature in features){
  PrintSeuratGraph(namecard = "CP5",sobject=mgAVM02,graphtype = "feature",feature = feature,cellnames = c("sample_description","Control_P5"))
}

for(feature in features){
  PrintSeuratGraph(namecard = "DP5",sobject=mgAVM02,graphtype = "feature",feature = feature,cellnames = c("sample_description","Deprived_P5"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "CP7",sobject=mgAVM02,graphtype = "feature",feature = feature,cellnames = c("sample_description","Control_P7"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "DP7",sobject=mgAVM02,graphtype = "feature",feature = feature,cellnames = c("sample_description","Deprived_P7"))
}

groups = c("age","sex","seurat_clusters")

for(group in groups){
  PrintSeuratGraph(namecard = "CP5",sobject=mgAVM02,graphtype = "dim",group = group, feature = group,cellnames = c("sample_description","Control_P5"))
}

for(group in groups){
  PrintSeuratGraph(namecard = "DP5",sobject=mgAVM02,graphtype = "dim",group = group, feature = group,cellnames = c("sample_description","Deprived_P5"))
}
for(group in groups){
  PrintSeuratGraph(namecard = "CP7",sobject=mgAVM02,graphtype = "dim",group = group, feature = group,cellnames = c("sample_description","Control_P7"))
}
for(group in groups){
  PrintSeuratGraph(namecard = "DP7",sobject=mgAVM02,graphtype = "dim",group = group, feature = group,cellnames = c("sample_description","Deprived_P7"))
}

for(feature in features){
  PrintSeuratGraph(namecard = "MGonlyvlnP5",sobject=mgAVM02,graphtype = "violin",group = "seurat_clusters",feature = feature,split = "sample_description",cellnames = c("age","P5"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "MGonlyvlnP7",sobject=mgAVM02,graphtype = "violin",group = "seurat_clusters",feature = feature,split = "sample_description",cellnames = c("age","P7"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "MGonlyvlnCtrl",sobject=mgAVM02,graphtype = "violin",group = "seurat_clusters",feature = feature,split = "sample_description",cellnames = c("condition","Control"))
}
for(feature in features){
  PrintSeuratGraph(namecard = "MGonlyvlnDep",sobject=mgAVM02,graphtype = "violin",group = "seurat_clusters",feature = feature,split = "sample_description",cellnames = c("condition","Deprived"))
}



r = table(mgAVM02$seurat_clusters,mgAVM02$sample_description)
Control_P5 = r[,1]*100/colSums(r)[1]
Control_P7 = r[,2]*100/colSums(r)[2]
Deprived_P5 = r[,3]*100/colSums(r)[3]
Deprived_P7 = r[,4]*100/colSums(r)[4]
relconpct = cbind(Control_P5,Deprived_P5,Control_P7,Deprived_P7)
relconpct

#notes: P5 - 4,5 P7 - 1,3 P5 Dep = 8 (at the expense of 0), P5 Ctrl - 5,4 
setEPS()
postscript("clustercomppctall.eps")
barplot(relconpct, main="Cluster composition by percent of sample",
        xlab="Cluster", ylab = "Number of Cells", ylim = c(0,100), col=c("darkblue","lightblue","red","orange","purple","lightgrey","darkgrey","magenta","yellow","black","red3","cornflowerblue"),axisnames = T,
        width = .2,xlim = c(0,2),legend = rownames(relconpct), space = 0.6,cex.names = 0.7,axis.lty = 1)
dev.off()

#Find out the relative proportion of male/female cells per cluster
ratio = table(mgAVM02$sex,mgAVM02$seurat_clusters)
Female = ratio[1,]/rowSums(ratio)[1]
Male = ratio[2,]/rowSums(ratio)[2] #normalize each row to 1000 cells per sex
ratio = rbind("Female" = Female*1000,"Male" = Male*1000)
rowSums(ratio)
#find out what percent of each cluster is made by each sex
Female = ratio[1,]*100/colSums(ratio)
Male = ratio[2,]*100/colSums(ratio)
ratio2 = rbind(Female,Male)
setEPS()
postscript("sexpct.eps")
barplot(ratio2, main="Cluster composition by sex",
        xlab="% of cluster", ylab = "Cluster", ylim = c(0,4), col=c("darkblue","lightblue"),axisnames = T,
        width = .2,xlim = c(0,100),legend = rownames(ratio2), space = 0.6,cex.names = 0.7,axis.lty = 1,horiz = T)
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

VlnPlot(mgAVM02,features = c("nCount_RNA","percent.mito"),group.by = "sample_description")


#VolcanoPlot

fc = read.csv("/Users/whippoorwill/Desktop/plots/allmarkers_vargenesMG_0215.csv",stringsAsFactors = F)
colnames(fc)[8] = "Gene"

fc = read.csv("/Users/whippoorwill/Desktop/plots/allmarkers_vargenesMG_0308_8.csv",stringsAsFactors = F)
colnames(fc)[1] = "Gene"
newlist = list()
#Split by cluster
for (i in c(0:10)){
  newlist[[i+1]] = fc[fc$cluster == i,]
}


#select a single cluster
cluster = 8

  fc = newlist[[cluster+1]]
  #Format data table
  fc = fc[!is.na(fc$avg_logFC),]
  colorkeysdown = fc$Gene[fc$avg_logFC < -log2(1.15) & fc$p_val_adj < 10e-25]
  colorkeysup = fc$Gene[fc$avg_logFC > log2(1.15) & fc$p_val_adj < 10e-25]
  
  a = fc$Gene[fc$avg_logFC < -log2(1.7) & fc$p_val_adj < 10e-25]
  length(a)
  b = fc$Gene[fc$avg_logFC > log2(1.7) & fc$p_val_adj < 10e-25]
  length(b)
  allcolors = rep("darkgrey",length(fc$Gene))
  names(allcolors) = fc$Gene
  allcolors[names(allcolors) %in% colorkeysdown] = "blue"
  names(allcolors)[allcolors == "blue"] = "d"
  allcolors[names(allcolors) %in% colorkeysup]= "red"
  names(allcolors)[allcolors == "red"] = "u"
  names(allcolors)[allcolors == "darkgrey"] = "-"
  
  setEPS()
  postscript(paste0("volcano_22620_",cluster,"down.eps"))
  EnhancedVolcano(fc,
                  lab = fc$Gene,
                  x = 'avg_logFC',
                  y = 'p_val_adj',
                  xlim = c(-2, 2),
                  title = paste0("LDAVM02_",cluster),
                  subtitle = "",
                  drawConnectors = F,
                  legendPosition = 'right',
                  legendVisible = F,
                  #shade = genestolabel,
                  pCutoff = 10e-25,
                  FCcutoff = log2(1.15),
                  selectLab = a,
                  transcriptPointSize = 1.5,
                  transcriptLabSize = 4.0,
                  col=c('black', 'black', 'black', 'red3'),
                  colCustom = allcolors,
                  gridlines.major = F,
                  gridlines.minor = F,
                  colAlpha = 1)
  dev.off()
  
  setEPS()
  postscript(paste0("volcano_22620_",cluster,"up.eps"))
  EnhancedVolcano(fc,
                  lab = fc$Gene,
                  x = 'avg_logFC',
                  y = 'p_val_adj',
                  xlim = c(-3, 4),
                  title = paste0("LDAVM02_",cluster),
                  subtitle = "",
                  drawConnectors = F,
                  legendPosition = 'right',
                  legendVisible = F,
                  #shade = genestolabel,
                  pCutoff = 10e-25,
                  FCcutoff = log2(1.15),
                  selectLab = b,
                  transcriptPointSize = 1.5,
                  transcriptLabSize = 4.0,
                  col=c('black', 'black', 'black', 'red3'),
                  colCustom = allcolors,
                  gridlines.major = F,
                  gridlines.minor = F,
                  colAlpha = 1)
  dev.off()
  
cluster = cluster+1
  