
version
packageStatus()


library(Seurat)
packageVersion("Seurat")
library(ggplot2); theme_set(theme_classic())

library(Matrix)

#Exclude clusters 2,9 and calculate markers
allclusters = c("0","1","3","5","7","8","11","13","16","17")
Idents(mg) = "seurat_clusters"
markersP5ctrl2 = FindMarkers(mg,ident.1 = "2",ident.2 = allclusters,test.use = "MAST")
markers9 = FindMarkers(mg,ident.1 = "9",ident.2 = allclusters,test.use = "MAST")
markers2_9 = cbind(markersP5ctrl2,markers9)
write.csv(markers9,"markers9vallbut2.csv")
write.csv(markersP5ctrl2,"markersP5ctrl2vallbut9.csv")


m2 = markersP5ctrl2[markersP5ctrl2$p_val_adj<0.01,]
m9 = markers9[markers9$p_val_adj<0.01,]

m2up = m2[m2$avg_logFC>0,]
m9up = m9[m9$avg_logFC>0,]
m2down = m2[m2$avg_logFC<0,]
m9down = m9[m9$avg_logFC<0,]

coup = m2up[rownames(m2up)%in% rownames(m9up),]
codown = m2down[rownames(m2down)%in% rownames(m9down),]
cross1 = m2up[rownames(m2up)%in% rownames(m9down),]
cross2 = m2down[rownames(m2down)%in% rownames(m9up),]
load(file="~/Desktop/LD_AVM02/Data/MG_all_Scran.RData")

MG_all


#For both of these, see if I can remove the legend so the plots are the same size
cellnames=experiment.aggregate@cell.names[experiment.aggregate@meta.data$orig.ident=="A1"]
plot2 = TSNEPlot(object=experiment.aggregate,cells.use=cellnames,pt.size=0.5,group.by = "region", do.label=F,plot.title="WT",no.legend = T)
setEPS()
postscript("WTtsneregion.eps")
plot2
dev.off()


cellnames1=experiment.aggregate@cell.names[experiment.aggregate@meta.data$SexID=="Female"&experiment.aggregate@meta.data$orig.ident=="A3"]
cellnames2=experiment.aggregate@cell.names[experiment.aggregate@meta.data$orig.ident=="A1"]
plot4 = TSNEPlot(object=experiment.aggregate,cells.use=c(cellnames1,cellnames2),pt.size=0.5,do.label=F,group.by = "orig.ident",plot.order = "A3", colors.use = c("darkgrey","red3"),plot.title="Mcherry Hi",no.legend. = T)
setEPS()
postscript("mCherrytsneall.eps")
plot4
dev.off()

plot4 = TSNEPlot(object=experiment.aggregate,cells.use=c(cellnames1),pt.size=0.5,do.label=F,group.by = "orig.ident", colors.use = c("red3"),plot.title="Mcherry Hi",no.legend. = T)

setEPS()
postscript("mCherrytsneallmCherry.eps")
plot4
dev.off()

plot4 = TSNEPlot(object=experiment.aggregate,cells.use=c(cellnames2),pt.size=0.5,do.label=F,group.by = "orig.ident", colors.use = c("darkgrey"),plot.title="Mcherry Hi",no.legend. = T)

setEPS()
postscript("mCherrytsneallWT.eps")
plot4
dev.off()


##TSNE plots that justify region labels
genes = c("Gad2","Satb2","Prox1","Socs2","Amigo2","Rgs14","Map3k15","Foxg1","Sv2b","Cdh24")
#DG marker: Prox1
#CA1 marker: Satb2
#IN marker: Gad2
#CA2 markers: amigo2, rgs14, map3k15
#mossy cells: foxg1, sv2b
#CA3: cdh24
for(gene in genes){
    setEPS()
    postscript(paste0(gene,"tsne.eps"))
    FeaturePlot(experiment.aggregate, gene, cols.use = c("lightgrey", "blue"))
    dev.off()
}
           

#Vln plots that justify region labels
genes = c("Gad2","Satb2","Prox1","Socs2","Amigo2","Rgs14","Map3k15","Foxg1","Sv2b","Cdh24")
i = 10
gene = genes[i]
setEPS()
postscript(paste0(gene,"vln.eps"))
VlnPlot(experiment.aggregate, gene, group.by = "region")
dev.off()




##Plot Il33 within the WT sample:
cellnames=experiment.aggregate@cell.names[experiment.aggregate@meta.data$orig.ident=="A1"]
FeaturePlot(object = experiment.aggregate, cells.use = cellnames, c("Il33"), overlay= F,cols.use= c("lightgrey","red"),)

setEPS()
postscript("Il33inWT.eps")
FeaturePlot(object = experiment.aggregate, cells.use = cellnames, c("Il33"), overlay= F,cols.use= c("lightgrey","red"))
dev.off()


#*************************************************************************************************
##Pull out DG cells and cluster them separately
cellnames1=experiment.aggregate@cell.names[experiment.aggregate@meta.data$region=="DG"&experiment.aggregate@meta.data$mouseIDfull %in% c("A1 Female","A1 Male","A3 Female")]
DG.subset<- SubsetData(
  object = experiment.aggregate,
  subset.raw=T,
  do.clean=T, #do.clean removes things like scaledata, tsne, pca so that you can recalculate them below
  cells.use = cellnames1)

##Renormalize using subset
DG.subset <- NormalizeData(
  object = DG.subset,
  normalization.method = "LogNormalize",
  scale.factor = 10000)

##Set up a method for filtering genes with low expression
FilterGenes <- 
  function (object, min.value=1, min.cells = 0, genes = NULL) {
    parameters.to.store <- as.list(environment(), all = TRUE)[names(formals("FilterGenes"))]
    object <- Seurat:::SetCalcParams(object = object, calculation = "FilterGenes", ... = parameters.to.store)
    genes.use <- rownames(object@data)
    
    if (!is.null(genes)) {
      genes.use <- intersect(genes.use, genes)
      object@data <- object@data[genes.use, ]
      return(object)
    } else if (min.cells > 0) {
      num.cells <- Matrix::rowSums(object@data > min.value)
      genes.use <- names(num.cells[which(num.cells >= min.cells)])
      object@data <- object@data[genes.use, ]
      return(object)
    } else {
      return(object)
    }
  }

#Filter genes - set min.value and min.cells as needed. This will remove any genes that don't have at least (min.value) counts in (min.cells) # of cells.
DG.subset <- FilterGenes(object = DG.subset, min.value = 1, min.cells = 5)
DG.subset

#Determine the desired number of variable genes
length(FindVariableGenes(
  object = DG.subset,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  x.low.cutoff = 0.1,
  x.high.cutoff = 5,
  y.cutoff = 0.1, do.plot=T)@var.genes) 

#Set variable genes
DG.subset <- FindVariableGenes(
  object = DG.subset,
  mean.function = ExpMean,
  dispersion.function = LogVMR,
  x.low.cutoff = 0.1,
  x.high.cutoff = 5,
  y.cutoff = 0.1, do.plot=T) #3147 genes

#find out how many of the variable genes are mitochondrial - a good indication that you are not stringent enough
grep("mt-",DG.subset@var.genes,value=TRUE) #5 when y.cutoff = 2

length(DG.subset@var.genes)

#Set a metadata column with meaningful names if necessary. Make sure that the replacement levels are in the same order as their corresponding levels from the original metadat column (show in the table command)
fulldescription<-DG.subset@meta.data$mouseIDfull
names<-rownames(DG.subset@meta.data)
names(fulldescription)<-names
table(fulldescription)
levels(fulldescription)<-c("Control Female","Control Male", "mCherry Female","ST2cKO Male","ST2cKO Enriched Female","ST2cKO Enriched Male")

DG.subset <- AddMetaData( 
  object = DG.subset,
  metadata = fulldescription,
  col.name = "fulldescription")

#Scale using variable genes, then run PCA and plot the results
DG.subset <- ScaleData(
  object = DG.subset,
  genes.use=DG.subset@var.genes)

DG.subset <- RunPCA(
  object = DG.subset,
  pc.genes = DG.subset@var.genes,
  do.print = TRUE,
  pcs.print = 1:5,
  genes.print = 5,
  pcs.compute = 100,
  maxit = 500)


#For all the following: add in a full box, and make the x and y dimensions equal
#For supplement: #Make a PCA plot separated by "Control F","Control M", "mCherry Female"
setEPS()
postscript("DGPCAplot.eps")
PCAPlot(object = DG.subset,group.by= "fulldescription",dim.1 = 1,dim.2 = 2,cols.use = c("darkgrey","lightgrey","red3"))
dev.off()

#For figure: Make pca plot showing control/Mcherry in grey/red
setEPS()
postscript("DGPCAplotWTmCherry.eps")
PCAPlot(object = DG.subset,group.by= "orig.ident",dim.1 = 1,dim.2 = 2,
        cols.use = c("darkgrey","red3"),no.legend = T,coord.fixed = F,vector.friendly = T,
        do.return = T,pt.size = 3) + xlim(-11,15) + ylim(-11,15)
dev.off()

# On Seurat 2:

mat <- DG.subset@scale.data
pca <- DG.subset@dr$pca

# Get the total variance:
total_variance <- sum(matrixStats::rowVars(mat))

eigValues = (pca@sdev)^2  ## EigenValues
varExplained = eigValues / total_variance


varExplained2 = eigValues / sum(eigValues)

PCElbowPlot(DG.subset,num.pc = 40)
#For figure: Make pca plot showing control/Mcherry in blue/red
setEPS()
postscript("DGPCAplotWTmCherryblue.eps")
PCAPlot(object = DG.subset,group.by= "orig.ident",dim.1 = 1,dim.2 = 2,cols.use = c("cornflowerblue","red3"), no.legend = T,coord.fixed = F,vector.friendly = T,
        do.return = T,pt.size = 3) + xlim(-11,15) + ylim(-11,15)
dev.off()

#Just wt, don't use red
setEPS()
postscript("DGPCAplotclustersWTonly.eps")
PCAPlot(
  object = DG.subset,
  group.by= "finalcluster",cells.use = DG.subset@cell.names[DG.subset@meta.data$fulldescription %in% c("Control Female","Control Male")],
  dim.1 = 1,
  dim.2 = 2,
  cols.use = c("darkgrey","cornflowerblue"),no.legend = T,coord.fixed = F,vector.friendly = T,
  do.return = T,pt.size = 3) + xlim(-11,15) + ylim(-11,15)
dev.off()

#Find clusters: will take a long time to run
DG.subset <- FindClusters(
  object = DG.subset, 
  reduction.type = "pca", 
  dims.use = use.pcs, 
  resolution = 1, #pick multiple resolutions; smaller = fewer clusters seq(0.5,4,0.5), #
  print.output = FALSE, 
  save.SNN = TRUE,
  force.recalc = TRUE
)
PrintFindClustersParams(object = DG.subset)

#This code tells you how many resolutions you calculated clusters for and prints out the number of clusters for each one
sapply(grep("^res",colnames(DG.subset@meta.data),value = TRUE),
       function(x) length(unique(DG.subset@meta.data[,x])))


#Assign a cluster ID to each cell; to color with higher resolution, change this id to the resolution you care about
DG.subset <- SetAllIdent(DG.subset, id = "res.0.5")

#Run t-sne. Set the seed to get reproducible plots.
set.seed(1)
DG.subset <- RunTSNE(
  object = DG.subset,
  reduction.use = "pca",
  dims.use = use.pcs,
  do.fast = TRUE)



#The following will allow you to calculate relative abundance of particular samples in particular clusters
t<-table(DG.subset@meta.data$fulldescription,DG.subset@meta.data$res.0.5)
c2=colSums(t)
r=100*t[3,]/c2
l=100*t[4,]/c2
z=100*colSums(t[5:6,])/c2
t<-rbind(t,"Total"=c2,"%mCherry"=r,"%ST2cKO"=l,"%Enriched"=z)
t[8,]<-sapply(t[8,],function(x)round(x,digits=2))

t[9,]<-sapply(t[9,],function(x)round(x,digits=2))
t[10,]<-sapply(t[10,],function(x)round(x,digits=2))
t

#T-sne plots for the DG subset specifically. each "Cellnames" simply selects the desired cells to be shown in the next t-sne plot.


save(object=DG.subset,file="DG.subset.RData")

#Find markers for each cluster
markersDG = FindAllMarkers(
  object = DG.subset, 
  only.pos = FALSE, 
  min.pct = 0.25, #gene must be present in 25% of the cells in the cluster
  logfc.threshold = 0.25
  )
dim(markersDG)
head(markersDG)
write.csv(markersDG,"markers_DGsubset.csv")
#clusters 2 and 4 have high IL-33 (therefore Xist)


DG.subset<-SetAllIdent(DG.subset,id="clustersampleID")



#VolcanoPlot

fc = read.csv("~/Desktop/LD_AVM02/spreadsheets/markersdep5micro11.csv",stringsAsFactors = F)
colnames(fc)[1] = "Gene"
fc = fc[!is.na(fc$avg_logFC),]
colorkeysdown = fc$Gene[fc$avg_logFC < -log2(1.3) & fc$p_val_adj < 10e-8]
colorkeysup = fc$Gene[fc$avg_logFC > log2(1.3) & fc$p_val_adj < 10e-8]
allcolors = rep("darkgrey",length(fc$Gene))
names(allcolors) = fc$Gene
allcolors[names(allcolors) %in% colorkeysdown] = "blue"
names(allcolors)[allcolors == "blue"] = "d"
allcolors[names(allcolors) %in% colorkeysup]= "red"
names(allcolors)[allcolors == "red"] = "u"
names(allcolors)[allcolors == "darkgrey"] = "-"

library(EnhancedVolcano)
EnhancedVolcano(fc,
lab = fc$Gene,
x = 'avg_logFC',
y = 'p_val_adj',
xlim = c(-5, 5),
title = 'MG P5 DEP CLUSTER',
subtitle = "",
drawConnectors = F,
legendPosition = 'right',
legendVisible = F,
#shade = genestolabel,
pCutoff = 10e-8,
FCcutoff = log2(1.2),
#selectLab = genestolabel,
transcriptPointSize = 1.5,
transcriptLabSize = 3.0,
col=c('black', 'black', 'black', 'red3'),
colCustom = allcolors,
gridlines.major = F,
gridlines.minor = F,
colAlpha = 1)


setEPS()
postscript("mgp5depvolc.eps")
EnhancedVolcano(fc,
                lab = fc$Gene,
                x = 'avg_logFC',
                y = 'p_val_adj',
                xlim = c(-5, 5),
                title = 'MG P5 DEP CLUSTER',
                subtitle = "",
                drawConnectors = F,
                legendPosition = 'right',
                legendVisible = F,
                #shade = genestolabel,
                pCutoff = 10e-8,
                FCcutoff = log2(1.2),
                #selectLab = genestolabel,
                transcriptPointSize = 1.5,
                transcriptLabSize = 3.0,
                col=c('black', 'black', 'black', 'red3'),
                colCustom = allcolors,
                gridlines.major = F,
                gridlines.minor = F,
                colAlpha = 1)
dev.off()

#Make a chart showing the relative contribution of mCherry population to each cluster
#Then maybe do this on the full sample for WT and mCherry only but with all regions
TSNEPlot(DG.subset,group.by = "fulldescription",vector.friendly = T, no.legend = F,colors.use = c("lightblue","blue","red3"))
TSNEPlot(DG.subset,group.by = "res.1",vector.friendly = T, no.legend = F,do.label = F)

relcon = table(droplevels(DG.subset@meta.data$fulldescription),DG.subset@meta.data$res.1)
relcon

DG.subset = SetAllIdent(DG.subset,id = "res.1")
DG.subset<- BuildClusterTree(
  DG.subset,
  pcs.use = use.pcs,
  do.reorder = F,
  reorder.numeric = F,
  do.plot=F)

PlotClusterTree(DG.subset, cex=1, use.edge.length=F,node.depth=2)
library(ape)

# Stacked Bar Plot with Colors and Legend
setEPS()
postscript("clustercomp.eps")
barplot(relcon, main="Cluster composition by sample",
        xlab="Cluster", ylab = "Number of Cells", ylim = c(0,1300), col=c("blue","lightblue","red"),
        legend = rownames(relcon))
dev.off()

Control_Female = relcon[1,]/colSums(relcon)
Control_Male = relcon[2,]/colSums(relcon)
mCherry_Female = relcon[3,]/colSums(relcon)
relconpct = rbind(Control_Female,Control_Male,mCherry_Female)
relconpct

setEPS()
postscript("clustercomppct.eps")
barplot(relconpct, main="Cluster composition by percent of sample",
        xlab="Cluster", ylab = "Number of Cells", ylim = c(0,2), col=c("darkblue","lightblue","red"),axis.lty = 1, 
        width = .5,xlim = c(0,5),legend = rownames(relcon), space = 0.5)
dev.off()

#contribution of mcherry to each cluster
m = experiment.aggregate@meta.data$mouseIDfull
names(m) = experiment.aggregate@cell.names
levels(m) = c("Control Female","Control Male","mCherry Female","ST2cKO Male","ST2cKO Enriched Female","ST2cKO Enriched Male")
experiment.aggregate = AddMetaData(experiment.aggregate, metadata = m, col.name = "fulldescription")
cellnames=experiment.aggregate@cell.names[experiment.aggregate@meta.data$mouseIDfull%in%c("A1 Female","A1 Male","A3 Female")]

sub<- SubsetData(
object = experiment.aggregate,
subset.raw=T,
do.clean=F, #do.clean removes things like scaledata, tsne, pca so that you can recalculate them below
cells.use = cellnames)


relcon = table(mg$fulldescription,mg$finalcluster)
relcon

# Stacked Bar Plot with Colors and Legend
setEPS()
postscript("clustercompall.eps")
barplot(relcon, main="Cluster composition by sample",
xlab="Cluster", ylab = "Number of Cells", ylim = c(0,3000), col=c("blue","lightblue","red"),
width = .2,xlim = c(0,6), space = 0.5,cex.names = 0.2,axis.lty = 1,
legend = rownames(relcon))
dev.off()

Control_Female = relcon[1,]/colSums(relcon)
Control_Male = relcon[2,]/colSums(relcon)
mCherry_Female = relcon[3,]/colSums(relcon)
relconpct = rbind(Control_Female,Control_Male,mCherry_Female)
relconpct

setEPS()
postscript("clustercomppctall.eps")
barplot(relconpct, main="Cluster composition by percent of sample",
xlab="Cluster", ylab = "Number of Cells", ylim = c(0,2), col=c("darkblue","lightblue","red"),axisnames = T,
width = .2,xlim = c(0,6),legend = rownames(relcon), space = 0.5,cex.names = 0.2,axis.lty = 1)
dev.off()
