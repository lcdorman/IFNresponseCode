
#associated pzfx file shows LFC versus internal control - within each experiment
hanson = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/Comparison with Cluster 8/hanson.csv",stringsAsFactors = F)
leah = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/Spreadsheets/allmarkers_vargenesMG_0501.csv",stringsAsFactors = F)
c8 = leah[leah$cluster == "8" & leah$p_val_adj < 0.005,]
c8 = c8[c8$avg_logFC>0.2,]

#colnames(hanson) = paste0(colnames(hanson),hanson[1,])
#hanson = hanson[,c(1,8,20:ncol(hanson))]

genes = c8$gene
hansongenes = hanson$MouseGene.Symbol
hansongenes = hansongenes[hansongenes %in% genes]

hanson8 = hanson[hanson$MouseGene.Symbol %in% hansongenes,]
rownames(hanson8) = hanson8$MouseGene.Symbol
hanson = hanson[,-1]

AD = grep("AD",colnames(hanson8))
AD = c(AD,grep("hMAPT",colnames(hanson8)))
AD = c(AD,grep("APP",colnames(hanson8)))
AD = c(AD,grep("5xFAD",colnames(hanson8)))
AD = c(AD,grep("PS2APP",colnames(hanson8)))
Trem2 = grep("Trem2",colnames(hanson8))
HD = grep("HD",colnames(hanson8))
FTD = grep("FTD",colnames(hanson8))
ALS = grep("ALS",colnames(hanson8))
ALS = c(ALS,grep("SOD1",colnames(hanson8)))
LPS = grep("LPS",colnames(hanson8))
injury = grep("transischemia",colnames(hanson8))
injury = c(injury,grep("cuprizone",colnames(hanson8)))

aging = grep("22M",colnames(hanson8))
aging = c(aging,grep("13months",colnames(hanson8)))


features = list("AD"=AD,"Trem2"=Trem2,"HD" = HD,"FTD" = FTD,"ALS" = ALS,"LPS" = LPS,"injury" = injury,"aging" = aging)

for (i in 1:length(features)){
  x = hanson8[,c(1,4,features[[i]])]
  pval = x[,c(1,2,grep("adjP",colnames(x)))]
  lfc = x[,c(1,2,grep("lfc",colnames(x)))]
  keep = apply(pval,1,function(a)min(a,na.rm = T))
  lfc = lfc[keep<0.01,]
  features[[i]] = lfc
  print(table(lfc$X.9Myeloid.Activation..Coarse.))
}



#Essentially - 21 out of 31 interferon module genes are identified in cluster 8



Trem2hanson = hanson8[,c(1,4,Trem2)]
HDhanson = hanson8[,c(1,4,HD)]
FTDhanson = hanson8[,c(1,4,FTD)]
ALShanson = hanson8[,c(1,4,ALS)]
LPShanson = hanson8[,c(1,4,LPS)]
injuryhanson = hanson8[,c(1,4,injury)]
aginghanson = hanson8[,c(1,4,aging)]


hammondIR2 = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/Comparison with Cluster 8/IR2.csv",stringsAsFactors = F)
hammondIR2 = hammondIR2[2:nrow(hammondIR2),1:6]
colnames(hammondIR2) = hammondIR2[1,]
hammondIR2 = hammondIR2[2:nrow(hammondIR2),]
for (x in 1:nrow(hammondIR2)){
  fc =  log2(as.numeric(hammondIR2[x,"Fold Change"]))
  hammondIR2[x,"Fold Change"] = fc
}


hammondoa3 = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/Comparison with Cluster 8/hammondoa3.csv",stringsAsFactors = F)
hammondoa3 = hammondoa3[2:nrow(hammondoa3),1:6]
colnames(hammondoa3) = hammondoa3[1,]
hammondoa3 = hammondoa3[2:nrow(hammondoa3),]
for (x in 1:nrow(hammondoa3)){
  fc =  log2(as.numeric(hammondoa3[x,"Fold Change"]))
  hammondoa3[x,"Fold Change"] = fc
}


kerenshauldam = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/Comparison with Cluster 8/kerenshauldam.csv",stringsAsFactors = F)
#kerenshauldam = kerenshauldam[kerenshauldam$DAM.FDR.p.value>5,]

hammondIR2 = hammondIR2[hammondIR2$Gene %in% hansongenes,]
hammondoa3 = hammondoa3[hammondoa3$Gene %in% hansongenes,]
kerenshauldam = kerenshauldam[kerenshauldam$Gene %in% hansongenes,]

x = as.data.frame(hansongenes)
x$kerenshauldam = rep(NA,39)
x$hammondoa3 = rep(NA,39)
x$hammondir2 = rep(NA,39)

y = match(x$hansongenes,kerenshauldam$Gene.name)
x$kerenshauldam = kerenshauldam[y,"Fold.change..DAM.to.homeostatic.microglia."]


y = match(x$hansongenes,hammondIR2$Gene)
x$hammondir2 = hammondIR2[y,"Fold Change"]

y = match(x$hansongenes,hammondoa3$Gene)
x$hammondoa3 = hammondoa3[y,"Fold Change"]
write.csv(x,file = "/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/Comparison with Cluster 8/singlecelladdition.csv")

t = as.data.frame(a)
t$a = as.character(t$a)
rownames(t) = a
for (i in 1:nrow(t)){
  if (t$a[i] %in% cluster8) {t$cluster8[i] = 1} else {t$cluster8[i] = 0}
  if (t$a[i] %in% hammond9) {t$hammond9[i] = 1}else {t$hammond9[i] = 0}
  if (t$a[i] %in% hammondoa3) {t$hammondoa3[i] = 1}else {t$hammondoa3[i] = 0}
  if (t$a[i] %in% injury) {t$injury[i] = 1}else {t$injury[i] = 0}
  if (t$a[i] %in% dam) {t$dam[i] = 1}else {t$dam[i] = 0}
  if (t$a[i] %in% dam2) {t$dam2[i] = 1}else {t$dam2[i] = 0}
}
t = t[,2:ncol(t)]
t$sums = rowSums(t)
t = t[t$sums>0,]
t = t[order(t$sums),]

#DAM comparison
#fc > 1.5 (lfc > .3), p< 10^-8
#dam signature: 
dam = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/DAM:PAM comparison/mmc3.csv",stringsAsFactors = F)
damgenes = dam$Gene.name
damgenes = damgenes[damgenes%in%rownames(mgAVM02)]
dam = dam[dam$Gene.name %in% damgenes,]
mgAVM02$damgene = PercentageFeatureSet(mgAVM02,features = damgenes)
setEPS()
postscript("~/Desktop/Plots/damgenevln.eps")
VlnPlot(mgAVM02,"damgene",group.by = "finalclusters",pt.size = 0)
dev.off()

#dam genes enrichment vs enrichment in cluster 4
cluster4 = read.csv("/Users/whippoorwill/Desktop/Sequencing/LD_AVM02/Spreadsheets/allmarkers_vargenesMG_0501 copy.csv",stringsAsFactors = F)
cluster4 = cluster4[cluster4$cluster == "4",]
dam$cluster4 = rep("NA",nrow(dam))
for (i in 1:nrow(dam)){
  gene = dam[i,"Gene.name"]
  if (gene %in% cluster4$gene){dam$cluster4[i] = cluster4[cluster4$gene == gene,"avg_logFC"]}
}

write.csv(dam,file = "~/Desktop/damcomparison.csv")
#PAM clusters

pamli = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/DAM:PAM comparison/li_cluster1.csv",stringsAsFactors = F)

pamh = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/DAM:PAM comparison/hammond_cluster4.csv",stringsAsFactors = F)


pamli =pamli$gene
pamh = pamh$Cluster.4.markers

genes = pamli[pamli%in%pamh]
genes
genes = genes[genes%in% rownames(mgAVM02)]

pamli = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/DAM:PAM comparison/li_cluster1.csv",stringsAsFactors = F)
pamh = read.csv("/Users/whippoorwill/Dropbox (Anna Molofsky Lab)/2020-Leah-barrelCortex/Manuscript data/Glia RNASeq/P5_P7 MG rnaseq/DAM:PAM comparison/hammond_cluster4.csv",stringsAsFactors = F)

pam = pamli[pamli$gene %in% genes,]
pam$h = rep("NA",nrow(pam))
pam$cluster4 = rep("NA",nrow(pam))
for (i in 1:nrow(pam)){
  gene = pam[i,"gene"]
  if (gene %in% pamh$Cluster.4.markers){pam$h[i] = pamh[pamh$Cluster.4.markers == gene,"lfc"]}
  if (gene %in% cluster4$gene){pam$cluster4[i] = cluster4[cluster4$gene == gene,"avg_logFC"]}
}
write.csv(pam,file = "~/Desktop/pam_cluster4.csv")
mgAVM02$pamgene = PercentageFeatureSet(mgAVM02,features = genes)
setEPS()
postscript("~/Desktop/Plots/pamgenevln.eps")
VlnPlot(mgAVM02,"pamgene",group.by = "finalclusters",pt.size = 0)
dev.off()

#Hammond cluster 9 (P540 + LPC)

c9 = read.csv("/Users/whippoorwill/Desktop/Paper/Dam:pam comparison/hammondcluster9.csv",stringsAsFactors = F)
c9 = c9$Gene
genes = c9[c9%in%rownames(mgAVM02)]
genes
exclude = grep("-",genes)
exclude = c(exclude,grep("^Rp",genes))
genes = genes[-exclude]
mgAVM02$c9gene = PercentageFeatureSet(mgAVM02,features = genes)
setEPS()
postscript("~/Desktop/Plots/c9genevln.eps")
VlnPlot(mgAVM02,"c9gene",group.by = "finalclusters",pt.size = 0)
dev.off()



#investigate differences when you choose certain genes for the eigengene
set.seed(NULL)
r1 = sample(pamgenes,size = 20,replace = F)
r1
mgAVM02 = PercentageFeatureSet(mgAVM02,features = r1, col.name = "r1")
VlnPlot(mgAVM02,"r1",group.by = "finalclusters",pt.size = 0.001)

all = GetAssayData(mgAVM02, slot = "data")
allpam = all[rownames(all) %in% genes,]
alldam = all[rownames(all) %in% damgenes,]
allpam4 = allpam[,colnames(allpam) %in% cells4]
alldam4 = alldam[,colnames(alldam) %in% cells4]

allpamelse = allpam[,!colnames(allpam) %in% cells4]
alldamelse = alldam[,!colnames(alldam) %in% cells4]

dammatrix = rbind("alldamelse" = Matrix::rowMeans(alldamelse),"alldam4" = Matrix::rowMeans(alldam4))
pammatrix = rbind("allpamelse" = Matrix::rowMeans(allpamelse),"allpam4" = Matrix::rowMeans(allpam4))

damratio = dammatrix[2,]/dammatrix[1,]
pamratio = pammatrix[2,]/pammatrix[1,]

damratio = damratio[order(damratio,decreasing = T)]
pamratio = pamratio[order(pamratio,decreasing = T)]

write.csv(damratio,file = "~/Desktop/damratio.csv")
write.csv(pamratio,file = "~/Desktop/pamratio.csv")