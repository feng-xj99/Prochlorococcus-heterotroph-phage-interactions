load(file = 'info.rdata')
gene_matrix<- read.csv("genes.TMM.EXPR.matrix.csv",header = TRUE,sep=",")
library(dplyr)
library(tibble)
gene_info<-read.csv("emapper.annotations.csv",header = TRUE,sep=",")
#gene_exp_1_name<-rownames_to_column(gene_exp_1,var='gene_id')
gene_matrix_info<-left_join(gene_matrix,gene_info,by=c("X"="query"))
gene_matrix_info_filter<-filter(gene_matrix_info,evalue!="")
gene_exp_matrix<-select(gene_matrix_info_filter, 1:22)%>%column_to_rownames(var='X')%>%select("P11","P12", "P13", "P21", "P22", "P23","P31", "P32" ,"P33", "P71", "P72", "P73", "P51", "P52", "P53", "P41", "P42", "P43", "P61", "P62", "P63")%>%filter()

gene_exp_matrix<-column_to_rownames(gene_matrix,var='X')
gene_exp_matrix_log<-log10(gene_exp_matrix+1)
#dim(total_matrix_t)
#BiocManager::install("WGCNA")
#BiocManager::install(c("AnnotationDbi","impute","GO.db","preprocessGore"))
library(WGCNA)
library(reshape2)

###检查缺失值太多的基因样本
gene_exp_matrix_log_t<-as.data.frame(t(gene_exp_matrix_log))
#gene_exp_log<-log10(gene_exp_1+1)
#gene_exp_1_t<-as.data.frame(t(gene_exp_log))
gsg_total=goodSamplesGenes(gene_exp_matrix_log_t,verbose = 3)
gsg_total$allOK
if(!gsg_total$allOK)
{if(sum(!gsg_total$goodGenes)>0)
  printFlush(paste("Removing genes",paste(names(gene_exp_matrix_log_t)[!gsg_total$goodGenes],collapse = ",")));
  if(sum(!gsg_total$goodSamples)>0)
    printFlush(paste("Removingsamples",paste(rownames(gene_exp_matrix_log_t)[!gsg_total$goodSamples],collapse = ",")));
  gene_exp_matrix_log_t=gene_exp_matrix_log_t[gsg_total$goodGenes,gsg_total$goodSamples]
}
#gene_exp_1_t
gene_exp_matrix_log_t<-gene_exp_matrix_log_t[1:21,]

#gene_exp_matrix_log_t <-select(gene_exp_matrix_log_t,-'R8W_16',-'R8W_26',-'R8W_20')

sampleTree=hclust(dist(gene_exp_matrix_log_t),method = "complete")
sizeGrWindow(12,9)
par(cex=0.6)
par(mar=c(0,4,2,0))
plot(sampleTree,main = "Sample clustering to detectoutliers",sub="",xlab = "",cex.lab=1.5, cex.axis=1.5,cex.main=2)

# Drawing Threshold Cut Lines
abline(h = 120,col="red"); 
# Determining clusters
clust = cutreeStatic(sampleTree, cutHeight =120, minSize = 10)

table(clust)	
keepSamples = (clust==1)	# Clust 1 into keepSamples
dat_total_matrix_t = gene_exp_matrix_log_t[keepSamples, ]		

nGenes =ncol(dat_total_matrix_t)
nSamples =nrow(dat_total_matrix_t)

# Setting the soft threshold tuning parameter range
powers =c(c(1:10),seq(from = 12, to=30,by=2))
# Network Analysis
sft = pickSoftThreshold(dat_total_matrix_t, powerVector = powers, verbose = 5)
# Drawing
sizeGrWindow(9, 5)	

par(mfrow =c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="SoftThreshold(power)",ylab="ScaleFreeTopologyModelFit,signedR^2",type="n",
     main =paste("Scaleindependence"));

text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.20,col="red")
# Mean Connectivity and SoftThreshold
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="SoftThreshold(power)",ylab="MeanConnectivity", type="n",
     main =paste("Meanconnectivity"))

text(sft$fitIndices[,1], sft$fitIndices[,5],labels=powers, cex=cex1,col="red")

sft$powerEstimate
net = blockwiseModules(dat_total_matrix_t,power= 8,	
                       TOMType ="unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase ="femaleMouseTOM",
                       maxBlockSize=6008,
                       verbose = 3)
table(net$colors)

# Visualization modules
sizeGrWindow(12, 9)

mergedColors = labels2colors(net$colors)

pdf("Fig.6A.pdf", width=12, height=6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Modulecolors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file="FemaleLiver-02-networkConstruction-auto.RData")

lnames =load(file="FemaleLiver-02-networkConstruction-auto.RData");
lnames


MEs0 = moduleEigengenes(dat_total_matrix_t, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
MEDiss = 1-cor(MEs); 
METree = hclust(as.dist(MEDiss), method = "complete"); 
pdf("Module correlation coefficient.pdf")
plotEigengeneNetworks(MEs, "Eigengene adjacency heatmap",
                      marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE,
                      xLabelsAngle = 90)
plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
abline(h=0.99, col="red")
dev.off()


##Here we extract the genes for each module.

gene.names<-data.frame(moduleLabels)
gene.names$color<-moduleColors
gene.names<-rownames_to_column(gene.names,var = "gene")
sort_gene.names<-gene.names[order(gene.names[,"color"]),]
write.csv(sort_gene.names,"Modules_gene.csv")

library(flashClust)
#hier1=flashClust(as.dist(MEDiss), method='complete')

#set the diagonal of the dissimilarity to NA
#diag(MEDiss) = NA
#pdf('Modules_heatmap.pdf')
#TOMplot(METree, hier1, as.character(mergedColors))
#TOMplot(plotTOM, net$dendrograms, mergedColors, main = "Network heatmap plot, all genes")
#write.csv(MEDiss,"dissTOM.txt")
library(tibble)
sample_info<-read.csv("samples.csv",header = TRUE,sep=",")
traitData<-select(sample_info,-no,-type,-lg_pro,-lg_total,-lg_bac,-lg_phage_C,-lg_phage_R,-lg_vbr)%>%column_to_rownames(var='name')
moduleTraitCor =cor(MEs, traitData, use ="p");	# 计算相关性系数
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)	# 计算P


pdf("Sample dendrogram and trait heatmap.pdf",width =12, height=8)
par(mar =c(6, 15, 3, 3))
traitColors = numbers2colors(traitData, signed = FALSE);

plotDendroAndColors(sampleTree, traitColors,
                    groupLabels =names(traitData),
                    main ="Sample dendrogram and trait heatmap")

dev.off()

pdf("Fig.6B.pdf",width =10, height=8)
par(mar =c(6, 15, 3, 3))

textMatrix <- paste(
  signif(moduleTraitCor, 2),
  "\n(",
  signif(moduleTraitPvalue, 1), ")",
  sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = names(traitData),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  zlim = c(-1,1),
  main = paste("Module-trait relationships"))
dev.off()

##Network visualization
#TOM = TOMsimilarityFromExpr(dat_total_matrix_t, power=12,networkType = 'signed',TOMType = 'signed')
load(net$TOMFiles[1], verbose=T)

## Loading objects:
##   TOM


TOM <- as.matrix(TOM)


dissTOM = 1-TOM
# Transform dissTOM with a power to make moderately strong 
# connections more visible in the heatmap

plotTOM = dissTOM^7

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA
# Call the plot function

pdf("HP.pdf",width =20, height=20)
par(mar =c(3, 3, 3, 3))
TOMplot(plotTOM, net$dendrograms, moduleColors,col=gplots::colorpanel(250,'red','orange','lemonchiffon'),
        main = "Network heatmap plot, all genes")
dev.off()
##turquoise
modules = c("turquoise")
probes = names (dat_total_matrix_t)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes [ inModule];
modTOM = TOM[ inModule,inModule] ;
dimnames (modTOM) = list (modProbes,modProbes)
cyt = exportNetworkToCytoscape (
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-",paste(modules,collapse="-"),
                   ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-",paste(modules,collapse="-"),
                   ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.785,
  nodeNames = modProbes ,
  nodeAttr = moduleColors[inModule]
)
##blue
modules = c("blue")
probes = names (dat_total_matrix_t)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes [ inModule];
modTOM = TOM[ inModule,inModule] ;
dimnames (modTOM) = list (modProbes,modProbes)
cyt = exportNetworkToCytoscape (
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-",paste(modules,collapse="-"),
                   ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-",paste(modules,collapse="-"),
                   ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.75,
  nodeNames = modProbes ,
  nodeAttr = moduleColors[inModule]
)
##yellow  
modules = c("yellow")
probes = names (dat_total_matrix_t)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes [ inModule];
modTOM = TOM[ inModule,inModule] ;
dimnames (modTOM) = list (modProbes,modProbes)
cyt = exportNetworkToCytoscape (
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-",paste(modules,collapse="-"),
                   ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-",paste(modules,collapse="-"),
                   ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.2,
  nodeNames = modProbes ,
  nodeAttr = moduleColors[inModule]
)
##brown
modules = c("brown")
probes = names (dat_total_matrix_t)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes [ inModule];
modTOM = TOM[ inModule,inModule] ;
dimnames (modTOM) = list (modProbes,modProbes)
cyt = exportNetworkToCytoscape (
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-",paste(modules,collapse="-"),
                   ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-",paste(modules,collapse="-"),
                   ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.3,
  nodeNames = modProbes ,
  nodeAttr = moduleColors[inModule]
)


ADJ=abs(cor(dat_total_matrix_t,use="p"))^6

Alldegrees =intramodularConnectivity(ADJ, moduleColors)
write.csv(Alldegrees, file = "intramodularConnectivity.csv")

# Plot the correlation scatterplot of GS and connectivity:

nSamples = nrow(dat_total_matrix_t)

traitData_names <- colnames(traitData)

for (j in traitData_names) {
  j_df <- as.data.frame(traitData[[j]])
  names(j_df) <- "[[j]]"
  modNames = substring(names(MEs), 3)
  # Calculate the P-value of MM
  geneModuleMembership = as.data.frame(cor(dat_total_matrix_t, MEs, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
  names(geneModuleMembership) = paste("MM", modNames, sep="")
  names(MMPvalue) = paste("p.MM", modNames, sep="")
  geneTraitSignificance = as.data.frame(cor(dat_total_matrix_t,j_df, use = "p"))
  GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
  names(geneTraitSignificance) = paste("GS.", j, sep="")
  names(GSPvalue) = paste("p.GS.",j, sep="")
  
  colorlevels=unique(moduleColors)
  items_to_remove <- c("grey")
  colorlevels <- setdiff(colorlevels, items_to_remove)
  pdf(paste0("GS vs. degree_",j,".pdf"),width = 14,height = 8)

  par(mfrow=c(3,4))
  par(mar = c(5,5,3,3))
  for (i in c(1:length(colorlevels))) 
  {
    whichmodule=colorlevels[[i]]; 
    restrict1 = (moduleColors==whichmodule);
    verboseScatterplot(Alldegrees$kWithin[restrict1], 
                       geneTraitSignificance[restrict1,1], 
                       col=moduleColors[restrict1],
                       main=whichmodule, 
                       xlab = "Connectivity", ylab = "Gene Significance", abline = TRUE)
  }

   dev.off()
}


####hub gene
hub<-read.csv("hub.csv",header = TRUE,sep = ",")
gene_info<-read.csv("emapper.annotations.csv",header = TRUE,sep = ",")
tpm_hub<-left_join(hub,gene_matrix,by=c("gene_id"="X"))%>%left_join(gene_info,by=c("gene_id"="query"))%>%select(gene_id,Preferred_name,Description,3:24)
%>%select(-affinity,-GO,Preferred_name,-KEGG_Pathway)

tpm_hub$name<-paste(tpm_hub$gene_id,"_",tpm_hub$Preferred_name)
library(reshape2)
sample_type<-select(sample_info,name,type,day)
sample_type$type1<-paste(sample_type$type,"_",sample_type$day)
tpm_hub_plot<-melt(tpm_hub,id.vars =c("gene_id","color","Description","name"))%>%left_join(sample_type,by=c("variable"="name"))%>%slice(61:1320)
library(Rmisc)
tpm_hub_plot$value <- as.numeric(tpm_hub_plot$value)
tpm_hub_plot_1_bubble<- summarySE(tpm_hub_plot, measurevar="value", groupvars=c("name","type1","Description","color"))
tpm_hub_plot_1_bubble$lg<-log10(tpm_hub_plot_1_bubble$value+1)
tpm_hub_plot_1_bubble$lg <- ifelse(tpm_hub_plot_1_bubble$lg == 0, NA, tpm_hub_plot_1_bubble$lg)

library(ggplot2)
tpm_hub_plot_1_bubble$type1<-factor(tpm_hub_plot_1_bubble$type1,levels=c("P _ 12", "PR _ 12", "PA _ 12", "PAC _ 12", "PAR _ 12",  "PC _ 42","PAC _ 42"))
p_tpm_hub_plot_1_bubble<-ggplot(tpm_hub_plot_1_bubble,aes(x=type1,y=name))+
  scale_fill_gradientn(colours = colorRampPalette(c("navy","white","firebrick3"))(100))+
  theme_bw()+geom_point(aes(size=`lg`, fill =`lg`),shape=21)+
  theme(panel.grid=element_line(colour ="#696969", size=0.2,linetype=3), axis.text.x =element_text(angle =45))+xlab(NULL) + ylab(NULL)+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.8, linetype="solid"))+
  facet_wrap(~color,nrow = 1,scales = 'free')

ggsave (p_tpm_hub_plot_1_bubble, filename="SF.12.pdf", width=20, height=6)


