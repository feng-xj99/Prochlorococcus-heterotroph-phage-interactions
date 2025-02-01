library(dplyr)

######Building the database
gene_info_1<-read.csv("emapper.annotations.csv",header = TRUE, sep=",")
gene_new<-dplyr::select(gene_info_1,GID=query,Gene_Symbol=Description,G0=GOs,KO=KEGG_ko,Pathway=KEGG_Pathway,gene_name=Preferred_name)
gene_info_new<-dplyr::select(gene_new,GID,Gene_Symbol,gene_name)%>%dplyr::filter(!Gene_Symbol=='-')

#gene_info_pro<-slice(gene_info_1,1:199)
gene_info<-dplyr::select(gene_info_1,query,Description,Preferred_name,GOs,KEGG_ko,KEGG_Pathway,KEGG_Module,CAZy)%>%filter(GOs!="-")
gene_info_pro<-filter(gene_infoGOs!="_")
write.table (gene_info,file ="gene_info.xls", sep ="\t", row.names =F)

#pathway2gene3<-read.table("pathway2gene.xls",header = TRUE,sep ="\t")
#pathway2name<-read.table("pathway2name.xls",header = TRUE,sep ="\t")
library(clusterProfiler)
pathway3gene<-read.table("pathway2gene.xls",header = TRUE,sep ="\t")
colnames(pathway3gene)<-c("GID","Pathway","ko")
pathway2gene<-dplyr::select(pathway3gene,Pathway,GID)
library(dplyr)
pathway<-read.csv("KEGG.profile.csv",header = TRUE, sep=",")
pathway2name<-left_join(pathway3gene,pathway,by=c("Pathway"="Pathway"))%>%dplyr::select(Pathway,Definition)
colnames(pathway2name)<-c("Pathway","Pathway_name")


library(magrittr)
###test
difgene_P1P5<-read.table("P1_vs_P5_DE.tsv",header = TRUE, sep="\t")%>%dplyr::select(gene_id,log2FoldChange)
colnames(difgene_P1P5)<-c("GID","P1_vs_P5")
gene_P1P5<-pull(difgene_P1P5,GID)

de_keggP1_vs_P2<-enricher(gene_P1P5,
                          TERM2GENE = pathway2gene,
                          TERM2NAME = pathway2name,
                          qvalueCutoff=0.05,
                          pvalueCutoff=0.05)

##cycle
library(clusterProfiler)
library(tidyverse)
list<-read.csv("list2.tsv")
list1<-list$list
items_to_remove <- c("P1_vs_P4") #"P2_vs_P4")
filtered_list <- setdiff(list1, items_to_remove)
plots_list3 <- list()
plots_list4 <- list()
library(gridExtra)
for (i in filtered_list) {
  dif_DE<-read.table(paste0(i, "_DE.tsv"),header = TRUE, sep="\t")%>%dplyr::select(gene_id,log2FoldChange)
  colnames(dif_DE) <- c("GID", paste0(i))
  difgene_DE<-pull(dif_DE,GID)
  DE_kegg<-enricher(difgene_DE,
                            TERM2GENE = pathway2gene,
                            TERM2NAME = pathway2name,
                            qvalueCutoff=0.05,
                            pvalueCutoff=0.05)
  DE_kegg_df <- as.data.frame(DE_kegg)
  write.table(DE_kegg_df, file = paste0("de_kegg_", i, "_df.xls"), sep = "\t", row.names = FALSE)
  P <- barplot(DE_kegg,showCategory = 10,main = i)
  plots_list3[[i]] <- P
  
  #表达量图
  colnames(dif_DE)<-c("GID","logFC")
  geneList<-dif_DE$logFC
  names(geneList)<-dif_DE$GID
  geneList<-sort(geneList,decreasing = TRUE)
  
  P1<-clusterProfiler::cnetplot(DE_kegg,
                                color.params = list(foldChange =geneList,edge = TRUE),
                                node_laber="category",
                                category = custom_color_function,
                                showCategory = 5,
                                circular=TRUE) 
  
  plots_list4[[i]] <- P1
  
}  
pdf("SF.9D-I.pdf",width =30, height=40)
big_plot <- grid.arrange(grobs = plots_list3, ncol = 3)
dev.off()  
pdf("kegg_enrich.pdf",width =30, height=40)
big_plot <- grid.arrange(grobs = plots_list4, ncol = 3)
dev.off() 
 


##补充blue/green模块富集
Modules_gene<-read.table('Modules_gene.csv',header = TRUE,sep = ",")
EXP_gene<-read.table('genes.TMM.EXPR.matrix.csv',header = TRUE,sep = ",")
Modules_EXP_info<-left_join(Modules_gene,EXP_gene,by=c('gene'='X'))%>%left_join(gene_info,by=c('gene'='query'))
write.table (Modules_EXP_info,file ="Modules_EXP_info.xls", sep ="\t", row.names =F)

blue_k<-filter(Modules_EXP_info,color=="blue")%>%select('gene','KEGG_ko')
brown_k<-filter(Modules_EXP_info,color=="brown")%>%select('gene','KEGG_ko')
turquoise_k<-filter(Modules_EXP_info,color=="turquoise")%>%select('gene','KEGG_ko')
yellow_k<-filter(Modules_EXP_info,color=="yellow")%>%select('gene','KEGG_ko')

blue_k<-pull(blue_k,gene)
brown_k<-pull(brown_k,gene)
turquoise_k<-pull(turquoise_k,gene)
yellow_k<-pull(yellow_k,gene)

pathway2gene<-read.table("pathway2gene.xls",header = TRUE,sep ="\t")
pathway2name<-read.table("pathway2name.xls",header = TRUE,sep ="\t")
colnames(pathway2gene)<-c("GID","Pathway","ko")
pathway2gene<-dplyr::select(pathway2gene,"Pathway","GID")
colnames(pathway2name)<-c("Pathway","Pathway_name")

pathway3gene<-read.table("pathway2gene.xls",header = TRUE,sep ="\t")
colnames(pathway3gene)<-c("GID","Pathway","ko")
pathway2gene<-dplyr::select(pathway3gene,Pathway,GID)
library(dplyr)
pathway<-read.csv("KEGG.profile.csv",header = TRUE, sep=",")
pathway2name<-left_join(pathway3gene,pathway,by=c("Pathway"="Pathway"))%>%dplyr::select(Pathway,Definition)
colnames(pathway2name)<-c("Pathway","Pathway_name")

library(enrichplot)

blue_k_kegg<-enricher(blue_k,
                         TERM2GENE = pathway2gene,
                         TERM2NAME = pathway2name,
                         qvalueCutoff=0.05,
                         pvalueCutoff=0.05)

brown_k_kegg<-enricher(brown_k,
                                TERM2GENE = pathway2gene,
                                TERM2NAME = pathway2name,
                                qvalueCutoff=0.05,
                                pvalueCutoff=0.05)

turquoise_k_kegg<-enricher(turquoise_k,
                                TERM2GENE = pathway2gene,
                                TERM2NAME = pathway2name,
                                qvalueCutoff=0.05,
                                pvalueCutoff=0.05)
yellow_k_kegg<-enricher(yellow_k,
                                TERM2GENE = pathway2gene,
                                TERM2NAME = pathway2name,
                                qvalueCutoff=0.05,
                                pvalueCutoff=0.05)
blue_k_df<-as.data.frame(blue_k_kegg)
write.table (blue_k_df,file ="blue_k_df.xls", sep ="\t", row.names =F)
library(tidyr)
blue_k_df_gene<-separate_rows(blue_k_df,geneID,sep='/',convert=F)%>%select(ID,Description,geneID)
blue_k_df_gene_des<-left_join(blue_k_df_gene,tpm_modules_node,by=c("geneID"="gene"))
write.table (blue_k_df_gene_des,file ="blue_k_df_gene_des.xls", sep ="\t", row.names =F)

turquoise_k_df<-as.data.frame(turquoise_k_kegg)
write.table (turquoise_k_df,file ="turquoise_k_df.xls", sep ="\t", row.names =F)

p_blue_k<-barplot(blue_k_kegg,showCategory = 10)
p_brown_k<-barplot(brown_k_kegg,showCategory = 10)
p_turquoise_k<-barplot(turquoise_k_kegg,showCategory = 10)
p_yellow_k<-barplot(yellow_k_kegg,showCategory = 10)



pdf("Fig.6CD.pdf",width =16, height=6)
p_blue_k+p_turquoise_k
dev.off()


