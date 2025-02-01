library(MicrobiotaProcess)
library(phyloseq)
library(tidyverse)
library(RColorBrewer)
setwd("/Users/fengxuejin/Desktop/app_16s")
otu <- "filter_table_rename.qza"
rep <- "rep-seqs_raname.qza"
tree <- "rooted-tree-no-mc.qza"
tax <- "taxonomy-rename.qza"
sample <- "metadata2.tsv"

ps_dada2 <- import_qiime2(otuqza=otu,taxaqza=tax,refseqqza=rep,mapfilename=sample,treeqza=tree)

ps <- ps_dada2

ps
observed <- estimate_richness(ps, measures = c('Observed'))
explore.df <- cbind(observed, sample_sums(ps),sample_data(ps)$Title)
colnames(explore.df) <- c('Observed', 'Sample_Sums',"Title")
observed_mean <- mean(explore.df$Observed)
sample_sum_mean <- mean(explore.df$Sample_Sums)
ggplot(data = explore.df, aes(x = Sample_Sums, y = Observed,color=Title)) + 
  geom_point() +
  geom_smooth(method="auto", se=TRUE, fullrange=FALSE, level=0.95, 
              inherit.aes = F, mapping = aes(Sample_Sums, Observed),
              data = explore.df)

#2_Rarefaction
scales_y <- list(
  ACE = scale_y_continuous(limits = c(50, 162), breaks = seq(50, 162, 30)),
  Chao1 = scale_y_continuous(limits = c(50, 160), breaks = seq(50, 160, 30)),
  Observe = scale_y_continuous(limits = c(50, 155), breaks = seq(50, 155, 30)),
  Shannon = scale_y_continuous(limits = c(1.7, 2.5), breaks = seq(1.7, 2.6, 0.25)),
  Simpson = scale_y_continuous(limits = c(0.10, 0.40), breaks = seq(0.10, 0.40, 0.1))
)

p_rare <- ggrarecurve(obj=ps_dada2, indexNames=c("Observe","Chao1","ACE","Shannon", "Simpson"),
                      chunks=800) +
  theme(legend.spacing.y=unit(0.02,"cm"),legend.text=element_text(size=4))+
  theme_bw()
ggsave (p_rare, filename="p_rare-new-1.pdf", width=12, height=4)

#3_Alpha_diversity in SF.7

ps_dada2@sam_data[["Title"]]<-factor(ps_dada2@sam_data[["Title"]],levels = c("P", "PC","PR", "PA","PAC", "PAR"))
alphaobj <- get_alphaindex(ps_dada2)
head(as.data.frame(alphaobj))
p_alpha <- ggbox(alphaobj, geom="violin", factorNames="Title",testmethod = "t.test",signifmap = TRUE) + 
  scale_fill_manual(values=c("#6BB658","#ECC847","#005493","#6BB658","#ECC847","#005493"))+
  theme(strip.background = element_rect(colour=NA, fill="grey"))

ggsave (p_alpha, filename="SF.7.pdf", width=14, height=4)

 #4_composition
phytax <- get_taxadf(obj=ps_dada2,taxlevel=6)
p <- phyloseq::otu_table(phytax) %>% as.data.frame()
write.table (p,file ="phylumt_6-new-1.xls", sep ="\t",col.names = NA)
phytax <- get_taxadf(obj=ps_dada2,taxlevel=2)
p <- phyloseq::otu_table(phytax) %>% as.data.frame()
write.table (p,file ="phylumt_2-new-1.xls", sep ="\t",col.names = NA)


#<1% to other
library(magrittr)
computed_persent <- function(path) {
  data <- path %>%
    read.delim(check.names = FALSE, row.names = 1)
  data2 <- data %>%
    mutate(sum = rowSums(.), persent = sum / sum(sum) * 100, 
           sum = NULL,) %>%
    rbind(filter(., persent < 0.1) %>% colSums()) %>%
    mutate(ASV_ID = c(data %>% rownames(), "others"))
  filter(data2[1:(nrow(data2) - 1),], persent > 0.1) %>%
    rbind(data2[nrow(data2),]) %>%
    select(ncol(.), 1:(ncol(.) - 2)) %>%
    set_rownames(seq_len(nrow(.))) %>%
    return()
}

path <- "phylumt_2-new-1-12.xls"
phylumt_6_other  <- computed_persent(path)
write.table (phylumt_6_other ,file ="phylumt_2_other-new-1-12.xls", sep ="\t", row.names =F) 


###total_genus in SF.5
library(reshape2)

phylumt_6_other <-read.table("phylumt_6_other-new-1.xls.xls",header = TRUE)%>% melt()
meta <- "metadata2.tsv" %>% read.delim()
phylumt_6_other<-left_join(phylumt_6_other,meta,by=c("variable"="sample.id"))
write.table (phylumt_6_other,file ="com_6-new-1.xls", sep ="\t", row.names =F) 

col=c("#6BB658","#EC5158","#FEC868","#ADBED2","#ECD59F","#F0C8DC","#E8CAA4","#93c47d","#F0E8E2","#ffd966","#EDBC7A","#FEE08B","#F46D43","#7FCDBB","#D3E7EE","#FFB6B9","#89A4C8","#F7E3AF","#A6C48A","#F3AB9D","#FFD3B6","#E5D67B","#E5EDCF","#B6D7A8","#F9B5AC","#A8D5E2","#F9F9C5","#7F7F7F")
phylumt_6_other$Title<-factor(phylumt_6_other$Title,levels = c("P","PC","PR","PA","PAC","PAR"))
phylumt_6_other$ASV_ID<-factor(phylumt_6_other$ASV_ID,levels = c("g__g__g__Prochlorococcus","g__g__g__Alteromonas","g__g__g__Psychroserpens","g__g__g__Marinobacter","g__g__g__Henriciella","g__g__g__Hwanghaeella","g__g__g__Maricaulis","g__g__g__Altererythrobacter","g__g__g__Amorphus","g__g__g__Aquicoccus","g__g__g__Lentilitoribacter","g__g__g__Limnobacter","g__g__g__Magnetospira","g__g__g__Methylophaga","g__g__g__Minwuia","g__g__g__Muricauda","g__g__g__Nitratireductor","g__g__g__Oryzicola","g__g__g__Parahaliea","g__g__g__Poriferisphaera","g__g__g__Rhodophyticola","g__g__g__Roseibium","g__g__g__Roseitalea","g__g__g__Roseovarius","g__g__g__Thalassospira","g__g__g__Tsuneonella","g__g__un_f__f__Erythrobacteraceae","others"))

phylumt_6_other_mean<- summarySE(phylumt_6_other, measurevar="value", groupvars=c("ASV_ID","Title","Time"))

p_com <- ggplot(phylumt_6_other, aes(x = factor(Time), y = value, fill = ASV_ID)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Title, scales = "free", space = "free_x") +
  labs(x = "", y = "Proportions") +
  labs(fill = "") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = col)+
  theme_bw()
ggsave (p_com, filename="p_com-time.pdf", width=15, height=5)

###12_genuns in SF.6A

phylumt_6_other_12 <-read.table("phylumt_6-new-1-12.xls",header = TRUE)%>% melt()
meta <- "metadata2.tsv" %>% read.delim()
phylumt_6_other_12<-left_join(phylumt_6_other_12,meta,by=c("variable"="sample.id"))
write.table (phylumt_6_other_12,file ="com_6-new-12.xls", sep ="\t", row.names =F) 

col=c("#6BB658","#EC5158","#FEC868","#ADBED2","#ECD59F","#F0C8DC","#E8CAA4","#93c47d","#F0E8E2","#ffd966","#EDBC7A","#F46D43","#7FCDBB","#D3E7EE","#FFB6B9","#89A4C8","#F7E3AF","#A6C48A","#F3AB9D","#FFD3B6","#E5D67B","#E5EDCF","#A8D5E2","#7F7F7F")
phylumt_6_other_12$Title<-factor(phylumt_6_other_12$Title,levels = c("P","PC","PR","PA","PAC","PAR"))
phylumt_6_other_12$Name<-factor(phylumt_6_other_12$Name,levels = c("P11_2","P12_2","P13_2","P21_2","P22_2","P23_2","P21_5","P22_5","P23_5","P41_2","P42_2","P43_2","P51_2","P52_2","P53_2","P31_2","P32_2","P33_2","P31_5","P32_5","P33_5","P61_2","P62_2","P63_2"))
phylumt_6_other_12$ASV_ID<-factor(phylumt_6_other_12$ASV_ID,levels = c("g__g__g__Prochlorococcus","g__g__g__Alteromonas","g__g__g__Psychroserpens","g__g__g__Marinobacter","g__g__g__Henriciella","g__g__g__Hwanghaeella","g__g__g__Maricaulis","g__g__g__Altererythrobacter","g__g__g__Amorphus","g__g__g__Aquicoccus","g__g__g__Lentilitoribacter","g__g__g__Magnetospira","g__g__g__Methylophaga","g__g__g__Minwuia","g__g__g__Muricauda","g__g__g__Nitratireductor","g__g__g__Oryzicola","g__g__g__Parahaliea","g__g__g__Poriferisphaera","g__g__g__Rhodophyticola","g__g__g__Roseibium","g__g__g__Roseitalea","g__g__g__Tsuneonella","others"))

p_com <- ggplot(phylumt_6_other_12, aes(x = factor(Name), y = value, fill = ASV_ID)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(~Title, scales = "free", space = "free_x") +
  labs(x = "", y = "Proportions") +
  labs(fill = "") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = col)+
  theme_bw()
ggsave (p_com, filename="SF.6A.pdf", width=10, height=5)

#Flattening
library(vegan)
library(reshape2)
phylumt_6<-read.csv("phylumt_6_other-new-1.xls.xls",header = TRUE,sep = "\t",row.names=1)
colSums(phylumt_6)
sum_phylumt_6<-colSums(phylumt_6)
phylumt_6_Flattening = as.data.frame(t(rrarefy(t(phylumt_6), min(sum_phylumt_6))))
colSums(phylumt_6_Flattening)
write.table (phylumt_6_Flattening,file ="phylumt_6_Flattening.xls", sep ="\t", col.names = NA) 
