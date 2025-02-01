setwd("/Users/xj/Desktop/Alter-pro-phage/fig_r")
library(Rmisc)
abundance<- read.csv("abundance.tsv",header = TRUE,sep = '\t')
abundance$Pro <- as.numeric(abundance$Pro)
#Pro<-na.omit(Pro) 整个表格中
abundance_pro<- subset(abundance,Pro!="NA")
Pro<- summarySE(abundance_pro, measurevar="Pro", groupvars=c("type","pro_type","time"))

abundance$Total <- as.numeric(abundance$Total)
abundance_Total<- subset(abundance,Total!="NA")
Total<- summarySE(abundance_Total, measurevar="Total", groupvars=c("type","pro_type","time"))

abundance$Bac <- as.numeric(abundance$Bac)
abundance_Bac<- subset(abundance,Bac!="NA")
Bac<- summarySE(abundance_Bac, measurevar="Bac", groupvars=c("type","pro_type","time"))

abundance$Pt <- as.numeric(abundance$Pt)
abundance_Pt<- subset(abundance,Pt!="NA")
Pt<- summarySE(abundance_Pt, measurevar="Pt", groupvars=c("type","pro_type","time"))

Pro_ID<-paste(Pro$pro_type,Pro$type,Pro$time)
Pro_ID<-gsub(" ", "_",Pro_ID,fixed = TRUE)
Pro$Pro_ID<-Pro_ID
Pro$parameter<-paste("pro")

Total_ID<-paste(Total$pro_type,Total$type,Total$time)
Total_ID<-gsub(" ", "_",Total_ID,fixed = TRUE)
Total$Total_ID<-Total_ID
Total$parameter<-paste("Total")

Bac_ID<-paste(Bac$pro_type,Bac$type,Bac$time)
Bac_ID<-gsub(" ", "_",Bac_ID,fixed = TRUE)
Bac$Bac_ID<-Bac_ID
Bac$parameter<-paste("Bac")

Pt_ID<-paste(Pt$pro_type,Pt$type,Pt$time)
Pt_ID<-gsub(" ", "_",Pt_ID,fixed = TRUE)
Pt$Pt_ID<-Pt_ID
Pt$parameter<-paste("pt")

colnames(Pt) <- c("type", "pro_type", "time", "N", "value","sd", "se","ci","ID","parameter")
colnames(Pro) <- c("type", "pro_type", "time", "N", "value","sd", "se","ci","ID","parameter")
colnames(Bac) <- c("type", "pro_type", "time", "N", "value","sd", "se","ci","ID","parameter")
colnames(Total) <- c("type", "pro_type", "time", "N", "value","sd", "se","ci","ID","parameter")

final<-rbind(Pro,Bac)%>%rbind(Total)%>%rbind(Pt)

final_ID_type<-paste(final$pro_type,final$type)
final_ID_type<-gsub(" ", "_",final_ID_type,fixed = TRUE)
final$ID_type<-final_ID_type



write.table(final,"final.tsv",quote = FALSE,row.names = FALSE, col.names = TRUE,sep = "\t")



final<- read.csv("final-1.tsv",header = TRUE,sep = "\t")

final_ID_parameter<-paste(final$pro_type,final$parameter)
final_ID_parameter<-gsub(" ", "_",final_ID_parameter,fixed = TRUE)
final$ID_parameter<-final_ID_parameter
final$parameter <- factor(final$parameter, levels=c("Total","pro","Bac","pt"))

final$parameter <- factor(final$parameter, levels=c("Total","pro","Bac","pt"))
final$type <- factor(final$type, levels=c("P", "PC","PR", "PA","PAC", "PAR"))
final$ID_parameter <- factor(final$ID_parameter, levels=c("NATL1A_Total","NATL1A_pro","NATL1A_Bac","NATL1A_pt"))

p_2<-ggplot(final, aes(x=time, y=value,fill=type,color=type,shape=type)) + 
  geom_line(size=0.75) +
  geom_point(size = 3,stroke =0.75)+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=0.5,color = "gray50") +
  scale_color_manual(values=c("#6BB658","#ECC847","#005493","#6BB658","#ECC847","#005493"))+
  scale_fill_manual(values=c("#6BB658","#ECC847","#005493","#6BB658","#ECC847","#005493"))+
  scale_shape_manual(values = c(1,2,6,16,17,25))+
  facet_wrap((~ID_parameter),nrow = 2,scales = "free")
ggsave (p_2, filename="Fig.1.pdf", width=8, height=6)


final$ID_type<-factor(final$ID_type, levels=c("NATL1A_P","NATL1A_PC","NATL1A_PR","NATL1A_PA","NATL1A_PAC","NATL1A_PAR"))
final$parameter1 <- factor(final$parameter1, levels=c("Total","pro","Bac","pt_c","pt_r"))
p_1<-ggplot(final, aes(x=time, y=value,fill=parameter1,color=parameter1,shape=type)) + 
  geom_line(size=0.75) +
  geom_point(size = 3,stroke =0.75)+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=1.5,size=0.5,color = "gray30") +
  scale_color_manual(values=c("#767171","#6BB658","#EC5158","#ECC847","#005493"))+
  scale_fill_manual(values=c("#767171","#6BB658","#EC5B4C","#ECC847","#005493"))+
  scale_shape_manual(values = c(1,2,6,16,17,25))+
  facet_wrap((~ID_type),nrow = 2,scale="free_y")
ggsave (p_1, filename="SF.3.pdf", width=15, height=8)


abundance<- read.csv("abundance.tsv",header = TRUE,sep = '\t')
abundance$vbr2 <- as.numeric(abundance$vbr2)
abundance_vbr2<- subset(abundance,vbr2!="NA")
vbr2<- summarySE(abundance_vbr2, measurevar="vbr2", groupvars=c("type","pro_type","time"))
vbr2$type <- factor(vbr2$type, levels=c("PC","PR", "PAC", "PAR"))
vbr2$type2 <- ifelse(vbr2$type %in% c("PC", "PAC"), "pt_c", 
                   ifelse(vbr2$type %in% c("PR", "PAR"), "pt_r", NA))
p_vbr2<-ggplot(vbr2, aes(x=time, y=vbr2,fill=type,color=type,shape=type)) + 
  geom_line(size=0.75) +
  geom_point(size = 3,stroke =0.75)+
  geom_errorbar(aes(ymin=vbr2-se, ymax=vbr2+se), width=1.5,size=0.5,color = "gray50") +
  scale_color_manual(values=c("#ECC847","#005493","#ECC847","#005493"))+
  scale_fill_manual(values=c("#ECC847","#005493","#ECC847","#005493"))+
  scale_shape_manual(values = c(2,6,17,25))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_wrap((~type2))
ggsave (p_vbr2, filename="Fig.2.pdf", width=8.5, height=4)

library(ggpubr)
library(reshape2)
cy_ssc_3<- read.csv("cy_ssc_final.tsv",header = TRUE,sep="\t")
cy_ssc_3$type2<-factor(cy_ssc_3$type2,levels = c("NATL1A_max","NATL1A_u1","NATL1A_PerCP.Cy5.5.A","NATL1A_SSC.A"))
cy_ssc_3$type2.1<-factor(cy_ssc_3$type2.1,levels = c("NATL1A_max","NATL1A_u1","NATL1A_PerCP.Cy5.5.A","NATL1A_SSC.A"))
my_comparisons <- list(c("P", "PC"),  c("PA", "PAC"), c("P", "PR"), c("PA", "PAR"),c("P", "PA"), c("PC", "PAC"),c("PR", "PAR"))
my_comparisons_1 <- list(c("PR", "PAR"), c("PC", "PAC"))
cy_ssc_3<- subset(cy_ssc_3,value!="NA")
cy_ssc_3$type1 <- factor(cy_ssc_3$type1, levels=c("P", "PC","PR", "PA","PAC","PAR"))
p_cy_ssc_3<-ggplot(cy_ssc_3,aes(x=type1,y=value,fill=type1,color=type1,shape=type1))+
  geom_point(aes(fill=type1,color=type1),size=2.5)+
  scale_fill_manual(values=c("#FFFFFF","#FFFFFF","#FFFFFF","#6BB658","#ECC847","#005493"))+
  scale_color_manual(values=c("#6BB658","#ECC847","#005493","#6BB658","#ECC847","#005493"))+
  scale_shape_manual(values = c(1,1,1,1,16,16,16,16))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  stat_compare_means(aes(label=..p.signif..),method = "t.test",comparisons=my_comparisons)+
  facet_wrap(~type2.1,scales = "free",nrow =2)
ggsave (p_cy_ssc_3, filename="SF.4.pdf", width=9, height=8)


### Reproducibility in SF.2 
library(Rmisc)
abundance<- read.csv("abundance-new-2.csv",header = TRUE)

abundance$Pro <- as.numeric(abundance$Pro)
#Pro<-na.omit(Pro) 
abundance_pro<- subset(abundance,Pro!="NA")
Pro<- summarySE(abundance_pro, measurevar="Pro", groupvars=c("type","pro_type","time"))
Pro$type <- factor(Pro$type, levels=c("P", "PC","PR", "PA","PAC", "PAR"))
p_pro<-ggplot(Pro, aes(x=time, y=Pro,fill=pro_type,color=pro_type,shape=pro_type)) + 
  geom_line(size=0.75) +
  geom_point(size = 3,stroke =0.75)+
  geom_errorbar(aes(ymin=Pro-se, ymax=Pro+se), width=2,size=0.5, color = "black") +
  scale_color_manual(values=c("#005493","#ECC847"))+
  scale_fill_manual(values=c("#005493","#ECC847"))+
  scale_shape_manual(values = c(1,2,6,16,17,25))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_wrap(~type, ncol = 1)
ggsave (p_pro, filename="SF.2-2.pdf", width=5, height=18)


abundance$Total <- as.numeric(abundance$Total)
abundance_Total<- subset(abundance,Total!="NA")
Total<- summarySE(abundance_Total, measurevar="Total", groupvars=c("type","pro_type","time"))
Total$type <- factor(Total$type, levels=c("P", "PC","PR", "PA","PAC", "PAR"))
p_Total<-ggplot(Total, aes(x=time, y=Total,fill=pro_type,color=pro_type,shape=pro_type)) + 
  geom_line(size=0.75) +
  geom_point(size = 3,stroke =0.75)+
  geom_errorbar(aes(ymin=Total-se, ymax=Total+se), width=2,size=0.5, color = "black") +
  scale_color_manual(values=c("#005493","#ECC847"))+
  scale_fill_manual(values=c("#005493","#ECC847"))+
  scale_shape_manual(values = c(2,1))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_wrap(~type, ncol = 1)
ggsave (p_Total, filename="SF.2-1.pdf", width=5, height=18)

abundance$Pt <- as.numeric(abundance$Pt)
abundance_Pt<- subset(abundance,Pt!="NA")
Pt<- summarySE(abundance_Pt, measurevar="Pt", groupvars=c("type","pro_type","time"))
Pt$type <- factor(Pt$type, levels=c("PC","PR", "PAC", "PAR"))
p_Pt<-ggplot(Pt, aes(x=time, y=Pt,fill=pro_type,color=pro_type,shape=pro_type)) + 
  geom_line(size=0.75) +
  geom_point(size = 3,stroke =0.75)+
  geom_errorbar(aes(ymin=Pt-se, ymax=Pt+se), width=2,size=0.5, color = "black") +
  scale_color_manual(values=c("#005493","#ECC847"))+
  scale_fill_manual(values=c("#005493","#ECC847"))+
  scale_shape_manual(values = c(2,1))+
  theme_bw()+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  facet_wrap(~type, ncol = 1)
ggsave (p_Pt, filename="SF.2-3.pdf", width=5, height=12)
