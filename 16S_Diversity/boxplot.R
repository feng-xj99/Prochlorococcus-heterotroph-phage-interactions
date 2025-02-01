library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)

different<-read.csv("different.csv",header = TRUE,sep = ",")
different$type<-factor(different$type,levels=c("PCP","PACPA","PARPA","PAP"))

p <- ggplot(different, aes(x = log2FoldChange,y = reorder(X, log2FoldChange),  fill = log2FoldChange > 0)) + 
  geom_bar(position = position_dodge(0.5), stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("FALSE" = "navy", "TRUE" = "firebrick3"),labels = c("log2FC < 0", "log2FC > 0")) +  # 设置蓝色和红色
  facet_wrap(~type, scales = "free", nrow = 2) + 
  theme_bw()
ggsave (p, filename="p_diff.pdf", width=10, height=7)
