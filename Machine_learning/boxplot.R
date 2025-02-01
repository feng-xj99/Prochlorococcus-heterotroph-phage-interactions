library(ggplot2)
library(ggpubr)
library(dplyr)
library(RColorBrewer)


data<-read.delim('PCP.txt',header=T,sep = '\t',stringsAsFactors = F, check.names = F)

# First, we need to convert the dataframe from wide to long format
data_long <- reshape2::melt(data, id.vars = c("name","group"))
data_long_filter<-filter(data_long,variable=="Marinobacter"|
                           variable=="Psychroserpens"|
                           variable=="Prochlorococcus"|
                           variable=="Maricaulis"|
                           variable=="others"|
                           variable=="Parahaliea"|
                           variable=="Limnobacter"|
                           variable=="Roseibium"|
                           variable=="Poriferisphaera"|
                           variable=="Aquicoccus"|
                           variable=="Minwuia"|
                           variable=="Magnetospira")

data_long_filter$value1<-data_long_filter$value/29324*100
data_long_filter$variable<-factor(data_long_filter$variable,levels = unique(data_long_filter$variable))
data_long_filter$variable<-factor(data_long_filter$variable,levels = c("Marinobacter","Psychroserpens","Prochlorococcus","Maricaulis","others","Parahaliea","Limnobacter","Roseibium","Poriferisphaera","Aquicoccus","Minwuia","Magnetospira"))

# Run the Wilcoxon test for each variable
stat.test1 <- compare_means(
  value1 ~ group, data = data_long_filter, group.by = "variable",
  method = "wilcox.test")

p <- ggplot(data_long_filter, aes(x = variable, y = value1)) +
  stat_boxplot(aes(fill = factor(group)), geom = 'errorbar', width = 0.8, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(fill = factor(group)), color = "black", position = position_dodge(width = 0.8), size = 1, shape = 21, alpha = 0.4) +
  geom_boxplot(aes(fill = factor(group)), position = position_dodge(width=0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("0" = "#6BB658", "1" = "#ECC847"), labels = c("No Cyanophage", "Cyanophage")) +
  theme_test()  +
  labs(y = "Relative Abundance(%)", x = "")  +
  theme(
    axis.text.x = element_text(size = 10, color = 'black', angle = 45, hjust = 1, face = "italic"),  # italics and bold
    axis.text.y = element_text(size = 10,  color = 'black', face="bold"),  # Normal and bold, note that it should be "bold" instead of "blod"
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.justification = "right",
    legend.position = "right") +
  stat_pvalue_manual(stat.test1, x = "variable", label = "p.signif", y.position = 80, position = position_dodge(width = 1),size=5)+
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) - 0.4, xend = as.numeric(variable) + 0.4, y = 78, yend = 78))+
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) - 0.4, xend = as.numeric(variable) - 0.4, y = 78, yend = 76)) +
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) + 0.4, xend = as.numeric(variable) + 0.4, y = 78, yend = 76))

ggsave('PCP.pdf',width=7,height=7, plot = p)

data<-read.delim('PACPA.txt',header=T,sep = '\t',stringsAsFactors = F, check.names = F)

# First, we need to convert the dataframe from wide to long format
data_long <- reshape2::melt(data, id.vars = c("name","group"))
data_long_filter<-filter(data_long,variable=="Psychroserpens"|
                           variable=="Prochlorococcus"|
                           variable=="Altererythrobacter"|
                           variable=="Henriciella"|
                           variable=="Tsuneonella"|
                           variable=="others")

data_long_filter$value1<-data_long_filter$value/29324*100
data_long_filter$variable<-factor(data_long_filter$variable,levels = unique(data_long_filter$variable))
data_long_filter$variable<-factor(data_long_filter$variable,levels = c("Psychroserpens","Prochlorococcus","Altererythrobacter","Henriciella","Tsuneonella","others"))
# Run the Wilcoxon test for each variable
stat.test1 <- compare_means(
  value1 ~ group, data = data_long_filter, group.by = "variable",
  method = "wilcox.test")

p <- ggplot(data_long_filter, aes(x = variable, y = value1)) +
  stat_boxplot(aes(fill = factor(group)), geom = 'errorbar', width = 0.8, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(fill = factor(group)), color = "black", position = position_dodge(width = 0.8), size = 1, shape = 21, alpha = 0.4) +
  geom_boxplot(aes(fill = factor(group)), position = position_dodge(width=0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("0" = "#6BB658", "1" = "#ECC847"), labels = c("No Cyanophage", "Cyanophage")) +
  theme_test()  +
  labs(y = "Relative Abundance(%)", x = "")  +
  theme(
    axis.text.x = element_text(size = 10, color = 'black', angle = 45, hjust = 1, face = "italic"),  # italics and bold
    axis.text.y = element_text(size = 10,  color = 'black', face="bold"),  # Normal and bold, note that it should be "bold" instead of "blod"
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.justification = "right",
    legend.position = "right") +
  stat_pvalue_manual(stat.test1, x = "variable", label = "p.signif", y.position = 80, position = position_dodge(width = 1),size=5)+
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) - 0.4, xend = as.numeric(variable) + 0.4, y = 78, yend = 78))+
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) - 0.4, xend = as.numeric(variable) - 0.4, y = 78, yend = 76)) +
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) + 0.4, xend = as.numeric(variable) + 0.4, y = 78, yend = 76))

ggsave('PACPA.pdf',width=7,height=7, plot = p)


data<-read.delim('PARPA.txt',header=T,sep = '\t',stringsAsFactors = F, check.names = F)

# First, we need to convert the dataframe from wide to long format
data_long <- reshape2::melt(data, id.vars = c("name","group"))
data_long_filter<-filter(data_long,variable=="Marinobacter"|
                           variable=="Tsuneonella"|
                           variable=="Thalassospira")

data_long_filter$value1<-data_long_filter$value/29324*100
data_long_filter$variable<-factor(data_long_filter$variable,levels = unique(data_long_filter$variable))
data_long_filter$variable<-factor(data_long_filter$variable,levels = c("Marinobacter","Tsuneonella","Thalassospira"))
# Run the Wilcoxon test for each variable
stat.test1 <- compare_means(
  value1 ~ group, data = data_long_filter, group.by = "variable",
  method = "wilcox.test")

p <- ggplot(data_long_filter, aes(x = variable, y = value1)) +
  stat_boxplot(aes(fill = factor(group)), geom = 'errorbar', width = 0.8, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(fill = factor(group)), color = "black", position = position_dodge(width = 0.8), size = 1, shape = 21, alpha = 0.4) +
  geom_boxplot(aes(fill = factor(group)), position = position_dodge(width=0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("0" = "#6BB658", "1" = "#005493"), labels = c("No Alterophage", "Alterophage")) +
  theme_test()  +
  labs(y = "Relative Abundance(%)", x = "")  +
  theme(
    axis.text.x = element_text(size = 10, color = 'black', angle = 45, hjust = 1, face = "italic"),  # italics and bold
    axis.text.y = element_text(size = 10,  color = 'black', face="bold"),  # Normal and bold, note that it should be "bold" instead of "blod"
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.justification = "right",
    legend.position = "right") +
  stat_pvalue_manual(stat.test1, x = "variable", label = "p.signif", y.position = 5, position = position_dodge(width = 1),size=5)+
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) - 0.4, xend = as.numeric(variable) + 0.4, y = 5, yend = 5))+
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) - 0.4, xend = as.numeric(variable) - 0.4, y = 5, yend = 4.8)) +
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) + 0.4, xend = as.numeric(variable) + 0.4, y = 5, yend = 4.8))

ggsave('PARPA.pdf',width=7,height=7, plot = p)





data<-read.delim('PAP.txt',header=T,sep = '\t',stringsAsFactors = F, check.names = F)

# First, we need to convert the dataframe from wide to long format
data_long <- reshape2::melt(data, id.vars = c("name","group"))
data_long_filter<-filter(data_long,variable=="Alteromonas"|
                           variable=="Marinobacter")

data_long_filter$value1<-data_long_filter$value/29324*100
data_long_filter$variable<-factor(data_long_filter$variable,levels = unique(data_long_filter$variable))
data_long_filter$variable<-factor(data_long_filter$variable,levels = c("Alteromonas","Marinobacter"))
# Run the Wilcoxon test for each variable
stat.test1 <- compare_means(
  value1 ~ group, data = data_long_filter, group.by = "variable",
  method = "wilcox.test")

p <- ggplot(data_long_filter, aes(x = variable, y = value1)) +
  stat_boxplot(aes(fill = factor(group)), geom = 'errorbar', width = 0.8, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(fill = factor(group)), color = "black", position = position_dodge(width = 0.8), size = 1, shape = 21, alpha = 0.4) +
  geom_boxplot(aes(fill = factor(group)), position = position_dodge(width=0.8), width = 0.6, outlier.shape = NA) +
  scale_fill_manual(values = c("0" = "#6BB658", "1" = "#767171"), labels = c("No Alteromonas", "Alteromonas")) +
  theme_test()  +
  labs(y = "Relative Abundance(%)", x = "")  +
  theme(
    axis.text.x = element_text(size = 10, color = 'black', angle = 45, hjust = 1, face = "italic"),  # italics and bold
    axis.text.y = element_text(size = 10,  color = 'black', face="bold"),  # Normal and bold, note that it should be "bold" instead of "blod"
    panel.border = element_rect(colour = "black", fill = NA, size = 1),
    legend.justification = "right",
    legend.position = "right") +
  stat_pvalue_manual(stat.test1, x = "variable", label = "p.signif", y.position = 41, position = position_dodge(width = 1),size=5)+
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) - 0.4, xend = as.numeric(variable) + 0.4, y = 40, yend = 40))+
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) - 0.4, xend = as.numeric(variable) - 0.4, y = 40, yend = 39)) +
  geom_segment(data = stat.test1, aes(x = as.numeric(variable) + 0.4, xend = as.numeric(variable) + 0.4, y = 40, yend = 39))+
  scale_y_continuous(limits = c(0, 50)) 
ggsave('PAP.pdf',width=7,height=7, plot = p)


