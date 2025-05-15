rm(list = ls())
options(stringsAsFactors = F)

library(ggplot2)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(dplyr)

##########################################BA_1#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_BA1_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('1873','1477','141','354','1438','1200')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)


theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black",  aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,700)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group A")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm'))+
  coord_cartesian(ylim = c(1000,4000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group B vs. Group A", face = "bold", color = "black", size =14))
dev.off()


##########################################BA_2#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_BA2_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('1681','1140','178','93','1239','1337')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,500)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group A2")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5))+
  coord_cartesian(ylim = c(1000,4000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group B vs. Group A2", face = "bold", color = "black", size =14))
dev.off()

##########################################BA_3#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_BA3_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('700','411','473','311','980','947')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,550)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group A3")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5))+
  coord_cartesian(ylim = c(600,2200)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group B vs. Group A3", face = "bold", color = "black", size =14))
dev.off()
##########################################D1C#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_D1C_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('2378','2694','113','36','1160','1238')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)


theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black",  aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,200)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group A")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm'))+
  coord_cartesian(ylim = c(1000,6000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group D vs. Group C", face = "bold", color = "black", size =14))
dev.off()


##########################################D2C#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_D2C_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('577','1135','15','99','1145','1438')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)


theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black",  aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,200)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group A")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm'))+
  coord_cartesian(ylim = c(500,3500)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group D2 vs. Group C", face = "bold", color = "black", size =14))
dev.off()

##########################################D3C#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_D3C_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('64','95','236','197','2323','2846')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)


theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black",  aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,600)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group A")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.title = element_blank(),
        legend.key.size = unit(0.3, 'cm'))+
  coord_cartesian(ylim = c(2000,6000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group D3 vs. Group C", face = "bold", color = "black", size =14))
dev.off()

##########################################cA_1#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_CA1_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('2540','2594','113','222','1411','1312')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,500)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group C vs. Group A")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5))+
  coord_cartesian(ylim = c(1000,6000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group C vs. Group A", face = "bold", color = "black", size =14))
dev.off()
##########################################cA_2#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_CA2_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('311','177','189','163','2122','1376')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,1000)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group C vs. Group A2")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5))+
  coord_cartesian(ylim = c(1200,4000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group C vs. Group A2", face = "bold", color = "black", size =14))
dev.off()
##########################################cA_3#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_CA3_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('327','283','960','225','2947','2586')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,2000)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group C vs. Group A3")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5))+
  coord_cartesian(ylim = c(2000,6000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group C vs. Group A3", face = "bold", color = "black", size =14))
dev.off()
##########################################BD_1#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_BD1_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('2665','2509','152','102','897','660')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,500)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group D")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5))+
  coord_cartesian(ylim = c(1000,6000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group B vs. Group D", face = "bold", color = "black", size =14))
dev.off()
##########################################BD_2#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_BD2_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('1894','1616','233','119','430','623')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,700)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group D2")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5))+
  coord_cartesian(ylim = c(1000,5000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group B vs. Group D2", face = "bold", color = "black", size =14))
dev.off()
##########################################BD_3#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_BD3_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('485','351','1293','395','725','809')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))


fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  scale_y_continuous(breaks=c(0,250,500,750,1000,1250,1500,1750,2000))+ 
  # geom_bar(stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group B vs. Group D3")+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(face="bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,2000)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

annotate_figure(fig1, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group B vs. Group D3", face = "bold", color = "black", size =14))
dev.off()
##########################################A1E#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_A1E_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('5072','4688','754','1733','2188','1777')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,800)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

fig2 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group A vs. Group E")+
  labs(x=NULL,y=NULL) +  
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
  theme(axis.text.y = element_text(face="bold",colour = "black", size = 9),
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        legend.title = element_blank(),
        plot.title = element_text(face="bold", hjust = 0.5))+
  coord_cartesian(ylim = c(1500,12000)) +  
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

figure<- ggarrange(fig2,fig1,heights=c(2/3, 1/3),ncol = 1, nrow = 2,  
                   common.legend = T,legend="top",align = "v")
annotate_figure(figure, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group A vs. Group E", face = "bold", color = "black", size =14))
dev.off()


##########################################A2E#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_A2E_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('267','409','65','193','265','145')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_bar(stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group A2 vs. Group E")+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(face="bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,800)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

annotate_figure(fig1, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group A2 vs. Group E", face = "bold", color = "black", size =14))
dev.off()

##########################################A3E#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_A3E_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('1025','1298','1219','1024','1020','1429')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_bar(stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group A3 vs. Group E")+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(face="bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,3000)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

annotate_figure(fig1, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group A3 vs. Group E", face = "bold", color = "black", size =14))
dev.off()

##########################################FE#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_FE_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('411','237','89','67','196','154')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_bar(stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group F vs. Group E")+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(face="bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,700)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

annotate_figure(fig1, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group F vs. Group E", face = "bold", color = "black", size =14))
dev.off()

##########################################EG#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_EG_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('227','206','80','28','993','700')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  # geom_bar(stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group F vs. Group E")+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(face="bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank())+
  # coord_cartesian(ylim = c(0,1800), breaks = seq(0, 1800,200)) + 
  scale_y_continuous(limits = c(0, 1800), breaks = seq(0, 1800, 200))+
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

annotate_figure(fig1, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group F vs. Group E", face = "bold", color = "black", size =14))
dev.off()


##########################################FG#######################################
setwd("C:/Users/lgs96/OneDrive - Virginia Tech/Desktop/figures")
pdf(file = "barplot_FG_key.pdf",width = 3,height = 5)
specie <- c(rep("Enhancers",2),rep("Promoters",2),rep("DEGs",2))
condition <- rep(c("increase" , "decrease"), 3)
value <- c('238','222','258','390','574','648')
value <- as.numeric(value)
data <- data.frame(specie,condition,value)

theme <- theme(panel.background = element_blank(),panel.border=element_rect(fill=NA,linetype = 1, color="black",size=0.5))

fig1 <- ggplot(data,aes(fill=condition, y=value, x=specie, label= value)) + 
  geom_bar(position = position_stack(reverse = TRUE), stat="identity", colour="black", size = 0.5, aes(fill = condition)) +
  scale_y_continuous(breaks=c(0,200,400,600,800,1000,1200))+ 
  # geom_bar(stat="identity", colour="black", size = 0.5, aes(fill = condition)) + 
  # geom_text(size = 3, position = position_stack(vjust = 0.5))+
  # ggtitle("Group F vs. Group G")+
  labs(x=NULL,y=NULL)+    
  scale_fill_manual(values = c("#0c84c6","#CF583D"))+
  theme(axis.text.x = element_text(face="bold", colour = "black", angle=45 ,hjust = 1, size =9), 
        axis.text.y = element_text(face="bold", colour = "black",size =9),
        legend.text = element_text(size=9), 
        legend.key = element_rect(colour = "black", size= 0.05),
        legend.key.size = unit(0.3, 'cm'),
        plot.title = element_text(face="bold", hjust = 0.5),
        legend.position = "top",
        legend.title = element_blank())+
  coord_cartesian(ylim = c(0,1300)) + 
  guides(fill = guide_legend(reverse = TRUE))+
  scale_x_discrete(limits = c("Enhancers", "Promoters", "DEGs"))+
  theme

annotate_figure(fig1, left = text_grob("Differential Expressed Peaks/Genes", face = "bold", color = "black", rot = 90 , size =12),
                top = text_grob("Group F vs. Group G", face = "bold", color = "black", size =14))
dev.off()
