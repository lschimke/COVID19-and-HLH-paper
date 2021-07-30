#Packages 
install.packages("ggpubr")
install.packages("ggplot2")
install.packages("ggExtra")

#Libraries
library(ggpubr)
library(ggplot2)
#add histogram:
library(ggExtra) 

#OLFM4 x MMP8
p <- ggplot(tableMLR, aes(x=LCN2,y=DEFA4,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="LCN2", y= "DEFA4") + 
  xlim(7,17)+
  ylim(5, 15)+
  scale_color_manual("Group:", values = c(CONTROL = "#9c9c9c", 
        COVID19_nonICU = "lightgreen", COVID19_ICU = "royalblue"), labels=c(CONTROL = "Control",COVID19_nonICU = "COVID-19_nonICU", COVID19_ICU ="COVID-19_ICU")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "#9c9c9c", 
        COVID19_nonICU = "lightgreen", COVID19_ICU = "royalblue"), labels=c(CONTROL = "Control",COVID19_nonICU = "COVID-19_nonICU", COVID19_ICU ="COVID-19_ICU"))  + 
  geom_smooth(method="lm", label.x = 9, se = FALSE) + theme_classic(base_size = 8) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 7, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)


#vizualize the plot
p









#ELANE x DEFA4
p <- ggplot(tableMLR, aes(x=ELANE,y=DEFA4,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="ELANE", y= "DEFA4") + 
  xlim(0,12)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 0, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p



#ARG1 x CTSG
p <- ggplot(tableMLR, aes(x=ARG1,y=CTSG,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="ARG1", y= "CTSG") + 
  xlim(5,15)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 5, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p




#ARG1 x IL1R1
p <- ggplot(tableMLR, aes(x=ARG1,y=IL1R1,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="ARG1", y= "IL1R1") + 
  xlim(5,15)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 5, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p


#IL10 x MMP8
p <- ggplot(tableMLR, aes(x=IL10,y=MMP8,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="IL10", y= "MMP8") + 
  xlim(0,10)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 0, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p


#CXCL1 x DEFA4
p <- ggplot(tableMLR, aes(x=CXCL1,y=DEFA4,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="CXCL1", y= "DEFA4") + 
  xlim(2,12)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 2, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p


#IL1R1 x LCN2
p <- ggplot(tableMLR, aes(x=IL1R1,y=LCN2,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="IL1R1", y= "LCN2") + 
  xlim(5,15)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 5, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p

#OLFM4 x CEACAM8
p <- ggplot(tableMLR, aes(x=OLFM4,y=CEACAM8,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="OLFM4", y= "CEACAM8") + 
  xlim(5,15)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 5, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p


#DOCK4 x SOD2
p <- ggplot(tableMLR, aes(x=IL1RAP,y=IL1R1,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="IL1RAP", y= "IL1R1") + 
  xlim(8,14)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 0, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p

#CEACAM8 x MMP8
p <- ggplot(tableMLR, aes(x=CEACAM8,y=MMP8,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="CEACAM8", y= "MMP8") + 
  xlim(5,15)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 5, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p

#MGAM x IL1R1
p <- ggplot(tableMLR, aes(x=MGAM,y=IL1R1,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="MGAM", y= "IL1R1") + 
  xlim(5,15)+
  ylim(0,15)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 12, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p

#PTX3 x CTSG
p <- ggplot(tableMLR, aes(x=RETN,y=IL10,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="RETN", y= "IL10") + 
  xlim(5,12)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 5, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p

#CEACAM1 x ELANE
p <- ggplot(tableMLR, aes(x=DEFA4,y=CXCL1,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="DEFA4", y= "CXCL1") + 
  xlim(5,15)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 5, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p

#MMP8 x CTSG
p <- ggplot(tableMLR, aes(x=MMP8,y=CTSG,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="MMP8", y= "CTSG") + 
  xlim(5,15)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 5, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p

#TNFAIP3 x IL10
p <- ggplot(tableMLR, aes(x=TNFAIP3,y=IL10,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="TNFAIP3", y= "IL10") + 
  xlim(7,14)+
  scale_color_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                          COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19")) + 
  scale_fill_manual("Group:", values = c(CONTROL = "goldenrod3", 
                                         COVID19 = "royalblue"), labels=c(CONTROL = "Control",COVID19 = "COVID-19"))  + 
  geom_smooth(method="lm", label.x = 15, se = FALSE) + theme_classic(base_size = 12) +
  theme(legend.position = "bottom") +
  #se a analise de residuo funcionar entao trocar o se = TRUE
  stat_cor(cor.coef.name = "rho", aes(color = Gene.Symbol), label.x = 7, label.y.npc = "top", size = 6, show.legend = FALSE) +
  theme(text = element_text(size = 30))  

#O boxplot mostra a mediana 

#Rho é para R de spearman
#e o p valeu represents the null hypothesis is = 0

#Add: you can choose among "density", "histogram", "boxplot","violin", "densigram"
p <- ggMarginal(p, type="boxplot", groupColour = TRUE, groupFill = TRUE)

#vizualize the plot
p



#Para analise do residuo. 
fit=lm(CEACAM8~ARG1+Gene.Symbol, data = tableMLR)
summary(fit)
plot(fit)

#Ver se o residuo esta adequado, se nao, nao usar o intervalo de confianca para a reta (a sombra ao redor da reta) 
#https://www.originlab.com/doc/origin-help/residual-plot-analysis

# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/

