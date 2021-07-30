#Packages 
install.packages("ggpubr")
install.packages("ggplot2")
install.packages("ggExtra")

#Libraries
library(ggpubr)
library(ggplot2)
#add histogram:
library(ggExtra) 

#ARG1 x CEACAM8
p <- ggplot(tableMLR, aes(x=ARG1,y=CEACAM8,color=Gene.Symbol, fill=Gene.Symbol)) +
  geom_point(alpha=0.95) + 
  labs(x="ARG1", y= "CEACAM8") + 
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

##############################

#ELANE x DEFA4 #change variables as for your interest
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

################################

#Para analise do residuo. 
fit=lm(CEACAM8~ARG1+Gene.Symbol, data = tableMLR)
summary(fit)
plot(fit)

#Ver se o residuo esta adequado, se nao, nao usar o intervalo de confianca para a reta (a sombra ao redor da reta) 
#https://www.originlab.com/doc/origin-help/residual-plot-analysis

# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/78-perfect-scatter-plots-with-correlation-and-marginal-histograms/

