#FROM: https://www.r-graph-gallery.com/115-study-correlations-with-a-correlogram.html

install.packages("corrgram")
install.packages("ggplot")
install.packages("psych")
install.packages("inlmisc")

library(corrgram)
library(readxl)  #to import and export tables 
library(ggplot2)
library(psych)
library(inlmisc)

Cor2 <- read_excel("/Users/iplace/Documents/R scripts/R scripts for COVID19 and HLH paper/Correlograms/Cor2.xlsx")# add your directory where your store your datafile table
View(Cor2)

graphic = corrgram(Cor2, order=TRUE, lower.panel=panel.fill, upper.panel=panel.pie, 
         text.panel=panel.txt, main="Gene expression", cex.labels = 1,
         col = colorRampPalette(c("white", "light green", "royalblue")))


AddGradientLegend(
  breaks = seq(-1,1,0.1),
  pal = colorRampPalette(c("white", "light green", "royalblue")),
  n = 3,
  labels = TRUE,
  strip.dim = c(1,8),
  loc = "right",
  inset = c(0.04,0.1))

#########################
#other graphic options to plot the correlogram

corrgram(Cor2, order=TRUE, lower.panel=panel.ellipse, upper.panel=panel.pts, 
         text.panel=panel.txt, main="Gene expression", 
         col = colorRampPalette(c("green", "white", "royalblue")))


corrgram(Cor2, order=TRUE, lower.panel=panel.fill, upper.panel=panel.ellipse, 
         text.panel=panel.txt, main="Gene expression", cex.labels = 2,
         col = colorRampPalette(c("green", "white", "royalblue")))

corrgram(Cor2, order=TRUE, lower.panel=panel.shade, upper.panel=panel.pts, 
         text.panel=panel.txt, main="Gene expression", cex.labels = 2, cl.lim = c(-1, 1), 
         col = colorRampPalette(c("green", "white", "royalblue")))

corrgram(Cor2, order=TRUE, lower.panel=panel.conf, upper.panel=panel.pts, 
         text.panel=panel.txt, main="Gene expression", cl.lim = c(-1, 1), 
         col = colorRampPalette(c("green", "white", "royalblue")))
                                                                                                                          


#Inserir scala do grafico?
#mudar de auadrado para bolas?
#adicionar histogram


# Quick display of two cabapilities of GGally, to assess the distribution and correlation of variables
library(GGally)

# From the help page:
data(flea)
ggpairs(flea, columns = 2:4, ggplot2::aes(colour=species))






