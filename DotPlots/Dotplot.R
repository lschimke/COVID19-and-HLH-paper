install.packages("tidyverse")
library(tidyverse)
library(viridis)
library(ggplot2)


df <- NEWDOTPLOt#insert your datatable 

df$Term <- factor(df$Term, levels = df$Term[order(df$Combined.Score)])

theme_dotplot <- theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(1.0)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank())
        

# create the plot
df <- df %>% filter(Adjusted.P.value < 0.05)
ggplot(df, aes(x = Combined.Score, y = Term, color=Adjusted.P.value)) +
  geom_point(mapping = aes(size = Combined.Score)) +
 # scale_x_continuous(limits = c(0, 70)) +
  theme_dotplot+
  xlab ("Combined.Score") +
  ylab("GO cathegory") +
  ggtitle("DotPlot by GO")+
  scale_color_gradient(low="green", high="royalblue")

