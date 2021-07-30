#FROM: https://www.datanovia.com/en/blog/how-to-create-a-beautiful-plots-in-r-with-summary-statistics-labels/

#BOXPLOTS with values:
# p values to your data 
install.packages("ggpubr")
install.packages("lemon")

library(ggpubr)
library(lemon)

table <- read.delim # file.txt of your datable
table$gene <- as.factor(table$gene)

#Choose display order

# option 1: manually added order for 2 groups
table$group <- factor(table$group, levels=c("CONTROL","COVID_ICU"))
#table's original order
table$gene <- factor(table$gene, levels=unique(table$gene))

# option 2: for 3 groups or more, just increase/include the group names as below
table$group <- factor(table$group, levels=c("CONTROL_nonICU","CONTROL_ICU", "COVID_nonICU","COVID_ICU"))
table$gene <- factor(table$gene, levels=unique(table$gene))

################################
# pallete can be either a list of colors (c("cornflowerblue", "darkgoldenrod", "brown")), or a pallete ("Dark2")
# color names at: www.datanovia.com/en/blog/awesome-list-of-657-r-color-names/

#OPTION 1: TWO GROUP COMPARISON 
#INSIDE BOX. FROM: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/#compare-two-paired-samples
ggboxplot(table, x = "group", y = "counts",
          color = "group", palette = "Dark2",
          add = c("jitter", "median_iqr"),
          facet.by = "gene", short.panel.labs = TRUE)+ 
  stat_compare_means(method= "t.test", label = "p.format")


#if you want to change colors see codes at: http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually

#OUTSIDE BOX
# Create a box plot
# ref https://www.datanovia.com/en/blog/how-to-create-a-beautiful-plots-in-r-with-summary-statistics-labels/
# pallete can be either a list of colors (c("cornflowerblue", "darkgoldenrod", "brown")), or a pallete ("Dark2")

ggboxplot(table, x = "gene", y = "counts",
           color = "group", palette = c("cornflowerblue", "darkgoldenrod", "brown"),
           add = c("jitter", "median_iqr"))+ 
   stat_compare_means(method= "t.test", label = "p.format", 
                      aes(group = group))
 
 
 # Change plot type to violin
 ggviolin(table, x = "gene", y = "counts", 
   add = c("jitter", "median_iqr"),
   color = "group", palette = "Dark2")+
stat_compare_means(method= "t.test", label = "p.format", aes(group = group))
 

 
 
 #OPTION 2: MULTIPLE GROUPS COMPARISON 
#FROM: http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/76-add-p-values-and-significance-levels-to-ggplots/#compare-two-paired-samples
 
 my_comparisons <- list( c("CONTROL_nonICU", "CONTROL_ICU"), c("CONTROL_nonICU", "COVID_nonICU"), c ("CONTROL_nonICU", "COVID_ICU"), c("CONTROL_ICU", "COVID_nonICU"), c("CONTROL_ICU", "COVID_ICU"), c("COVID_nonICU", "COVID_ICU"))
 
 
 ggboxplot(table, x = "group", y = "counts",
           color = "group", palette = c("#ff7373", "#eaac43", "#38c3e1", "#1657e1"),
           add = c("jitter", "median_iqr"))+ 
    
  scale_y_continuous(breaks=c(-5, 0, 5, 10, 15, 20, 25, 30, 35), name="counts")+
    expand_limits(y = c(-5, 25), size=3)+
    facet_rep_wrap(~gene, strip.position="bottom", ncol=5)+ #if you have more than 2 lines, use ( facet_rep_wrap(~gene, strip.position="bottom") ) instead of facet_Wrap
    stat_compare_means(method= "t.test", label = "p.format", label.y = c(18, 21, 24, 27, 30, 33), size = 3, #location of statistical bars and font size
                       aes(group = group), comparisons = my_comparisons)+
        theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_line(),
          strip.background = element_rect(color="white", fill="white"),
          strip.placement = "outside",
          panel.border = element_blank(), 
          panel.background = element_blank(), 
          panel.grid = element_blank(), 
          panel.spacing.x = unit(2,"line"))
   
 
 
