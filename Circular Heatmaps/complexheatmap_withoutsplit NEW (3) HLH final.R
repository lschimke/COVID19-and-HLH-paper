#Install and access packages
install.packages("BiocManager")
library(BiocManager)
BiocManager::install("ComplexHeatmap")
install.packages("circlize")
library(circlize)
library(ComplexHeatmap)


#before import your dataset run the heatmap normalize (i.e. by log2, Z score, inputation) the data. You can use perseus.

#PLOT RETANGULAR HEATMAP
Heatmap(DEGs, 
        name="Gene Expression", 
        row_title = "Genes",
        row_title_side = "right",
        row_title_gp = gpar(fontsize = 14, fontface = "bold"),
        row_names_gp = gpar(fontsize= 8),
        row_title_rot = 270,
        row_dend_side = "left", 
        row_dend_width = unit(1, "cm"),
        column_title = "Samples",
        column_title_side = "bottom",  
        column_title_gp = gpar(fontsize = 14, fontface = "bold"),
        column_names_gp = gpar(fontsize= 8),
        column_dend_side = "top",
        column_dend_height = unit(1, "cm"))


#CIRCULAR HEATMAP WITH BOXPLOT
circos.clear()

DEGs = as.matrix(DEGs)

column_od = hclust(dist(t(DEGs)), method = "single")$order

#Ajust here according the expression (counts) range
col_fun1 = colorRamp2(c(-0, 2.15, 4.0), c("white", "green", "royal blue"))

circos.par(gap.after = 20) #gap t the end to include sample names
circos.heatmap(DEGs[,column_od], col= col_fun1, 
               rownames.cex = 1.5, #gene size
               rownames.font = par("font"),
               rownames.side = "outside",
               cluster = FALSE,
               track.height=0.4)

circos.track(ylim = range(DEGs), panel.fun = function(x, y) {
        m = DEGs[CELL_META$subset, 1:9, drop = FALSE]
        m = m[CELL_META$row_order, , drop = FALSE]
        n = nrow(m)
        
        # circos.boxplot is applied on matrix columns, so here we transpose it.
        circos.boxplot(t(m), pos = 1:n - 0.5, pch = 16, cex = 0.3)
        circos.lines(CELL_META$cell.xlim, c(0, 0), lty = 2, col = "grey")
}, cell.padding = c(0, 1, 2, 4))


circos.track(track.index = 2, panel.fun = function(x, y) {
        if(CELL_META$sector.numeric.index == 1) { # the last sector
                cn = colnames(DEGs[,column_od])
                n = length(cn)
                circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 
                            1:n - 0.5, cn, 
                            cex = 0.4, adj = c(0, 0.5), facing = "inside")
        }
}, bg.border = TRUE)

lgd = Legend(title = "-Log10 p.adj", col_fun = col_fun1)
grid.draw(lgd)

circos.clear()


 
