## Performs hierarchical clustering and draw heatmap

## Load libraries
library(plotly)
library(ggplot2)
library(RColorBrewer)
library(heatmaply)
#library(processx)

## Read input matrix with 1row and 1col as headers
matrix <- read.delim("matrix.csv", row.names=1)
rownames(matrix) <- paste("<i>", rownames(matrix), "</i>")
rownames(matrix) <- paste("<b>", rownames(matrix), "</b>")

origin = matrix[, c("Origin")]

#Set color palette
#To show the available palettes use
#display.brewer.all()
ByPal <- colorRampPalette(brewer.pal(6,"RdBu"))

#Run heatmaply on a dataset where the first row is the column header, the first column lists 
#the species and the second column contains the groups by origin of isolation

heatmap = heatmaply(matrix[,-1], fontsize_row = 6, plot_method = "plotly", colorbar_len = 0.2, row_side_palette= ByPal, 
          row_side_colors=data.frame("Origin" = origin, check.names=FALSE))

#To run orca you need to have orca installed and in the path
orca(heatmap, "Heatmap.svg")