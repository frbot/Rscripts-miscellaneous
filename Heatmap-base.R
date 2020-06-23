## Performs hierarchical clustering and draw heatmap

## Load libraries
library(plotly)
library(ggplot2)
library(RColorBrewer)
library(heatmaply)



## Read input matrix with 1row and 1col as headers

matrix <- read.delim("matrix.txt", row.names=1)

#rownames(matrix) <- paste("<b>", rownames(matrix), "</b>")
#colnames(matrix) <- paste("<b>", colnames(matrix), "</b>")

heatmaply(matrix, fontsize_row = 8, plot_method = "plotly", colorbar_len = 0.1, colors = PRGn)