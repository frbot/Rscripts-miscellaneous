#R script which uses dendextend to make entangled trees
#Francesca Bottacini
#f.bottacini@umail.ucc.ie


# load required packages
library("ape")
library("dendextend")
#library("dendextendRcpp")

# read in tree files
tree1 <- read.tree("tree1.nwk")
tree2 <- read.tree("tree2.nwk")

# if your trees are not ultrametric, use this and apply denextend command on ctree1 and ctree2
# if your trees are already ultrametric, use tree1 and tree2, you can skip this block
# lambda=1 indicates clock like diversification, a rough ultrametricization method
ctree1 <- chronos(tree1, lambda=1)
ctree2 <- chronos(tree2, lambda=1)

# get measure of similary based on Robinson-Roulds distance (= symmetric distance)
entanglement(ctree1,ctree2)

# create entanglement plot, pdf dimensions are in inches. lwd is tangle line width, cex the taxon label size
#pdf("Dendextend_entanglement.pdf", width=7, height=7)
#tanglegram(ctree1,ctree2, lwd=2, cex_main = 0.3, common_subtrees_color_branches = TRUE)
#dev.off()

## Colour the tree
labels = tree1 %>% labels
labels <- as.data.frame(labels)

meta <- read.delim("tree-metadata.txt")

labels2 <- merge(labels, meta, by.x="labels", sort=F)

#Colour based on origin of isolation
#cols <- as.character(labels2$Colour)

tanglegram(ctree1,ctree2, lwd=2, cex_main = 0.3, margin_inner = 6, color_lines = cols, highlight_branches_lwd = FALSE)

