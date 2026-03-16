# Open FigTree.tre in FigTree software and export as NEXUS format
# Save the exported file as FigTree_NEXUS.tre

# Multiply branch lengths by 100 and force the tree to be ultrametric
library(ape)
library("phytools")
tree <- read.nexus("FigTree_NEXUS.tre", force.multi = TRUE)
tree2 <- tree[[1]]
tree2$edge.length <- tree2$edge.length * 100
tree3 <- force.ultrametric(tree2)
is.ultrametric(tree3)
write.tree(tree3, "Prunus.nw.ultra.tre")
