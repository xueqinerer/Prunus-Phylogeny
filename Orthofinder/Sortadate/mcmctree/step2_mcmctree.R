library("MCMCtreeR")
phy <- readMCMCtree("FigTree.tre", from.file = TRUE)
MCMC.tree.plot(phy, analysis.type = "MCMCtree", cex.tips = 0.6,  add.abs.time=FALSE,
               time.correction = 90, plot.type = "phylogram", lwd.bar = 2, cex.age = 0.4,
               scale.res = c("Eon", "Period", "Epoch", "Age"), node.method = "bar", col.age = "navy", cex.labels = 0.5, relative.height = 0.07, 
               no.margin = TRUE, label.offset = 1, ladderize.tree = FALSE)