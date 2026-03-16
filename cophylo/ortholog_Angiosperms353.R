setwd("C:/Users/301/Desktop/李属文章数据/tree/ortholg_vs_353/Coalescence")
library(ape)
library(phytools)
list.files()
orto_1to1 <- read.tree("./Prunus_Astral_BS10_treeshrink_1to1.tre")
orto_1to1 $edge.length <- NULL
plotTree(orto_1to1 ,use.edge.length=FALSE)
nodelabels()
orto_1to1_122 <- ape::rotate(orto_1to1, 122)
orto_1to1_123 <- ape::rotate(orto_1to1_122, 123)
orto_1to1_133 <- ape::rotate(orto_1to1_123, 133)
plot(orto_1to1_133)

MI <- read.tree("./Prunus_Astral_BS10_treeshrink_MI.tre")
MI$edge.length <- NULL
plotTree(MI,use.edge.length=FALSE)
nodelabels()
MI_122 <- ape::rotate(MI, 122)
MI_123 <- ape::rotate(MI_122, 123)
MI_133 <- ape::rotate(MI_123, 133)
plot(MI_133)

MO <- read.tree("./Prunus_Astral_BS10_treeshrink_MO.tre")
MO$edge.length <- NULL
plotTree(MO,use.edge.length=FALSE)
nodelabels()
MO_122 <- ape::rotate(MO, 122)
MO_123 <- ape::rotate(MO_122, 123)
MO_133 <- ape::rotate(MO_123, 133)
plot(MO_133)

Ang_353<- read.tree("Prunus_353_Astral_BS10.tre")
Ang_353$edge.length <- NULL
plotTree(Ang_353,use.edge.length=FALSE)
nodelabels()
Ang_353_163 <- ape::rotate(Ang_353, 163)
plot(Ang_353_163)

obj11 <- cophylo(Ang_353, MO_133, rotate = TRUE)
obj12 <- cophylo(Ang_353, MI_133,rotate = TRUE)
obj13 <- cophylo(Ang_353, orto_1to1_133, rotate = TRUE)

pdf("Angiosperms353&Ortholog_MO_MI_1to1.pdf", height = 50, width = 34)
par(mfrow=c(2,2))
plot.cophylo(obj11, link.type = "curved", 
             link.lty = "solid", fsize=2, link.lwd=1, tip.lty=2, tip.lwd=1.0, lwd=1,
             tip.len=0, pts=TRUE, link.col=make.transparent("blue", 0.9))
nodelabels.cophylo(text=obj11$trees[[1]]$node.label,frame="none",adj=c(1,-0.4),cex=1.5, which="left", col="black")
nodelabels.cophylo(text=obj11$trees[[2]]$node.label,frame="none",adj=c(-0.3,-0.4),cex=1.5, which="right", col="black")
mtext("(a)",side=3,line=-4, at=par("usr")[1]+0.05*diff(par("usr")[1:2]), cex=2)
mtext("Angiosperms353", cex=3, side=2, line=-3, col="red")
mtext("MO(5477 gene)", cex=3, side=4, line=-1.5, col="red")

plot.cophylo(obj12, link.type = "curved", 
             link.lty = "solid", fsize=2, link.lwd=1, tip.lty=2, tip.lwd=1.0, lwd=1,
             tip.len=0, pts=TRUE, link.col=make.transparent("blue", 0.9))
nodelabels.cophylo(text=obj12$trees[[1]]$node.label,frame="none",adj=c(1,-0.4),cex=1.5, which="left", col="black")
nodelabels.cophylo(text=obj12$trees[[2]]$node.label,frame="none",adj=c(-0.3,-0.4),cex=1.5, which="right", col="black")
mtext("(b)",side=3,line=-4, at=par("usr")[1]+0.05*diff(par("usr")[1:2]), cex=2)
mtext("Angiosperms353", cex=3, side=2, line=-3, col="red")
mtext("MI(8316 gene)", cex=3, side=4, line=-1.5, col="red")

plot.cophylo(obj13, link.type = "curved", 
             link.lty = "solid", fsize=2, link.lwd=1, tip.lty=2, tip.lwd=1.0, lwd=1,
             tip.len=0, pts=TRUE, link.col=make.transparent("blue", 0.9))
nodelabels.cophylo(text=obj13$trees[[1]]$node.label,frame="none",adj=c(1,-0.4),cex=1.5, which="left", col="black")
nodelabels.cophylo(text=obj13$trees[[2]]$node.label,frame="none",adj=c(-0.3,-0.4),cex=1.5, which="right", col="black")
mtext("(c)",side=3,line=-4, at=par("usr")[1]+0.05*diff(par("usr")[1:2]), cex=2)
mtext("Angiosperms353", cex=3, side=2, line=-3, col="red")
mtext("1to1(3736 gene)", cex=3, side=4, line=-1.5, col="red")

dev.off()


