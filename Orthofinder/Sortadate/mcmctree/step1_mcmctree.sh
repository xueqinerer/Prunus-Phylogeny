# Step 1: Open FigTree.tre in FigTree software and export as NEXUS format
# Save the exported file as FigTree_NEXUS.tre

# Step 2: Reroot the MO tree and adjust topology if needed
# (adjust node rotations using the R script below)
library(ape)
library(phytools)

orto_MO <- read.tree("Prunus_Astral_BS10_treeshrink_MO.tre")
orto_MO$edge.length <- NULL
plotTree(orto_MO, use.edge.length=FALSE)
nodelabels()
orto_MO_122 <- ape::rotate(orto_MO, 122)
orto_MO_123 <- ape::rotate(orto_MO_122, 123)
orto_MO_133 <- ape::rotate(orto_MO_123, 133)
plot(orto_MO_133)
orto_MO_109 <- ape::rotate(orto_MO_133, 109)

plot(orto_MO_109)

write.tree(orto_MO_109, "MO.tre")

nw_topology -I MO.tre > MO_NEW.tre
# After confirming topology and rerooting, update MO_NEW.tre
# Rename abbreviated taxon names to full names:
# sed -i 's/Prunus_arborea_var/Prunus_arborea_var._montana/g;s/Prunus_davidiana_var/Prunus_davidiana_var._potaninii/g;s/Prunus_salicina_var/Prunus_salicina_var._mandshurica/g;s/glandulifolia/maackiiNCBI/g' MO_NEW.tre

python /path/to/software/AMAS/AMAS.py convert -d dna -f fasta -i *.fasta -u phylip
cat *phy > ../Prunus_mcmctree/Prunus.phylip
cd ../Prunus_mcmctree
mkdir Hessian
mcmctree mcmctree.ctl
mkdir approx1
mcmctree mcmctree.ctl
# usedata = 2    * 0: no data; 1: seq like; 2: normal approximation; 3: out.BV (in.BV)
