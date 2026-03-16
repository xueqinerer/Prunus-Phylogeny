# Step 1: Open FigTree.tre in FigTree software and export as NEXUS format
# Save the exported file as FigTree_NEXUS.tre

# Note: remove the extra variable prefix (single quotes around variable names)
# Remove outgroup taxon
pxrmt -t Prunus.nw.ultra.tre -n Lyonothamnus_floribundus > Prunus.mcmctree.dated_no_outgroup.tre
sed -i 's/Prunus_arborea_var/Prunus_arborea_var._montana/g;s/Prunus_salicina_var/Prunus_salicina_var._mandshurica/g;s/Prunus_davidiana_var/Prunus_davidiana_var._potaninii/g' Prunus.mcmctree.dated_no_outgroup.tre


library("BAMMtools")
Prunus <- read.tree("./Prunus.mcmctree.dated_no_outgroup.tre")
#setBAMMpriors(phy = Prunus, outfile = NULL)
setBAMMpriors(phy = Prunus, total.taxa = 352, outfile = NULL)
#expectedNumberOfShifts        #lambdaInitPrior       #lambdaShiftPrior
            #1.00000000             #2.43831651             #0.01826394
          # muInitPrior
           # 2.43831651
