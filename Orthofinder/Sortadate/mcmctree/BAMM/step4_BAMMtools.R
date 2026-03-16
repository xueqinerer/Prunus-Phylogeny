library("BAMMtools")
Prunus <- read.tree("./Prunus.mcmctree.dated_no_outgroup.tre")
#setBAMMpriors(phy = Prunus, outfile = NULL)
setBAMMpriors(phy = Prunus,total.taxa = 352, outfile = NULL)
#expectedNumberOfShifts        #lambdaInitPrior       #lambdaShiftPrior
            #1.00000000             #2.43831651             #0.01826394
          # muInitPrior
           # 2.43831651
