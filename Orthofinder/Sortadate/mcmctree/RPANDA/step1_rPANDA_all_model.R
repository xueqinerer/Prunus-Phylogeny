# Required R packages
library(RPANDA)
library(picante)
library(pspline)
options(rgl.useNULL = TRUE)
library(qpcR)

# Load custom RPANDA functions
source("./Diversification/fit_bd.R")
source("./Diversification/fit_env_bd.R")
source("./Diversification/likelihood_bd.R")
source("./Diversification/Phi.R")
source("./Diversification/Psi.R")
source("./Diversification/integrate.R")

# Command-line argument: tree index
args <- commandArgs(trailingOnly = TRUE)
i <- as.numeric(args[1])

# Variables
#missing.lineages <- 0
names <- c("Clematis_tree")

# Tree path
tree_path <- "/path/to/BAMM/Prunus.mcmctree.dated_no_outgroup.tre"
Clematis_tree <- read.tree(tree_path)

# Confirmation
cat("Successfully read tree from:\n", tree_path, "\n")

# Prepare for model
phyloi <- Clematis_tree
tot_time <- max(node.age(phyloi)$ages)
f <- 83 / 343
cond <- "crown"
posteriors <- 1

# Environmental data
CO2_data0 <- read.table("../data/500kyrCO2_combined_62.75Ma.txt", header = TRUE)
CO2_data1 <- as.vector(scale(CO2_data0[,2]))
CO2_data <- cbind(CO2_data0[,1], CO2_data1)

Temp_data0 <- read.csv("../data/temperature-0.5ma-FABIEN-approx-63.csv", header = TRUE)
Temp_data1 <- as.vector(scale(Temp_data0[,2]))
Temp_data <- cbind(Temp_data0[,1], Temp_data1)

cat("Running model for tree", i, "\n")
#####################################################
###### Pure birth / constant Birth-death         ####
#####################################################
# BCST (Pure birth)
print("BCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){0}
lamb_par<-c(0.1)
mu_par<-c()
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
treei_BCST<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCST)

# BCST DCST (constant Birth-death)
print("BCST DCST")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(treei_BCST$lamb_par[1])
mu_par<-c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCSTDCST<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCSTDCST)

#####################################################
###### CO2Dependence (exponential variation) ######
#####################################################
# BCO2Var EXPO
print("BCO2Var EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}
#lamb_par<-c(0.1,0.01)
lamb_par<-c(0.03317581,0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
treei_BCO2Var_EXPO<-fit_env_bd(phyloi,CO2_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCO2Var_EXPO)

# BCO2Var DCST EXPO
print("BCO2Var DCST EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(treei_BCO2Var_EXPO$lamb_par[1]),treei_BCO2Var_EXPO$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCO2VarDCST_EXPO<-fit_env_bd(phyloi,CO2_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCO2VarDCST_EXPO)
  
  
# BCST DCO2Var EXPO
print("BCST DCO2Var EXPO")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.001)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCSTDCO2Var_EXPO<-fit_env_bd(phyloi,CO2_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCSTDCO2Var_EXPO)
  
  
# BCO2Var DCO2Var EXPO
print("BCO2Var DCO2Var EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(abs(treei_BCO2VarDCST_EXPO$lamb_par[1]),treei_BCO2VarDCST_EXPO$lamb_par[2])
mu_par<-c(treei_BCSTDCO2Var_EXPO$mu_par[1], treei_BCSTDCO2Var_EXPO$mu_par[2])
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCO2VarDCO2Var_EXPO<-fit_env_bd(phyloi,CO2_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCO2VarDCO2Var_EXPO)
  
################################################
###### CO2Dependence (linear variation) ######
################################################
# print(i)
# BCO2Var LIN
print("BCO2Var LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.01,0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
treei_BCO2Var_LIN<-fit_env_bd(phyloi,CO2_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCO2Var_LIN)
 
# BCO2Var DCST LIN
print("BCO2Var DCST LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(treei_BCO2Var_LIN$lamb_par[1]),treei_BCO2Var_LIN$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCO2VarDCST_LIN<-fit_env_bd(phyloi,CO2_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCO2VarDCST_LIN)
  
  
# BCST DCO2Var LIN
print("BCST DCO2Var LIN")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]+y[2]*x}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
#mu_par<-c(0.01,0.001)
mu_par<-c(0.01,0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCSTDCO2Var_LIN<-fit_env_bd(phyloi,CO2_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCSTDCO2Var_LIN)
  
  
# BCO2Var DCO2Var LIN
print("BCO2Var DCO2Var LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]+y[2]*x}
lamb_par<-c(abs(treei_BCO2VarDCST_LIN$lamb_par[1]),0)
mu_par<-c(0.05,0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCO2VarDCO2Var_LIN<-fit_env_bd(phyloi,CO2_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCO2VarDCO2Var_LIN)
  
############# CO2_results ###########################################
results<-matrix(NA,19,9)
colnames(results)<-c("Models","Parameters","logL","AICc","Lambda","Alpha","Mu","Beta","Akaike_w")
#Models
results[,1]<-c("BCSTDCST",
  "BCO2VarDCST_EXPO",
  "BCSTDCO2Var_EXPO",
  "BCO2VarDCO2Var_EXPO",
  "BCO2VarDCST_LIN",
  "BCSTDCO2Var_LIN",
  "BCO2VarDCO2Var_LIN",
  "BTempVarDCST_EXPO",
  "BCSTDTempVar_EXPO",
  "BTempVarDTempVar_EXPO",
  "BTempVarDCST_LIN",
  "BCSTDTempVar_LIN",
  "BTempVarDTempVar_LIN",
  "BTimeVarDCST_EXPO",
  "BCSTDTimeVar_EXPO",
  "BTimeVarDTimeVar_EXPO",
  "BTimeVarDCST_LIN",
  "BCSTDTimeVar_LIN",
  "BTimeVarDTimeVar_LIN" )
#Parameters
  results[1,2]<-2
  results[2,2]<-3
  results[3,2]<-3
  results[4,2]<-4
  results[5,2]<-3
  results[6,2]<-3
  results[7,2]<-4

#logL
  results[1,3]<-round(treei_BCSTDCST$LH,3)
  results[2,3]<-round(treei_BCO2VarDCST_EXPO$LH,3)
  results[3,3]<-round(treei_BCSTDCO2Var_EXPO$LH,3)
  results[4,3]<-round(treei_BCO2VarDCO2Var_EXPO$LH,3)
  results[5,3]<-round(treei_BCO2VarDCST_LIN$LH,3)
  results[6,3]<-round(treei_BCSTDCO2Var_LIN$LH,3)
  results[7,3]<-round(treei_BCO2VarDCO2Var_LIN$LH,3)
  
  
#AICc             
  results[1,4]<-round(treei_BCSTDCST$aicc,3)
  results[2,4]<-round(treei_BCO2VarDCST_EXPO$aicc,3)
  results[3,4]<-round(treei_BCSTDCO2Var_EXPO$aicc,3)
  results[4,4]<-round(treei_BCO2VarDCO2Var_EXPO$aicc,3)
  results[5,4]<-round(treei_BCO2VarDCST_LIN$aicc,3)
  results[6,4]<-round(treei_BCSTDCO2Var_LIN$aicc,3)
  results[7,4]<-round(treei_BCO2VarDCO2Var_LIN$aicc,3)
  
  
#Lambda0          
  results[1,5]<-round((treei_BCSTDCST$lamb_par[1]),3)
  results[2,5]<-round((treei_BCO2VarDCST_EXPO$lamb_par[1]),3)
  results[3,5]<-round((treei_BCSTDCO2Var_EXPO$lamb_par[1]),3)
  results[4,5]<-round((treei_BCO2VarDCO2Var_EXPO$lamb_par[1]),3)
  results[5,5]<-round((treei_BCO2VarDCST_LIN$lamb_par[1]),3)
  results[6,5]<-round((treei_BCSTDCO2Var_LIN$lamb_par[1]),3)
  results[7,5]<-round((treei_BCO2VarDCO2Var_LIN$lamb_par[1]),3)
  
  
#Alpha CO2     
  
  results[2,6]<-round(treei_BCO2VarDCST_EXPO$lamb_par[2],4)
  results[4,6]<-round(treei_BCO2VarDCO2Var_EXPO$lamb_par[2],4)
  results[5,6]<-round(treei_BCO2VarDCST_LIN$lamb_par[2],4)
  results[7,6]<-round(treei_BCO2VarDCO2Var_LIN$lamb_par[2],4)
  
  
#Mu0              
  results[1,7]<-round((treei_BCSTDCST$mu_par[1]),3)
  results[2,7]<-round((treei_BCO2VarDCST_EXPO$mu_par[1]),3)
  results[3,7]<-round((treei_BCSTDCO2Var_EXPO$mu_par[1]),3)
  results[4,7]<-round((treei_BCO2VarDCO2Var_EXPO$mu_par[1]),3)
  results[5,7]<-round((treei_BCO2VarDCST_LIN$mu_par[1]),3)
  results[6,7]<-round((treei_BCSTDCO2Var_LIN$mu_par[1]),3)
  results[7,7]<-round((treei_BCO2VarDCO2Var_LIN$mu_par[1]),3)
  
  
#Beta CO2       
  results[3,8]<-round(treei_BCSTDCO2Var_EXPO$mu_par[2],4)
  results[4,8]<-round(treei_BCO2VarDCO2Var_EXPO$mu_par[2],4)
  results[6,8]<-round(treei_BCSTDCO2Var_LIN$mu_par[2],4)
  results[7,8]<-round(treei_BCO2VarDCO2Var_LIN$mu_par[2],4)
  



#####################################################
###### TempDependence (exponential variation) ######
#####################################################
# BTempVar EXPO
print("BTempVar EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}
#lamb_par<-c(0.1,0.01)
lamb_par<-c(0.03317581,0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
treei_BTempVar_EXPO<-fit_env_bd(phyloi,Temp_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTempVar_EXPO)

# BTempVar DCST EXPO
print("BTempVar DCST EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(treei_BTempVar_EXPO$lamb_par[1]),treei_BTempVar_EXPO$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BTempVarDCST_EXPO<-fit_env_bd(phyloi,Temp_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTempVarDCST_EXPO)
  
  
# BCST DTempVar EXPO
print("BCST DTempVar EXPO")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.001)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCSTDTempVar_EXPO<-fit_env_bd(phyloi,Temp_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCSTDTempVar_EXPO)
  
  
# BTempVar DTempVar EXPO
print("BTempVar DTempVar EXPO")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(abs(treei_BTempVarDCST_EXPO$lamb_par[1]),treei_BTempVarDCST_EXPO$lamb_par[2])
mu_par<-c(treei_BCSTDTempVar_EXPO$mu_par[1], treei_BCSTDTempVar_EXPO$mu_par[2])
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BTempVarDTempVar_EXPO<-fit_env_bd(phyloi,Temp_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTempVarDTempVar_EXPO)
  
################################################
###### TempDependence (linear variation) ######
################################################
# print(i)
# BTempVar LIN
print("BTempVar LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){0}
lamb_par<-c(0.01,0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
treei_BTempVar_LIN<-fit_env_bd(phyloi,Temp_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTempVar_LIN)
 
# BTempVar DCST LIN
print("BTempVar DCST LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(treei_BTempVar_LIN$lamb_par[1]),treei_BTempVar_LIN$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BTempVarDCST_LIN<-fit_env_bd(phyloi,Temp_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTempVarDCST_LIN)
  
  
# BCST DTempVar LIN
print("BCST DTempVar LIN")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]+y[2]*x}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
#mu_par<-c(0.01,0.001)
mu_par<-c(0.01,0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCSTDTempVar_LIN<-fit_env_bd(phyloi,Temp_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCSTDTempVar_LIN)
  
  
# BTempVar DTempVar LIN
print("BTempVar DTempVar LIN")
f.lamb<-function(t,x,y){y[1]+y[2]*x}
f.mu<-function(t,x,y){y[1]+y[2]*x}
lamb_par<-c(abs(treei_BTempVarDCST_LIN$lamb_par[1]),0)
mu_par<-c(0.05,0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BTempVarDTempVar_LIN<-fit_env_bd(phyloi,Temp_data,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTempVarDTempVar_LIN)
  
#Parameters
  results[8,2]<-3
  results[9,2]<-3
  results[10,2]<-4
  results[11,2]<-3
  results[12,2]<-3
  results[13,2]<-4
  
#logL
  results[8,3]<-round(treei_BTempVarDCST_EXPO$LH,3)
  results[9,3]<-round(treei_BCSTDTempVar_EXPO$LH,3)
  results[10,3]<-round(treei_BTempVarDTempVar_EXPO$LH,3)
  results[11,3]<-round(treei_BTempVarDCST_LIN$LH,3)
  results[12,3]<-round(treei_BCSTDTempVar_LIN$LH,3)
  results[13,3]<-round(treei_BTempVarDTempVar_LIN$LH,3)
  
#AICc             
  results[8,4]<-round(treei_BTempVarDCST_EXPO$aicc,3)
  results[9,4]<-round(treei_BCSTDTempVar_EXPO$aicc,3)
  results[10,4]<-round(treei_BTempVarDTempVar_EXPO$aicc,3)
  results[11,4]<-round(treei_BTempVarDCST_LIN$aicc,3)
  results[12,4]<-round(treei_BCSTDTempVar_LIN$aicc,3)
  results[13,4]<-round(treei_BTempVarDTempVar_LIN$aicc,3)
  
#Lambda0          
  results[8,5]<-round((treei_BTempVarDCST_EXPO$lamb_par[1]),3)
  results[9,5]<-round((treei_BCSTDTempVar_EXPO$lamb_par[1]),3)
  results[10,5]<-round((treei_BTempVarDTempVar_EXPO$lamb_par[1]),3)
  results[11,5]<-round((treei_BTempVarDCST_LIN$lamb_par[1]),3)
  results[12,5]<-round((treei_BCSTDTempVar_LIN$lamb_par[1]),3)
  results[13,5]<-round((treei_BTempVarDTempVar_LIN$lamb_par[1]),3)
  
#Alpha Temp     
  
  results[8,6]<-round(treei_BTempVarDCST_EXPO$lamb_par[2],4)
  results[10,6]<-round(treei_BTempVarDTempVar_EXPO$lamb_par[2],4)
  results[11,6]<-round(treei_BTempVarDCST_LIN$lamb_par[2],4)
  results[13,6]<-round(treei_BTempVarDTempVar_LIN$lamb_par[2],4)
  
  
#Mu0              
  results[8,7]<-round((treei_BTempVarDCST_EXPO$mu_par[1]),3)
  results[9,7]<-round((treei_BCSTDTempVar_EXPO$mu_par[1]),3)
  results[10,7]<-round((treei_BTempVarDTempVar_EXPO$mu_par[1]),3)
  results[11,7]<-round((treei_BTempVarDCST_LIN$mu_par[1]),3)
  results[12,7]<-round((treei_BCSTDTempVar_LIN$mu_par[1]),3)
  results[13,7]<-round((treei_BTempVarDTempVar_LIN$mu_par[1]),3)
  
  
#Beta Temp       
  results[9,8]<-round(treei_BCSTDTempVar_EXPO$mu_par[2],4)
  results[10,8]<-round(treei_BTempVarDTempVar_EXPO$mu_par[2],4)
  results[12,8]<-round(treei_BCSTDTempVar_LIN$mu_par[2],4)
  results[13,8]<-round(treei_BTempVarDTempVar_LIN$mu_par[2],4)
  
#####################################################
###### Time Dependence (exponential variation) ######
#####################################################
# BTimeVar EXPO
print("BTimeVar EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){0}
lamb_par<-c(0.1,0.01)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=T

	treei_BTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVar_EXPO)
	

# BTimeVar DCST EXPO
print("BTimeVar DCST EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]}
lamb_par<-c(treei_BTimeVar_EXPO$lamb_par[1],treei_BTimeVar_EXPO$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=F

	treei_BTimeVarDCST_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDCST_EXPO)


# BCST DTimeVar EXPO
print("BCST DTimeVar EXPO")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.01)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=T; fix.mu=F

	treei_BCSTDTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BCSTDTimeVar_EXPO)


# BTimeVar DTimeVar EXPO
print("BTimeVar DTimeVar EXPO")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(treei_BTimeVarDCST_EXPO$lamb_par[1],0.01)
mu_par<-c(0.05,0.01)
cst.lamb=F; cst.mu=F; expo.lamb=T; expo.mu=T; fix.mu=F

	treei_BTimeVarDTimeVar_EXPO<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
	print(treei_BTimeVarDTimeVar_EXPO)


################################################
###### Time Dependence (linear variation) ######
################################################

# BTimeVar LIN
print("BTimeVar LIN")
f.lamb<-function(x,y){y[1]+(y[2]*x)}
f.mu<-function(x,y){0}
lamb_par<-c(0.01,0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T
treei_BTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTimeVar_LIN)
	
	
# BTimeVar DCST LIN
print("BTimeVar DCST LIN")
f.lamb<-function(x,y){y[1]+(y[2]*x)}
f.mu<-function(x,y){y[1]}
lamb_par<-c(abs(treei_BTimeVar_LIN$lamb_par[1]),treei_BTimeVar_LIN$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BTimeVarDCST_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTimeVarDCST_LIN)
	
	
# BCST DTimeVar LIN
print("BCST DTimeVar LIN")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]+(y[2]*x)}
lamb_par<-c(treei_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BCSTDTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BCSTDTimeVar_LIN)
	
	
# BTimeVar DTimeVar LIN
print("BTimeVar DTimeVar LIN")
f.lamb<-function(x,y){y[1]+(y[2]*x)}
f.mu<-function(x,y){y[1]+(y[2]*x)}
lamb_par<-c(abs(treei_BTimeVarDCST_LIN$lamb_par[1]),0)
mu_par<-c(0.05,0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F
treei_BTimeVarDTimeVar_LIN<-fit_bd(phyloi,tot_time,f.lamb,f.mu,lamb_par,mu_par,f=f,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond=cond)
print(treei_BTimeVarDTimeVar_LIN)
	
###################################
###########results############
###################################
#Parameters
results[14,2]<-3
results[15,2]<-3
results[16,2]<-4
results[17,2]<-3
results[18,2]<-3
results[19,2]<-4

	#logL

	results[14,3]<-round(treei_BTimeVarDCST_EXPO$LH,3)
	results[15,3]<-round(treei_BCSTDTimeVar_EXPO$LH,3)
	results[16,3]<-round(treei_BTimeVarDTimeVar_EXPO$LH,3)
	results[17,3]<-round(treei_BTimeVarDCST_LIN$LH,3)
	results[18,3]<-round(treei_BCSTDTimeVar_LIN$LH,3)
	results[19,3]<-round(treei_BTimeVarDTimeVar_LIN$LH,3)
	
	#AICc
	results[14,4]<-round(treei_BTimeVarDCST_EXPO$aicc,3)
	results[15,4]<-round(treei_BCSTDTimeVar_EXPO$aicc,3)
	results[16,4]<-round(treei_BTimeVarDTimeVar_EXPO$aicc,3)
	results[17,4]<-round(treei_BTimeVarDCST_LIN$aicc,3)
	results[18,4]<-round(treei_BCSTDTimeVar_LIN$aicc,3)
	results[19,4]<-round(treei_BTimeVarDTimeVar_LIN$aicc,3)
	
	
	#Lambda0
	results[14,5]<-round(abs(treei_BTimeVarDCST_EXPO$lamb_par[1]),3)
	results[15,5]<-round(abs(treei_BCSTDTimeVar_EXPO$lamb_par[1]),3)
	results[16,5]<-round(abs(treei_BTimeVarDTimeVar_EXPO$lamb_par[1]),3)
	results[17,5]<-round(abs(treei_BTimeVarDCST_LIN$lamb_par[1]),3)
	results[18,5]<-round(abs(treei_BCSTDTimeVar_LIN$lamb_par[1]),3)
	results[19,5]<-round(abs(treei_BTimeVarDTimeVar_LIN$lamb_par[1]),3)
	
	#Alpha Time

	results[14,6]<-round(treei_BTimeVarDCST_EXPO$lamb_par[2],4)
	results[16,6]<-round(treei_BTimeVarDTimeVar_EXPO$lamb_par[2],4)
	results[17,6]<-round(treei_BTimeVarDCST_LIN$lamb_par[2],4)
	results[19,6]<-round(treei_BTimeVarDTimeVar_LIN$lamb_par[2],4)
	
	#Mu0
	results[14,7]<-round(abs(treei_BTimeVarDCST_EXPO$mu_par[1]),3)
	results[15,7]<-round(abs(treei_BCSTDTimeVar_EXPO$mu_par[1]),3)
	results[16,7]<-round(abs(treei_BTimeVarDTimeVar_EXPO$mu_par[1]),3)
	results[17,7]<-round(abs(treei_BTimeVarDCST_LIN$mu_par[1]),3)
	results[18,7]<-round(abs(treei_BCSTDTimeVar_LIN$mu_par[1]),3)
	results[19,7]<-round(abs(treei_BTimeVarDTimeVar_LIN$mu_par[1]),3)
	
	#Beta Time
	results[15,8]<-round(treei_BCSTDTimeVar_EXPO$mu_par[2],4)
	results[16,8]<-round(treei_BTimeVarDTimeVar_EXPO$mu_par[2],4)
	results[18,8]<-round(treei_BCSTDTimeVar_LIN$mu_par[2],4)
	results[19,8]<-round(treei_BTimeVarDTimeVar_LIN$mu_par[2],4)
	


#Akaike weights
  all_AICc<-c(results[,4])
  results[1,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[1],3)
  results[2,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[2],3)
  results[3,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[3],3)
  results[4,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[4],3)
  results[5,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[5],3)
  results[6,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[6],3)
  results[7,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[7],3)
  results[8,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[8],3)
  results[9,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[9],3)
  results[10,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[10],3)
  results[11,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[11],3)
  results[12,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[12],3)
  results[13,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[13],3)
  results[14,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[14],3)
  results[15,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[15],3)
  results[16,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[16],3)
  results[17,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[17],3)
  results[18,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[18],3)
  results[19,9]<-round(akaike.weights(as.numeric(all_AICc))$weights[19],3)

best_index <- which.max(as.numeric(results[, "Akaike_w"]))
best_model <- results[best_index, "Models"]
best_weight <- results[best_index, "Akaike_w"]

cat("Best model based on Akaike weight:", best_model, "\n")
cat("Akaike weight of the best model:", best_weight, "\n")
resi_all<-list("Clade"=names,
                 "Clade_age"=tot_time,
                 "Taxon_sampling"=83,
                 "Clade_size" = 343, 
                 "Sampling_fraction"=f,
                 "best_model"=best_model,
                 "BCSTDCST"=treei_BCSTDCST, 
                 "BCO2VarDCST_EXPO"=treei_BCO2VarDCST_EXPO,
                 "BCSTDCO2Var_EXPO"=treei_BCSTDCO2Var_EXPO,
                 "BCO2VarDCO2Var_EXPO"=treei_BCO2VarDCO2Var_EXPO,
                 "BCO2VarDCST_LIN"=treei_BCO2VarDCST_LIN,
                 "BCSTDCO2Var_LIN"=treei_BCSTDCO2Var_LIN,
                 "BCO2VarDCO2Var_LIN"=treei_BCO2VarDCO2Var_LIN,
                 "BTempVarDCST_EXPO"=treei_BTempVarDCST_EXPO,
                 "BCSTDTempVar_EXPO"=treei_BCSTDTempVar_EXPO,
                 "BTempVarDTempVar_EXPO"=treei_BTempVarDTempVar_EXPO,
                 "BTempVarDCST_LIN"=treei_BTempVarDCST_LIN,
                 "BCSTDTempVar_LIN"=treei_BCSTDTempVar_LIN,
                 "BTempVarDTempVar_LIN"=treei_BTempVarDTempVar_LIN,
                 "BTimeVarDCST_EXPO"=treei_BTimeVarDCST_EXPO,
                 "BCSTDTimeVar_EXPO"=treei_BCSTDTimeVar_EXPO,
                 "BTimeVarDTimeVar_EXPO"=treei_BTimeVarDTimeVar_EXPO,
                 "BTimeVarDCST_LIN"=treei_BTimeVarDCST_LIN,
                 "BCSTDTimeVar_LIN"=treei_BCSTDTimeVar_LIN,
                 "BTimeVarDTimeVar_LIN"=treei_BTimeVarDTimeVar_LIN)
  
Clematis_res_all<-c(resi_all,list(resi_all))
finalClematis_all <- results
final_table_all <- as.data.frame(results)

output_dir <- file.path("results", paste0("Clematis_tact", i))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


write.csv(final_table_all, file = file.path(output_dir, paste0("results.tact", i, "_all_models.csv")))
save(final_table_all,     file = file.path(output_dir, paste0("final_table.tact", i, "_all_models.Rdata")))
save(finalClematis_all,   file = file.path(output_dir, paste0("final.tact", i, "_all_models.Rdata")))
save(Clematis_res_all,    file = file.path(output_dir, paste0("tact", i, "_all_res_list.Rdata")))
