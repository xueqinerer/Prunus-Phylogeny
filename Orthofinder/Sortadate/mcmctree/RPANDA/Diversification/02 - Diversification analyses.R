#Required R packages and source the R codes and functions
library(picante)
library(pspline)
library(qpcR)
library(TreePar)

source("fit_bd.R"); source("fit_env_bd.R")
source("likelihood_bd.R")
source("Phi.R"); source("Psi.R")
source("integrate.R")

Aromobatidae <-read.tree("Aromobatidae.tre")
Centrolenidae <-read.tree("Centrolenidae.tre")
Dendrobatidae <-read.tree("Dendrobatidae.tre")
Hemiphractidae <-read.tree("Hemiphractidae.tre")
Leptodactylidae <-read.tree("Leptodactylidae.tre")
Liolaemidae <-read.tree("Liolaemidae.tre")
Liolaemus <-read.tree("Liolaemus.tre")

phylo <- list(Aromobatidae, Centrolenidae, Dendrobatidae, Hemiphractidae, Leptodactylidae, Liolaemidae, Liolaemus)
names <- list("Aromobatidae", "Centrolenidae", "Dendrobatidae", "Hemiphractidae", "Leptodactylidae", "Liolaemidae", "Liolaemus")
missing.lineages = c(9,32,61,32,37,63,74)

final<-list(); all_res<-list()

##
## Fit the models 
##

for (i in 1:length(phylo)){

	tree <- phylo[[i]]
	branch.times <- getx(tree)
	total.time <- max(node.age(tree)$ages)
	total.diversity <- Ntip(tree)+missing.lineages[i]
	sampling.fraction <- Ntip(tree)/(Ntip(tree)+missing.lineages[i])

cat(i, "tree out of ",length(phylo),": This is the clade",names[[i]],"that originated", round(total.time,2), "Myrs ago", "and includes", total.diversity, "extant species", "\n")

##
## Null models (for model testing purpose)
##
cat("Constant-rate diversification models", "\n")

# BCST (pure birth)
cat("BCST: pure birth model", "\n")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){0}
lamb_par<-c(0.1)
mu_par<-c()
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

	model_BCST<-fit_bd(tree,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BCST)

# BCST DCST (constant birth-death)
cat("BCST DCST: constant birth-death", "\n")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]}
lamb_par<-c(model_BCST$lamb_par[1])
mu_par<-c(0.01)
cst.lamb=T; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

	model_BCSTDCST<-fit_bd(tree,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BCSTDCST)


##
## Time Dependence
##
cat("Time-dependent diversification models", "\n")

# BTimeVar EXPO
cat("BTimeVar EXPO: speciation varies exponentially through time without extinction", "\n")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){0}
lamb_par<-c(model_BCST$lamb_par[1], 0.01)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=T

	model_BTimeVar_EXPO<-fit_bd(tree,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BTimeVar_EXPO)

# BTimeVar DCST EXPO
cat("BTimeVar DCST EXPO: speciation varies exponentially through time with constant extinction", "\n")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]}
lamb_par<-c(model_BTimeVar_EXPO$lamb_par[1],model_BTimeVar_EXPO$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=T; expo.mu=F; fix.mu=F

	model_BTimeVarDCST_EXPO<-fit_bd(tree,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BTimeVarDCST_EXPO)

# BCST DTimeVar EXPO
cat("BCST DTimeVar EXPO: constant speciation and extinction varies exponentially through time", "\n")
f.lamb<-function(x,y){y[1]}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(model_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0.001)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=T; fix.mu=F

	model_BCSTDTimeVar_EXPO<-fit_bd(tree,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BCSTDTimeVar_EXPO)

# BTimeVar DTimeVar EXPO
cat("BTimeVar DTimeVar EXPO: both speciation and extinction vary exponentially through time", "\n")
f.lamb<-function(x,y){y[1]*exp(y[2]*x)}
f.mu<-function(x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(model_BTimeVarDCST_EXPO$lamb_par[1],0.01)
mu_par<-c(0.01,0.001)
cst.lamb=F; cst.mu=F; expo.lamb=T; expo.mu=T; fix.mu=F

	model_BTimeVarDTimeVar_EXPO<-fit_bd(tree,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BTimeVarDTimeVar_EXPO)


##
## Temperature Dependence 
##
cat("Temperature-dependent diversification models", "\n")
env_data <- read.table("../Environmental data/merged_veizer_westerhold_Ts.txt", h=T)
df <- smooth.spline(x= env_data[,1], env_data[,2])$df; print(df)

# BTempVar EXPO
cat("BTempVar EXPO: speciation varies exponentially with temperature without extinction", "\n")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(model_BTimeVar_EXPO$lamb_par[1], 0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

	model_BTempVar_EXPO<-fit_env_bd(tree,env_data, df=80,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BTempVar_EXPO)

# BTempVar DCST EXPO
cat("BTempVar DCST EXPO: speciation varies exponentially with temperature and extinction is constant", "\n")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(model_BTempVar_EXPO$lamb_par[1]),model_BTempVar_EXPO$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

	model_BTempVarDCST_EXPO<-fit_env_bd(tree,env_data, df=80,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BTempVarDCST_EXPO)

# BCST DTempVar EXPO
cat("BCST DTempVar EXPO: constant speciation and extinction varies exponentially with temperature", "\n")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(model_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	model_BCSTDTempVar_EXPO<-fit_env_bd(tree,env_data, df=80,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BCSTDTempVar_EXPO)

# BTempVar DTempVar EXPO
cat("BTempVar DTempVar EXPO: both speciation and extinction vary exponentially with temperature", "\n")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(abs(model_BTempVarDCST_EXPO$lamb_par[1]),model_BTempVarDCST_EXPO$lamb_par[2])
mu_par<-c(0.01,0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	model_BTempVarDTempVar_EXPO<-fit_env_bd(tree,env_data, df=80,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BTempVarDTempVar_EXPO)


##
## Andean Dependence
##
cat("Andean-uplift-dependent diversification models", "\n")

if (i==1){env_data <- read.table("../Environmental data/Andes_mean_elevations_no_basins_NORTHCENTRAL.txt", h=T)}
if (i==2){env_data <- read.table("../Environmental data/Andes_mean_elevations_no_basins_NORTHCENTRAL.txt", h=T)}
if (i==3){env_data <- read.table("../Environmental data/Andes_mean_elevations_no_basins_NORTHCENTRAL.txt", h=T)}
if (i==4){env_data <- read.table("../Environmental data/Andes_mean_elevations_no_basins_NORTHCENTRAL.txt", h=T)}
if (i==5){env_data <- read.table("../Environmental data/Andes_mean_elevations_no_basins_ALL.txt", h=T)}
if (i==6){env_data <- read.table("../Environmental data/Andes_mean_elevations_no_basins_CENTRALSOUTH.txt", h=T)}
if (i==7){env_data <- read.table("../Environmental data/Andes_mean_elevations_no_basins_CENTRALSOUTH.txt", h=T)}

df <- smooth.spline(x= env_data[,1], env_data[,2])$df; print(df)

# BAndesVar EXPO
cat("BAndesVar EXPO: speciation varies exponentially with Andes without extinction", "\n")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){0}
lamb_par<-c(abs(model_BTimeVar_EXPO$lamb_par[1]),0)
mu_par<-c()
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=T

	model_BAndesVar_EXPO<-fit_env_bd(tree,env_data, df=80,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BAndesVar_EXPO)

# BAndesVar DCST EXPO
cat("BAndesVar DCST EXPO: speciation varies exponentially with Andes and extinction is constant", "\n")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]}
lamb_par<-c(abs(model_BAndesVar_EXPO$lamb_par[1]),model_BAndesVar_EXPO$lamb_par[2])
mu_par<-c(0.01)
cst.lamb=F; cst.mu=T; expo.lamb=F; expo.mu=F; fix.mu=F

	model_BAndesVarDCST_EXPO<-fit_env_bd(tree,env_data, df=80,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BAndesVarDCST_EXPO)

# BCST DAndesVar EXPO
cat("BCST DAndesVar EXPO: constant speciation and extinction varies exponentially with Andes", "\n")
f.lamb<-function(t,x,y){y[1]}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(model_BCSTDCST$lamb_par[1])
mu_par<-c(0.01,0)
cst.lamb=T; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	model_BCSTDAndesVar_EXPO<-fit_env_bd(tree,env_data, df=80,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BCSTDAndesVar_EXPO)

# BAndesVar DAndesVar EXPO
cat("BAndesVar DAndesVar EXPO: both speciation and extinction vary exponentially with Andes", "\n")
f.lamb<-function(t,x,y){y[1]*exp(y[2]*x)}
f.mu<-function(t,x,y){y[1]*exp(y[2]*x)}
lamb_par<-c(abs(model_BAndesVarDCST_EXPO$lamb_par[1]),model_BAndesVarDCST_EXPO$lamb_par[2])
mu_par<-c(0.01,0)
cst.lamb=F; cst.mu=F; expo.lamb=F; expo.mu=F; fix.mu=F

	model_BAndesVarDAndesVar_EXPO<-fit_env_bd(tree,env_data, df=80,total.time,f.lamb,f.mu,lamb_par,mu_par,f=sampling.fraction,cst.lamb=cst.lamb,cst.mu=cst.mu,expo.lamb=expo.lamb,expo.mu=expo.mu,fix.mu=fix.mu,cond="crown")
	print(model_BAndesVarDAndesVar_EXPO)


##
## RESULTS
##

	results<-matrix(NA,14,10)
	colnames(results)<-c("Models","Rate_variation","NP","logL","AICc","Akaike_w","Lambda","Alpha","Mu","Beta")

#Models
	results[,1]<-c("BCST","BCSTDCST",
	"BTimeVar","BTimeVarDCST","BCSTDTimeVar","BTimeVarDTimeVar",
	"BTempVar","BTempVarDCST","BCSTDTempVar","BTempVarDTempVar",
	"BAndesVar","BAndesVarDCST","BCSTDAndesVar","BAndesVarDAndesVar")

#Rate variation
	results[,2]<-c("constant","constant",
	"exponential","exponential","exponential","exponential",
	"exponential","exponential","exponential","exponential",
	"exponential","exponential","exponential","exponential")

#Number of parameters
	results[1,3]<-1
	results[2,3]<-2
	results[3,3]<-2
	results[4,3]<-3
	results[5,3]<-3
	results[6,3]<-4
	results[7,3]<-2
	results[8,3]<-3
	results[9,3]<-3
	results[10,3]<-4
	results[11,3]<-2
	results[12,3]<-3
	results[13,3]<-3
	results[14,3]<-4

#log-likelihood
	results[1,4]<-round(model_BCST$LH,3)
	results[2,4]<-round(model_BCSTDCST$LH,3)
	results[3,4]<-round(model_BTimeVar_EXPO$LH,3)
	results[4,4]<-round(model_BTimeVarDCST_EXPO$LH,3)
	results[5,4]<-round(model_BCSTDTimeVar_EXPO$LH,3)
	results[6,4]<-round(model_BTimeVarDTimeVar_EXPO$LH,3)
	results[7,4]<-round(model_BTempVar_EXPO$LH,3)
	results[8,4]<-round(model_BTempVarDCST_EXPO$LH,3)
	results[9,4]<-round(model_BCSTDTempVar_EXPO$LH,3)
	results[10,4]<-round(model_BTempVarDTempVar_EXPO$LH,3)
	results[11,4]<-round(model_BAndesVar_EXPO$LH,3)
	results[12,4]<-round(model_BAndesVarDCST_EXPO$LH,3)
	results[13,4]<-round(model_BCSTDAndesVar_EXPO$LH,3)
	results[14,4]<-round(model_BAndesVarDAndesVar_EXPO$LH,3)

#AICc
	results[1,5]<-round(model_BCST$aicc,3)
	results[2,5]<-round(model_BCSTDCST$aicc,3)
	results[3,5]<-round(model_BTimeVar_EXPO$aicc,3)
	results[4,5]<-round(model_BTimeVarDCST_EXPO$aicc,3)
	results[5,5]<-round(model_BCSTDTimeVar_EXPO$aicc,3)
	results[6,5]<-round(model_BTimeVarDTimeVar_EXPO$aicc,3)
	results[7,5]<-round(model_BTempVar_EXPO$aicc,3)
	results[8,5]<-round(model_BTempVarDCST_EXPO$aicc,3)
	results[9,5]<-round(model_BCSTDTempVar_EXPO$aicc,3)
	results[10,5]<-round(model_BTempVarDTempVar_EXPO$aicc,3)
	results[11,5]<-round(model_BAndesVar_EXPO$aicc,3)
	results[12,5]<-round(model_BAndesVarDCST_EXPO$aicc,3)
	results[13,5]<-round(model_BCSTDAndesVar_EXPO$aicc,3)
	results[14,5]<-round(model_BAndesVarDAndesVar_EXPO$aicc,3)

#Akaike weights
	all_AICc<-c(results[,5])
	results[1,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[1],3)
	results[2,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[2],3)
	results[3,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[3],3)
	results[4,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[4],3)
	results[5,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[5],3)
	results[6,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[6],3)
	results[7,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[7],3)
	results[8,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[8],3)
	results[9,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[9],3)
	results[10,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[10],3)
	results[11,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[11],3)
	results[12,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[12],3)
	results[13,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[13],3)
	results[14,6]<-round(akaike.weights(as.numeric(all_AICc))$weights[14],3)
	
#Speciation (lambda) parameter
	results[1,7]<-round(abs(model_BCST$lamb_par[1]),4)
	results[2,7]<-round(abs(model_BCSTDCST$lamb_par[1]),4)
	results[3,7]<-round(abs(model_BTimeVar_EXPO$lamb_par[1]),4)
	results[4,7]<-round(abs(model_BTimeVarDCST_EXPO$lamb_par[1]),4)
	results[5,7]<-round(abs(model_BCSTDTimeVar_EXPO$lamb_par[1]),4)
	results[6,7]<-round(abs(model_BTimeVarDTimeVar_EXPO$lamb_par[1]),4)
	results[7,7]<-round(abs(model_BTempVar_EXPO$lamb_par[1]),4)
	results[8,7]<-round(abs(model_BTempVarDCST_EXPO$lamb_par[1]),4)
	results[9,7]<-round(abs(model_BCSTDTempVar_EXPO$lamb_par[1]),4)
	results[10,7]<-round(abs(model_BTempVarDTempVar_EXPO$lamb_par[1]),4)
	results[11,7]<-round(abs(model_BAndesVar_EXPO$lamb_par[1]),4)
	results[12,7]<-round(abs(model_BAndesVarDCST_EXPO$lamb_par[1]),4)
	results[13,7]<-round(abs(model_BCSTDAndesVar_EXPO$lamb_par[1]),4)
	results[14,7]<-round(abs(model_BAndesVarDAndesVar_EXPO$lamb_par[1]),4)
	
#Speciation variation (alpha) parameter
	results[3,8]<-round(model_BTimeVar_EXPO$lamb_par[2],4)
	results[4,8]<-round(model_BTimeVarDCST_EXPO$lamb_par[2],4)
	results[6,8]<-round(model_BTimeVarDTimeVar_EXPO$lamb_par[2],4)
	results[7,8]<-round(model_BTempVar_EXPO$lamb_par[2],4)
	results[8,8]<-round(model_BTempVarDCST_EXPO$lamb_par[2],4)
	results[10,8]<-round(model_BTempVarDTempVar_EXPO$lamb_par[2],4)
	results[11,8]<-round(model_BAndesVar_EXPO$lamb_par[2],4)
	results[12,8]<-round(model_BAndesVarDCST_EXPO$lamb_par[2],4)
	results[14,8]<-round(model_BAndesVarDAndesVar_EXPO$lamb_par[2],4)
	
#Extinction (mu) parameter
	results[2,9]<-round(abs(model_BCSTDCST$mu_par[1]),4)
	results[4,9]<-round(abs(model_BTimeVarDCST_EXPO$mu_par[1]),4)
	results[5,9]<-round(abs(model_BCSTDTimeVar_EXPO$mu_par[1]),4)
	results[6,9]<-round(abs(model_BTimeVarDTimeVar_EXPO$mu_par[1]),4)
	results[8,9]<-round(abs(model_BTempVarDCST_EXPO$mu_par[1]),4)
	results[9,9]<-round(abs(model_BCSTDTempVar_EXPO$mu_par[1]),4)
	results[10,9]<-round(abs(model_BTempVarDTempVar_EXPO$mu_par[1]),4)
	results[12,9]<-round(abs(model_BAndesVarDCST_EXPO$mu_par[1]),4)
	results[13,9]<-round(abs(model_BCSTDAndesVar_EXPO$mu_par[1]),4)
	results[14,9]<-round(abs(model_BAndesVarDAndesVar_EXPO$mu_par[1]),4)
	
#Extinction variation (beta) parameter
	results[5,10]<-round(model_BCSTDTimeVar_EXPO$mu_par[2],4)
	results[6,10]<-round(model_BTimeVarDTimeVar_EXPO$mu_par[2],4)
	results[9,10]<-round(model_BCSTDTempVar_EXPO$mu_par[2],4)
	results[10,10]<-round(model_BTempVarDTempVar_EXPO$mu_par[2],4)
	results[13,10]<-round(model_BCSTDAndesVar_EXPO$mu_par[2],4)
	results[14,10]<-round(model_BAndesVarDAndesVar_EXPO$mu_par[2],4)
	
	print(results)

	all_res[[i]]<-list("Clade_name"=names[[i]],"Clade_age"=total.time,"Taxon_sampling"=Ntip(tree),"Sampling_fraction"=sampling.fraction,
	"BCST"=model_BCST,"BCSTDCST"=model_BCSTDCST,
	"BTimeVar_EXPO"=model_BTimeVar_EXPO,"BTimeVarDCST_EXPO"=model_BTimeVarDCST_EXPO,"BCSTDTimeVar_EXPO"=model_BCSTDTimeVar_EXPO,"BTimeVarDTimeVar_EXPO"=model_BTimeVarDTimeVar_EXPO,
	"BTempVar_EXPO"=model_BTempVar_EXPO,"BTempVarDCST_EXPO"=model_BTempVarDCST_EXPO,"BCSTDTempVar_EXPO"=model_BCSTDTempVar_EXPO,"BTempVarDTempVar_EXPO"=model_BTempVarDTempVar_EXPO,
	"BAndesVar_EXPO"=model_BAndesVar_EXPO,"BAndesVarDCST_EXPO"=model_BAndesVarDCST_EXPO,"BCSTDAndesVar_EXPO"=model_BCSTDAndesVar_EXPO,"BAndesVarDAndesVar_EXPO"=model_BAndesVarDAndesVar_EXPO)
	
	final[[i]]<-results
}

	write.table(final,file="Results_Andean_clades_PDDM.txt",quote=FALSE,sep="\t",row.names=FALSE)
	save(all_res,file="Results_Andean_clades_PDDM.Rdata")
