setwd("/path/to/RPANDA")

list.files()

library(pspline)
library(ggplot2)
library(picante)
library(scales)
library(RPANDA)


# Load RPANDA results from previous run
load("tactNA_all_res_list.Rdata")

###############################
## Plots for the best models ##
###############################

## Geological time scale background
add_geochrono <- function(Y1,Y2){
  polygon(-c(145,145,100.5,100.5),     c(Y1,Y2,Y2,Y1), col = "#A0C96D", lwd = 0.5) # Lower Cretaceous
  polygon(-c(100.5,100.5,66,66),       c(Y1,Y2,Y2,Y1), col = "#BAD25F", lwd = 0.5) # Upper Cretaceous
  polygon(-c(66,66,56,56),             c(Y1,Y2,Y2,Y1), col = "#F8B77D", lwd = 0.5) # Paleocene
  polygon(-c(56,56,33.9,33.9),         c(Y1,Y2,Y2,Y1), col = "#FAC18A", lwd = 0.5) # Eocene
  polygon(-c(33.9,33.9,23.03,23.03),   c(Y1,Y2,Y2,Y1), col = "#FBCC98", lwd = 0.5) # Oligocene
  polygon(-c(23.03,23.03,5.33,5.33),   c(Y1,Y2,Y2,Y1), col = "#FFED00", lwd = 0.5) # Miocene
  polygon(-c(5.33,5.33,2.58,2.58),     c(Y1,Y2,Y2,Y1), col = "#FFF7B2", lwd = 0.5) # Pliocene
  polygon(-c(2.58,2.58,0.0117,0.0117), c(Y1,Y2,Y2,Y1), col = "#FFF1C4", lwd = 0.5) # Pleistocene
  polygon(-c(0.0117,0.0117,0,0),       c(Y1,Y2,Y2,Y1), col = "#FEF6F2", lwd = 0.5) # Holocene
}

## Plot colors
colors <- c("chartreuse3","dodgerblue","orange","darkgray")

##################
## Prunus       ##
##################

pdf(file="Diversification dynamics of Prunus.pdf", height=12, width=8)

par(mfrow=c(3,1), mar=c(3,3,0.5,0.5))
resi_all <- Clematis_res_all
age <- resi_all$Clade_age

# --- Panel 1: Time-dependent diversification (BTimeVarDTimeVar_LIN) ---
speciation_par1 <- resi_all$BTimeVarDTimeVar_LIN$lamb_par[1]
speciation_par2 <- resi_all$BTimeVarDTimeVar_LIN$lamb_par[2]

# Define speciation function (linear time dependence)
speciation <- function(time) {
  speciation_par1 + speciation_par2 * time
}

# Define extinction function (linear time dependence)
extinction_par1 <- resi_all$BTimeVarDTimeVar_LIN$mu_par[1]
extinction_par2 <- resi_all$BTimeVarDTimeVar_LIN$mu_par[2]
extinction <- function(time) {
  extinction_par1 + extinction_par2 * time
}

plot(age, age, type='n',
     ylim=c(0, round(max(speciation(c(seq(age, 0, by=-1), 0)),
                         extinction(c(seq(age, 0, by=-1), 0)), na.rm=TRUE), 3)),
     xlim=c(-age, 0), axes=F, ylab='', xlab='',
     cex.axis=0.5, cex.lab=0.5)

rect(-66, 0, -age, 6.7, col="whitesmoke", border=NA)
rect(-56, 0, -66, 6.7, col="white", border=NA)
rect(-56, 0, -33.9, 6.7, col="whitesmoke", border=NA)
rect(-33.9, 0, -23.03, 6.7, col="white", border=NA)
rect(-23.03, 0, -5.33, 6.7, col="whitesmoke", border=NA)
rect(-5.33, 0, -2.58, 6.7, col="white", border=NA)
rect(-2.58, 0, 0, 6.7, col="whitesmoke", border=NA)

lines(c(seq(-age, 0, by=1), 0), speciation(c(seq(age, 0, by=-1), 0)),
      type="l", col=alpha(colors[3], 0.7), lwd=2)
lines(c(seq(-age, 0, by=1), 0), extinction(c(seq(age, 0, by=-1), 0)),
      type="l", col=alpha(colors[4], 0.7), lwd=2, lty="dashed")

abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -age), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

axis(side=1, at=seq(-round(age, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)

max_speciation <- round(max(speciation(c(seq(age, 0, by=-1), 0)), na.rm=TRUE), 3)
max_extinction <- round(max(extinction(c(seq(age, 0, by=-1), 0)), na.rm=TRUE), 3)
max_y <- max(max_speciation, max_extinction)

axis(side=2, at=pretty(c(0, max_y)), col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)
legend("topright", bty="n", c("Speciation rate", "Extinction rate"), lwd="2",
       col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7)), lty=c(1, 2), cex=0.8)

add_geochrono(0, -0.2)
mtext("BTimeVarDTimeVar_LIN  AICweight = 0.767", side=3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side=1, line=1.5, cex=0.8)

# --- Panel 2: Temperature-dependent diversification (BTempVarDTempVar_LIN) ---
Temperature <- read.csv("./data/temperature-0.5ma-FABIEN-approx-63.csv", header=T)
Temperature_data1 <- as.vector(scale(Temperature[,2]))
Temperature_data <- cbind(Temperature[,1], Temperature_data1)
colnames(Temperature_data) <- c("Age", "Temperature")
Temperature <- as.data.frame(Temperature_data)

Tm_spline <- sm.spline(Temperature[,1], Temperature[,2])
Tm_time <- function(x) { predict(Tm_spline, x) }

speciation_par1 <- resi_all$BTempVarDTempVar_LIN$lamb_par[1]
speciation_par2 <- resi_all$BTempVarDTempVar_LIN$lamb_par[2]
speciation <- function(x) { abs(speciation_par1 + speciation_par2 * Tm_time(x)) }
speciation.min <- round(min(speciation(Temperature[,"Age"])),3)
speciation.max <- round(max(speciation(Temperature[,"Age"])),3)

extinction_par1 <- resi_all$BTempVarDTempVar_LIN$mu_par[1]
extinction_par2 <- resi_all$BTempVarDTempVar_LIN$mu_par[2]
extinction <- function(x) { abs(extinction_par1 + extinction_par2 * Tm_time(x)) }
extinction.min <- round(min(extinction(Temperature[,"Age"])),3)
extinction.max <- round(max(extinction(Temperature[,"Age"])),3)

plot(age, age, type='n', ylim=c(0, max(speciation.max, extinction.max)),
     xlim=c(-age,0), axes=F, ylab='', xlab='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)

rect(-66,0,-max(Temperature$Age),1, col="whitesmoke", border=NA)
rect(-56,0,-66,1, col="white", border=NA)
rect(-56,0,-33.9,1, col="whitesmoke", border=NA)
rect(-33.9,0,-23.03,1, col="white", border=NA)
rect(-23.03,0,-5.33,1, col="whitesmoke", border=NA)
rect(-5.33,0,-2.58,1, col="white", border=NA)
rect(-2.58,0,0,1, col="whitesmoke", border=NA)

lines(-Temperature[,"Age"], speciation(Temperature[,"Age"]), type="l", col=alpha(colors[3], 0.7), lwd=2)
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -max(Temperature$Age)), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

axis(side=1, seq(-round(age, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, seq(0, speciation.max, by=0.1), col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)

add_geochrono(0, -1)
lines(-Temperature[,"Age"], extinction(Temperature[,"Age"]), type="l", col=alpha(colors[4], 0.7), lwd="2", lty="dashed")

legend("topright", bty="n", c("Speciation rate", "Extinction rate"), lwd=2,
       col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7)), lty=c(1, 2), cex=0.8)
mtext("BTempVarDTempVar_LIN  AICweight = 0.165", side=3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side=1, line=1.5, cex=0.8)


# --- Panel 3: CO2-dependent diversification (BCO2VarDCST_EXPO) ---
CO2 <- read.table("./data/500kyrCO2_combined_62.75Ma.txt", h=T)
CO2_data1 <- as.vector(scale(CO2[,2]))
CO2_data <- cbind(CO2[,1], CO2_data1)
colnames(CO2_data) <- c("Age", "CO2")
CO2 <- as.data.frame(CO2_data)
CO2_spline <- sm.spline(CO2[,1], CO2[,2])
CO2_time <- function(x){ predict(CO2_spline, x) }

speciation_par1 <- resi_all$BCO2VarDCST_EXPO$lamb_par[1]
speciation_par2 <- resi_all$BCO2VarDCST_EXPO$lamb_par[2]
speciation <- function(x){ abs(speciation_par1 * exp(speciation_par2 * CO2_time(x))) }

plot(age, age, type='n', ylim=c(0, round(max(speciation(CO2[,1])),1)),
     xlim=c(-age,0), axes=F, ylab='', xlab='', main='', bty='n', las=1, cex.axis=0.8, cex.lab=0.5)

rect(-66,-0.001,-max(CO2$Age),2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)
lines(-CO2[,1], speciation(CO2[,1]), type="l", col=alpha(colors[3], 0.7), lwd="2")
abline(v=c(-2.58, -5.33, -max(CO2$Age)), col="black", lty="dotted", lwd="1")
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

axis(side=1, seq(-round(age+1, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, seq(0, round(max(speciation(CO2[,"Age"])),1), by=0.1), col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)
add_geochrono(0, -1)
extinction_par1 <- abs(resi_all$BCO2VarDCST_EXPO$lamb_par[1])
extinction_par2 <- 0
lines(c(seq(-age,0,by=1),0), extinction_par1*exp(extinction_par2*c(seq(age,0,by=-1),0)),
      type="l", col=alpha(colors[4], 0.7), lwd="2", lty="dashed")
axis(side=2, seq(0, speciation.max, by=0.1), col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)

legend("topright", bty="n", c("Speciation rate", "Extinction rate"), lwd="3",
       col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7)), lty=c(1, 2), cex=0.8)
mtext("BCO2VarDCST_EXPO  AICweight = 0.022", side=3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side=1, line=1.7, cex=0.8)

dev.off()



## Net diversification rate plots

pdf(file="Net Diversification dynamics of Prunus.pdf", height=12, width=8)

par(mfrow=c(3,1), mar=c(3,3,0.5,0.5))

speciation_par1 <- resi_all$BTimeVarDTimeVar_LIN$lamb_par[1]
speciation_par2 <- resi_all$BTimeVarDTimeVar_LIN$lamb_par[2]

# Define speciation function (linear)
speciation <- function(time) {
  speciation_par1 + speciation_par2 * time
}

extinction_par1 <- resi_all$BTimeVarDTimeVar_LIN$mu_par[1]
extinction_par2 <- resi_all$BTimeVarDTimeVar_LIN$mu_par[2]

# Define extinction function (linear)
extinction <- function(time) {
  extinction_par1 + extinction_par2 * time
}

# Define net diversification function
net_diversification <- function(time) {
  speciation(time) - extinction(time)
}

# Plot
plot(age, age, type='n',
     ylim=c(0, round(max(speciation(c(seq(age, 0, by=-1), 0)),
                         extinction(c(seq(age, 0, by=-1), 0)),
                         net_diversification(c(seq(age, 0, by=-1), 0)), na.rm=TRUE), 3)),
     xlim=c(-age, 0), axes=F, ylab='', xlab='',
     cex.axis=0.5, cex.lab=0.5)

# Geological period background
rect(-66, 0, -age, 6.7, col="whitesmoke", border=NA)
rect(-56, 0, -66, 6.7, col="white", border=NA)
rect(-56, 0, -33.9, 6.7, col="whitesmoke", border=NA)
rect(-33.9, 0, -23.03, 6.7, col="white", border=NA)
rect(-23.03, 0, -5.33, 6.7, col="whitesmoke", border=NA)
rect(-5.33, 0, -2.58, 6.7, col="white", border=NA)
rect(-2.58, 0, 0, 6.7, col="whitesmoke", border=NA)

# Plot speciation, extinction, and net diversification
lines(c(seq(-age, 0, by=1), 0), speciation(c(seq(age, 0, by=-1), 0)),
      type="l", col=alpha(colors[3], 0.7), lwd=2)
lines(c(seq(-age, 0, by=1), 0), extinction(c(seq(age, 0, by=-1), 0)),
      type="l", col=alpha(colors[4], 0.7), lwd=2, lty="dashed")
lines(c(seq(-age, 0, by=1), 0), net_diversification(c(seq(age, 0, by=-1), 0)),
      type="l", col=alpha(colors[2], 0.7), lwd=2, lty="dashed")

# Geological period boundaries
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -age), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

axis(side=1, at=seq(-round(age, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, at=pretty(c(0, max_y)), col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)

legend("topright", bty="n",
       c("Speciation rate", "Extinction rate", "Net diversification rate"),
       lwd=2,
       col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7), alpha(colors[2], 0.7)),
       lty=c(1, 2, 2), cex=0.8)

add_geochrono(0, -0.2)
mtext("BTimeVarDTimeVar_LIN  AICweight = 0.767", side=3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side=1, line=1.5, cex=0.8)


# Temperature-dependent diversification (BTempVarDTempVar_LIN)
Temperature <- read.csv("./data/temperature-0.5ma-FABIEN-approx-63.csv", header=T)
Temperature_data1 <- as.vector(scale(Temperature[,2]))
Temperature_data <- cbind(Temperature[,1], Temperature_data1)
colnames(Temperature_data) <- c("Age", "Temperature")
Temperature <- as.data.frame(Temperature_data)

Tm_spline <- sm.spline(Temperature[,1], Temperature[,2])
Tm_time <- function(x) { predict(Tm_spline, x) }

speciation_par1 <- resi_all$BTempVarDTempVar_LIN$lamb_par[1]
speciation_par2 <- resi_all$BTempVarDTempVar_LIN$lamb_par[2]
speciation <- function(x) { abs(speciation_par1 + speciation_par2 * Tm_time(x)) }
speciation.min <- round(min(speciation(Temperature[,"Age"])),3)
speciation.max <- round(max(speciation(Temperature[,"Age"])),3)

extinction_par1 <- resi_all$BTempVarDTempVar_LIN$mu_par[1]
extinction_par2 <- resi_all$BTempVarDTempVar_LIN$mu_par[2]
extinction <- function(x) { abs(extinction_par1 + extinction_par2 * Tm_time(x)) }
extinction.min <- round(min(extinction(Temperature[,"Age"])),3)
extinction.max <- round(max(extinction(Temperature[,"Age"])),3)

# Net diversification (set negative values to 0)
net_diversification <- function(x) {
  nd <- speciation(x) - extinction(x)
  nd[nd < 0] <- 0
  return(nd)
}

plot(age, age, type='n',
     ylim=c(0, max(speciation.max, extinction.max, net_div.max)),
     xlim=c(-age,0), axes=F, ylab='', xlab='',
     bty='n', las=1, cex.axis=0.5, cex.lab=0.5)

rect(-66,0,-max(Temperature$Age),1, col="whitesmoke", border=NA)
rect(-56,0,-66,1, col="white", border=NA)
rect(-56,0,-33.9,1, col="whitesmoke", border=NA)
rect(-33.9,0,-23.03,1, col="white", border=NA)
rect(-23.03,0,-5.33,1, col="whitesmoke", border=NA)
rect(-5.33,0,-2.58,1, col="white", border=NA)
rect(-2.58,0,0,1, col="whitesmoke", border=NA)

lines(-Temperature[,"Age"], speciation(Temperature[,"Age"]), type="l", col=alpha(colors[3], 0.7), lwd=2)
lines(-Temperature[,"Age"], extinction(Temperature[,"Age"]), type="l", col=alpha(colors[4], 0.7), lwd=2, lty="dashed")
lines(-Temperature[,"Age"], net_diversification(Temperature[,"Age"]), type="l", col=alpha(colors[2], 0.7), lwd=2, lty="dashed")

abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -max(Temperature$Age)), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

axis(side=1, seq(-round(age, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, seq(0, max(speciation.max, extinction.max, net_div.max), by=0.1),
     col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)

add_geochrono(0, -0.1)
legend("topright", bty="n",
       c("Speciation rate", "Extinction rate", "Net diversification rate"),
       lwd=2,
       col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7), alpha(colors[2], 0.7)),
       lty=c(1, 2, 2), cex=0.8)
mtext("BTempVarDTempVar_LIN  AICweight = 0.165", side=3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side=1, line=1.5, cex=0.8)

# CO2-dependent diversification (BCO2VarDCST_EXPO)
CO2 <- read.table("./data/500kyrCO2_combined_62.75Ma.txt", h=T)
CO2_data1 <- as.vector(scale(CO2[,2]))
CO2_data <- cbind(CO2[,1], CO2_data1)
colnames(CO2_data) <- c("Age", "CO2")
CO2 <- as.data.frame(CO2_data)

CO2_spline <- sm.spline(CO2[,1], CO2[,2])
CO2_time <- function(x) { predict(CO2_spline, x) }

speciation_par1 <- resi_all$BCO2VarDCST_EXPO$lamb_par[1]
speciation_par2 <- resi_all$BCO2VarDCST_EXPO$lamb_par[2]
speciation <- function(x) { abs(speciation_par1 * exp(speciation_par2 * CO2_time(x))) }
speciation.max <- round(max(speciation(CO2[,1])), 1)

extinction_par1 <- abs(resi_all$BCO2VarDCST_EXPO$mu_par[1])
extinction_par2 <- 0
extinction <- function(x) { extinction_par1 * exp(extinction_par2 * CO2_time(x)) }
extinction.max <- round(max(extinction(CO2[,1])), 1)

net_diversification <- function(x) { speciation(x) - extinction(x) }
net_div.max <- round(max(net_diversification(CO2[,1])), 1)
net_div.min <- round(min(net_diversification(CO2[,1])), 1)

plot(age, age, type='n',
     ylim=c(0, max(speciation.max, extinction.max)),
     xlim=c(-age, 0), axes=F, ylab='', xlab='',
     main='', bty='n', las=1, cex.axis=0.8, cex.lab=0.5)

rect(-66, -0.001, -max(CO2$Age), max(speciation.max, extinction.max, net_div.max), col="whitesmoke", border=NA)
rect(-56, -0.001, -66, max(speciation.max, extinction.max, net_div.max), col="white", border=NA)
rect(-56, -0.001, -33.9, max(speciation.max, extinction.max, net_div.max), col="whitesmoke", border=NA)
rect(-33.9, -0.001, -23.03, max(speciation.max, extinction.max, net_div.max), col="white", border=NA)
rect(-23.03, -0.001, -5.33, max(speciation.max, extinction.max, net_div.max), col="whitesmoke", border=NA)
rect(-5.33, -0.001, -2.58, max(speciation.max, extinction.max, net_div.max), col="white", border=NA)
rect(-2.58, -0.001, 0, max(speciation.max, extinction.max, net_div.max), col="whitesmoke", border=NA)

lines(-CO2[,1], speciation(CO2[,1]), type="l", col=alpha(colors[3], 0.7), lwd=2)
lines(-CO2[,1], extinction(CO2[,1]), type="l", alpha(colors[4], 0.7), lwd=2, lty="dashed")
lines(-CO2[,1], net_diversification(CO2[,1]), type="l", col=alpha(colors[2], 0.7), lwd=2, lty="dashed")

abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -max(CO2$Age)), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

axis(side=1, seq(-round(age+1, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, seq(0, max(speciation.max, extinction.max, net_div.max), by=0.1),
     col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)

add_geochrono(0, -0.2)

legend("topright", bty="n",
       c("Speciation rate", "Extinction rate", "Net diversification rate"),
       lwd=2,
       col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7), alpha(colors[2], 0.7)),
       lty=c(1, 2, 2), cex=0.8)

mtext("BCO2VarDCST_EXPO  AICweight = 0.022", side=3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side=1, line=1.7, cex=0.8)
dev.off()
