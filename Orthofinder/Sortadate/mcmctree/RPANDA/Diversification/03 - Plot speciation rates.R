##
## Plot Speciation rates
##

library(picante); library(pspline); library(scales)

## Fit of the models
attach('./Results_CO2_Fab_PDDM.Rdata')

## Temperature data
Temperature <- read.table("./data/temperature.txt", h=T)

## Geological time scale
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

## Plot the rates through time

colors<-c("red","chartreuse3","dodgerblue","orange","purple","brown")

##################
## Aromobatidae ##
##################

pdf(file="Diversification dynamics of Aromobatidae.pdf", height=3.8, width=4)
i=1
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5))

age<-all_res[[i]]$Clade_age

# Environment-independent diversification
speciation_par1<-all_res[[i]]$BTimeVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BTimeVar_EXPO$lamb_par[2]

plot(age,age, type='n', ylim=c(0, round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),1)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(c(seq(-age,0,by=1),0), speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0)), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),1), by=0.02), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

# Environment-dependent diversification
Andes <- read.table("../Environmental data/Andes_mean_elevations_no_basins_NORTHCENTRAL.txt", h=T)
Andes_subset<-subset(Andes, Age<age)
Andes_spline<-sm.spline(Andes_subset[,"Age"], Andes_subset[,"Altitude"])
Andes_time<-function(x){predict(Andes_spline,x)}

speciation_par1<-all_res[[i]]$BAndesVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BAndesVar_EXPO$lamb_par[2]
speciation<-function(x){speciation_par1*exp(speciation_par2*Andes_time(x))}

plot(age,age, type='n', ylim=c(0, round(max(speciation(Andes_subset[,"Age"])),1)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(-Andes_subset[,"Age"], speciation(Andes_subset[,"Age"]), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, round(max(speciation(Andes_subset[,"Age"])),1), by=0.02), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -1)

dev.off()
#dev.print(pdf, file="Diversification dynamics of Aromobatidae.pdf")


###################
## Centrolenidae ##
###################

pdf(file="Diversification dynamics of Centrolenidae.pdf", height=3.8, width=4)
i=2
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5))

age<-all_res[[i]]$Clade_age

# Environment-independent diversification
speciation_par1<-all_res[[i]]$BTimeVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BTimeVar_EXPO$lamb_par[2]

plot(age,age, type='n', ylim=c(0, round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(c(seq(-age,0,by=1),0), speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0)), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2), by=0.05), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

# Environment-dependent diversification
Temp_subset<-subset(Temperature, Age<age)
Temp_spline<-sm.spline(Temp_subset[,"Age"], Temp_subset[,"Temperature"], df=66)
Temp_time<-function(x){predict(Temp_spline,x)}

speciation_par1<-all_res[[i]]$BTempVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BTempVar_EXPO$lamb_par[2]
speciation<-function(x){speciation_par1*exp(speciation_par2*Temp_time(x))}

plot(age,age, type='n', ylim=c(0, round(max(speciation(Temp_subset[,"Age"])),2)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(-Temp_subset[,"Age"], speciation(Temp_subset[,"Age"]), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, round(max(speciation(Temp_subset[,"Age"])),2), by=0.02), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

dev.off()
#dev.print(pdf, file="Diversification dynamics of Centrolenidae.pdf")


###################
## Dendrobatidae ##
###################

pdf(file="Diversification dynamics of Dendrobatidae.pdf", height=3.8, width=4)
i=3
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5))

age<-all_res[[i]]$Clade_age

# Environment-independent diversification
speciation_par1<-all_res[[i]]$BCST$lamb_par[1]
speciation_par2<-0

plot(age,age, type='n', ylim=c(0, (round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)+0.02)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(c(seq(-age,0,by=1),0), speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0)), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, (round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)+0.02), by=0.05), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

# Environment-dependent diversification
Temp_subset<-subset(Temperature, Age<age)
Temp_spline<-sm.spline(Temp_subset[,"Age"], Temp_subset[,"Temperature"], df=66)
Temp_time<-function(x){predict(Temp_spline,x)}

speciation_par1<-all_res[[i]]$BTempVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BTempVar_EXPO$lamb_par[2]
speciation<-function(x){speciation_par1*exp(speciation_par2*Temp_time(x))}

plot(age,age, type='n', ylim=c(0, round(max(speciation(Temp_subset[,"Age"])),2)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(-Temp_subset[,"Age"], speciation(Temp_subset[,"Age"]), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, round(max(speciation(Temp_subset[,"Age"])),2), by=0.02), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

dev.off()
#dev.print(pdf, file="Diversification dynamics of Dendrobatidae.pdf")


####################
## Hemiphractidae ##
####################

pdf(file="Diversification dynamics of Hemiphractidae.pdf", height=3.8, width=4)
i=4
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5))

age<-all_res[[i]]$Clade_age

# Environment-independent diversification
speciation_par1<-all_res[[i]]$BCST$lamb_par[1]
speciation_par2<-0

plot(age,age, type='n', ylim=c(0, (round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)+0.02)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(c(seq(-age,0,by=1),0), speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0)), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, (round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)+0.02), by=0.01), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

# Environment-dependent diversification
Andes <- read.table("../Environmental data/Andes_mean_elevations_no_basins_NORTHCENTRAL.txt", h=T)
Andes_subset<-subset(Andes, Age<age)
Andes_spline<-sm.spline(Andes_subset[,"Age"], Andes_subset[,"Altitude"])
Andes_time<-function(x){predict(Andes_spline,x)}

speciation_par1<-all_res[[i]]$BAndesVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BAndesVar_EXPO$lamb_par[2]
speciation<-function(x){speciation_par1*exp(speciation_par2*Andes_time(x))}

plot(age,age, type='n', ylim=c(0, round(max(speciation(Andes_subset[,"Age"])),1)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(-Andes_subset[,"Age"], speciation(Andes_subset[,"Age"]), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, round(max(speciation(Andes_subset[,"Age"])),1), by=0.02), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -1)

dev.off()
#dev.print(pdf, file="Diversification dynamics of Hemiphractidae.pdf")


#####################
## Leptodactylidae ##
#####################

pdf(file="Diversification dynamics of Leptodactylidae.pdf", height=3.8, width=4)
i=5
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5))

age<-all_res[[i]]$Clade_age

# Environment-independent diversification
speciation_par1<-all_res[[i]]$BCST$lamb_par[1]
speciation_par2<-0

plot(age,age, type='n', ylim=c(0, (round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)+0.02)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(c(seq(-age,0,by=1),0), speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0)), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, (round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)+0.02), by=0.01), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

# Environment-dependent diversification
Andes <- read.table("../Environmental data/Andes_mean_elevations_no_basins_ALL.txt", h=T)
Andes_subset<-subset(Andes, Age<age)
Andes_spline<-sm.spline(Andes_subset[,"Age"], Andes_subset[,"Altitude"])
Andes_time<-function(x){predict(Andes_spline,x)}

speciation_par1<-all_res[[i]]$BAndesVarDAndesVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BAndesVarDAndesVar_EXPO$lamb_par[2]
speciation<-function(x){speciation_par1*exp(speciation_par2*Andes_time(x))}

plot(age,age, type='n', ylim=c(0, round(max(speciation(Andes_subset[,"Age"])),2)+0.03), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(-Andes_subset[,"Age"], speciation(Andes_subset[,"Age"]), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, round(max(speciation(Andes_subset[,"Age"])),2)+0.03, by=0.05), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

extinction_par1<-all_res[[i]]$BAndesVarDAndesVar_EXPO$mu_par[1]
extinction_par2<-all_res[[i]]$BAndesVarDAndesVar_EXPO$mu_par[2]
extinction<-function(x){extinction_par1*exp(extinction_par2*Andes_time(x))}
lines(-Andes_subset[,"Age"], extinction(Andes_subset[,"Age"]), type="l", col=alpha(colors[i], 0.7),lwd="4", lty="dashed")

legend("top", bty="n", c("Speciation rate", "Extinction rate"),lwd="2",col=c(alpha(colors[i], 0.7), alpha(colors[i], 0.7)), lty = c(1, 2), cex=0.6)

dev.off()
#dev.print(pdf, file="Diversification dynamics of Leptodactylidae.pdf")


#################
## Liolaemidae ##
#################

pdf(file="Diversification dynamics of Liolaemidae.pdf", height=3.8, width=4)
i=6
par(mfrow=c(2,1), mar=c(2,2,0.5,0.5))

age<-all_res[[i]]$Clade_age

# Environment-independent diversification
speciation_par1<-all_res[[i]]$BTimeVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BTimeVar_EXPO$lamb_par[2]

plot(age,age, type='n', ylim=c(0, (round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)+0.05)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(c(seq(-age,0,by=1),0), speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0)), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, (round(max(speciation_par1*exp(speciation_par2*c(seq(age,0,by=-1),0))),2)+0.05), by=0.1), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

# Environment-dependent diversification
Andes <- read.table("../Elevation data/Andes_mean_elevations_no_basins_CENTRALSOUTH.txt", h=T) 
Andes_subset<-subset(Andes, Age<age)
Andes_spline<-sm.spline(Andes_subset[,"Age"], Andes_subset[,"Altitude"])
Andes_time<-function(x){predict(Andes_spline,x)}

speciation_par1<-all_res[[i]]$BAndesVar_EXPO$lamb_par[1]
speciation_par2<-all_res[[i]]$BAndesVar_EXPO$lamb_par[2]
speciation<-function(x){speciation_par1*exp(speciation_par2*Andes_time(x))}

plot(age,age, type='n', ylim=c(0, (round(max(speciation(Andes_subset[,"Age"])),2)+0.06)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.5, cex.lab=0.5)
rect(-66,-0.001,-100.5,2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)

lines(-Andes_subset[,"Age"], speciation(Andes_subset[,"Age"]), type="l", col=alpha(colors[i], 0.7),lwd="4")
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-round(age+10, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.5, padj=-2)
axis(side=2, seq(0, (round(max(speciation(Andes_subset[,"Age"])),2)+0.06), by=0.1), col="black",col.axis="black",las=2, cex.axis=0.5, hadj=0.5)
add_geochrono(0, -0.2)

dev.off()
#dev.print(pdf, file="Diversification dynamics of Liolaemidae.pdf")



##
## Plot phylogenetic trees
##

library(picante); library(strap)

##################
## Aromobatidae ##
##################
tree <-read.tree("Aromobatidae.tre")
tree$root.time <- max(node.age(tree)$ages)

pdf(file="Phylogeny of Aromobatidae.pdf", height=3.8, width=4)
geoscalePhylo(tree=ladderize(tree,right=T), show.tip.label=F, units=c("Period","Epoch"), boxes="Epoch", vers="ICS2015", erotate=90, arotate=0, cex.tip=0.3, cex.age=0.6, cex.ts=0.4, label.offset=0.5, x.lim=c(-1, round(tree$root.time, -1)), lwd=1, width=1, quat.rm=T)
dev.off()

###################
## Centrolenidae ##
###################

tree <-read.tree("Centrolenidae.tre")
tree$root.time <- max(node.age(tree)$ages)

pdf(file="Phylogeny of Centrolenidae.pdf", height=3.8, width=4)
geoscalePhylo(tree=ladderize(tree,right=T), show.tip.label=F, units=c("Period","Epoch"), boxes="Epoch", vers="ICS2015", erotate=90, arotate=0, cex.tip=0.3, cex.age=0.6, cex.ts=0.4, label.offset=0.5, x.lim=c(-1, (round(tree$root.time)+2)), lwd=1, width=1, quat.rm=T)
dev.off()

###################
## Dendrobatidae ##
###################

tree <-read.tree("Dendrobatidae.tre")
tree$root.time <- max(node.age(tree)$ages)

pdf(file="Phylogeny of Dendrobatidae.pdf", height=3.8, width=4)
geoscalePhylo(tree=ladderize(tree,right=T), show.tip.label=F, units=c("Period","Epoch"), boxes="Epoch", vers="ICS2015", erotate=90, arotate=0, cex.tip=0.3, cex.age=0.6, cex.ts=0.4, label.offset=0.5, x.lim=c(-1, round(tree$root.time, -1)), lwd=1, width=1, quat.rm=T)
dev.off()

####################
## Hemiphractidae ##
####################

tree <-read.tree("Hemiphractidae.tre")
tree$root.time <- max(node.age(tree)$ages)

pdf(file="Phylogeny of Hemiphractidae.pdf", height=3.8, width=4)
geoscalePhylo(tree=ladderize(tree,right=T), show.tip.label=F, units=c("Period","Epoch"), boxes="Epoch", vers="ICS2015", erotate=90, arotate=0, cex.tip=0.3, cex.age=0.6, cex.ts=0.4, label.offset=0.5, x.lim=c(-1, round(tree$root.time, -1)), lwd=1, width=1, quat.rm=T)
dev.off()

#####################
## Leptodactylidae ##
#####################

tree <-read.tree("Leptodactylidae.tre")
tree$root.time <- max(node.age(tree)$ages)

pdf(file="Phylogeny of Leptodactylidae.pdf", height=3.8, width=4)
geoscalePhylo(tree=ladderize(tree,right=T), show.tip.label=F, units=c("Period","Epoch"), boxes="Epoch", vers="ICS2015", erotate=90, arotate=0, cex.tip=0.3, cex.age=0.6, cex.ts=0.4, label.offset=0.5, x.lim=c(-1, round(tree$root.time, -1)), lwd=1, width=1, quat.rm=T)
dev.off()

#################
## Liolaemidae ##
#################

tree <-read.tree("Liolaemidae.tre")
tree$root.time <- max(node.age(tree)$ages)

pdf(file="Phylogeny of Liolaemidae.pdf", height=3.8, width=4)
geoscalePhylo(tree=ladderize(tree,right=T), show.tip.label=F, units=c("Period","Epoch"), boxes="Epoch", vers="ICS2015", erotate=90, arotate=0, cex.tip=0.3, cex.age=0.6, cex.ts=0.4, label.offset=0.5, x.lim=c(-1, round(tree$root.time, -1)), lwd=1, width=1, quat.rm=T)
dev.off()
