## Plot of the uplift dynamic of the Andes
library(pspline); library(scales)

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

## Temperatures
Temperature<-read.table("../Environmental data/merged_veizer_westerhold_Ts.txt",header=T)

#par(mfrow=c(1,2), mar=c(4, 5, 1, 1))
par(mar=c(4, 5, 1, 1))
age=(0:(80-1))* -1

colors<-c("dodgerblue","blue")

plot(age,age,type = 'n', ylim = c(0, 35), xlim = c(-80,0), axes=F, ylab = '', xlab = 'Time (Ma)',main='', bty='n', las=1)
rect(-100.5,-100,-105, 5000, col="white", border=NA)
rect(-66,-100,-100.5, 5000, col="whitesmoke", border=NA)
rect(-56,-100,-66, 5000, col="white", border=NA)
rect(-56,-100,-33.9, 5000, col="whitesmoke", border=NA)
rect(-33.9,-100,-23.03, 5000, col="white", border=NA)
rect(-23.03,-100,-5.33, 5000, col="whitesmoke", border=NA)
rect(-5.33,-100,-2.58, 5000, col="white", border=NA)
rect(-2.58,-100,0, 5000, col="whitesmoke", border=NA)

abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-90, 0, by=10), col="black",col.axis="black",las=1)
axis(side=2, seq(0, 35, by=5), col="blue",col.axis="blue",las=2)
mtext("Global sea-surface temperatures (°C)",side=2,col="blue",line=4)

points(-Temperature[,"Age"], Temperature[,"Temperature"], pch=19, lwd=0.5, col=alpha(colors[1], 0.15))
temp.spl<-smooth.spline(-Temperature[,"Age"], Temperature[,"Temperature"])
lines(temp.spl, col=alpha(colors[2], 0.7), lwd=3)

add_geochrono(0, -10)
dev.print(pdf, "Fig. 3 - Global paleotemperatures.pdf")

## Andes
Andes<-read.table("../Environmental data/Andes_mean_elevations_no_basins_ALL.txt", h=T)
North<-read.table("../Environmental data/Andes_mean_elevations_no_basins_NORTHCENTRAL.txt", h=T)
South<-read.table("../Environmental data/Andes_mean_elevations_no_basins_CENTRALSOUTH.txt", h=T)

#par(new = TRUE)

par(mar=c(4, 5, 1, 1))
age=(0:(80-1))* -1

colors<-c("red","red", "chartreuse3", "dodgerblue")
#plot(age,age,type = 'n', ylim = c(-100, 4500), xlim = c(-80,0), axes=F, ylab = '', xlab = '',main='', bty='n', las=1)

plot(age,age,type = 'n', ylim = c(-100, 4500), xlim = c(-80,0), axes=F, ylab = '', xlab = 'Time (Ma)',main='', bty='n', las=1)
rect(-100.5,-100,-105, 5000, col="white", border=NA)
rect(-66,-100,-100.5, 5000, col="whitesmoke", border=NA)
rect(-56,-100,-66, 5000, col="white", border=NA)
rect(-56,-100,-33.9, 5000, col="whitesmoke", border=NA)
rect(-33.9,-100,-23.03, 5000, col="white", border=NA)
rect(-23.03,-100,-5.33, 5000, col="whitesmoke", border=NA)
rect(-5.33,-100,-2.58, 5000, col="white", border=NA)
rect(-2.58,-100,0, 5000, col="whitesmoke", border=NA)

abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -100.5),col="black",lty="dotted",lwd="1")
abline(v=c(-66),col="red",lty="dashed",lwd="2")

axis(side=1, seq(-90, 0, by=10), col="black",col.axis="black",las=1)
axis(side=2, seq(-100, 4500, by=200), col="black",col.axis="black",las=2)
mtext("Andean elevation (m)",side=2,col="black",line=4)

points(-Andes[,"Age"], Andes[,"Altitude"], pch=17, lwd=0.5, col=alpha(colors[1], 0.2))
temp.spl<-smooth.spline(-Andes[,"Age"], Andes[,"Altitude"], df=80)
lines(temp.spl, col=alpha(colors[2], 0.7), lwd=4)
#abline(h=0)

temp.spl<-smooth.spline(-North[,"Age"], North[,"Altitude"], df=80)
lines(temp.spl, col=alpha(colors[3], 0.7), lwd=4)

temp.spl<-smooth.spline(-South[,"Age"], South[,"Altitude"], df=80)
lines(temp.spl, col=alpha(colors[4], 0.7), lwd=4)

legend("topleft", c("All Andes", "Northern and Central Andes", "Southern and Central Andes"), lwd=3, col=c("red", "chartreuse3", "dodgerblue"), bty="n")
add_geochrono(-100,-300)

dev.print(pdf, "Fig. 2 - Uplift history of the Andes.pdf")
