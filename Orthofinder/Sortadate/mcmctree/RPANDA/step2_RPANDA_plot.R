setwd("C:/Users/301/Desktop/李属文章数据/BAMM_35_data_5/RPANDA")

list.files()

library(pspline)
library(ggplot2)
library(picante)
library(scales)
library(RPANDA)


#/data/xueqin/Project/Prunus/Prunus_353_easy353/add_Prunus_2024.7.24/08_27_delete_lower_50_FNA/Prunus_0.12/Prunus_gene_trimal/single_gene_iqtree/Prunus_treeshrink/Prunus_sortadate_3,1,2/tree_sortadate_35/Prunus_mcmctree/approx01/BAMM/Diverstification/RPANDA/rpanda3/refs/final_all_model_RPANDA_2025.6.12/final.tactNA_all_models.Rdata
###############################
## Plots for the best models ##
###############################
load("tactNA_all_res_list.Rdata")
#attach('final.tactNA_all_models.Rdata')
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

colors<-c("chartreuse3","dodgerblue","orange","darkgray")

##################
## Fabaceae ##
##################

pdf(file="Diversification dynamics of Prunus.pdf", height=12, width=8)

#par(mfrow=c(3,1), mar=c(2,2,0.5,0.5))
par(mfrow=c(3,1), mar=c(3,3,0.5,0.5))
resi_all <- Clematis_res_all  # 你前面用了 list(resi_all)，所以需要用 [[1]]
age <- resi_all$Clade_age  

# Environment-independent diversification
# BTimeVarDTimeVar_LIN
#plot_fit_bd(resi_all$BTimeVarDTimeVar_LIN,resi_all$Clade_age)

#?plot_fit_bd

speciation_par1 <- resi_all$BTimeVarDTimeVar_LIN$lamb_par[1]
speciation_par2 <- resi_all$BTimeVarDTimeVar_LIN$lamb_par[2]

# Define the speciation function based on time (linear relationship)
speciation <- function(time) {
  speciation_par1 + speciation_par2 * time
}

# Define the extinction function based on time (linear relationship)
extinction_par1 <- resi_all$BTimeVarDTimeVar_LIN$mu_par[1]
extinction_par2 <- resi_all$BTimeVarDTimeVar_LIN$mu_par[2]
extinction <- function(time) {
  extinction_par1 + extinction_par2 * time
}

# Plot
plot(age, age, type='n', 
     ylim=c(0, round(max(speciation(c(seq(age, 0, by=-1), 0)), 
                         extinction(c(seq(age, 0, by=-1), 0)), na.rm = TRUE), 3)), 
     xlim=c(-age, 0), axes=F, ylab='', xlab='', 
     #main=paste0(names, "--Tm"), bty='n', las=1, 
     cex.axis=0.5, cex.lab=0.5)

rect(-66, 0, -age, 6.7, col="whitesmoke", border=NA)
rect(-56, 0, -66, 6.7, col="white", border=NA)
rect(-56, 0, -33.9, 6.7, col="whitesmoke", border=NA)
rect(-33.9, 0, -23.03, 6.7, col="white", border=NA)
rect(-23.03, 0, -5.33, 6.7, col="whitesmoke", border=NA)
rect(-5.33, 0, -2.58, 6.7, col="white", border=NA)
rect(-2.58, 0, 0, 6.7, col="whitesmoke", border=NA)

# Plot speciation and extinction rates
lines(c(seq(-age, 0, by=1), 0), speciation(c(seq(age, 0, by=-1), 0)), 
      type="l", col=alpha(colors[3], 0.7), lwd=2)
lines(c(seq(-age, 0, by=1), 0), extinction(c(seq(age, 0, by=-1), 0)), 
      type="l", col=alpha(colors[4], 0.7), lwd=2, lty="dashed")

# Add vertical lines for geological periods
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -age), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

# Add axes
axis(side=1, at=seq(-round(age, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)

# Calculate the maximum value for the y-axis
max_speciation <- round(max(speciation(c(seq(age, 0, by=-1), 0)), na.rm = TRUE), 3)
max_extinction <- round(max(extinction(c(seq(age, 0, by=-1), 0)), na.rm = TRUE), 3)
max_y <- max(max_speciation, max_extinction)

# Add y-axis
axis(side=2, at=pretty(c(0, max_y)), col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)
legend("topright", bty="n", c("Speciation rate", "Extinction rate"),lwd="2",col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7)), lty = c(1, 2), cex=0.8)

# Add geological time scale
add_geochrono(0, -0.2)
# Add title
mtext("BTimeVarDTimeVar_LIN  AICweight = 0.767", side = 3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side = 1, line = 1.5, cex=0.8)

# temperature-dependent diversification
#BTempVarDTempVar_LIN
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

### Extinction
extinction_par1 <- resi_all$BTempVarDTempVar_LIN$mu_par[1]
extinction_par2 <- resi_all$BTempVarDTempVar_LIN$mu_par[2]
extinction <- function(x) { abs(extinction_par1 + extinction_par2 * Tm_time(x)) }
extinction.min <- round(min(extinction(Temperature[,"Age"])),3)
extinction.max <- round(max(extinction(Temperature[,"Age"])),3)

# Plot
plot(age, age, type='n', ylim=c(0, max(speciation.max, extinction.max)), xlim=c(-age,0), axes=F, ylab='', xlab='', 
     #main=paste0(names, "--Tm"), 
     bty='n', las=1, cex.axis=0.5, cex.lab=0.5)

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

# Add X-axis
axis(side=1, seq(-round(age, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, seq(0, speciation.max, by=0.1), col="black",col.axis="black",las=2, cex.axis=0.8, hadj=0.5)
# Add Y-axis with specific tick marks
# y_ticks <- seq(0, round(max(speciation(Temperature[,"Age"])), 3), by=0.05)  # Define specific tick marks
# axis(side=2, at=y_ticks, col="black", col.axis="black", las=2, cex.axis=0.5, hadj=0.7)

# Add geological time scale
add_geochrono(0, -1)
lines(-Temperature[,"Age"], extinction(Temperature[,"Age"]), type="l", col=alpha(colors[4], 0.7), lwd="2", lty="dashed")
#segments(-max(Temperature[,"Age"]), extinction(Temperature[,"Age"]), 0, extinction(Temperature[,"Age"]), col=scales::alpha("black", 0.7), lwd=2, lty="dashed")

# Add legend
legend("topright", bty="n", c("Speciation rate", "Extinction rate"), lwd=2, col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7)), lty = c(1, 2), cex=0.8)

# Add text
mtext("BTempVarDTempVar_LIN  AICweight = 0.165", side = 3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side = 1, line = 1.5, cex=0.8)


# CO2-dependent diversification
# BCO2VarDCST_EXPO
CO2 <- read.table("./data/500kyrCO2_combined_62.75Ma.txt", h=T)
# Environmental data
#CO2_data0 <- read.csv("../data/ea49_age.csv", header = TRUE)
CO2_data1 <- as.vector(scale(CO2[,2]))
CO2_data <- cbind(CO2[,1], CO2_data1)
colnames(CO2_data) <- c("Age", "CO2")
CO2 <- as.data.frame(CO2_data)
# Andes_subset<-subset(Andes, Age<age)
CO2_spline<-sm.spline(CO2[,1], CO2[,2])
CO2_time<-function(x){predict(CO2_spline,x)}

# speciation
speciation_par1 <- resi_all$BCO2VarDCST_EXPO$lamb_par[1]
speciation_par2 <- resi_all$BCO2VarDCST_EXPO$lamb_par[2]
#speciation <- function(x) { speciation_par1 * exp(speciation_par2 * CO2_time(x)) }
speciation<-function(x){abs(speciation_par1*exp(speciation_par2*CO2_time(x)))}

plot(age,age, type='n', ylim=c(0, round(max(speciation(CO2[,1])),1)), xlim=c(-age,0), axes=F, ylab='', xlab='',main='', bty='n', las=1, cex.axis=0.8, cex.lab=0.5)


# 背景色块保持不变
rect(-66,-0.001,-max(CO2$Age),2, col="whitesmoke", border=NA)
rect(-56,-0.001,-66,2, col="white", border=NA)
rect(-56,-0.001,-33.9,2, col="whitesmoke", border=NA)
rect(-33.9,-0.001,-23.03,2, col="white", border=NA)
rect(-23.03,-0.001,-5.33,2, col="whitesmoke", border=NA)
rect(-5.33,-0.001,-2.58,2, col="white", border=NA)
rect(-2.58,-0.001,0,2, col="whitesmoke", border=NA)
lines(-CO2[,1], speciation(CO2[,1]), type="l", col=alpha(colors[3], 0.7),lwd="2")
abline(v=c(-2.58, -5.33, -max(CO2$Age)),col="black",lty="dotted",lwd="1")
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

axis(side=1, seq(-round(age+1, -1), 0, by=10), col="black",col.axis="black",las=1, cex.axis=0.8,padj=-1)
axis(side=2, seq(0, round(max(speciation(CO2[,"Age"])),1), by=0.1), col="black",col.axis="black",las=2, cex.axis=0.8, hadj=0.5)
add_geochrono(0, -1)
extinction_par1<-abs(resi_all$BCO2VarDCST_EXPO$lamb_par[1])
extinction_par2<-0
lines(c(seq(-age,0,by=1),0), extinction_par1*exp(extinction_par2*c(seq(age,0,by=-1),0)), type="l", col=alpha(colors[4], 0.7),lwd="2", lty="dashed")
axis(side=2, seq(0, speciation.max, by=0.1), col="black",col.axis="black",las=2, cex.axis=0.8, hadj=0.5)
#lines(-CO2[,"Age"], extinction(CO2[,"Age"]), type="l", col=alpha(colors[4], 0.7),lwd="2", lty="dashed")

legend("topright", bty="n", c("Speciation rate", "Extinction rate"),lwd="3",col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7)), lty = c(1, 2), cex=0.8)
# 标题
mtext("BCO2VarDCST_EXPO  AICweight = 0.022", 
      side=3, at=-40, line=-3, cex=0.8)
# mtext("Speciation rate (events/lineage/Myr)", 
#       side=2, line=2, cex=0.6)  # 调整line参数使标题不重叠
mtext("Time (Myr)", side = 1, line = 1.7, cex=0.8)

dev.off() 



##add the net diversification

pdf(file="Net Diversification dynamics of Prunus.pdf", height=12, width=8)

#par(mfrow=c(3,1), mar=c(2,2,0.5,0.5))
par(mfrow=c(3,1), mar=c(3,3,0.5,0.5))

speciation_par1 <- resi_all$BTimeVarDTimeVar_LIN$lamb_par[1]
speciation_par2 <- resi_all$BTimeVarDTimeVar_LIN$lamb_par[2]

# 定义 speciation 函数（线性）
speciation <- function(time) {
  speciation_par1 + speciation_par2 * time
}

extinction_par1 <- resi_all$BTimeVarDTimeVar_LIN$mu_par[1]
extinction_par2 <- resi_all$BTimeVarDTimeVar_LIN$mu_par[2]

# 定义 extinction 函数（线性）
extinction <- function(time) {
  extinction_par1 + extinction_par2 * time
}

# 定义 net diversification 函数
net_diversification <- function(time) {
  speciation(time) - extinction(time)
}

# 绘图
plot(age, age, type='n', 
     ylim=c(0, round(max(speciation(c(seq(age, 0, by=-1), 0)), 
                         extinction(c(seq(age, 0, by=-1), 0)), 
                         net_diversification(c(seq(age, 0, by=-1), 0)), na.rm = TRUE), 3)), 
     xlim=c(-age, 0), axes=F, ylab='', xlab='', 
     cex.axis=0.5, cex.lab=0.5)

# 背景色块（地质时期）
rect(-66, 0, -age, 6.7, col="whitesmoke", border=NA)
rect(-56, 0, -66, 6.7, col="white", border=NA)
rect(-56, 0, -33.9, 6.7, col="whitesmoke", border=NA)
rect(-33.9, 0, -23.03, 6.7, col="white", border=NA)
rect(-23.03, 0, -5.33, 6.7, col="whitesmoke", border=NA)
rect(-5.33, 0, -2.58, 6.7, col="white", border=NA)
rect(-2.58, 0, 0, 6.7, col="whitesmoke", border=NA)

# 绘制 speciation、extinction 和 net diversification
lines(c(seq(-age, 0, by=1), 0), speciation(c(seq(age, 0, by=-1), 0)), 
      type="l", col=alpha(colors[3], 0.7), lwd=2)
lines(c(seq(-age, 0, by=1), 0), extinction(c(seq(age, 0, by=-1), 0)), 
      type="l", col=alpha(colors[4], 0.7), lwd=2, lty="dashed")
lines(c(seq(-age, 0, by=1), 0), net_diversification(c(seq(age, 0, by=-1), 0)), 
      type="l", col=alpha(colors[2], 0.7), lwd=2, lty="dashed")

# 地质时期分界线
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -age), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

# 坐标轴
axis(side=1, at=seq(-round(age, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, at=pretty(c(0, max_y)), col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)

# 图例（更新，包含 net diversification）
legend("topright", bty="n", 
       c("Speciation rate", "Extinction rate", "Net diversification rate"),
       lwd=2, 
       col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7), alpha(colors[2], 0.7)), 
       lty=c(1, 2, 2), 
       cex=0.8)

# 添加地质年代
add_geochrono(0, -0.2)

# 标题
mtext("BTimeVarDTimeVar_LIN  AICweight = 0.767", side = 3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side = 1, line = 1.5, cex=0.8)


# temperature-dependent diversification
#BTempVarDTempVar_LIN
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

### Extinction
extinction_par1 <- resi_all$BTempVarDTempVar_LIN$mu_par[1]
extinction_par2 <- resi_all$BTempVarDTempVar_LIN$mu_par[2]
extinction <- function(x) { abs(extinction_par1 + extinction_par2 * Tm_time(x)) }
extinction.min <- round(min(extinction(Temperature[,"Age"])),3)
extinction.max <- round(max(extinction(Temperature[,"Age"])),3)

### Net Diversification (新增)
net_diversification <- function(x) {
  nd <- speciation(x) - extinction(x)
  nd[nd < 0] <- 0  # 将负值设为0
  # 或者使用 nd[nd < 0] <- 0  # 如果不想显示0值线
  return(nd)
}

# Plot
plot(age, age, type='n', 
     ylim=c(0, max(speciation.max, extinction.max, net_div.max)),  # y轴从0开始
     xlim=c(-age,0), axes=F, ylab='', xlab='', 
     bty='n', las=1, cex.axis=0.5, cex.lab=0.5)

# 背景色块（地质时期）
rect(-66,0,-max(Temperature$Age),1, col="whitesmoke", border=NA)
rect(-56,0,-66,1, col="white", border=NA)
rect(-56,0,-33.9,1, col="whitesmoke", border=NA)
rect(-33.9,0,-23.03,1, col="white", border=NA)
rect(-23.03,0,-5.33,1, col="whitesmoke", border=NA)
rect(-5.33,0,-2.58,1, col="white", border=NA)
rect(-2.58,0,0,1, col="whitesmoke", border=NA)

# 绘制 speciation、extinction 和 net diversification
lines(-Temperature[,"Age"], speciation(Temperature[,"Age"]), type="l", col=alpha(colors[3], 0.7), lwd=2)
lines(-Temperature[,"Age"], extinction(Temperature[,"Age"]), type="l", col=alpha(colors[4], 0.7), lwd=2, lty="dashed")
lines(-Temperature[,"Age"], net_diversification(Temperature[,"Age"]), type="l", col=alpha(colors[2], 0.7), lwd=2, lty="dashed")
# 地质时期分界线
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -max(Temperature$Age)), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

# 坐标轴
axis(side=1, seq(-round(age, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, seq(0, max(speciation.max, extinction.max, net_div.max), by=0.1),  # 调整 y 轴刻度
     col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)

# 添加地质年代
add_geochrono(0, -0.1)
# 图例（更新，包含 net diversification）
legend("topright", bty="n", 
       c("Speciation rate", "Extinction rate", "Net diversification rate"),
       lwd=2, 
       col=c(alpha(colors[3], 0.7), alpha(colors[4], 0.7), alpha(colors[2], 0.7)), 
       lty=c(1, 2, 2), 
       cex=0.8)

# 标题
mtext("BTempVarDTempVar_LIN  AICweight = 0.165", side = 3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side = 1, line = 1.5, cex=0.8)

# CO2-dependent diversification
# BCO2VarDCST_EXPO
CO2 <- read.table("./data/500kyrCO2_combined_62.75Ma.txt", h=T)
# Environmental data
CO2_data1 <- as.vector(scale(CO2[,2]))
CO2_data <- cbind(CO2[,1], CO2_data1)
colnames(CO2_data) <- c("Age", "CO2")
CO2 <- as.data.frame(CO2_data)

CO2_spline <- sm.spline(CO2[,1], CO2[,2])
CO2_time <- function(x) { predict(CO2_spline, x) }

# Speciation
speciation_par1 <- resi_all$BCO2VarDCST_EXPO$lamb_par[1]
speciation_par2 <- resi_all$BCO2VarDCST_EXPO$lamb_par[2]
speciation <- function(x) { abs(speciation_par1 * exp(speciation_par2 * CO2_time(x))) }
speciation.max <- round(max(speciation(CO2[,1])), 1)

# Extinction
extinction_par1 <- abs(resi_all$BCO2VarDCST_EXPO$mu_par[1])  # 注意这里应该是mu_par不是lamb_par
extinction_par2 <- 0
extinction <- function(x) { extinction_par1 * exp(extinction_par2 * CO2_time(x)) }
extinction.max <- round(max(extinction(CO2[,1])), 1)

# Net Diversification (新增)
net_diversification <- function(x) {
  speciation(x) - extinction(x)
}
net_div.max <- round(max(net_diversification(CO2[,1])), 1)
net_div.min <- round(min(net_diversification(CO2[,1])), 1)

# Plot
plot(age, age, type='n', 
     ylim=c(0, max(speciation.max, extinction.max)),  # 调整y轴范围
     xlim=c(-age, 0), axes=F, ylab='', xlab='', 
     main='', bty='n', las=1, cex.axis=0.8, cex.lab=0.5)

# 背景色块（地质时期）
rect(-66, -0.001, -max(CO2$Age), max(speciation.max, extinction.max, net_div.max), col="whitesmoke", border=NA)
rect(-56, -0.001, -66, max(speciation.max, extinction.max, net_div.max), col="white", border=NA)
rect(-56, -0.001, -33.9, max(speciation.max, extinction.max, net_div.max), col="whitesmoke", border=NA)
rect(-33.9, -0.001, -23.03, max(speciation.max, extinction.max, net_div.max), col="white", border=NA)
rect(-23.03, -0.001, -5.33, max(speciation.max, extinction.max, net_div.max), col="whitesmoke", border=NA)
rect(-5.33, -0.001, -2.58, max(speciation.max, extinction.max, net_div.max), col="white", border=NA)
rect(-2.58, -0.001, 0, max(speciation.max, extinction.max, net_div.max), col="whitesmoke", border=NA)

# 绘制三条曲线
lines(-CO2[,1], speciation(CO2[,1]), type="l", col=alpha(colors[3], 0.7), lwd=2)
lines(-CO2[,1], extinction(CO2[,1]), type="l", alpha(colors[4], 0.7), lwd=2, lty="dashed")
lines(-CO2[,1], net_diversification(CO2[,1]), type="l", col=alpha(colors[2], 0.7), lwd=2, lty="dashed")  # 新增净多样化速率

# 地质时期分界线
abline(v=c(-2.58, -5.33, -23.03, -33.9, -56, -max(CO2$Age)), col="black", lty="dotted", lwd=1)
abline(v=c(-55.52), col="red", lty="dashed", lwd=2)

# 坐标轴
axis(side=1, seq(-round(age+1, -1), 0, by=10), col="black", col.axis="black", las=1, cex.axis=0.8, padj=-1)
axis(side=2, seq(0, max(speciation.max, extinction.max, net_div.max), by=0.1),  # 调整y轴刻度
     col="black", col.axis="black", las=2, cex.axis=0.8, hadj=0.5)

# 添加地质年代
add_geochrono(0, -0.2)  

legend("topright", bty="n", 
       c("Speciation rate", "Extinction rate", "Net diversification rate"),
       lwd=2, 
       col=c(alpha(colors[3], 0.7),alpha(colors[4], 0.7), alpha(colors[2], 0.7)), 
       lty=c(1, 2, 2), 
       cex=0.8)

mtext("BCO2VarDCST_EXPO  AICweight = 0.022", side=3, at=-40, line=-3, cex=0.8)
mtext("Time (Myr)", side=1, line=1.7, cex=0.8)
dev.off()