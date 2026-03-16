setwd("C:/Users/301/Desktop/李属文章数据/BAMM_MO_35,1,3,2")
library("BAMMtools")
#Analysis of rate shifts
library("ape")
# BiocManager::install("treeio")
library(treeio)
# BiocManager::install("ggtree")
library(ggtree)
# install.packages("ggplot2")
library(ggplot2)
library("scales")
library(phytools)
library("coda")
list.files()

#Assessing MCMC convergence
mcmcout <- read.csv("BAMM_Prunus_Prunus_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.25 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
#check the effective sample sizes 
library("coda")
logLik <- effectiveSize(postburn$logLik) # calculates autocorrelation function
N_shift <- effectiveSize(postburn$N_shift) #effective sample size on N-shifts
ESSample <- cbind.data.frame("Prunus", logLik, N_shift)
#at least 200


tree <- read.tree("Prunus.mcmctree.dated_no_outgroup.tre")
#tree1 <- ladderize(tree)
plotTree(tree, use.edge.length=FALSE)
nodelabels()
#错误于UseMethod("offspring"): 
#"offspring"没有适用于"NULL"目标对象的方法
#这种报错，一般是安装包冲突了，指定安装包和函数就行
# tree_97<- ape::rotate(tree,97)
# plotTree(tree_97)
# 
# tree_170<- ape::rotate(tree_97,170)
# plotTree(tree_170)
# tree_188<- ape::rotate(tree_170,188)
# plotTree(tree_188)

#write.tree(tree_188,"Prunus.treePL.dated_no_outgroup_change_orientation.tre")
edata <- getEventData(tree, eventdata = "BAMM_Prunus_Prunus_event_data.txt", burnin = 0.25) 
#How many rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)

shift_probs <- summary(edata)
# Analyzed 9001 posterior samples

saveRDS(edata, file=paste("Prunus_edata.rds", sep=""))
edata<- readRDS("./Prunus_edata.rds")
TR <- getTipRates(edata, returnNetDiv = FALSE, statistic = "median")
write.csv(TR,"Prunus_BAMM_TipRates.csv")
file <- paste0("Prunus_BAMM_TipRates.csv", sep="")
cat("Tip_label,Rate\n", file=file)
write.table(TR$lambda.avg, sep=",", file=file, col.names = FALSE, quote=FALSE, append = TRUE)
rtt <- getRateThroughTimeMatrix(edata)
Mean.lamda <- apply(rtt$lambda, 2, quantile,c(0.5))
Mean.mu <- apply(rtt$mu, 2, quantile, c(0.5))
mean.netdiv <- Mean.lamda - Mean.mu
Mean.Rate.Matrix <- cbind.data.frame(rtt$times, Mean.lamda, Mean.mu, mean.netdiv)
write.csv(Mean.Rate.Matrix, paste("Prunus_Mean_Rate_Matrix.csv", sep=""), row.names = FALSE, quote=FALSE)
saveRDS(rtt, paste("Prunus_RateThroughTimeMatrix.rds", sep=""))

# 加载已保存的速率矩阵
Mean.Rate.Matrix <- read.csv("Prunus_Mean_Rate_Matrix.csv")

# 计算整个李属的平均物种形成率（用中位数，抗极端值）
#average_lambda <- median(Mean.Rate.Matrix$Mean.lamda)

# 也可以用算术平均值（根据需要选择）
 average_lambda <- mean(Mean.Rate.Matrix$Mean.lamda)

# 输出结果
cat("李属的平均物种形成率（λ）：", round(average_lambda, 5), "\n")
#李属的平均物种形成率（λ）： 0.09796  

# 计算整个李属的平均净多样化速率（推荐用中位数，抗极端值）
#average_netdiv <- median(Mean.Rate.Matrix$mean.netdiv)

# 也可以用算术平均值（根据需求选择）
 average_netdiv <- mean(Mean.Rate.Matrix$mean.netdiv)

# 输出结果
cat("李属的平均净多样化速率：", round(average_netdiv, 5), "\n")
#李属的平均净多样化速率： 0.08483 


# 加载已有的速率矩阵（rtt）
rtt <- readRDS("Prunus_RateThroughTimeMatrix.rds")

# ------------------------------
# 1. 计算每个时间点的中位数速率（减少极端值影响）
# ------------------------------
# 每个时间点的物种形成率（lambda）中位数
lambda_per_time <- apply(rtt$lambda, 2, median)  # 按列（时间点）取中位数
# 每个时间点的灭绝率（mu）中位数
mu_per_time <- apply(rtt$mu, 2, median)
# 每个时间点的净多样化率（lambda - mu）中位数
netdiv_per_time <- lambda_per_time - mu_per_time


# ------------------------------
# 2. 计算整体平均速率（Mean）和标准误（SE）
# ------------------------------
# 平均物种形成率及其标准误
Mean_Speciation <- mean(lambda_per_time)
SE_Speciation <- sd(lambda_per_time) / sqrt(length(lambda_per_time))  # SE = 标准差 / 样本量开方

# 平均灭绝率及其标准误
Mean_Extinction <- mean(mu_per_time)
SE_Extinction <- sd(mu_per_time) / sqrt(length(mu_per_time))

# 平均净多样化率及其标准误
Mean_NetDiv <- mean(netdiv_per_time)
SE_NetDiv <- sd(netdiv_per_time) / sqrt(length(netdiv_per_time))


# ------------------------------
# 3. 输出结果（整理成表格）
# ------------------------------
result <- data.frame(
  Metric = c("Mean Speciation", "SE(Speciation)", 
             "Mean Extinction", "SE(Extinction)", 
             "Mean NetDiv", "SE(NetDiv)"),
  Value = c(Mean_Speciation, SE_Speciation,
            Mean_Extinction, SE_Extinction,
            Mean_NetDiv, SE_NetDiv)
)

# 保留5位小数
result$Value <- round(result$Value, 5)

# 打印结果
print(result)

# 保存为CSV（可选）
write.csv(result, "Prunus_Mean_Rates_with_SE.csv", row.names = FALSE)



#Compute credible set of shift configurations for more information:
#See ?credibleShiftSet and ?getBestShiftConfiguration
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
#[1] 3
summary(css)
#Distinct shift configurations in credible set:  3


#绘制不同f值下shift点的位置
pdf("./credibleShift_.pdf",width= 24,height = 16)
plot.credibleshiftset(css, shiftColor = "red", pal="temperature", lwd=1.5)
#Omitted 0 plotshttp://127.0.0.1:21755/graphics/plot_zoom_png?width=1536&height=824
dev.off()


#绘制最佳shift点
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
#x <-plot.bammdata(edata,lwd=2,pal="temperature") 
pdf("bestShift_.pdf",width= 16,height = 24)
x <- plot.bammdata(edata, lwd=1,labels =TRUE,cex=1,spex="netdiv")
addBAMMshifts(best, col="black",cex=2.5)
addBAMMlegend(x,location = c(90, 92, 60, 90),side=4)
#title(sub="3000w + tact_monophyly + 0.25burnin",cex.sub=3, line=-1)
dev.off()


coll <- c("black", "cyan4", "deeppink",  "blue", "darkorange")
inter.cc <- "gray70"
# aa <- outlier(rtt$lambda-rtt$mu)

ratemax <- max(fivenum(rtt$lambda)[5], fivenum(rtt$mu)[5])
ratemin <- min(fivenum(rtt$lambda)[1], fivenum(rtt$mu)[1])
# plotRateThroughTime(edata, ratetype = "speciation", avgCol="deeppink",intervalCol="deeppink",ylim = c(ratemin, ratemax), axis.labels =FALSE)
# plotRateThroughTime(edata, ratetype = "netdiv", avgCol="black", intervalCol = "gray70",add=TRUE)
# plotRateThroughTime(edata, ratetype = "extinction",  avgCol ="blue",intervalCol = "blue", add=TRUE)
# legend(max(edata$times), ratemax, legend=c( "Speciation", "Extinction", "Net Diver."), border=FALSE, merge=TRUE, seg.len=0.8, bty = "n", col=c("red","black", "blue"), cex=0.8, lty=c(1,1,1), lwd=1)
cc <- read.table("./netdivRate.interpolation/env_data/500kyrCO2_combined_62.75Ma.txt",header = T)
smoothingSpline_co2 <-  smooth.spline(cc$CO2 ~ cc$Age , spar=0.35)
InfTemp <- read.csv("./netdivRate.interpolation/env_data/temperature-0.5ma-FABIEN-approx-63.csv", header = TRUE)
smoothingSpline_tm <-  smooth.spline(InfTemp$Temperature ~ InfTemp$Age, spar=0.35)
# 取log2
smoothingSpline_co2_log2 <- smoothingSpline_co2
smoothingSpline_co2_log2$y <- log2(smoothingSpline_co2$y)

#上下排版的bamm多样化速率图
pdf("best_diversification_shift_speciation_tm_co2_final.pdf", height = 62, width = 42)

# 先控制画图范围
par(fig=c(0,1,0.48,1), new=FALSE, mgp=c(0,2,2), oma=c(8,15,3,6))

# 第一个plot
plot(smoothingSpline_tm, col='#FF0000', axes=F, xlim=c(max(rtt$times),0), type="l", lty=5, lwd=5, xlab="", ylab="")
axis(2,line=0, cex.axis=5,lwd=6,lwd.ticks = 6,tck=-0.008)
mtext("Temperature (°C)", line=8, side=2, cex=6, col="black")
title(main="Best Diversification Shifts (speciation)", font.main=1, cex.main=5)
text(x=-1.2, y=28.5,
     labels=expression("speciation rate (species Ma"^{-1} * "capita"^{-1} * ")"),
     adj=c(1,1), srt=90, cex=3.5)

par(new=TRUE)
plot(1, type="n", xlim=c(max(rtt$times),0), ylim=range(c(smoothingSpline_tm$y, 0)), axes=FALSE, xlab="", ylab="")
legend(x=45, y=28, legend=c("Temperature", "CO2"), col=c("#c04851", "black"), lty=5, lwd=5, cex=5, bty="n")

# 控制画图范围
par(fig=c(0,1,0,0.52), new=TRUE, mgp=c(0,2,2), oma=c(8,15,3,6))

# 第二个plot
par(fig=c(0,1,0,0.52),new=TRUE)
plot(smoothingSpline_co2_log2, col='black', axes=F, xlim=c(max(rtt$times),0), type="l", lty=5, lwd=5, xlab="", ylab="")
axis(2, line=0, cex.axis=5,at=seq(8, 10, by=1),labels=round(seq(8, 10, by=1), 1),lwd=6,lwd.ticks = 6,tck=-0.008)
axis(1, line=0, cex.axis=5, padj=0.8,lwd=6,lwd.ticks = 6,tck=-0.008)
mtext(expression(log[2](CO[2]~ppm)), line=8, side=2, cex=6, col="black")
mtext("Time before present (Myr)",line=9,side=1,cex=6)

par(fig=c(0,1,0,1),new=T)
x <- plot.bammdata(best, lwd=2.5, spex="netdiv",pal = "Spectral",vtheta=5,rbf=0.0000001,breaksmethod = "jenks")
par(fig=c(0,1,0,1),new=T)
addBAMMshifts(best, cex=5,par.reset=FALSE)
addBAMMlegend(x,location = c(65, 68, 40, 60),side=4,cex.axis = 1.5)
dev.off()

pdf("bestShift_addBAMMlegend.pdf",width= 16,height = 24)
x <- plot.bammdata(best, lwd=2.5,labels =TRUE, spex="s",pal = "Spectral",vtheta=5,rbf=0.0000001,breaksmethod = "jenks")
addBAMMshifts(best, cex=5,par.reset=FALSE)
addBAMMlegend(x,location = c(90, 95, 60, 90),side=4, cex=3)
dev.off()


#bamm曲线
# ===========================
# BAMM 曲线绘制（新版配色 + 灰色置信区间）
# ===========================

colors <- c(
  "Net Diversification" = "#758963",  # olive green
  "extinction" = "#234862",           # dark blue
  "speciation" = "#87171c"            # dark red
)

inter.cc <- "gray70"

ratemax <- max(fivenum(rtt$lambda)[5], fivenum(rtt$mu)[5])
ratemin <- min(fivenum(rtt$lambda)[1], fivenum(rtt$mu)[1])

pdf("./bamm_Prunus_final_peise.pdf", height = 12, width = 16)
par(new = TRUE, mgp = c(0, 0.6, 0))

# ---- 1. 物种形成速率 ----
plotRateThroughTime(
  rtt,
  ratetype = "speciation",
  axis = FALSE,
  useMedian = TRUE,
  lwd = 2,
  avgCol = colors["speciation"],
  intervalCol = inter.cc,      # 灰色置信区间
  ylim = c(ratemin, ratemax),
  smooth = TRUE,
  cex.lab = 1,
  xline = 3,
  yline = 2.8,
  cex.axis = 1
)

# ---- 2. 净多样化速率 ----
plotRateThroughTime(
  rtt,
  ratetype = "netdiv",
  axis.labels = FALSE,
  useMedian = TRUE,
  lwd = 2,
  avgCol = colors["Net Diversification"],
  intervalCol = inter.cc,      # 灰色置信区间
  smooth = TRUE,
  cex.lab = 1.3,
  xline = 2.2,
  yline = 3.2,
  cex.axis = 1,
  add = TRUE
)

# ---- 3. 灭绝速率 ----
plotRateThroughTime(
  rtt,
  ratetype = "extinction",
  axis.labels = FALSE,
  useMedian = TRUE,
  lwd = 2,
  avgCol = colors["extinction"],
  intervalCol = inter.cc,      # 灰色置信区间
  smooth = TRUE,
  cex.lab = 1,
  xline = 2,
  yline = 2.8,
  cex.axis = 1,
  add = TRUE
)

# ---- 坐标轴标签 ----
mtext(expression("Diversification rate (Myr"^"-1"*")"), line = 2.5, side = 2, cex = 1.3)
mtext("Time before present (Myr)", line = 2.5, side = 1, cex = 1.3)

# ---- 图例 ----
legend(
  max(rtt$times), ratemax,
  legend = c("Speciation", "Net Diversification", "Extinction"),
  border = FALSE,
  merge = TRUE,
  seg.len = 0.8,
  bty = "n",
  col = c(colors["speciation"], colors["Net Diversification"], colors["extinction"]),
  cex = 0.8,
  lty = c(1, 1, 5),
  lwd = 2
)

dev.off()

# 1. 先把整条时间轴翻转成“距根”
time_from_root <- max(rtt$times) - rtt$times

# 2. 找到 speciation 中位数峰值对应的“距根”时间
lambda_median <- apply(rtt$lambda, 2, median)
peak_pos      <- which.max(lambda_median)
peak_rate     <- lambda_median[peak_pos]
peak_time_root <- time_from_root[peak_pos]   # 这就是你图上看到的 50 Myr

cat("speciation 中位数最大值 =", peak_rate,
    "\n对应距根时间 =", peak_time_root, "Myr\n")
# speciation 中位数最大值 = 0.1320683 
# 对应距根时间 = 55.39561 Myr


## ================================================================
## 0. 前置计算（时间翻转 + 三曲线峰值）
## ================================================================
lambda_median  <- apply(rtt$lambda, 2, median)
mu_median      <- apply(rtt$mu,     2, median)
net_median     <- lambda_median - mu_median
time_from_root <- max(rtt$times) - rtt$times   # 距根时间

pos_lam <- which.max(lambda_median)
pos_mu  <- which.max(mu_median)
pos_net <- which.max(net_median)

x_peak  <- time_from_root[c(pos_lam, pos_net, pos_mu)]   # 横坐标（距根）
y_peak  <- c(lambda_median[pos_lam],
             net_median[pos_net],
             mu_median[pos_mu])

## 颜色与图例顺序一致
col_peak <- c(colors["speciation"],
              colors["Net Diversification"],
              colors["extinction"])

## ================================================================
## 1. 打开 PDF
## ================================================================
pdf("./bamm_Prunus_final_peise_peak.pdf", height = 12, width = 16)
par(new = TRUE, mgp = c(0, 0.6, 0))

## ================================================================
## 2. 画三条曲线（你的原代码，未改动）
## ================================================================
plotRateThroughTime(
  rtt,
  ratetype = "speciation",
  axis = FALSE,
  useMedian = TRUE,
  lwd = 2,
  avgCol = colors["speciation"],
  intervalCol = inter.cc,
  ylim = c(ratemin, ratemax),
  smooth = TRUE,
  cex.lab = 1,
  xline = 3,
  yline = 2.8,
  cex.axis = 1
)

plotRateThroughTime(
  rtt,
  ratetype = "netdiv",
  axis.labels = FALSE,
  useMedian = TRUE,
  lwd = 2,
  avgCol = colors["Net Diversification"],
  intervalCol = inter.cc,
  smooth = TRUE,
  cex.lab = 1.3,
  xline = 2.2,
  yline = 3.2,
  cex.axis = 1,
  add = TRUE
)

plotRateThroughTime(
  rtt,
  ratetype = "extinction",
  axis.labels = FALSE,
  useMedian = TRUE,
  lwd = 2,
  avgCol = colors["extinction"],
  intervalCol = inter.cc,
  smooth = TRUE,
  cex.lab = 1,
  xline = 2,
  yline = 2.8,
  cex.axis = 1,
  add = TRUE
)

## ================================================================
## 3. 坐标轴 & 图例（你的原代码）
## ================================================================
mtext(expression("Diversification rate (Myr"^"-1"*")"), line = 2.5, side = 2, cex = 1.3)
mtext("Time before present (Myr)", line = 2.5, side = 1, cex = 1.3)

legend(
  max(rtt$times), ratemax,
  legend = c("Speciation", "Net Diversification", "Extinction"),
  border = FALSE, merge = TRUE, seg.len = 0.8, bty = "n",
  col = col_peak, cex = 0.8, lty = c(1, 1, 5), lwd = 2
)

## ================================================================
## 4. 标记峰值点（实心圆 + 文字）
## ================================================================
for(i in 1:3){
  points(x_peak[i], y_peak[i], pch = 19, col = col_peak[i], cex = 1.4)
  text(x_peak[i], y_peak[i] + 0.005,
       labels = sprintf("%.3f_%.1f Myr", y_peak[i], x_peak[i]),
       cex = 0.75, col = "black", pos = 3)
}

## ================================================================
## 5. 关闭 PDF
## ================================================================
dev.off()



#add tem + co2
#add tem + co2
pdf("./bamm_Prunus_temp_co2.pdf", height = 12, width = 16)

# 1. 画温度变化曲线（右轴）
plot(smoothingSpline_tm, col = "#c04851", axes = FALSE,
     xlim = c(max(rtt$times), 0), type = "l", lty = 5, lwd = 2,
     xlab = "", ylab = "")
axis(4, line = -1, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.2, 0))
mtext("Temperature (°C)", line = -2.5, side = 4, cex = 1.1)

# 2. 叠加 CO? 曲线（左轴）
par(new = TRUE)
plot(smoothingSpline_co2_log2, col = "black", axes = FALSE,
     xlim = c(max(rtt$times), 0), type = "l", lty = 5, lwd = 2,
     xlab = "", ylab = "")
axis(2, line = 1.5, tck = -0.01, cex.axis = 0.9, mgp = c(1, 0.2, 0))
mtext(expression(log[2](CO[2]~ppm)), line = 1.5, side = 2, cex = 1.1)

# 3. 叠加 speciation 速率曲线
par(new = TRUE)
plotRateThroughTime(rtt, ratetype = "speciation", axis = FALSE, useMedian = TRUE,
                    lwd = 2, avgCol = "red", intervalCol = "red",
                    ylim = c(ratemin, ratemax), smooth = TRUE,
                    cex.lab = 1, xline = 3, yline = 2.8, cex.axis = 1)

# 4. 叠加 netdiv 曲线
plotRateThroughTime(rtt, ratetype = "netdiv", axis.labels = FALSE, useMedian = TRUE,
                    lwd = 2, intervalCol = "gray70", avgCol = "black",
                    smooth = TRUE, cex.lab = 1.3, xline = 2.2, yline = 3.2,
                    cex.axis = 1, add = TRUE)

# 5. 叠加 extinction 曲线
plotRateThroughTime(rtt, ratetype = "extinction", axis.labels = FALSE, useMedian = TRUE,
                    lwd = 2, avgCol = "blue", intervalCol = "blue",
                    smooth = TRUE, cex.lab = 1, xline = 2, yline = 2.8,
                    cex.axis = 1, add = TRUE)

# 6. 添加主 Y 轴标签和 X 轴标签
mtext(expression("Diversification rate (Myr"^"-1"*")"), line = 1.5, side = 2, cex = 1.3)
mtext("Time before present (Myr)", line = 2.5, side = 1, cex = 1.3)

# 7. 添加物种多样化图例
legend("topright", inset = 0.02,
       legend = c("Speciation", "Net Diver.", "Extinction"),
       col = c("red", "black", "blue"), lty = c(1, 1, 5), lwd = 2,
       cex = 1, bty = "n")

# 8. 添加温度与CO2图例
legend("topleft", inset = 0.02,
       legend = c("Temperature", "CO2"),
       col = c("#c04851", "black"), lty = 5, lwd = 2,
       cex = 1, bty = "n")

dev.off()