setwd("/path/to/BAMM")
library("BAMMtools")
# Analysis of rate shifts
library("ape")
library(treeio)
library(ggtree)
library(ggplot2)
library("scales")
library(phytools)
library("coda")
list.files()

# Assess MCMC convergence
mcmcout <- read.csv("BAMM_Prunus_Prunus_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)
burnstart <- floor(0.25 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]
# Check the effective sample sizes
library("coda")
logLik <- effectiveSize(postburn$logLik)
N_shift <- effectiveSize(postburn$N_shift) # effective sample size on N-shifts
ESSample <- cbind.data.frame("Prunus", logLik, N_shift)
# ESS should be at least 200

tree <- read.tree("Prunus.mcmctree.dated_no_outgroup.tre")
plotTree(tree, use.edge.length=FALSE)
nodelabels()
# Note: if "offspring" UseMethod error occurs, it usually indicates package conflicts.
# Fix by specifying the package and function explicitly, e.g., ape::rotate()

edata <- getEventData(tree, eventdata = "BAMM_Prunus_Prunus_event_data.txt", burnin = 0.25)
# How many rate shifts
post_probs <- table(postburn$N_shifts) / nrow(postburn)
names(post_probs)

shift_probs <- summary(edata)
# Analyzed 9001 posterior samples

saveRDS(edata, file=paste("Prunus_edata.rds", sep=""))
edata <- readRDS("./Prunus_edata.rds")
TR <- getTipRates(edata, returnNetDiv = FALSE, statistic = "median")
write.csv(TR, "Prunus_BAMM_TipRates.csv")
file <- paste0("Prunus_BAMM_TipRates.csv", sep="")
cat("Tip_label,Rate\n", file=file)
write.table(TR$lambda.avg, sep=",", file=file, col.names = FALSE, quote=FALSE, append = TRUE)
rtt <- getRateThroughTimeMatrix(edata)
Mean.lamda <- apply(rtt$lambda, 2, quantile, c(0.5))
Mean.mu <- apply(rtt$mu, 2, quantile, c(0.5))
mean.netdiv <- Mean.lamda - Mean.mu
Mean.Rate.Matrix <- cbind.data.frame(rtt$times, Mean.lamda, Mean.mu, mean.netdiv)
write.csv(Mean.Rate.Matrix, paste("Prunus_Mean_Rate_Matrix.csv", sep=""), row.names = FALSE, quote=FALSE)
saveRDS(rtt, paste("Prunus_RateThroughTimeMatrix.rds", sep=""))

# Load the saved rate matrix
Mean.Rate.Matrix <- read.csv("Prunus_Mean_Rate_Matrix.csv")

# Calculate the mean speciation rate of Prunus
average_lambda <- mean(Mean.Rate.Matrix$Mean.lamda)
cat("Mean speciation rate (lambda) of Prunus:", round(average_lambda, 5), "\n")

# Calculate the mean net diversification rate of Prunus
average_netdiv <- mean(Mean.Rate.Matrix$mean.netdiv)
cat("Mean net diversification rate of Prunus:", round(average_netdiv, 5), "\n")


# Load the saved rate matrix (rtt)
rtt <- readRDS("Prunus_RateThroughTimeMatrix.rds")

# 1. Calculate median rate at each time point (reduces influence of outliers)
lambda_per_time <- apply(rtt$lambda, 2, median)  # column-wise median (per time point)
mu_per_time <- apply(rtt$mu, 2, median)
netdiv_per_time <- lambda_per_time - mu_per_time

# 2. Calculate overall mean and standard error (SE)
Mean_Speciation <- mean(lambda_per_time)
SE_Speciation <- sd(lambda_per_time) / sqrt(length(lambda_per_time))
Mean_Extinction <- mean(mu_per_time)
SE_Extinction <- sd(mu_per_time) / sqrt(length(mu_per_time))
Mean_NetDiv <- mean(netdiv_per_time)
SE_NetDiv <- sd(netdiv_per_time) / sqrt(length(netdiv_per_time))

# 3. Output results as a table
result <- data.frame(
  Metric = c("Mean Speciation", "SE(Speciation)",
             "Mean Extinction", "SE(Extinction)",
             "Mean NetDiv", "SE(NetDiv)"),
  Value = c(Mean_Speciation, SE_Speciation,
            Mean_Extinction, SE_Extinction,
            Mean_NetDiv, SE_NetDiv)
)
result$Value <- round(result$Value, 5)
print(result)
write.csv(result, "Prunus_Mean_Rates_with_SE.csv", row.names = FALSE)


# Compute credible set of shift configurations
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95)
css$number.distinct
summary(css)
# Distinct shift configurations in credible set: 3

# Plot shift locations under different f values
pdf("./credibleShift_.pdf", width=24, height=16)
plot.credibleshiftset(css, shiftColor = "red", pal="temperature", lwd=1.5)
dev.off()

# Plot the best shift configuration
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1, threshold=5)
pdf("bestShift_.pdf", width=16, height=24)
x <- plot.bammdata(edata, lwd=1, labels=TRUE, cex=1, spex="netdiv")
addBAMMshifts(best, col="black", cex=2.5)
addBAMMlegend(x, location=c(90, 92, 60, 90), side=4)
dev.off()


coll <- c("black", "cyan4", "deeppink",  "blue", "darkorange")
inter.cc <- "gray70"

ratemax <- max(fivenum(rtt$lambda)[5], fivenum(rtt$mu)[5])
ratemin <- min(fivenum(rtt$lambda)[1], fivenum(rtt$mu)[1])

cc <- read.table("./netdivRate.interpolation/env_data/500kyrCO2_combined_62.75Ma.txt", header=T)
smoothingSpline_co2 <-  smooth.spline(cc$CO2 ~ cc$Age, spar=0.35)
InfTemp <- read.csv("./netdivRate.interpolation/env_data/temperature-0.5ma-FABIEN-approx-63.csv", header=TRUE)
smoothingSpline_tm <-  smooth.spline(InfTemp$Temperature ~ InfTemp$Age, spar=0.35)
# Log2 transform CO2
smoothingSpline_co2_log2 <- smoothingSpline_co2
smoothingSpline_co2_log2$y <- log2(smoothingSpline_co2$y)

# Vertically stacked BAMM diversification rate plot with temperature and CO2
pdf("best_diversification_shift_speciation_tm_co2_final.pdf", height=62, width=42)

# Control plot area for upper panel
par(fig=c(0,1,0.48,1), new=FALSE, mgp=c(0,2,2), oma=c(8,15,3,6))

# Upper panel: temperature curve
plot(smoothingSpline_tm, col='#FF0000', axes=F, xlim=c(max(rtt$times),0), type="l", lty=5, lwd=5, xlab="", ylab="")
axis(2, line=0, cex.axis=5, lwd=6, lwd.ticks=6, tck=-0.008)
mtext("Temperature (\u00b0C)", line=8, side=2, cex=6, col="black")
title(main="Best Diversification Shifts (speciation)", font.main=1, cex.main=5)
text(x=-1.2, y=28.5,
     labels=expression("speciation rate (species Ma"^{-1} * "capita"^{-1} * ")"),
     adj=c(1,1), srt=90, cex=3.5)

par(new=TRUE)
plot(1, type="n", xlim=c(max(rtt$times),0), ylim=range(c(smoothingSpline_tm$y, 0)), axes=FALSE, xlab="", ylab="")
legend(x=45, y=28, legend=c("Temperature", "CO2"), col=c("#c04851", "black"), lty=5, lwd=5, cex=5, bty="n")

# Control plot area for lower panel
par(fig=c(0,1,0,0.52), new=TRUE, mgp=c(0,2,2), oma=c(8,15,3,6))

# Lower panel: CO2 curve
par(fig=c(0,1,0,0.52), new=TRUE)
plot(smoothingSpline_co2_log2, col='black', axes=F, xlim=c(max(rtt$times),0), type="l", lty=5, lwd=5, xlab="", ylab="")
axis(2, line=0, cex.axis=5, at=seq(8, 10, by=1), labels=round(seq(8, 10, by=1), 1), lwd=6, lwd.ticks=6, tck=-0.008)
axis(1, line=0, cex.axis=5, padj=0.8, lwd=6, lwd.ticks=6, tck=-0.008)
mtext(expression(log[2](CO[2]~ppm)), line=8, side=2, cex=6, col="black")
mtext("Time before present (Myr)", line=9, side=1, cex=6)

par(fig=c(0,1,0,1), new=T)
x <- plot.bammdata(best, lwd=2.5, spex="netdiv", pal="Spectral", vtheta=5, rbf=0.0000001, breaksmethod="jenks")
par(fig=c(0,1,0,1), new=T)
addBAMMshifts(best, cex=5, par.reset=FALSE)
addBAMMlegend(x, location=c(65, 68, 40, 60), side=4, cex.axis=1.5)
dev.off()

pdf("bestShift_addBAMMlegend.pdf", width=16, height=24)
x <- plot.bammdata(best, lwd=2.5, labels=TRUE, spex="s", pal="Spectral", vtheta=5, rbf=0.0000001, breaksmethod="jenks")
addBAMMshifts(best, cex=5, par.reset=FALSE)
addBAMMlegend(x, location=c(90, 95, 60, 90), side=4, cex=3)
dev.off()


# BAMM rate-through-time curves (custom color scheme + gray confidence intervals)
colors <- c(
  "Net Diversification" = "#758963",  # olive green
  "extinction" = "#234862",           # dark blue
  "speciation" = "#87171c"            # dark red
)

inter.cc <- "gray70"
ratemax <- max(fivenum(rtt$lambda)[5], fivenum(rtt$mu)[5])
ratemin <- min(fivenum(rtt$lambda)[1], fivenum(rtt$mu)[1])

pdf("./bamm_Prunus_final_peise.pdf", height=12, width=16)
par(new=TRUE, mgp=c(0, 0.6, 0))

# Speciation rate
plotRateThroughTime(
  rtt, ratetype="speciation", axis=FALSE, useMedian=TRUE, lwd=2,
  avgCol=colors["speciation"], intervalCol=inter.cc,
  ylim=c(ratemin, ratemax), smooth=TRUE, cex.lab=1, xline=3, yline=2.8, cex.axis=1
)

# Net diversification rate
plotRateThroughTime(
  rtt, ratetype="netdiv", axis.labels=FALSE, useMedian=TRUE, lwd=2,
  avgCol=colors["Net Diversification"], intervalCol=inter.cc,
  smooth=TRUE, cex.lab=1.3, xline=2.2, yline=3.2, cex.axis=1, add=TRUE
)

# Extinction rate
plotRateThroughTime(
  rtt, ratetype="extinction", axis.labels=FALSE, useMedian=TRUE, lwd=2,
  avgCol=colors["extinction"], intervalCol=inter.cc,
  smooth=TRUE, cex.lab=1, xline=2, yline=2.8, cex.axis=1, add=TRUE
)

mtext(expression("Diversification rate (Myr"^"-1"*")"), line=2.5, side=2, cex=1.3)
mtext("Time before present (Myr)", line=2.5, side=1, cex=1.3)

legend(
  max(rtt$times), ratemax,
  legend=c("Speciation", "Net Diversification", "Extinction"),
  border=FALSE, merge=TRUE, seg.len=0.8, bty="n",
  col=c(colors["speciation"], colors["Net Diversification"], colors["extinction"]),
  cex=0.8, lty=c(1, 1, 5), lwd=2
)
dev.off()

# Find peak speciation rate and corresponding time from root
time_from_root <- max(rtt$times) - rtt$times
lambda_median  <- apply(rtt$lambda, 2, median)
mu_median      <- apply(rtt$mu,     2, median)
net_median     <- lambda_median - mu_median

pos_lam <- which.max(lambda_median)
pos_mu  <- which.max(mu_median)
pos_net <- which.max(net_median)

x_peak  <- time_from_root[c(pos_lam, pos_net, pos_mu)]
y_peak  <- c(lambda_median[pos_lam], net_median[pos_net], mu_median[pos_mu])
col_peak <- c(colors["speciation"], colors["Net Diversification"], colors["extinction"])

cat("Peak median speciation rate =", lambda_median[pos_lam],
    "\nCorresponding time from root =", time_from_root[pos_lam], "Myr\n")

pdf("./bamm_Prunus_final_peise_peak.pdf", height=12, width=16)
par(new=TRUE, mgp=c(0, 0.6, 0))

plotRateThroughTime(rtt, ratetype="speciation", axis=FALSE, useMedian=TRUE, lwd=2,
  avgCol=colors["speciation"], intervalCol=inter.cc, ylim=c(ratemin, ratemax),
  smooth=TRUE, cex.lab=1, xline=3, yline=2.8, cex.axis=1)
plotRateThroughTime(rtt, ratetype="netdiv", axis.labels=FALSE, useMedian=TRUE, lwd=2,
  avgCol=colors["Net Diversification"], intervalCol=inter.cc, smooth=TRUE,
  cex.lab=1.3, xline=2.2, yline=3.2, cex.axis=1, add=TRUE)
plotRateThroughTime(rtt, ratetype="extinction", axis.labels=FALSE, useMedian=TRUE, lwd=2,
  avgCol=colors["extinction"], intervalCol=inter.cc, smooth=TRUE,
  cex.lab=1, xline=2, yline=2.8, cex.axis=1, add=TRUE)

mtext(expression("Diversification rate (Myr"^"-1"*")"), line=2.5, side=2, cex=1.3)
mtext("Time before present (Myr)", line=2.5, side=1, cex=1.3)
legend(max(rtt$times), ratemax, legend=c("Speciation", "Net Diversification", "Extinction"),
  border=FALSE, merge=TRUE, seg.len=0.8, bty="n", col=col_peak, cex=0.8, lty=c(1,1,5), lwd=2)

# Mark peak points (filled circles + labels)
for(i in 1:3){
  points(x_peak[i], y_peak[i], pch=19, col=col_peak[i], cex=1.4)
  text(x_peak[i], y_peak[i] + 0.005,
    labels=sprintf("%.3f_%.1f Myr", y_peak[i], x_peak[i]),
    cex=0.75, col="black", pos=3)
}
dev.off()


# Add temperature and CO2 overlay
pdf("./bamm_Prunus_temp_co2.pdf", height=12, width=16)

# 1. Plot temperature curve (right axis)
plot(smoothingSpline_tm, col="#c04851", axes=FALSE,
     xlim=c(max(rtt$times), 0), type="l", lty=5, lwd=2, xlab="", ylab="")
axis(4, line=-1, tck=-0.01, cex.axis=0.9, mgp=c(1, 0.2, 0))
mtext("Temperature (\u00b0C)", line=-2.5, side=4, cex=1.1)

# 2. Overlay CO2 curve (left axis)
par(new=TRUE)
plot(smoothingSpline_co2_log2, col="black", axes=FALSE,
     xlim=c(max(rtt$times), 0), type="l", lty=5, lwd=2, xlab="", ylab="")
axis(2, line=1.5, tck=-0.01, cex.axis=0.9, mgp=c(1, 0.2, 0))
mtext(expression(log[2](CO[2]~ppm)), line=1.5, side=2, cex=1.1)

# 3. Overlay speciation rate
par(new=TRUE)
plotRateThroughTime(rtt, ratetype="speciation", axis=FALSE, useMedian=TRUE,
  lwd=2, avgCol="red", intervalCol="red", ylim=c(ratemin, ratemax),
  smooth=TRUE, cex.lab=1, xline=3, yline=2.8, cex.axis=1)

# 4. Overlay net diversification
plotRateThroughTime(rtt, ratetype="netdiv", axis.labels=FALSE, useMedian=TRUE,
  lwd=2, intervalCol="gray70", avgCol="black", smooth=TRUE,
  cex.lab=1.3, xline=2.2, yline=3.2, cex.axis=1, add=TRUE)

# 5. Overlay extinction
plotRateThroughTime(rtt, ratetype="extinction", axis.labels=FALSE, useMedian=TRUE,
  lwd=2, avgCol="blue", intervalCol="blue", smooth=TRUE,
  cex.lab=1, xline=2, yline=2.8, cex.axis=1, add=TRUE)

mtext(expression("Diversification rate (Myr"^"-1"*")"), line=1.5, side=2, cex=1.3)
mtext("Time before present (Myr)", line=2.5, side=1, cex=1.3)

legend("topright", inset=0.02,
  legend=c("Speciation", "Net Diver.", "Extinction"),
  col=c("red", "black", "blue"), lty=c(1,1,5), lwd=2, cex=1, bty="n")

legend("topleft", inset=0.02,
  legend=c("Temperature", "CO2"),
  col=c("#c04851", "black"), lty=5, lwd=2, cex=1, bty="n")

dev.off()
