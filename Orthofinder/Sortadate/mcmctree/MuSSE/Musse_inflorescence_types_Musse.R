setwd("/path/to/MuSSE")

rm(list=ls())
library(ape)
library(diversitree)
library(phytools)
library(geiger)

phy <- read.tree("../../Prunus.mcmctree.dated_no_outgroup.tre")

X <- read.csv("inflorescence_types1.csv", row.names = 1)

trait <- X$classify
names(trait) <- row.names(X)
comparison <- name.check(phy=phy, data=trait)
comparison
comparison <- name.check(phy=phy, data=trait)

if (is.list(comparison)) {
  # If there are unmatched names, drop tips not in data
  phy <- drop.tip(phy, comparison$tree_not_data)
} else {
  # If "OK" is returned, tree and data match perfectly
  message("Tree and data match perfectly, no tips dropped.")
}

samplingf <- c(0.22, 0.49, 0.125)
p <- starting.point.musse(phy, k=3)
p

# Starting parameter values:
# lambda1     lambda2     lambda3         mu1         mu2
# 0.045258787 0.045258787 0.045258787 0.000000000 0.000000000
# mu3         q12         q13         q21         q23
# 0.000000000 0.009051757 0.009051757 0.009051757 0.009051757
# q31         q32
# 0.009051757 0.009051757


make_musse_models <- function(lik, k) {
  # Null model: all rates constrained equal
  null.formulas <- c(
    lapply(2:k, function(i) as.formula(paste0("lambda", i, " ~ lambda1"))),
    lapply(2:k, function(i) as.formula(paste0("mu", i, " ~ mu1")))
  )
  lik.null <- do.call(constrain, c(list(lik), null.formulas))

  # Free lambda: allow different speciation rates, constrain extinction
  lambda.formulas <- lapply(2:k, function(i) as.formula(paste0("mu", i, " ~ mu1")))
  lik.lambda <- do.call(constrain, c(list(lik), lambda.formulas))

  # Free mu: allow different extinction rates, constrain speciation
  mu.formulas <- lapply(2:k, function(i) as.formula(paste0("lambda", i, " ~ lambda1")))
  lik.mu <- do.call(constrain, c(list(lik), mu.formulas))

  list(null = lik.null, lambda = lik.lambda, mu = lik.mu)
}


trait <- trait + 1
table(trait)

k <- 3
lik.musse <- make.musse(phy, trait, k, sampling.f = samplingf)

models <- make_musse_models(lik.musse, k)

fit.null   <- find.mle(models$null,   p[argnames(models$null)])
fit.full   <- find.mle(lik.musse,     p[argnames(lik.musse)])
fit.lambda <- find.mle(models$lambda, p[argnames(models$lambda)])
fit.mu     <- find.mle(models$mu,     p[argnames(models$mu)])


AnovaResults <- anova(fit.null,
                      all.different = fit.full,
                      free.lambda   = fit.lambda,
                      free.mu       = fit.mu)
AnovaResults

# Df   lnLik    AIC  ChiSq Pr(>|Chi|)
# minimal        8 -386.85 789.70
# all.different 12 -382.44 788.88 8.8243    0.06564 .
# free.lambda   10 -382.51 785.02 8.6882    0.01298 *
# free.mu       10 -383.68 787.35 6.3534    0.04172 *
#
# Model selection:
# AIC minimum indicates best-fit model (lower is better).
# free.lambda has lowest AIC -> preferred model
# free.lambda: allows different speciation (lambda) per state, same extinction (mu)
# free.mu: allows different extinction (mu) per state, same speciation (lambda)
# Speciation differences are more significant than extinction differences

write.csv(AnovaResults, "Musse_modeltest_infloresence_types1.csv")
aicw(setNames(AnovaResults$AIC, rownames(AnovaResults)))

# Calculate AIC weights
aic_table <- aicw(setNames(AnovaResults$AIC, rownames(AnovaResults)))

# Extract delta and w columns
delta_w <- aic_table[, c("delta", "w")]

# Merge with AnovaResults
AnovaResults_full <- cbind(AnovaResults, delta_w)

# Write to CSV
write.csv(AnovaResults_full, "Musse_modeltest_infloresence_types1_with_delta_w.csv", row.names = TRUE)


#############
#############################
# 1. Set prior
#############################
prior <- make.prior.exponential(1/2)

#############################
# 2. Ensure parameter order matches model
#############################
start_par <- fit.lambda$par
start_par <- start_par[argnames(models$lambda)]

#############################
# 3. Short MCMC to estimate step widths
#############################
prelim <- mcmc(models$lambda, start_par, nsteps=1000, prior=prior, w=0.1, print.every=100)

# Extract parameter columns (exclude generation column)
params_mat <- prelim[, argnames(models$lambda)]

# Compute step width as 5%-95% quantile range per parameter
w <- apply(params_mat, 2, function(x) diff(quantile(x, c(0.05, 0.95))))

# Check length
if(length(w) != length(argnames(models$lambda))) stop("Step width length does not match parameter count")

#############################
# 4. Long MCMC
#############################
nsteps_long <- 5000
mcmc.fit.full <- mcmc(models$lambda, start_par, nsteps=nsteps_long, prior=prior, w=w, print.every=100)

#############################
# 5. Discard burn-in
#############################
burnin <- 500
mcmc.fit.full <- mcmc.fit.full[(burnin+1):nrow(mcmc.fit.full), ]

#############################
# 6. Extract lambda, mu, net diversification
#############################
lambda <- mcmc.fit.full[, grep("lambda", colnames(mcmc.fit.full))]
mu     <- mcmc.fit.full[, grep("mu", colnames(mcmc.fit.full))]
net.div <- lambda - mu
colnames(net.div) <- paste0("lambda(", 1:ncol(lambda), ")")

#############################
# 7. Save to CSV
#############################
write.csv(lambda, "lambda_samples_infloresence_types_.csv", row.names = FALSE)
write.csv(mu, "mu_samples_infloresence_types_.csv", row.names = FALSE)
write.csv(net.div, "netdiv_samples_infloresence_types_.csv", row.names = FALSE)
save(lambda, mu, net.div, w, mcmc.fit.full, file = "Musse_results_infloresence_types_.RData")
load("Musse_results_infloresence_types_.RData")

#############################
# 8. Plot
pdf("Musse_inflorescence_types.pdf", width=12, height=5)

k <- ncol(lambda)
colors <- setNames(c('#99DD99', '#90CBFB', '#FDC187', "#234862")[1:k], 1:k)

# Two-panel layout
par(mfrow = c(1, 2))

# Left panel: speciation rate (lambda) density curves
profiles.plot(lambda,
              xlab = "Speciation rate (\u03bb)",
              ylab = "Probability density",
              legend.pos = "topright",
              legend.text = paste0("State ", 1:k),
              col.line = colors,
              lty = 1)

# Right panel: net diversification rate density curves
profiles.plot(net.div,
              xlab = "Net diversification rate (r = \u03bb \u2212 \u03bc)",
              ylab = "Probability density",
              legend.pos = "topright",
              col.line = setNames(colors, colnames(net.div)),
              lty = 1)

dev.off()

pdf("Musse_inflorescence_types_axis_separate.pdf", width=12, height=5)

k <- ncol(lambda)
colors <- setNames(c('#99DD99', '#90CBFB', '#FDC187', "#234862")[1:k], 1:k)

# Custom legend labels
legend_labels <- c("0 = Solitary", "1 = Corymbose", "2 = Racemose")

par(mfrow = c(1, 2))

# Left panel: speciation rate
profiles.plot(lambda,
              xlab = "", ylab = "",
              col.line = colors,
              lty = 1,
              axes = FALSE, frame.plot = FALSE)
legend("topright",
       legend = c("0 = Solitary", "1 = Corymbose", "2 = Racemose"),
       fill = colors,
       border = "black",
       bty = "n",
       cex = 1.1,
       pt.cex = 1.5)

box(bty = "n")
axis(1, line = 0, lwd = 2, tck = -0.02)
axis(2, line = 0, lwd = 2, tck = -0.02)
mtext("Speciation rate", side = 1, line = 2.5, cex = 1.3)
mtext("Probability density", side = 2, line = 2.5, cex = 1.3)

# Right panel: net diversification rate
profiles.plot(net.div,
              xlab = "", ylab = "",
              col.line = setNames(colors, colnames(net.div)),
              lty = 1,
              axes = FALSE, frame.plot = FALSE)
legend("topright",
       legend = c("0 = Solitary", "1 = Corymbose", "2 = Racemose"),
       fill = colors,
       border = "black",
       bty = "n",
       cex = 1.1,
       pt.cex = 1.5)
box(bty = "n")
axis(1, line = 0, lwd = 2, tck = -0.02)
axis(2, line = 0, lwd = 2, tck = -0.02)
mtext("Net diversification rate", side = 1, line = 2.5, cex = 1.3)
mtext("Probability density", side = 2, line = 2.5, cex = 1.3)
dev.off()
