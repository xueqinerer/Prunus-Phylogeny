setwd("C:/Users/301/Desktop/李属文章数据/BAMM_MO_35,1,3,2/Ancestral_traits/data/Inflorescence_2025.12.1/")


rm(list=ls())
library(ape)
library (diversitree)
library(phytools)
library(geiger)

phy<-read.tree("../../Prunus.mcmctree.dated_no_outgroup.tre")

X <- read.csv("inflorescence_types1.csv", row.names = 1)

trait <- X$classify
names(trait)<- row.names(X)
comparison <- name.check(phy=phy,data=trait)
comparison
comparison <- name.check(phy=phy, data=trait)

if (is.list(comparison)) {
  # 如果有不匹配的名字，删除树里没有的
  phy <- drop.tip(phy, comparison$tree_not_data)
} else {
  # 如果返回 "OK"，说明不用处理
  message("Tree and data match perfectly, no tips dropped.")
}

samplingf<-c(0.22, 0.49, 0.125)
p <- starting.point.musse(phy, k=3)
p

# lambda1     lambda2     lambda3         mu1         mu2 
# 0.045258787 0.045258787 0.045258787 0.000000000 0.000000000 
# mu3         q12         q13         q21         q23 
# 0.000000000 0.009051757 0.009051757 0.009051757 0.009051757 
# q31         q32 
# 0.009051757 0.009051757 



make_musse_models <- function(lik, k) {
  # 全部约束相同
  null.formulas <- c(
    lapply(2:k, function(i) as.formula(paste0("lambda", i, " ~ lambda1"))),
    lapply(2:k, function(i) as.formula(paste0("mu", i, " ~ mu1")))
  )
  lik.null <- do.call(constrain, c(list(lik), null.formulas))
  
  # extinction 相同
  lambda.formulas <- lapply(2:k, function(i) as.formula(paste0("mu", i, " ~ mu1")))
  lik.lambda <- do.call(constrain, c(list(lik), lambda.formulas))
  
  # speciation 相同
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
#   free.mu       10 -383.68 787.35 6.3534    0.04172 *
#   ---
#   Signif. codes:  
#   0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 模型选择原则
# AIC 最小优先
# AIC = 2 × Df - 2 × lnLik
# 越小越好，权衡拟合度和参数数量。
# 这里 free.lambda 的 AIC = 801.58 最小 → 是首选模型。
# 统计显著性（卡方检验）
# 与最简模型比较：
# all.different vs minimal：p < 0.001 → 显著改善拟合
# free.lambda vs minimal：p < 0.001 → 显著改善拟合
# free.mu vs minimal：p = 0.0059 → 也显著，但改善程度不如 free.lambda
# 与 AIC 结合来看，free.lambda 更好。

# 生物学解释
# free.lambda：允许不同状态有不同 speciation (λ)，但 extinction (μ) 相同
# free.mu：允许不同状态有不同 extinction (μ)，但 speciation (λ) 相同
# 数据显示 speciation 的差异比 extinction 更显著 → 也与 AIC 结果一致。

write.csv(AnovaResults,"Musse_modeltest_infloresence_types1.csv")
aicw(setNames(AnovaResults$AIC,rownames(AnovaResults)))
# fit    delta          w
# minimal       789.7050 4.688240 0.06181293
# all.different 788.8807 3.863906 0.09334287
# free.lambda   785.0167 0.000000 0.64434379
# free.mu       787.3516 2.334832 0.20050041
# 计算 AIC 权重
aic_table <- aicw(setNames(AnovaResults$AIC, rownames(AnovaResults)))

# aic_table 的列名是: fit, delta, w
# 我们只取 delta 和 w
delta_w <- aic_table[, c("delta", "w")]

# 合并到原来的 AnovaResults
AnovaResults_full <- cbind(AnovaResults, delta_w)

# 写入 CSV
write.csv(AnovaResults_full, "Musse_modeltest_infloresence_types1_with_delta_w.csv", row.names = TRUE)






#############
#############################
# 1️⃣ 设置 prior
#############################
prior <- make.prior.exponential(1/2)

#############################
# 2️⃣ 确保参数顺序匹配
#############################
start_par <- fit.lambda$par
start_par <- start_par[argnames(models$lambda)]  # 按模型参数顺序排列

#############################
# 3️⃣ 短 MCMC 预估步长 w
#############################
prelim <- mcmc(models$lambda, start_par, nsteps=1000, prior=prior, w=0.1, print.every=100)

# 提取参数列（不包括 generation 列）
params_mat <- prelim[, argnames(models$lambda)]

# 每列计算 5% 和 95% 分位数差作为步长
w <- apply(params_mat, 2, function(x) diff(quantile(x, c(0.05, 0.95))))

# 检查长度
if(length(w) != length(argnames(models$lambda))) stop("步长 w 长度与参数数量不匹配")

#############################
# 4️⃣ 长 MCMC
#############################
nsteps_long <- 5000
mcmc.fit.full <- mcmc(models$lambda, start_par, nsteps=nsteps_long, prior=prior, w=w, print.every=100)

#############################
# 5️⃣ 丢弃 burn-in
#############################
burnin <- 500
mcmc.fit.full <- mcmc.fit.full[(burnin+1):nrow(mcmc.fit.full), ]

#############################
# 6️⃣ 提取 λ、μ、净多样化率
#############################
lambda <- mcmc.fit.full[, grep("lambda", colnames(mcmc.fit.full))]
mu     <- mcmc.fit.full[, grep("mu", colnames(mcmc.fit.full))]
net.div <- lambda - mu
colnames(net.div) <- paste0("lambda(", 1:ncol(lambda), ")")

#############################
# 7️⃣ 输出 CSV 文件
#############################
write.csv(lambda, "lambda_samples_infloresence_types_.csv", row.names = FALSE)
write.csv(mu, "mu_samples_infloresence_types_.csv", row.names = FALSE)
write.csv(net.div, "netdiv_samples_infloresence_types_.csv", row.names = FALSE)
save(lambda, mu, net.div, w, mcmc.fit.full, file = "Musse_results_infloresence_types_.RData")
load("Musse_results_infloresence_types_.RData") 
#############################
# 8️⃣ 绘图

pdf("Musse_inflorescence_types.pdf", width = 12, height = 5)  # 横向更宽

k <- ncol(lambda)
colors <- setNames(c('#99DD99', '#90CBFB', '#FDC187', "#234862")[1:k], 1:k)

# 设置一行两列布局
par(mfrow = c(1, 2))

# 左图：物种形成率 λ 曲线
profiles.plot(lambda, 
              xlab = "Speciation rate (λ)", 
              ylab = "Probability density",
              legend.pos = "topright",
              legend.text = paste0("State ", 1:k),
              col.line = colors,
              lty = 1)

# 右图：净多样化率 r 曲线
profiles.plot(net.div, 
              xlab = "Net diversification rate (r = λ − μ)", 
              ylab = "Probability density",
              legend.pos = "topright",
              col.line = setNames(colors, colnames(net.div)),
              lty = 1)

dev.off()

pdf("Musse_inflorescence_types_axis_separate.pdf", width = 12, height = 5)

# 设置颜色
k <- ncol(lambda)
colors <- setNames(c('#99DD99', '#90CBFB', '#FDC187', "#234862")[1:k], 1:k)

# 自定义图例标签
legend_labels <- c("0 = Solitary", "1 = Corymbose", "2 = Racemose")

# 设置一行两列布局
par(mfrow = c(1, 2))

# -----------------------
# 左图：物种形成率 λ 曲线
# -----------------------
profiles.plot(lambda,
              xlab = "", ylab = "",
              #legend.pos = "topright",
              #legend.text = legend_labels,   # ✅ 使用自定义图例
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

# 添加坐标轴与标签
box(bty = "n")
axis(1, line = 0, lwd = 2, tck = -0.02)
axis(2, line = 0, lwd = 2, tck = -0.02)
mtext("Speciation rate", side = 1, line = 2.5, cex = 1.3)
mtext("Probability density", side = 2, line = 2.5, cex = 1.3)

# -----------------------
# 右图：净多样化率 r 曲线
# -----------------------
profiles.plot(net.div,
              xlab = "", ylab = "",
              #legend.pos = "topright",
              #legend.text = legend_labels,   # ✅ 同样使用自定义图例
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



