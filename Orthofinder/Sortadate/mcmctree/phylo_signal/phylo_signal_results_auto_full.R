# =========================================================
# 自动化离散性状系统发育信号与祖先状态重建脚本（最终版）
# 功能：
# - 自动选择最佳模型 (ER/SYM/ARD)
# - 离散性状自动处理，多态与二态均支持
# - 二态性状计算 D 值
# - 增加连续性指标：Pagel's Lambda, Blomberg's K
# - Q 网络图与 ACE 图使用完全一致的颜色
# =========================================================

# ========================== 加载包 ==========================
required_pkgs <- c("ape","phytools","geiger","caper","igraph","ggraph","ggplot2","viridis","grid","RColorBrewer")
#to_install <- required_pkgs[!required_pkgs %in% installed.packages()[,"Package"]]
#if(length(to_install)) install.packages(to_install)
lapply(required_pkgs, library, character.only=TRUE)

# ========================== 基本设置 ==========================
setwd("C:/Users/301/Desktop/李属文章数据/BAMM_MO_35,1,3,2/Ancestral_traits/data/csv/more_characteristic_2025.10.30")
treefile <- "../../../Prunus.mcmctree.dated_no_outgroup.tre"
csv_files <- list.files(path = "./", pattern = "\\.csv$", full.names = TRUE)
outdir <- "phylo_signal_results_auto_full"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

set.seed(123)
nsim <- 1000

# 统一调色板（你的指定 8 色）
base_palette <- c('#99DD99', '#90CBFB',"#a5aaa3", '#FDC187', "#015932", "#FC8D62", "#234862", "#c5cc95")

# ========================== 辅助函数 ==========================
wrap_tip <- function(x, width=20){
  sapply(x, function(name) paste(strwrap(name, width=width), collapse="\n"))
}

calc_AICc <- function(logLik, k, n){
  AIC <- -2*logLik + 2*k
  if(n - k - 1 > 0){
    AICc <- AIC + (2 * k * (k + 1)) / (n - k - 1)
  } else {
    AICc <- NA
  }
  return(AICc)
}

# ========================== 读树 ==========================
if(!file.exists(treefile)) stop("Tree file not found: ", treefile)
tree <- read.tree(treefile)

# ========================== 汇总表 ==========================
phylo_summary <- data.frame(
  Trait = character(),
  Type = character(),
  Best_Model = character(),
  LogLik = numeric(),
  AICc = numeric(),
  D = numeric(),
  D_P1 = numeric(),
  D_P0 = numeric(),
  Lambda = numeric(),
  Lambda_P = numeric(),
  K = numeric(),
  K_P = numeric(),
  Valid = logical(),
  Reason = character(),
  stringsAsFactors = FALSE
)

# ========================== 主循环 ==========================
for(file in csv_files){
  trait_name <- gsub(".csv","",basename(file))
  cat("Processing:", trait_name, "\n")
  
  # ---------------- 读取 CSV ----------------
  data <- tryCatch(read.csv(file, row.names=1, check.names = FALSE, stringsAsFactors = FALSE),
                   error=function(e) NULL)
  if(is.null(data)){
    phylo_summary <- rbind(phylo_summary, data.frame(
      Trait = trait_name, Type="Discrete", Best_Model=NA, LogLik=NA, AICc=NA,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="Failed to read CSV"
    ))
    next
  }
  
  # ---------------- 数据清理 ----------------
  data_clean <- data[!is.na(data[,1]) & !(data[,1] %in% c("?", "NA", "")), , drop=FALSE]
  data_clean[,1] <- trimws(as.character(data_clean[,1]))
  
  if(nrow(data_clean) < 10){
    phylo_summary <- rbind(phylo_summary, data.frame(
      Trait = trait_name, Type="Discrete", Best_Model=NA, LogLik=NA, AICc=NA,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="Too few species (<10)"
    ))
    next
  }
  
  # ---------------- 保留树上存在的物种 ----------------
  tip_keep <- intersect(tree$tip.label, rownames(data_clean))
  if(length(tip_keep) < 10){
    phylo_summary <- rbind(phylo_summary, data.frame(
      Trait = trait_name, Type="Discrete", Best_Model=NA, LogLik=NA, AICc=NA,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="Too few overlapping tips"
    ))
    next
  }
  
  tree_pruned <- keep.tip(tree, tip_keep)
  data_clean <- data_clean[tree_pruned$tip.label, , drop=FALSE]
  tree_pruned$tip.label <- rownames(data_clean)
  
  # ---------------- trait factor ----------------
  trait_factor <- as.factor(as.character(data_clean[,1]))
  names(trait_factor) <- rownames(data_clean)
  states <- sort(unique(trait_factor))
  m_states <- length(states)
  
  # ---------------- 颜色 ----------------
  if (m_states <= length(base_palette)) {
    state_colors <- sample(base_palette, m_states)
  } else {
    state_colors <- rep(base_palette, length.out=m_states)
  }
  names(state_colors) <- states
  tip_colors <- state_colors[as.character(trait_factor)]
  
  # ---------------- 模型拟合 ----------------
  models_try <- if(m_states == 2) c("ER","ARD") else c("ER","SYM","ARD")
  fit_results <- list()
  n_tips <- length(tree_pruned$tip.label)
  
  for(mod in models_try){
    f <- try(fitMk(tree_pruned, x=trait_factor, model=mod, method="ML"), silent=TRUE)
    if(inherits(f, "try-error")) next
    k_par <- if(mod=="ER") 1 else if(mod=="SYM") m_states*(m_states-1)/2 else m_states*(m_states-1)
    logL <- f$logLik
    AICc_val <- calc_AICc(logLik=logL, k=k_par, n=n_tips)
    fit_results[[mod]] <- list(f=f, logLik=logL, AICc=AICc_val)
  }
  
  if(length(fit_results)==0){
    phylo_summary <- rbind(phylo_summary, data.frame(
      Trait = trait_name, Type="Discrete", Best_Model=NA, LogLik=NA, AICc=NA,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="All model fits failed"
    ))
    next
  }
  
  best_model <- names(which.min(sapply(fit_results, `[[`, "AICc")))
  best_fit <- fit_results[[best_model]]$f
  best_logL <- fit_results[[best_model]]$logLik
  best_AICc <- fit_results[[best_model]]$AICc
  
  # ---------------- ACE 绘图 ----------------
  ace_res <- tryCatch(ace(trait_factor, tree_pruned, type="discrete", model=best_model, method="ML"), error=function(e) NULL)
  pdf(file.path(outdir, paste0(trait_name, "_ACE_", best_model, ".pdf")), width=10, height=14)
  try({
    plot(tree_pruned, show.tip.label=FALSE, no.margin=TRUE, direction="rightwards")
    title(main=paste("Ancestral state reconstruction:", trait_name, "(", best_model, ")"))
    if(!is.null(ace_res)){
      nodelabels(pie=ace_res$lik.anc, piecol=state_colors[colnames(ace_res$lik.anc)], cex=0.35, frame="none")
    }
    tiplabels(pch=21, bg=tip_colors, cex=1.2)
    legend("topright", legend=names(state_colors), fill=state_colors, bty="n", title="States")
  }, silent=TRUE)
  dev.off()
  
  # ---------------- make.simmap + Q 网络 ----------------
  sim_ok <- TRUE
  trees_sim <- tryCatch(make.simmap(tree_pruned, trait_factor, model=best_model, nsim=nsim, pi="estimated"),
                        error=function(e){ sim_ok <<- FALSE; NULL })
  if(!sim_ok || is.null(trees_sim)){
    phylo_summary <- rbind(phylo_summary, data.frame(
      Trait=trait_name, Type="Discrete", Best_Model=best_model, LogLik=best_logL, AICc=best_AICc,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="make.simmap failed"
    ))
    next
  }
  
  sim_summary <- summary(trees_sim, plot=FALSE)
  Q_list <- lapply(trees_sim, function(x) x$Q)
  Q_mean <- apply(simplify2array(Q_list), 1:2, mean)
  write.csv(Q_mean, file.path(outdir, paste0(trait_name,"_Q_matrix_mean.csv")), row.names=TRUE)
  
  transition_df <- as.data.frame(as.table(sim_summary$count))
  transition_df <- transition_df[transition_df$Freq>0,]
  if(nrow(transition_df)>0){
    colnames(transition_df) <- c("From","To","Average_Count")
    write.csv(transition_df, file.path(outdir, paste0(trait_name,"_transitions.csv")), row.names=FALSE)
  }
  
  # ---- Q 网络图颜色与 ACE 一致 ----
  edges <- which(Q_mean > 0 & row(Q_mean) != col(Q_mean), arr.ind=TRUE)
  if(nrow(edges) > 0){
    df_edges <- data.frame(
      from = rownames(Q_mean)[edges[,1]],
      to   = colnames(Q_mean)[edges[,2]],
      rate = Q_mean[edges]
    )
    g <- graph_from_data_frame(df_edges, directed=TRUE, vertices=data.frame(name=states))
    
    p_net <- ggraph(g, layout='circle') +
      geom_edge_fan(aes(width=rate, alpha=rate, color=rate),
                    arrow=arrow(type="closed", length=unit(3,'mm')),
                    end_cap=circle(6,'mm')) +
      geom_node_point(aes(color=name), size=16) +
      geom_node_text(aes(label=name), size=5, fontface="bold") +
      scale_color_manual(values=state_colors) +
      scale_edge_width(range=c(0.3,3)) +
      scale_edge_alpha(range=c(0.3,1)) +
      scale_edge_color_gradient(low='#FDC160', high="#FC8D62") +
      theme_void() +
      labs(title=paste("Transition rates among states:", trait_name))
    
    ggsave(file.path(outdir, paste0(trait_name,"_Q_network.pdf")), plot=p_net, width=10, height=8)
  }
  
  # ---------------- 二态性状 D 值 ----------------
  D_val <- D_p1 <- D_p0 <- NA
  if(m_states==2){
    comp_df <- data.frame(Species=names(trait_factor), Trait=trait_factor)
    comp_data <- tryCatch(comparative.data(tree_pruned, comp_df, names.col="Species", vcv=TRUE, warn.dropped=FALSE),
                          error=function(e) NULL)
    if(!is.null(comp_data)){
      D_res <- tryCatch(phylo.d(comp_data, binvar="Trait"), error=function(e) NULL)
      if(!is.null(D_res)){
        D_val <- D_res$DEst
        D_p1 <- D_res$Pval1
        D_p0 <- D_res$Pval0
      }
    }
  }
  
  # ---------------- Lambda & K ----------------
  lambda_val <- lambda_p <- K_val <- K_p <- NA
  numeric_trait <- as.numeric(trait_factor)
  names(numeric_trait) <- names(trait_factor)
  
  # Lambda
  lambda_res <- tryCatch({
    phylosig(tree_pruned, numeric_trait, method="lambda", test=TRUE)
  }, error=function(e) NULL)
  
  if(!is.null(lambda_res)){
    lambda_val <- lambda_res$lambda
    lambda_p   <- lambda_res$P
  }
  
  # Blomberg's K
  K_res <- tryCatch({
    phylosig(tree_pruned, numeric_trait, method="K", test=TRUE)
  }, error=function(e) NULL)
  
  if(!is.null(K_res)){
    K_val <- K_res$K
    K_p   <- K_res$P
  }
  
  # ---------------- 汇总 ----------------
  phylo_summary <- rbind(phylo_summary, data.frame(
    Trait=trait_name, Type="Discrete", Best_Model=best_model, LogLik=best_logL, AICc=best_AICc,
    D=D_val, D_P1=D_p1, D_P0=D_p0,
    Lambda=lambda_val, Lambda_P=lambda_p,
    K=K_val, K_P=K_p,
    Valid=TRUE, Reason="Success"
  ))
}

# ========================== 保存汇总表 ==========================
write.csv(phylo_summary, file.path(outdir,"phylo_signal_summary_auto_full.csv"), row.names=FALSE)
cat("✅ All done. Results saved in:", normalizePath(outdir), "\n")
