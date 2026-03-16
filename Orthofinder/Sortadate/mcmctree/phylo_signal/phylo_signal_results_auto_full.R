# =========================================================
# Automated phylogenetic signal and ancestral state reconstruction script
# Features:
# - Automatic best-model selection (ER/SYM/ARD)
# - Handles discrete traits with any number of states
# - Calculates D statistic for binary traits
# - Calculates continuous phylogenetic signal metrics: Pagel's Lambda, Blomberg's K
# - Q-network and ACE plots use consistent color schemes
# =========================================================

# ========================== Load packages ==========================
required_pkgs <- c("ape","phytools","geiger","caper","igraph","ggraph","ggplot2","viridis","grid","RColorBrewer")
#to_install <- required_pkgs[!required_pkgs %in% installed.packages()[,"Package"]]
#if(length(to_install)) install.packages(to_install)
lapply(required_pkgs, library, character.only=TRUE)

# ========================== Basic settings ==========================
setwd("/path/to/ancestral_traits/data/csv")
treefile <- "../../../Prunus.mcmctree.dated_no_outgroup.tre"
csv_files <- list.files(path = "./", pattern = "\\.csv$", full.names = TRUE)
outdir <- "phylo_signal_results_auto_full"
if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

set.seed(123)
nsim <- 1000

# Unified color palette (8 colors)
base_palette <- c('#99DD99', '#90CBFB',"#a5aaa3", '#FDC187', "#015932", "#FC8D62", "#234862", "#c5cc95")

# ========================== Helper functions ==========================
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

# ========================== Read tree ==========================
if(!file.exists(treefile)) stop("Tree file not found: ", treefile)
tree <- read.tree(treefile)

# ========================== Summary table ==========================
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

# ========================== Main loop ==========================
for(file in csv_files){
  trait_name <- gsub(".csv","",basename(file))
  cat("Processing:", trait_name, "\n")

  # ---------------- Read CSV ----------------
  data <- tryCatch(read.csv(file, row.names=1, check.names=FALSE, stringsAsFactors=FALSE),
                   error=function(e) NULL)
  if(is.null(data)){
    phylo_summary <- rbind(phylo_summary, data.frame(
      Trait=trait_name, Type="Discrete", Best_Model=NA, LogLik=NA, AICc=NA,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="Failed to read CSV"
    ))
    next
  }

  # ---------------- Data cleaning ----------------
  data_clean <- data[!is.na(data[,1]) & !(data[,1] %in% c("?", "NA", "")), , drop=FALSE]
  data_clean[,1] <- trimws(as.character(data_clean[,1]))

  if(nrow(data_clean) < 10){
    phylo_summary <- rbind(phylo_summary, data.frame(
      Trait=trait_name, Type="Discrete", Best_Model=NA, LogLik=NA, AICc=NA,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="Too few species (<10)"
    ))
    next
  }

  # ---------------- Keep only species in tree ----------------
  tip_keep <- intersect(tree$tip.label, rownames(data_clean))
  if(length(tip_keep) < 10){
    phylo_summary <- rbind(phylo_summary, data.frame(
      Trait=trait_name, Type="Discrete", Best_Model=NA, LogLik=NA, AICc=NA,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="Too few overlapping tips"
    ))
    next
  }

  tree_pruned <- keep.tip(tree, tip_keep)
  data_clean <- data_clean[tree_pruned$tip.label, , drop=FALSE]
  tree_pruned$tip.label <- rownames(data_clean)

  # ---------------- Convert to factor ----------------
  trait_factor <- as.factor(as.character(data_clean[,1]))
  names(trait_factor) <- rownames(data_clean)
  states <- sort(unique(trait_factor))
  m_states <- length(states)

  # ---------------- Assign colors ----------------
  if (m_states <= length(base_palette)) {
    state_colors <- sample(base_palette, m_states)
  } else {
    state_colors <- rep(base_palette, length.out=m_states)
  }
  names(state_colors) <- states
  tip_colors <- state_colors[as.character(trait_factor)]

  # ---------------- Model fitting ----------------
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
      Trait=trait_name, Type="Discrete", Best_Model=NA, LogLik=NA, AICc=NA,
      D=NA, D_P1=NA, D_P0=NA, Lambda=NA, Lambda_P=NA, K=NA, K_P=NA,
      Valid=FALSE, Reason="All model fits failed"
    ))
    next
  }

  best_model <- names(which.min(sapply(fit_results, `[[`, "AICc")))
  best_fit <- fit_results[[best_model]]$f
  best_logL <- fit_results[[best_model]]$logLik
  best_AICc <- fit_results[[best_model]]$AICc

  # ---------------- ACE plot ----------------
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

  # ---------------- make.simmap + Q network ----------------
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

  # ---- Q network (colors consistent with ACE plot) ----
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

  # ---------------- Binary trait D statistic ----------------
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

  # ---------------- Pagel's Lambda and Blomberg's K ----------------
  lambda_val <- lambda_p <- K_val <- K_p <- NA
  numeric_trait <- as.numeric(trait_factor)
  names(numeric_trait) <- names(trait_factor)

  # Pagel's Lambda
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

  # ---------------- Append to summary ----------------
  phylo_summary <- rbind(phylo_summary, data.frame(
    Trait=trait_name, Type="Discrete", Best_Model=best_model, LogLik=best_logL, AICc=best_AICc,
    D=D_val, D_P1=D_p1, D_P0=D_p0,
    Lambda=lambda_val, Lambda_P=lambda_p,
    K=K_val, K_P=K_p,
    Valid=TRUE, Reason="Success"
  ))
}

# ========================== Save summary table ==========================
write.csv(phylo_summary, file.path(outdir,"phylo_signal_summary_auto_full.csv"), row.names=FALSE)
cat("All done. Results saved in:", normalizePath(outdir), "\n")
