library(hisse)
library(ape)

# ========= Basic parameters =========
tree_path <- "/path/to/mcmctree/Prunus.mcmctree.dated_no_outgroup.tre"
out_dir <- "/path/to/mcmctree/Misse"
dir.create(out_dir, showWarnings = FALSE)

# ========= Read tree =========
tree <- read.tree(tree_path)
sampled_species <- length(tree$tip.label)
total_species <- 352   # Total species count in Prunus
f.vec <- 0.2357954545

# ========= Run MiSSE models (1-5 hidden states) =========
for (model.number in 1:5) {
  outfile <- file.path(out_dir, paste0("misse.", model.number, ".Rsave"))

  if (file.exists(outfile)) {
    cat("Model", model.number, "already exists, skipping.\n")
    next
  }

  misse.fit <- NA
  cat("Running MiSSE model", model.number, "...\n")
  if (model.number == 1) {
    misse.fit <- MiSSE(phy = tree, f = f.vec, turnover = c(1), eps = c(0))
  } else if (model.number == 2) {
    misse.fit <- MiSSE(phy = tree, f = f.vec, turnover = c(1,2), eps = c(0,0))
  } else if (model.number == 3) {
    misse.fit <- MiSSE(phy = tree, f = f.vec, turnover = c(1,2,3), eps = c(0,0,0))
  } else if (model.number == 4) {
    misse.fit <- MiSSE(phy = tree, f = f.vec, turnover = c(1,2,3,4), eps = c(0,0,0,0))
  } else if (model.number == 5) {
    misse.fit <- MiSSE(phy = tree, f = f.vec, turnover = c(1,2,3,4,5), eps = c(0,0,0,0,0))
  }

  save(tree, misse.fit, file = outfile)
  cat("Completed: model", model.number, "\n")
}

# ========= Summarize results =========
files <- list.files(path = out_dir, pattern = "misse.*\\.Rsave$", full.names = TRUE)
results <- list()

for (file in files) {
  load(file)
  split_name <- strsplit(basename(file), "\\.")[[1]]
  model.number <- split_name[2]

  if (exists("misse.fit") && !is.null(misse.fit)) {
    df <- data.frame(
      model.number = model.number,
      loglik = misse.fit$loglik,
      AIC = misse.fit$AIC,
      AICc = misse.fit$AICc,
      solution = I(list(misse.fit$solution))
    )
    results[[file]] <- df
  }
}

if (length(results) > 0) {
  results_all <- do.call(rbind, results)
  write.csv(results_all,
            file = file.path(out_dir, "misse.all_final_results.csv"),
            row.names = FALSE)
  saveRDS(results_all,
          file = file.path(out_dir, "misse.all_final_results.rds"))
  cat("All results summarized.\n")
}


# Visualization
library(hisse)
library(ape)
setwd("/path/to/mcmctree/Misse")
load("misse.1.Rsave")

phy <- read.tree("Prunus.mcmctree.dated_no_outgroup.tre")
turnover <- c(1,2)
eps <- c(1,1)
two.rate <- MiSSE(phy, f=0.2357954545, turnover=turnover, eps=eps)

two.rate.recon <- MarginReconMiSSE(phy=phy, f=0.2357954545, hidden.states=2,
pars=two.rate$solution, AIC=two.rate$AIC)
save(phy, two.rate.recon, file = "two.rate.recon.Rsave")

load("misse.vignette.Rsave")
class(two.rate.recon)

pdf("two.rate.recon.pdf", height=12, width=8)
plot.misse.states(two.rate.recon, rate.param="net.div", show.tip.label=TRUE, type="phylogram",
                  fsize=0.8, legend="none")
dev.off()
