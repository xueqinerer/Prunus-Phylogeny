# =========================================================
# HiSSE automated pipeline: per-trait CSV + model status + RData + summary CSV + log
# =========================================================

# ==========================
# 1. Load packages
# ==========================
library(hisse)
library(geiger)
library(dplyr)
library(ape)

# ==========================
# 2. Read tree
# ==========================
tree <- read.tree("../Prunus.mcmctree.dated_no_outgroup.tre")

# ==========================
# 3. Get CSV file list
# ==========================
csv_files <- list.files("traits", pattern="\\.csv$", full.names = TRUE)
if(length(csv_files) == 0){
  stop("No CSV files found in 'traits/' folder!")
}

# ==========================
# 4. Model names
# ==========================
model_names <- c("BiSSE_full", "BiSSE_null", "HiSSE_full", "HiSSE_cid2_null")

# ==========================
# 5. AIC weight calculation function
# ==========================
calc_aic_weights <- function(aic_vec){
  delta_aic <- aic_vec - min(aic_vec)
  weights <- exp(-0.5 * delta_aic) / sum(exp(-0.5 * delta_aic))
  return(list(delta_aic = delta_aic, weights = weights))
}

# ==========================
# 6. HiSSE model running function
# ==========================
RunModel <- function(Trait, Trait_name, tree){
  species_with_data <- names(Trait)

  # Drop tips without trait data
  phy <- drop.tip(tree, setdiff(tree$tip.label, species_with_data))

  if(length(phy$tip.label) == 0){
    cat("[WARN] No matching species for trait:", Trait_name, "\n")
    return(NULL)
  }

  # treedata matching
  tree_trait <- tryCatch({
    treedata(phy, Trait, sort=TRUE, warnings=FALSE)
  }, error=function(e){
    cat("[ERROR] treedata failed for trait:", Trait_name, "\n")
    return(NULL)
  })

  if(is.null(tree_trait)) return(NULL)

  phy <- tree_trait$phy
  sim.dat <- cbind(Species=names(Trait), Trait)
  sampling.f <- c(0.5, 0.5)  # Sampling fraction - adjust as needed

  # BiSSE transition matrix
  trans.rates.bisse <- TransMatMaker.old(hidden.states=FALSE)
  trans.rates.bisse <- ParEqual(trans.rates.bisse, c(1,2))

  # HiSSE transition matrix
  trans.rates.hisse <- TransMatMaker.old(hidden.states=TRUE)
  trans.rates.hisse <- ParDrop(trans.rates.hisse, c(3,5,8,10))
  trans.rates.hisse <- ParEqual(trans.rates.hisse, c(1,2,1,3,1,4,1,5,1,6,1,7,1,8))

  trait_model_results <- list()

  for(model.number in 1:4){
    cat("  Running model", model.number, "(", model_names[model.number], ") for trait:", Trait_name, "\n")
    hisse.fit <- NA
    try({
      if(model.number==1){
        hisse.fit <- hisse.old(phy, sim.dat, f=sampling.f, hidden.states=FALSE,
                               turnover.anc = c(1,2,0,0), eps.anc = c(1,1,0,0),
                               trans.rate = trans.rates.bisse, output.type="raw")
      }
      if(model.number==2){
        hisse.fit <- hisse.old(phy, sim.dat, f=sampling.f, hidden.states=FALSE,
                               turnover.anc=c(1,1,0,0), eps.anc=c(1,1,0,0),
                               trans.rate=trans.rates.bisse, output.type="raw")
      }
      if(model.number==3){
        hisse.fit <- hisse.old(phy, sim.dat, f=sampling.f, hidden.states=TRUE,
                               turnover.anc=c(1,2,3,4), eps.anc=c(1,1,1,1),
                               trans.rate=trans.rates.hisse, output.type="raw")
      }
      if(model.number==4){
        hisse.fit <- hisse.old(phy, sim.dat, f=sampling.f, hidden.states=TRUE,
                               turnover.anc=c(1,1,2,2), eps.anc=c(1,1,1,1),
                               trans.rate=trans.rates.hisse, output.type="raw")
      }
    }, silent=TRUE)

    if(!is.null(hisse.fit) && !all(is.na(hisse.fit))){
      # Save RData
      save(phy, hisse.fit, Trait_name, file=paste0(Trait_name, "_trait_model", model.number, ".RData"))

      # Extract parameters
      params <- hisse.fit$solution
      np <- length(params)

      # Build result dataframe
      result_df <- data.frame(
        Trait_name   = Trait_name,
        model.number = model.number,
        model.name   = model_names[model.number],
        hidden.states = hisse.fit$hidden.states,
        condition.on.survival = hisse.fit$condition.on.survival,
        loglik       = hisse.fit$loglik,
        AIC          = hisse.fit$AIC,
        AICc         = hisse.fit$AICc,
        np           = np,
        solution = I(list(hisse.fit$solution)),
        stringsAsFactors = FALSE
      )
      trait_model_results[[length(trait_model_results)+1]] <- result_df
    } else {
      cat("[WARN] Model", model.number, "failed for trait:", Trait_name, "\n")
    }
  }

  if(length(trait_model_results) > 0){
    df <- do.call(rbind, trait_model_results)

    # Calculate delta AIC and AIC weights
    aic_info <- calc_aic_weights(df$AIC)
    df$deltaAIC <- round(aic_info$delta_aic, 3)
    df$AIC_weight <- round(aic_info$weights, 3)

    return(df)
  } else {
    cat("[WARN] No successful models for trait:", Trait_name, "\n")
    return(NULL)
  }
}

# ==========================
# 7. Loop over all CSV files
# ==========================
all_trait_results <- list()

for(csv_file in csv_files){
  cat("Processing CSV:", csv_file, "\n")

  traits <- read.csv(csv_file, stringsAsFactors = FALSE)
  traits[traits == "?"] <- NA
  trait_names <- colnames(traits)[-1]

  # Build trait_list
  trait_list <- list()
  for(i in 2:ncol(traits)){
    trait <- as.numeric(traits[,i])
    names(trait) <- traits[,1]
    trait <- trait[!is.na(trait)]
    trait_list[[i-1]] <- trait
  }

  csv_trait_results <- list()
  for(i in 1:length(trait_list)){
    trait_binary <- ifelse(trait_list[[i]] > 0, 1, 0)
    res <- RunModel(trait_binary, trait_names[i], tree)
    if(!is.null(res)){
      csv_trait_results[[trait_names[i]]] <- res
    }
  }

  # Save per-CSV results
  if(length(csv_trait_results) > 0){
    df_csv <- do.call(rbind, csv_trait_results)
    csv_basename <- tools::file_path_sans_ext(basename(csv_file))
    write.csv(df_csv, paste0(csv_basename, "_trait_hisse.csv"), row.names = FALSE)
    all_trait_results[[csv_basename]] <- df_csv
  }
}

# ==========================
# 8. Merge all results and save summary CSV
# ==========================
if(length(all_trait_results) > 0){
  final_results <- do.call(rbind, all_trait_results)
  write.csv(final_results, "all_traits_hisse_summary.csv", row.names = FALSE)
  cat("All CSVs processed! Individual and summary CSV saved.\n")
} else {
  cat("[WARN] No results generated from any CSV.\n")
}
