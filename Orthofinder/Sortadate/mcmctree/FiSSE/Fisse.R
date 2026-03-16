setwd("/path/to/FiSSE")

#source("./traitDependent_functions.R")
library(ape)
# Read trait table
traits <- read.csv("Prunus_all_traits.csv", header = TRUE)
head(traits)
# Process ploidy column
traits$ploidy[traits$ploidy == "?"] <- NA        # Replace ? with NA
traits$ploidy <- as.numeric(traits$ploidy)       # Convert to numeric

# Check state distribution
table(traits$ploidy, useNA = "ifany")

# Filter out NA values
ploidy_df <- traits[!is.na(traits$ploidy), c("taxon_name", "ploidy")]

# Load the molecular phylogeny
library(ape)
tree <- read.tree("./Prunus.mcmctree.dated_no_outgroup.tre")

# Keep only species present in the tree
ploidy_df <- ploidy_df[ploidy_df$taxon_name %in% tree$tip.label, ]

# Check that each state has at least 2 species
state_counts <- table(ploidy_df$ploidy)
state_counts
if(any(state_counts < 2)) {
  warning("A state has fewer than 2 species; FiSSE may not run")
}

# Build the named vector for FiSSE
ploidy_vec <- ploidy_df$ploidy
names(ploidy_vec) <- ploidy_df$taxon_name

ploidy_vec
# Original full tree
tree <- read.tree("../Prunus.mcmctree.dated_no_outgroup.tre")

# Keep only species present in both tree and trait vector
common_species <- intersect(tree$tip.label, names(ploidy_vec))
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))

# Subset trait vector to matching species
ploidy_vec <- ploidy_vec[common_species]

length(tree_pruned$tip.label)   # Check species count after pruning
length(ploidy_vec)              # Should match tree_pruned

library(phytools)

is.ultrametric(tree_pruned)  # Check if tree is ultrametric
# If FALSE, force ultrametric
tree_pruned <- force.ultrametric(tree_pruned, method="extend")

# Run FiSSE binary test (reps = 1000 permutations)
set.seed(123)
fisse_result <- FISSE.binary(tree_pruned, ploidy_vec, reps = 1000)

# View results
fisse_result

write.csv(as.data.frame(fisse_result), "FISSE_ploidy_results.csv", row.names = FALSE)



####################################################
library(phytools)
library(ape)
library(geiger)
library(caper)
library(igraph)
library(ggraph)
library(ggplot2)
library(dplyr)
library(diversitree)

setwd("/path/to/FiSSE")

# Read trait table
traits <- read.csv("Prunus_all_traits.csv", header = TRUE)

# Load the molecular phylogeny
tree <- read.tree("../Prunus.mcmctree.dated_no_outgroup.tre")

# Ensure ultrametric
if(!is.ultrametric(tree)) tree <- force.ultrametric(tree, method="extend")

# Trait columns to analyze
trait_cols <- c("ploidy", "lifeform", "perianth")

# Create empty results data frame
all_results <- data.frame()

# Loop over trait columns
for(trait_name in trait_cols){

  cat("\n Processing trait:", trait_name, "\n")

  # Extract trait vector and ensure numeric
  trait_vec <- as.numeric(as.character(traits[[trait_name]]))
  names(trait_vec) <- traits$taxon_name

  # Replace invalid values with NA
  trait_vec[!(trait_vec %in% c(0,1))] <- NA
  trait_vec <- trait_vec[!is.na(trait_vec)]

  # Keep only species present in the tree
  common_species <- intersect(tree$tip.label, names(trait_vec))
  if(length(common_species) == 0){
    cat("No species from trait found in the tree. Skipping.\n")
    next
  }
  trait_vec <- trait_vec[common_species]
  tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))

  # Check that each state has at least 2 species
  state_counts <- table(trait_vec)
  if(length(state_counts) < 2 || any(state_counts < 2)){
    cat("Not enough species in each state. Skipping.\n")
    next
  }

  # Run FiSSE
  set.seed(123)
  fisse_result <- FISSE.binary(tree_pruned, trait_vec, reps = 1000)

  # Convert to data frame and add trait column
  res_df <- as.data.frame(fisse_result)
  res_df$trait <- trait_name

  # Append to summary table
  all_results <- rbind(all_results, res_df)

  cat("FiSSE finished for trait", trait_name, "\n")
}

# Save summary results
write.csv(all_results, "FISSE_all_traits_results.csv", row.names = FALSE)
cat("All FiSSE results saved to FISSE_all_traits_results.csv\n")
