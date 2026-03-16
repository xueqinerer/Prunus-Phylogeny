# A script to calculate DR statistics from Jetz et al. (2012)
rm(list=ls())
library("ape")
library("phytools")

# Input tree file from bash stream
treefile <- commandArgs(TRUE)

# Ensure treefile is provided
if (length(treefile) == 0) {
  stop("No tree file provided as input.")
}

# Define function from Jetz et al. (2012), also see Harvey et al. (2016)
DR_statistic <- function(tree, return.mean = FALSE) {
  rootnode <- length(tree$tip.label) + 1
  sprates <- numeric(length(tree$tip.label))
  for (i in 1:length(sprates)) {
    node <- i
    index <- 1
    qx <- 0
    while (node != rootnode) {
      el <- tree$edge.length[tree$edge[, 2] == node]
      node <- tree$edge[, 1][tree$edge[, 2] == node]
      qx <- qx + el * (1 / 2^(index - 1))
      index <- index + 1
    }
    sprates[i] <- 1 / qx
  }
  if (return.mean) {
    return(mean(sprates))
  } else {
    names(sprates) <- tree$tip.label
    return(sprates)
  }
}

# Function to adjust the tree to make it ultrametric if necessary
node_adjust_ultrametric <- function(phy, node.diff.log) {
  N <- Ntip(phy)
  tre_node_adjust <- reorder(phy, "postorder")
  e1 <- tre_node_adjust$edge[, 1]
  e2 <- tre_node_adjust$edge[, 2]
  EL <- tre_node_adjust$edge.length

  ages <- numeric(N + tre_node_adjust$Nnode)

  for (ii in seq_along(EL)) {
    if (ages[e1[ii]] == 0) {
      ages[e1[ii]] <- ages[e2[ii]] + EL[ii]
    } else {
      recorded_age <- ages[e1[ii]]
      new_age <- ages[e2[ii]] + EL[ii]
      if (recorded_age != new_age) {
        if (node.diff.log) {
          cat(sprintf("node %i age %.6f != %.6f\n", e1[ii], recorded_age, new_age))
        }
        EL[ii] <- recorded_age - ages[e2[ii]]
      }
    }
  }
  tre_node_adjust$edge.length <- EL
  return(tre_node_adjust)
}

# Read tree
tree <- read.tree("/path/to/BAMM/Prunus.mcmctree.dated_no_outgroup.tre")

# Generate the output file name
out.name <- gsub(".tre", ".DR.csv", basename(treefile))

# Define output directory and path
output_dir <- "./results/"
output_path <- file.path(output_dir, out.name)

# Print the output path for debugging
print(paste("Output will be saved to:", output_path))

# Ensure directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Calculate DR statistics
if (is.ultrametric(tree)) {
  DR <- DR_statistic(tree)
} else {
  tree1 <- node_adjust_ultrametric(tree, node.diff.log = FALSE)
  DR <- DR_statistic(tree1)
}

# Write results to file
write.csv(data.frame(Species = names(DR), DR = DR), output_path, row.names = FALSE, quote = FALSE)
