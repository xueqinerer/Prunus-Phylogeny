library("ape")
suppressWarnings(suppressMessages(library("phytools")))
source(paste0("/path/to/software/mad/", "mad.R", sep="")) # Load the MAD rerooting function

# Get list of all ".treefile" files in the gene tree directory
file_list <- list.files("/path/to/Prunus_gene_tre/gene_tree", pattern = "\.treefile$", full.names = TRUE)

# Loop through each file
for (file in file_list) {
  # Read the tree file
  g6825 <- read.tree(file)

  # Check if "Lyonothamnus_floribundus" is absent (trees lacking outgroup need MAD rerooting)
  if (!"Lyonothamnus_floribundus" %in% g6825$tip.label) {
    # Apply MAD (Minimum Ancestor Deviation) rerooting
    tmp <- mad(g6825)

    # Ladderize the tree
    Tg6689.mad <- ladderize(read.tree(text = tmp))

    # Generate output file path
    output_file <- paste0(file, ".tre")
    # Write tree to file
    write.tree(Tg6689.mad, file = output_file)

    # Print the processed file name
    cat("Processed file:", file, "\n")
  }
}
