#R
 setwd("./")
 clads <- "CLaDS_out.Rdata"
 load(clads)
 tips = sort(CladsOutput$tree$tip.label)
 stats <- data.frame(gelm = NA, length = NA, run = clads)
 res <- matrix(NA, nrow = length(tips), ncol = 1, dimnames = list(tips))
 d <- CladsOutput
 cladsrates <- d$lambdatip_map
 names(cladsrates) <- d$tree$tip.label
 res[, 1] <- cladsrates[tips]
 stats[1, "gelm"] <- d$gelm[2]
 stats[1, "length"] <- length(d$rtt_chains[[1]])
 write.csv(stats, "TACT2_clads.convergence_stats.csv", row.names = FALSE)
 res2 <- res[grep("xzx", dimnames(res)[[1]], invert = T), ]
 write.table(res2, "best_tree_100_imputations_TACT2_CLaDS_rates.txt",
             row.names = T, col.names = F)
