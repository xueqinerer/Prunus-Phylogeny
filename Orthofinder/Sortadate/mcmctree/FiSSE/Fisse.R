
setwd("C:/Users/301/Desktop/李属文章数据/BAMM_MO_35,1,3,2/Fisse")

#source("./traitDependent_functions.R")
library(ape)
# 读取性状表
traits <- read.csv("Prunus_all_traits.csv", header = TRUE)
head(traits)
# 处理 ploidy 列
traits$ploidy[traits$ploidy == "?"] <- NA        # 将 ? 转为 NA
traits$ploidy <- as.numeric(traits$ploidy)       # 转成数值型

# 检查状态分布
table(traits$ploidy, useNA = "ifany")

# 过滤掉 NA
ploidy_df <- traits[!is.na(traits$ploidy), c("taxon_name", "ploidy")]

# 载入分子树
library(ape)
tree <- read.tree("./Prunus.mcmctree.dated_no_outgroup.tre")

# 保留树上存在的物种
ploidy_df <- ploidy_df[ploidy_df$taxon_name %in% tree$tip.label, ]

# 检查是否每个状态至少有2个物种
state_counts <- table(ploidy_df$ploidy)
state_counts
if(any(state_counts < 2)) {
  warning("某个状态物种数 < 2，FiSSE 可能无法运行")
}

# 构建 FiSSE 可用向量
ploidy_vec <- ploidy_df$ploidy
names(ploidy_vec) <- ploidy_df$taxon_name

# 现在 ploidy_vec 就可以直接用于 runAnalyses 或 FISSE.binary
ploidy_vec
# 原始树
tree <- read.tree("../Prunus.mcmctree.dated_no_outgroup.tre")

# 保留 tree 上存在的物种
common_species <- intersect(tree$tip.label, names(ploidy_vec))
tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))

# 对 trait 向量也只保留这些物种
ploidy_vec <- ploidy_vec[common_species]

length(tree_pruned$tip.label)   # 检查剪枝后物种数
length(ploidy_vec)              # 与 tree_pruned 一致

library(phytools)

is.ultrametric(tree_pruned)  # TRUE / FALSE
# 如果 FALSE，可以用
tree_pruned <- force.ultrametric(tree_pruned, method="extend")

# FISSE.binary(tree, trait, reps = 1000)
set.seed(123)   # 保证可重复
fisse_result <- FISSE.binary(tree_pruned, ploidy_vec, reps = 1000)

# 查看结果
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

setwd("C:/Users/301/Desktop/李属文章数据/BAMM_MO_35,1,3,2/Fisse/")

# 读取性状表
traits <- read.csv("Prunus_all_traits.csv", header = TRUE)

# 载入分子树
tree <- read.tree("../Prunus.mcmctree.dated_no_outgroup.tre")

# 确保 ultrametric
if(!is.ultrametric(tree)) tree <- force.ultrametric(tree, method="extend")

# 要分析的性状列
trait_cols <- c("ploidy", "lifeform", "perianth")

# 创建空的结果数据框
all_results <- data.frame()

# 遍历性状列
for(trait_name in trait_cols){
  
  cat("\n➡️  Processing trait:", trait_name, "\n")
  
  # 提取向量，并确保 numeric
  trait_vec <- as.numeric(as.character(traits[[trait_name]]))
  names(trait_vec) <- traits$taxon_name
  
  # 替换不合法值为 NA
  trait_vec[!(trait_vec %in% c(0,1))] <- NA
  trait_vec <- trait_vec[!is.na(trait_vec)]
  
  # 保留树上存在的物种
  common_species <- intersect(tree$tip.label, names(trait_vec))
  if(length(common_species) == 0){
    cat("⚠️ No species from trait in the tree. Skipping.\n")
    next
  }
  trait_vec <- trait_vec[common_species]
  tree_pruned <- drop.tip(tree, setdiff(tree$tip.label, common_species))
  
  # 检查每个状态至少有2个物种
  state_counts <- table(trait_vec)
  if(length(state_counts) < 2 || any(state_counts < 2)){
    cat("⚠️ Not enough species in each state. Skipping.\n")
    next
  }
  
  # 运行 FiSSE
  set.seed(123)
  fisse_result <- FISSE.binary(tree_pruned, trait_vec, reps = 1000)
  
  # 转换为数据框并添加 trait 列
  res_df <- as.data.frame(fisse_result)
  res_df$trait <- trait_name
  
  # 将结果加入总表
  all_results <- rbind(all_results, res_df)
  
  cat("✅  FiSSE finished for trait", trait_name, "\n")
}

# 保存汇总结果
write.csv(all_results, "FISSE_all_traits_results.csv", row.names = FALSE)
cat("🎉 All FiSSE results saved to FISSE_all_traits_results.csv\n")



