library("ape")
suppressWarnings(suppressMessages(library("phytools")))
source(paste0("/data/xueqin/software/mad/","mad.R", sep="")) #source "mad" function

# 获取以".tre"结尾的文件列表
file_list <- list.files("/data/xueqin/Project/Prunus/Hybpiper_genome_orthofinder/Orthofinder_only_Prunus_9_namelist_08_27_delete_lower_50/Prunus_paralog/1old_treeshrink0.2_final/26_MO_final_alignments_trimmed/Prunus_treeshrink/Prunus_treeshrink_trimal/Prunus_gene_tre/gene_tree", pattern = "\\.treefile$", full.names = TRUE)

# 循环处理每个文件
for (file in file_list) {
  # 读取文件
  g6825 <- read.tree(file)
  
  # 检查是否存在"Lyonothamnus_floribundus"
  if (!"Lyonothamnus_floribundus" %in% g6825$tip.label) {
    # 计算MAD
    tmp <- mad(g6825)
    
    # 解析为树对象并重新排序
    Tg6689.mad <- ladderize(read.tree(text = tmp))
    
    # 构造输出文件路径
    output_file <- paste0(file,".tre")
    # 写入文件
    write.tree(Tg6689.mad, file = output_file)
         
    # 打印处理完成的文件名
    cat("Processed file:", file, "\n")
  }
}