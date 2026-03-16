# 加载库
library(ggplot2)
library(reshape2)

# 准备数据：来自你提供的表格
df <- data.frame(
  Tree = c("1to1 vs MI", "1to1 vs MO", "MI vs MO"),
  MatchingSplit = c(31, 22, 11),
  RF = c(9, 7, 4)
)

# 转为长格式
df_long <- melt(df, id.vars = "Tree", variable.name = "Metric", value.name = "Distance")

# 绘图
p <- ggplot(df_long, aes(x = Tree, y = Distance, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           color = "black", width = 0.7) +
  labs(title = "Tree Distance Comparison (MatchingSplit vs RF)",
       y = "Distance (lower means more similar)",
       x = "Tree Pair") +
  scale_fill_manual(values = c("MatchingSplit" = "#299e8c", "RF" = "#f2a361")) +
  theme_minimal(base_size = 20) +
  theme(
    legend.position = c(0.9, 0.9),
    text = element_text(size = 20),
    legend.title = element_blank(),
    axis.text.x = element_text(size = 16,vjust = 5),
    axis.text.y = element_text(size = 16)
  )

# 输出图像
ggsave("tree_distance_comparison.png", p, width = 5, height = 4, dpi = 300)

pdf("tree_distance_comparison.pdf", height = 12, width = 10)
print(p)
dev.off()