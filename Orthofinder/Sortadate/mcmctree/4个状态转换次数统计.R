setwd("C:/Users/301/Desktop/李属文章数据/BAMM_MO_35,1,3,2/Ancestral_traits/data/csv/")


library(phytools)
# -------------------------
# 1. 花序类型
# ------------------------
phy <- read.tree("../../Prunus.mcmctree.dated_no_outgroup.tre")
dat <- read.csv("../Inflorescence_2025.12.1/inflorescence_types1.csv", stringsAsFactors = FALSE)

dat_clean <- dat[!is.na(dat$classify) & dat$classify != "?", ]

# 只保留树上存在的物种
dat_clean <- dat_clean[dat_clean$taxon_name %in% phy$tip.label, ]

# 转为 factor
states_factor <- setNames(as.factor(dat_clean$classify), dat_clean$taxon_name)

# 3️⃣ 裁剪树，去掉没有性状的物种
tips_to_drop <- setdiff(phy$tip.label, names(states_factor))
phy_pruned <- drop.tip(phy, tips_to_drop)

# 确认匹配
all(names(states_factor) %in% phy_pruned$tip.label)

# 4️⃣ 运行 make.simmap
simmap <- make.simmap(phy_pruned, states_factor, model = "ER", nsim = 1000)

# 5️⃣ 提取转换次数
simmap_summary <- summary(simmap)
simmap_summary

# 1000 trees with a mapped discrete character with states:
#   0, 1, 2 
# 
# trees have 3.501 changes between states on average
# 
# changes are of the following types:
#   0,1   0,2   1,0   1,2   2,0   2,1
# x->y 0.152 0.085 1.155 0.382 0.852 0.875
# 
# mean total time spent in each state is:
#   0           1           2    total
# raw  534.7753052 575.9809750 966.1918350 2076.948
# prop   0.2574813   0.2773208   0.4651979    1.000

# 查看 summary 对象里都有什么
names(simmap_summary)



# -------------------------
# 2. 花瓣有无
# ------------------------

phy <- read.tree("../../Prunus.mcmctree.dated_no_outgroup.tre")
dat <- read.csv("petal_types.csv", stringsAsFactors = FALSE)

dat_clean <- dat[!is.na(dat$classify) & dat$classify != "?", ]

# 只保留树上存在的物种
dat_clean <- dat_clean[dat_clean$taxon_name %in% phy$tip.label, ]

# 转为 factor
states_factor <- setNames(as.factor(dat_clean$classify), dat_clean$taxon_name)

# 3️⃣ 裁剪树，去掉没有性状的物种
tips_to_drop <- setdiff(phy$tip.label, names(states_factor))
phy_pruned <- drop.tip(phy, tips_to_drop)

# 确认匹配
all(names(states_factor) %in% phy_pruned$tip.label)

# 4️⃣ 运行 make.simmap
simmap <- make.simmap(phy_pruned, states_factor, model = "ER", nsim = 1000)

# 5️⃣ 提取转换次数
simmap_summary <- summary(simmap)
simmap_summary
# 
# 1000 trees with a mapped discrete character with states:
#   0, 1 
# 
# trees have 2.157 changes between states on average
# 
# changes are of the following types:
#   0,1   1,0
# x->y 0.031 2.126
# 
# mean total time spent in each state is:
#   0            1    total
# raw  280.1501548 1796.7979604 2076.948
# prop   0.1348855    0.8651145    1.000


# -------------------------
# 3. 倍性
# ------------------------

phy <- read.tree("../../Prunus.mcmctree.dated_no_outgroup.tre")
dat <- read.csv("chromosome_type.csv", stringsAsFactors = FALSE)

dat_clean <- dat[!is.na(dat$classify) & dat$classify != "?", ]

# 只保留树上存在的物种
dat_clean <- dat_clean[dat_clean$taxon_name %in% phy$tip.label, ]

# 转为 factor
states_factor <- setNames(as.factor(dat_clean$classify), dat_clean$taxon_name)

# 3️⃣ 裁剪树，去掉没有性状的物种
tips_to_drop <- setdiff(phy$tip.label, names(states_factor))
phy_pruned <- drop.tip(phy, tips_to_drop)

# 确认匹配
all(names(states_factor) %in% phy_pruned$tip.label)

# 4️⃣ 运行 make.simmap
simmap <- make.simmap(phy_pruned, states_factor, model = "ER", nsim = 1000)

# 5️⃣ 提取转换次数
simmap_summary <- summary(simmap)
simmap_summary

# 1000 trees with a mapped discrete character with states:
#   0, 1 
# 
# trees have 10.345 changes between states on average
# 
# changes are of the following types:
#   0,1  1,0
# x->y 7.325 3.02
# 
# mean total time spent in each state is:
#   0           1   total
# raw  688.9315470 753.1889023 1442.12
# prop   0.4777212   0.5222788    1.00


# -------------------------
# 4.生活型（常绿和落叶）
# ------------------------

phy <- read.tree("../../Prunus.mcmctree.dated_no_outgroup.tre")
dat <- read.csv("lifeform_types.csv", stringsAsFactors = FALSE)

dat_clean <- dat[!is.na(dat$classify) & dat$classify != "?", ]

# 只保留树上存在的物种
dat_clean <- dat_clean[dat_clean$taxon_name %in% phy$tip.label, ]

# 转为 factor
states_factor <- setNames(as.factor(dat_clean$classify), dat_clean$taxon_name)

# 3️⃣ 裁剪树，去掉没有性状的物种
tips_to_drop <- setdiff(phy$tip.label, names(states_factor))
phy_pruned <- drop.tip(phy, tips_to_drop)

# 确认匹配
all(names(states_factor) %in% phy_pruned$tip.label)

# 4️⃣ 运行 make.simmap
simmap <- make.simmap(phy_pruned, states_factor, model = "ARD", nsim = 1000)

# 5️⃣ 提取转换次数
simmap_summary <- summary(simmap)
simmap_summary

# 1000 trees with a mapped discrete character with states:
#   0, 1 
# 
# trees have 2.083 changes between states on average
# 
# changes are of the following types:
#   0,1   1,0
# x->y 1.655 0.428
# 
# mean total time spent in each state is:
#   0           1    total
# raw  1782.3887491 294.5593661 2076.948
# prop    0.8581768   0.1418232    1.000






