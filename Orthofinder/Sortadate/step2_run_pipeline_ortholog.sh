#!/bin/bash

# 确认传入参数
if [ $# -ne 1 ]; then
    echo "Usage: $0 <number>"
    exit 1
fi

num=$1

# 路径定义，方便后面调用
base_path="./"
sortadate_path="${base_path}/sortadate_${num}"

# 运行主要流程
python /data/xiongtao/scripts/SortaDate/src/get_good_genes.py comb --max $num --order 3,1,2 --outf result_${num}.txt

cat result_${num}.txt | sed '1d' | cut -f1 -d' ' > tree_name_${num}.txt

sed 's/.treefile.tre//g' tree_name_${num}.txt> gene_sortadate_name_${num}.txt

mkdir -p ${sortadate_path}

# parallel 复制文件
parallel cp ../Prunus_gene_trimal/{} ${sortadate_path}/{} < gene_sortadate_name_${num}.txt

# 进入对应目录
cd ${sortadate_path}

# 创建Prunus_gene_trimal目录并移动fasta文件
mkdir -p Prunus_gene_trimal
mv *.fasta ./Prunus_gene_trimal/

# 生成supermatrix
pxcat -s ./Prunus_gene_trimal/*.fasta -p Prunus_sp_partition.txt -o Prunus_sp_supermatrix.fasta
#sed -i 's/Prunus_arborea_var/Prunus_arborea_var._montana/g;s/Prunus_davidiana_var/Prunus_davidiana_var._potaninii/g;s/Prunus_salicina_var/Prunus_salicina_var._mandshurica/g' Prunus_sp_supermatrix.fasta
# 进入Prunus_species_tree目录，准备RAxML分析
mkdir -p Prunus_species_tree
cd Prunus_species_tree/

iqtree2 -s ../Prunus_sp_supermatrix.fasta -p ../Prunus_sp_partition.txt -m MFP -g ../../Prunus_iqtree_BS10_species_re.tre -T 5 --prefix Prunus_Astral_iqtree
nw_reroot Prunus_Astral_iqtree.treefile Lyonothamnus_floribundus >Prunus_Astral_species_br_re.tre
