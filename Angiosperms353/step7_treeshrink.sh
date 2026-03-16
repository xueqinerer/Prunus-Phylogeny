#2.不删除外类群是长枝的情况(no_delete_outgroup_long_branch)
#基因
cp ../Prunus_gene_best_rt/RAxML_bestTree* ./
for i in *.tre; do tt=$(echo $i|sed 's/RAxML_bestTree.//g;s/.mafft.aln.tri.fasta.ML.tre/.tree/g;s/g//g'); mv $i $tt; done
for i in *.tree; do   pxrmt -t $i -n Lyonothamnus_floribundus > "${i%.tree}_modified.tree"; done
#手动删除$i的信息
for i in *.tree; do tt=$(echo $i|sed 's/_modified.tree/.tree/g;s/g//g'); mv $i $tt; done
#序列
cp ../Prunus_gene_trimal/*.fasta ./
#改名
for i in *.fasta; do mv $i $(basename $i .mafft.aln.tri.fasta).fasta; done
while read -r Line; do     mkdir -p $Line    ;     mv "${Line}.tree" "${Line}.fasta" "$Line"; done<../Angio353_gene_list.txt
#修改文件夹内的基因名和序列名，分别为input.tree和input.fasta
for name in *;do cd $name; mv $name.fasta input.fasta; mv $name.tree input.tree; cd ../; done
#treeshrink -q 0.2
nohup /home/xueqin/miniconda3/envs/ParaGone/bin/run_treeshrink.py -i ./ -t input.tree -q 0.2 -a input.fasta -m per-species -b 20 > ./input.tree.treeshrinklog.txt &
mkdir original_dir
#修改输出的fasta名称为基因名.output.fasta
for name in *; do     if [ -d "$name" ]; then         cd "$name";         mv output.fasta "$name".output.fasta;         cd ..;     fi; done
#拷贝序列
mkdir Prunus_treeshrink_trimal
#统计过滤情况：all_genes_info.txt
#bash stastics_output.sh
# 创建一个新的文件用于保存输出结果
output_file="all_genes_info.txt"
> "$output_file"  # 清空文件（如果文件已存在）

# 遍历所有基因文件夹
for dir in */; do
    # 检查是否是一个目录，并且该目录下有 output.txt 文件
    if [ -d "$dir" ] && [ -f "$dir/output.txt" ]; then
        # 提取基因名（文件夹的名字）
        gene_name=$(basename "$dir")
        
        # 读取 output.txt 中的内容
        output_content=$(cat "$dir/output.txt" | tr '\n' ';')  # 使用分号分隔行
        
        # 写入新的文件
        echo "$gene_name: $output_content" >> "$output_file"
    fi
done


#2.去掉外类群后将all_genes_info.txt更名为delete_outgroup.txt，删除这个名单的物种和序列
bash remove_delete_outgroup.txt.sh
for name in *; do     if [ -d "$name" ]; then         cd "$name";         mv filtered.fasta "$name".filtered.fasta;         cd ..;     fi; done
for name in *; do     if [ -d "$name" ]; then         cd "$name";         cp "$name".filtered.fasta ../Prunus_treeshrink_trimal;         cd ..;     fi; done
for name in *; do     if [ -d "$name" ]; then         cd "$name";         cp "$name".filtered.fasta ../Prunus_treeshrink_trimal;         cd ..;     fi; done
cd Prunus_treeshrink_trimal
#合并矩阵
pxcat -s ./Prunus_treeshrink_trimal/*.fasta -p Prunus_sp_partition_treeshrink_0.2.txt -o Prunus_sp_supermatrix_treeshrink_0.2.fasta

bash step4_iqtree2.sh
bash step5_reroot.sh
Rscript step6_reroot.R
bash step7_astral.sh

