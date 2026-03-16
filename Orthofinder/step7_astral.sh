mv *.tre ../Prunus_gene_best_rt/

#建物种树
cd ../Prunus_gene_best_rt/
cd ../
pxcat -s ./Prunus_gene_trimal/*.fasta -p Prunus_sp_partition.txt -o Prunus_sp_supermatrix.fasta

mkdir Prunus_species_tree
cd ./Prunus_gene_best_rt/
cat *.tre> ../Prunus_species_tree/Prunus_sp_gene.tre

cd ../Prunus_species_tree/
nw_ed Prunus_sp_gene.tre 'i & b<=10' o > Prunus_353-BS10.treefile
#astral-pro3构建物种树
/data/xueqin/software/ASTER/bin/astral-pro3 -i Prunus_353-BS10.treefile -o Prunus_sp_gene_astral-BS10-pro3.tre --root Lyonothamnus_floribundus
#astral构建物种树
java -jar /home/xueqin/data/software/ASTRAL-master/Astral/astral.5.7.8.jar -i Prunus_353-BS10.treefile --outgroup Lyonothamnus_floribundus -o Prunus_Astral_species.tre

#置根
nw_reroot -l Prunus_sp_gene_astral-BS10-pro3.tre Lyonothamnus_floribundus >Prunus-BS10-pro3_iqtree_BS10_species_re.tre
