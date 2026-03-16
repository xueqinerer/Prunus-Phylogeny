mv *.tre ../Prunus_gene_best_rt/

# Move all rerooted gene trees to the Prunus_gene_best_rt directory
cd ../Prunus_gene_best_rt/
cd ../
pxcat -s ./Prunus_gene_trimal/*.fasta -p Prunus_sp_partition.txt -o Prunus_sp_supermatrix.fasta

mkdir Prunus_species_tree
cd ./Prunus_gene_best_rt/
cat *.tre > ../Prunus_species_tree/Prunus_sp_gene.tre

cd ../Prunus_species_tree/
# Collapse gene tree nodes with bootstrap support <= 10
nw_ed Prunus_sp_gene.tre 'i & b<=10' o > Prunus_353-BS10.treefile

# Species tree inference using ASTRAL-pro3
/path/to/software/ASTER/bin/astral-pro3 -i Prunus_353-BS10.treefile -o Prunus_sp_gene_astral-BS10-pro3.tre --root Lyonothamnus_floribundus

# Species tree inference using ASTRAL
java -jar /path/to/software/ASTRAL/astral.5.7.8.jar -i Prunus_353-BS10.treefile --outgroup Lyonothamnus_floribundus -o Prunus_Astral_species.tre

# Reroot the species tree with Lyonothamnus_floribundus as outgroup
nw_reroot -l Prunus_sp_gene_astral-BS10-pro3.tre Lyonothamnus_floribundus > Prunus-BS10-pro3_iqtree_BS10_species_re.tre
