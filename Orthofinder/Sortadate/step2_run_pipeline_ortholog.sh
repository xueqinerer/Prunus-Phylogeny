#!/bin/bash

# Ensure an argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <number>"
    exit 1
fi

num=$1

# Set paths based on the input parameter
base_path="./"
sortadate_path="${base_path}/sortadate_${num}"

# Select high-quality genes for dating
python /path/to/scripts/SortaDate/src/get_good_genes.py comb --max $num --order 3,1,2 --outf result_${num}.txt

cat result_${num}.txt | sed '1d' | cut -f1 -d' ' > tree_name_${num}.txt

sed 's/.treefile.tre//g' tree_name_${num}.txt > gene_sortadate_name_${num}.txt

mkdir -p ${sortadate_path}

# Copy selected gene alignment files in parallel
parallel cp ../Prunus_gene_trimal/{} ${sortadate_path}/{} < gene_sortadate_name_${num}.txt

# Enter the target directory
cd ${sortadate_path}

# Move fasta files from Prunus_gene_trimal directory
mkdir -p Prunus_gene_trimal
mv *.fasta ./Prunus_gene_trimal/

# Generate supermatrix from selected gene alignments
pxcat -s ./Prunus_gene_trimal/*.fasta -p Prunus_sp_partition.txt -o Prunus_sp_supermatrix.fasta

# Create Prunus_species_tree directory for IQTree2 constrained analysis
mkdir -p Prunus_species_tree
cd Prunus_species_tree/

iqtree2 -s ../Prunus_sp_supermatrix.fasta -p ../Prunus_sp_partition.txt -m MFP -g ../../Prunus_iqtree_BS10_species_re.tre -T 5 --prefix Prunus_Astral_iqtree
nw_reroot Prunus_Astral_iqtree.treefile Lyonothamnus_floribundus > Prunus_Astral_species_br_re.tre
