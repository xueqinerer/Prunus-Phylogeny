#!/bin/bash

# Step 2: Remove samples with long branches in outgroup (no_delete_outgroup_long_branch)

# Copy gene trees from previous step
cp ../Prunus_gene_best_rt/RAxML_bestTree* ./
for i in *.tre; do tt=$(echo $i|sed 's/RAxML_bestTree.//g;s/.mafft.aln.tri.fasta.ML.tre/.tree/g;s/g//g'); mv $i $tt; done
for i in *.tree; do pxrmt -t $i -n Lyonothamnus_floribundus > "${i%.tree}_modified.tree"; done

# Rename modified tree files
for i in *.tree; do tt=$(echo $i|sed 's/_modified.tree/.tree/g;s/g//g'); mv $i $tt; done

# Copy alignment files
cp ../Prunus_gene_trimal/*.fasta ./

# Rename alignment files
for i in *.fasta; do mv $i $(basename $i .mafft.aln.tri.fasta).fasta; done
while read -r Line; do mkdir -p $Line; mv "${Line}.tree" "${Line}.fasta" "$Line"; done < ../Angio353_gene_list.txt

# Rename gene trees and alignments within each directory to input.tree and input.fasta
for name in *; do cd $name; mv $name.fasta input.fasta; mv $name.tree input.tree; cd ../; done

# Run TreeShrink with q-value 0.2 to remove outlier long branches
nohup /path/to/envs/ParaGone/bin/run_treeshrink.py -i ./ -t input.tree -q 0.2 -a input.fasta -m per-species -b 20 > ./input.tree.treeshrinklog.txt &
mkdir original_dir

# Rename output fasta files to species-named .output.fasta
for name in *; do
    if [ -d "$name" ]; then
        cd "$name"
        mv output.fasta "$name".output.fasta
        cd ..
    fi
done

# Create output directory for filtered alignments
mkdir Prunus_treeshrink_trimal

# Summarize filtering results across all gene folders into all_genes_info.txt
output_file="all_genes_info.txt"
> "$output_file"  # Clear existing content if file exists

# Iterate through all gene folders
for dir in */; do
    # Check if it's a directory containing an output.txt file
    if [ -d "$dir" ] && [ -f "$dir/output.txt" ]; then
        # Get the gene folder name
        gene_name=$(basename "$dir")

        # Get content from output.txt, using semicolons as separator
        output_content=$(cat "$dir/output.txt" | tr '\n' ';')

        # Write to the output file
        echo "$gene_name: $output_content" >> "$output_file"
    fi
done


# Step 2: Remove outgroup sequences marked in delete_outgroup.txt,
# then remove related sequences and symbols
bash remove_delete_outgroup.txt.sh
for name in *; do
    if [ -d "$name" ]; then
        cd "$name"
        mv filtered.fasta "$name".filtered.fasta
        cd ..
    fi
done
for name in *; do
    if [ -d "$name" ]; then
        cd "$name"
        cp "$name".filtered.fasta ../Prunus_treeshrink_trimal
        cd ..
    fi
done
cd Prunus_treeshrink_trimal

# Merge all filtered alignments into a supermatrix
pxcat -s ./Prunus_treeshrink_trimal/*.fasta -p Prunus_sp_partition_treeshrink_0.2.txt -o Prunus_sp_supermatrix_treeshrink_0.2.fasta

