cd ../
mkdir Prunus_gene_tre
cp ./Prunus_gene_tre_raw/*.treefile ./Prunus_gene_tre
cd Prunus_gene_tre
#置根
mkdir ../Prunus_gene_best_rt #（当前在gene_alns/Prunus_gene_tre/gene_tree目录下）

#!/bin/bash

  for file in *.treefile; do
  if grep -q "Lyonothamnus_floribundus" "$file"; then
    # Reroot the tree using Lyonothamnus_floribundus as the outgroup
    nw_reroot "$file" Lyonothamnus_floribundus > "$file.tre"
    
    # Move the rerooted tree to the Prunus_gene_best_rt directory
    # and rename it to best.rt.tre
    mv "$file.tre" ../Prunus_gene_best_rt
    
    # Remove the original tree file
    rm "$file"
  fi
done
