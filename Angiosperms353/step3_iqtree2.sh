cd Prunus_gene_trimal

for i in *.fasta;do
iqtree2 -s $i -m MFP --alrt 1000 -B 1000 -T 20 --prefix ../Prunus_gene_tre_raw/$i.tre
done
