mkdir FNA
mv *.FNA ./FNA
mkdir ./Prunus_gene_mafft
cd FNA
for gene in ./*.FNA
do
    # Remove extra information from sequence headers
    sed -i 's/_[0-9]*_recovered_.*$//g' "$gene"

    # Perform multiple sequence alignment and save output to file
    mafft --auto --thread 15 "$gene" > "${gene%.*}.mafft"
done

 mv ./*.mafft ../Prunus_gene_mafft


cd ../

mkdir ./Prunus_gene_trimal
for gene in ./Prunus_gene_mafft/*.mafft; do
    trimal -in "$gene" -out "${gene}.tri.fasta" -automated1
    mv "${gene}.tri.fasta" ./Prunus_gene_trimal
done
