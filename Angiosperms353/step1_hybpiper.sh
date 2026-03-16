while read name; 
do 
hybpiper assemble -t_dna con_353_finall.fa -r ${name}_R1.fq ${name}_R2.fq --prefix $name --diamond --hybpiper_output /data/xueqin/Project/Prunus/Prunus_353_easy353/add_Prunus_2024.7.24;
done < namelist_08_27_delete_lower_50.txt

hybpiper retrieve_sequences -t_dna con_353_finall.fa dna --sample_names namelist_08_27_delete_lower_50.txt
