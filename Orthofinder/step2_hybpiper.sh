while read name; 
do 
hybpiper assemble -t_aa con_orthofinder_only_Prunus_9_genome.fa -r ${name}_R1.fq ${name}_R2.fq --prefix $name --diamond --hybpiper_output /data/xueqin/Project/Prunus/Hybpiper_genome_orthofinder/Orthofinder_only_Prunus_9_namelist_08_27_delete_lower_50/FNA/FNA;
done < namelist_08_27_delete_lower_50.txt


**hybpiper paralog_retriever**
hybpiper paralog_retriever namelist_08_27_delete_lower_50.txt -t_dna con_orthofinder_only_Prunus_9_genome.fa
#paralogs_all
#paralogs_no_chimeras
#组装产生多个contigs时则会存在
#HybPiper检测到多个包含长编码序列的contigs时-默认情况下至少75%的基因序列
#paralog_report.tsv：如果> 1，则为潜在的旁系基因
