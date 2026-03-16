while read name;
do
hybpiper assemble -t_aa con_orthofinder_only_Prunus_9_genome.fa -r ${name}_R1.fq ${name}_R2.fq --prefix $name --diamond --hybpiper_output /path/to/output/FNA;
done < namelist_08_27_delete_lower_50.txt


# Run HybPiper paralog_retriever to extract paralogous sequences
hybpiper paralog_retriever namelist_08_27_delete_lower_50.txt -t_dna con_orthofinder_only_Prunus_9_genome.fa
# Output directories:
#   paralogs_all          - all detected paralog sequences
#   paralogs_no_chimeras  - paralogs with chimeric sequences removed
# Note: A warning is issued when multiple contigs are assembled for a single locus
# HybPiper flags loci where detected paralogs exceed 75% of the contig threshold
# In paralog_report.tsv, genes with coverage > 1 are considered potential paralogs
