mkdir Prunus_sortadate_1,3,2
cd ../Prunus_sortadate_1,3,2

# Get the root-to-tip variance with:
python /path/to/scripts/SortaDate/src/get_var_length.py ./Prunus_gene_best_rt/ --flend .tre --outf var --outg Lyonothamnus_floribundus
# Tree_reroot: reference species tree file used for subsequent rerooting (note the file paths)
# Lyonothamnus_floribundus: outgroup species used for rerooting

# Get the bipartition support with:
python /path/to/scripts/SortaDate/src/get_bp_genetrees.py ./Prunus_gene_best_rt/ Prunus_iqtree_BS10_species_re.tre --flend .tre --outf bp
# Combine the results from these two runs with:
python /path/to/scripts/SortaDate/src/combine_results.py var bp --outf comb

nohup bash run_pipeline_ortholog.sh 35 &
nohup bash run_pipeline_ortholog.sh 50 &
