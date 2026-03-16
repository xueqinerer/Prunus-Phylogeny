mkdir Prunus_sortadate_1,3,2
cd ../Prunus_sortadate_1,3,2
#Get the root-to-tip variance with:
python /data/xiongtao/scripts/SortaDate/src/get_var_length.py ./Prunus_gene_best_rt/ --flend .tre --outf var --outg Lyonothamnus_floribundus
# Tree_reroot：为包含所有置根后的基因树的文件夹(注意是最优树)
# Barbeya_oleoides：为使用的外类群的名字

#Get the bipartition support with:
python /data/xiongtao/scripts/SortaDate/src/get_bp_genetrees.py ./Prunus_gene_best_rt/ Prunus_iqtree_BS10_species_re.tre --flend .tre --outf bp
#Combine the results from these two runs with 
python /data/xiongtao/scripts/SortaDate/src/combine_results.py var bp --outf comb

nohup bash run_pipeline_ortholog.sh 35 &
nohup bash run_pipeline_ortholog.sh 50 &