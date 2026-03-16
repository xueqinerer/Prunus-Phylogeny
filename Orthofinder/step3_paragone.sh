paragone check_and_align ./paralogs_all --internal_outgroup Lyonothamnus_floribundus --pool 1 --threads 20
paragone alignment_to_tree 04_alignments_trimmed_cleaned --pool 1 --threads 20
paragone qc_trees_and_extract_fasta 04_alignments_trimmed_cleaned --treeshrink_q_value 0.20 --cut_deep_paralogs_internal_branch_length_cutoff 0.04
paragone align_selected_and_tree 04_alignments_trimmed_cleaned --pool 1 --threads 20
paragone prune_paralogs --mo
paragone final_alignments --mo --pool 1 --threads 20 --keep_intermediate_files
