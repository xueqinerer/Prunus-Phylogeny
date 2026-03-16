# Run CLaDS (Cladogenetic and Anagenetic Diversification Shifts) analysis
# using the PANDA Julia package
# Input: dated phylogenetic tree without outgroup
# Sampling fraction: f = 83/343 = 0.241982507

# Run CLaDS inference in Julia
julia
 using PANDA
 using Serialization

 my_tree = load_tree("../BAMM/Prunus.mcmctree.dated_no_outgroup.tre")
 output = infer_ClaDS(my_tree, print_state = 100, f = 0.241982507)
 open("output_data.jls", "w") do io
           serialize(io, output)
       end

# Export tip rates for use in R
Julia
 save_ClaDS_in_R(output, "CLaDS_out.Rdata")


 output = open(deserialize, "output_data.jls")
 # Plot phylogram with diversification rates mapped onto branches
 plot_CladsOutput(output, options = "type = 'phylogram'")
 # Plot Diversity Through Time (DTT) - dark line: LTT; thin line: MCMC samples; bold: 95% CI; dashed: mean
 plot_CladsOutput(output, method = "DTT")
 # Plot Mean Rate Through Time (RTT) - thin line: MCMC samples; bold: 95% CI; dashed: mean
 plot_CladsOutput(output, method = "RTT")

# Export tip rates for use in R
Julia
 save_ClaDS_in_R(output, "CLaDS_out.Rdata")
