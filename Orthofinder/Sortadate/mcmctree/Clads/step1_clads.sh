#Prunus.treePL.dated_no_outgroup_change_orientation.tre是Prunus.treePL.dated_no_outgroup.tre旋转分支后的结果（总状-伞房-单花）
  
# 运行：
julia
 using PANDA 
 using Serialization
 
 my_tree = load_tree("../BAMM/Prunus.mcmctree.dated_no_outgroup.tre")
 output = infer_ClaDS(my_tree, print_state = 100, f = 0.241982507)
 open("output_data.jls", "w") do io
           serialize(io, output)
       end

# 保存与计算tip rates
Julia
 save_ClaDS_in_R(output, "CLaDS_out.Rdata")


 output = open(deserialize, "output_data.jls")
 plot_CladsOutput(output, options = "type = 'phylogram'")  #绘出多样化速率在时间树上的树形图
 plot_CladsOutput(output, method = "DTT")  #绘出多样化速率随时间的变化图（Diversity through time plot）。黑色加粗的线表示LTT图，蓝色细线表示独立的MCMC迭代，蓝色粗线表示95%置信区间，绿色虚线表示点估计
 plot_CladsOutput(output, method = "RTT")  #绘出平均速率随时间的变化图（Mean rate through time plot）。蓝色细线表示独立的MCMC迭代，蓝色粗线表示95%置信区间，绿色虚线表示点估计

# 保存与计算tip rates
Julia
 save_ClaDS_in_R(output, "CLaDS_out.Rdata")

