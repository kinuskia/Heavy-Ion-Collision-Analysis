import matplotlib.pyplot as plt
import numpy as np 



#modes = [0, 1, 2, 3, 4, 5]
#modes = [0, 1]
percentiles = np.array([0, 20])
counter_fig = 0

#colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']


for p in range(0, len(percentiles)):
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.figure(figsize=(10,7.3))
	plt.rcParams.update({'font.size': 23})
	plt.rcParams['axes.titlepad'] = 10
	centrality = str(percentiles[p]) + '-' + str(percentiles[p]+1)
	r_trento, W_trento = np.loadtxt("output/weight_functions_"+centrality+".txt", unpack=True)
	r_magma, W_magma = np.loadtxt("Saclay_simplified/output/"+centrality+"/weight_functions_magma_"+centrality+".txt", unpack=True)
	plt.xlabel("$r$ (fm)")
	plt.ylabel("$W(r)$ (1/fm$^2$)")
	plt.plot(r_trento, W_trento, linestyle = "dashed", color='darkblue', linewidth=3, label="TrENTo")
	plt.plot(r_magma, W_magma, linestyle = "dotted", color='#db6e0d', linewidth=4, label="CGC")
	plt.legend(loc='best')
	centrality_class = centrality + '%'
	plt.title("$W(r)$, "+centrality_class)
	filename = "plots/W_comparison_" + str(percentiles[p]) + "-" + str(percentiles[p]+1) + ".pdf"
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close(counter_fig)




