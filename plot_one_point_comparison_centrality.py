import matplotlib.pyplot as plt
import numpy as np 


modes = [2, 4]
radials = [1,2, 3, 4, 5]

figs, axs = plt.subplots(len(modes), len(radials), figsize=(8.0,3.0), sharex="col", sharey="row")
colors = ['darkblue', '#db6e0d', '#2ca02c', '#d62728', '#9467bd']


counter_fig = 1
clm = np.loadtxt("output/clm.txt")
for l in radials:
	for m in range(0, len(modes)):
		mode = modes[m]
		# compute G_1^(mode) for each model as a function of centrality
		percentiles = np.arange(95)
	
		G_Trento = np.zeros(len(percentiles))
		G_Trento_err = np.zeros(len(percentiles))
		G_magma = np.zeros(len(percentiles))

		for p in range(0,len(percentiles)):
			centrality_class = str(percentiles[p]) + "-" + str(percentiles[p]+1)

			## import profiles

			filename_trento_one = 'output/one_point_' + centrality_class +'.txt'
			filename_trento_one_err = 'output/one_point_' + centrality_class +'_error.txt'
			filename_CGC_one = 'Saclay_simplified/output/'+centrality_class+'/background_coeffs_CGC_' + centrality_class +'.txt'
	

			one_points = np.loadtxt(filename_trento_one)
			one_points_err = np.loadtxt(filename_trento_one_err)
			one_points_CGC = np.loadtxt(filename_CGC_one)
	


			## compute correlators
			G_Trento[p] = one_points[mode, l-1]
			G_Trento_err[p] = one_points_err[mode, l-1]
			G_magma[p] = one_points_CGC[mode, l-1]



		# plot
		G_plus = G_Trento + G_Trento_err
		G_minus = G_Trento - G_Trento_err
		G_plus_smooth = np.zeros(len(G_plus))
		G_minus_smooth = np.zeros(len(G_plus))
		size = 3
		for i in range(0, len(G_plus)):
			G_plus_smooth[i] = 0
			G_minus_smooth[i] = 0
			if ((i < size) | (i > len(G_plus)-1-size)):
				G_plus_smooth[i] = G_plus[i]
				G_minus_smooth[i] = G_minus[i]
				continue
			for j in range(i-size, i+size+1):
				G_plus_smooth[i] += G_plus[j]/(2.*size+1)
				G_minus_smooth[i] += G_minus[j]/(2.*size+1)

		axs[m][l-1].fill_between(percentiles, G_plus_smooth, G_minus_smooth ,color = colors[0])
		axs[m][l-1].plot(percentiles, G_Trento, color = colors[0], label="TrENTo")
		axs[m][l-1].plot(percentiles, G_magma, color = colors[1], label="CGC")
		axs[m][l-1].axhline(y=0, color ="black", linestyle="dashed", linewidth=1.0, zorder = 1)
		axs[m][l-1].set_xticks([0,20,40,60,80])
		# if (l==1):
		# 	axs[mode][l-1].set_ylim(ymin=-0.02, ymax=l1max)
		# elif (l==2):
		# 	axs[mode][l-1].set_ylim(ymin=-0.005, ymax=l2max)
		# else:
		# 	axs[mode][l-1].set_ylim(ymin=-0.005, ymax=l3max)
		#axs[l-1][mode].set_xticks([10,30,50,70])
		#axs[l-1][mode].annotate(" \n$|G_"+str(l)+"^{("+str(mode)+")}|$",xy=(5,0.05))
		#axs[l-1][mode].set_yscale("log")
		if (mode == modes[-1]):
			axs[m][l-1].set_xlabel("Centrality (%)")
		#plt.title("$G_"+str(l)+"^{("+str(mode)+")}$")


#plt.legend(loc="best")


# annotate
for l in radials:
	for m in range(0, len(modes)):
		mode = modes[m]
		ymin, ymax = axs[m][l-1].get_ylim()
		xmin, xmax = axs[m][l-1].get_xlim()
		if (mode == 2):
			axs[m][l-1].annotate(" \n$\\bar{\\epsilon}_"+str(l)+"^{("+str(mode)+")}$",xy=(xmax*0.02,ymin*0.9))
		if (mode == 4):
			axs[m][l-1].annotate(" \n$\\bar{\\epsilon}_"+str(l)+"^{("+str(mode)+")}$",xy=(xmax*0.02,ymax*0.7))


handles, labels = axs[-1][-1].get_legend_handles_labels()
leg = figs.legend(handles, labels, loc='upper center', ncol=4)
counter = 0
for line in leg.get_lines():
	if (counter == 0):
		line.set_linewidth(5.0)
	counter = counter + 1
plt.savefig("plots/one_centrality.pdf", format="pdf", bbox_inches="tight")
		










