import matplotlib.pyplot as plt
import numpy as np 


modes = [0, 1, 2, 3, 4]
radials = [1,2,3]

figs, axs = plt.subplots(len(radials), len(modes), figsize=(8.0,4.0), sharex="col", sharey="row")
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

counter_fig = 1
clm = np.loadtxt("output/clm.txt")
for l in radials:
	for mode in modes:
		# compute G_1^(mode) for each model as a function of centrality
		percentiles = np.arange(95)
		G_IPSM = np.zeros(len(percentiles))
		G_Trento = np.zeros(len(percentiles))
		G_LargeNc = np.zeros(len(percentiles))
		G_magma = np.zeros(len(percentiles))

		for p in range(0,len(percentiles)):
			centrality_class = str(percentiles[p]) + "-" + str(percentiles[p]+1)

			## import profiles

			filename_trento = "output/two_point_random_" + centrality_class + "_m" + str(mode) + "_real.txt"
			filename_trento_one = 'output/one_point_' + centrality_class +'.txt'
			filename_Saclay_simple = "Saclay_simplified/output/"+centrality_class+"/two_point_random_connected_m_" + str(mode) + ".txt"
			filename_Saclay = "Saclay/output/"+centrality_class+"/Nr41/Nm64/m5.0e-3/two_point_random_connected_m_" + str(mode) + ".txt"
			filename_fit = "../IPSM-Fit/output/" + centrality_class + ".txt"
			One_sw2_N, sn2_N_N2 = np.loadtxt(filename_fit, unpack=True)

			trento = np.loadtxt(filename_trento)
			one_points = np.loadtxt(filename_trento_one)
			saclay_simple = np.loadtxt(filename_Saclay_simple)
			saclay = np.loadtxt(filename_Saclay)

			# for i in range(0, len(saclay[:,0])):
			# 	for j in range(0, len(saclay[:,0])):
			# 		if ((i==0)&(j==0)&(mode==0)):
			# 			saclay_simple[i][j] += 1./np.pi/np.pi*sn2_N_N2
			# 			saclay[i][j] += 1./np.pi/np.pi*sn2_N_N2
			# 		else:
			# 			saclay_simple[i][j] += (1.+ sn2_N_N2 )*one_points[mode, i]*one_points[mode, j]
			# 			saclay[i][j] += (1. + sn2_N_N2 )*one_points[mode, i]*one_points[mode, j]

			## compute correlators
			if (mode == 0):
				if (l == 1):	
					G_IPSM[p] = 1.0/np.pi/np.pi*(One_sw2_N +sn2_N_N2)
					G_Trento[p] = trento[0,0] -one_points[0,0]*one_points[0,0]
				else:
					G_IPSM[p] = 1.0/2./np.pi**2/clm[l-2, mode]*One_sw2_N
					G_Trento[p] = trento[l-1,l-1]
			else:
				G_IPSM[p] = (-1.0)**mode/2./np.pi**2/clm[l-1, mode]*One_sw2_N + (1+sn2_N_N2)*one_points[mode, l-1]*one_points[mode, l-1]
				G_Trento[p] = trento[l-1,l-1]
			
			G_LargeNc[p] = saclay[l-1,l-1]
			G_magma[p] = saclay_simple[l-1,l-1]

		# plot
		axs[l-1][mode].plot(percentiles, np.abs(G_IPSM), color = colors[0], label="IPSM")
		axs[l-1][mode].plot(percentiles, np.abs(G_Trento),color = colors[1], label="TrENTo")
		axs[l-1][mode].plot(percentiles, np.abs(G_LargeNc),color = colors[2], label="CGC large-$N_c$")
		axs[l-1][mode].plot(percentiles, np.abs(G_magma),color = colors[3], label="magma")
		axs[l-1][mode].set_xticks([0,20,40,60,80])
		if (l==1):
			axs[l-1][mode].set_ylim(ymin=-0.02, ymax=0.10)
		elif (l==2):
			axs[l-1][mode].set_ylim(ymin=-0.005, ymax=0.05)
		else:
			axs[l-1][mode].set_ylim(ymin=-0.005, ymax=0.05)
		#axs[l-1][mode].set_xticks([10,30,50,70])
		#axs[l-1][mode].annotate(" \n$|G_"+str(l)+"^{("+str(mode)+")}|$",xy=(5,0.05))
		#axs[l-1][mode].set_yscale("log")
		if (l == radials[-1]):
			axs[l-1][mode].set_xlabel("Centrality %")
		#plt.title("$G_"+str(l)+"^{("+str(mode)+")}$")


#plt.legend(loc="best")


# annotate
for l in radials:
	for mode in modes:
		ymin, ymax = axs[l-1][mode].get_ylim()
		xmin, xmax = axs[l-1][mode].get_xlim()
		axs[l-1][mode].annotate(" \n$\\left|G_"+str(l)+"^{("+str(mode)+")}\\right|$",xy=(xmax*0.06,ymax*0.8))


handles, labels = axs[-1][-1].get_legend_handles_labels()
figs.legend(handles, labels, loc='upper center', ncol=4)
plt.savefig("plots/G_centrality.pdf", format="pdf", bbox_inches="tight")
		










