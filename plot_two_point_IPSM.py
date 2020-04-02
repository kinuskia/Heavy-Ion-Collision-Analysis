import matplotlib.pyplot as plt
import numpy as np 

#percentiles = np.arange(21)
percentiles = np.array([0, 20])
modes = [0,1,2,3,4]
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

figs, axs = plt.subplots(len(modes), len(percentiles), figsize=(8.0,7.5),sharex="col", sharey="row")
for pp in range(0, len(percentiles)):
	p = percentiles[pp]


	
	#modes = [0, 1]

	centrality_class =  str(p) + '-' + str(p+1)



	clm = np.loadtxt("output/clm.txt")
	lMax = 7

	trento_gauge = 1

	

	# modulus plots
	for mode in modes:
		#plt.figure(counter_fig)
		#plt.rcParams.update({'font.size': 20})
		#from matplotlib.ticker import MaxNLocator
		#ax = plt.figure(counter_fig).gca()
		#ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		#ax.yaxis.set_major_locator(MaxNLocator(integer=True))
		#plt.figure(figsize=(5.0,1.5))

		filename_trento = "output/two_point_random_" + centrality_class + "_m" + str(mode) + "_real.txt"
		filename_trento_error = "output/two_point_random_" + centrality_class + "_m" + str(mode) + "_real_error.txt"
		filename_trento_one = 'output/one_point_' + centrality_class +'.txt'

		filename_Saclay_simple = "Saclay_simplified/output/"+centrality_class+"/two_point_random_connected_m_" + str(mode) + ".txt"
		filename_Saclay = "Saclay/output/"+centrality_class+"/Nr41/Nm64/m5.0e-3/two_point_random_connected_m_" + str(mode) + ".txt"

		trento = np.loadtxt(filename_trento)
		trento_error = np.loadtxt(filename_trento_error)
		trento_one_ml = np.loadtxt(filename_trento_one)
		saclay_simple = np.loadtxt(filename_Saclay_simple)
		saclay = np.loadtxt(filename_Saclay)




		# plot IPSM

		filename_fit = "../IPSM-Fit/output/" + centrality_class + ".txt"
		One_sw2_N, sn2_N_N2 = np.loadtxt(filename_fit, unpack=True)


		source_one_point = 'output/one_point_' + centrality_class + '.txt'
		one_points = np.loadtxt(source_one_point)
		l = np.zeros(lMax)
		y = np.zeros(lMax)
		y[0] = 1.0/np.pi/np.pi*(One_sw2_N +sn2_N_N2)
		for i in range(0, lMax):
			l[i] = i+1
			if mode == 0:
				if i == 0:
					continue
				else:
					y[i] = 1.0/2./np.pi**2/clm[i-1, mode]*One_sw2_N
			else:
				y[i] = (-1.0)**mode/2./np.pi**2/clm[i, mode]*One_sw2_N + (1+sn2_N_N2)*one_points[mode, i]*one_points[mode, i]

					
		axs[mode][pp].scatter(l, y, label="IPSM", s=100, color = colors[0], marker= "+")
		if (mode == modes[-1]):
			axs[mode][pp].set_xlabel("l")
		if (pp == percentiles[0]):
			axs[mode][pp].set_ylabel("$G_l^{(" + str(mode)  + ")}$")

		if (mode == modes[0]):
			axs[mode][pp].set_title(centrality_class+"%")




		# # add geometry part
		# source_one_point = 'output/one_point_'+centrality_class+'.txt'
		# one_points = np.loadtxt(source_one_point)
		# for i in range(0, lMax):
		# 	for j in range(0, lMax):
		# 		if ((i==0)&(j==0)&(mode==0)):
		# 			saclay_simple[i][j] += 1./np.pi/np.pi*sn2_N_N2
		# 			saclay[i][j] += 1./np.pi/np.pi*sn2_N_N2
		# 		else:
		# 			saclay_simple[i][j] += (1.+ sn2_N_N2 )*one_points[mode, i]*one_points[mode, j]
		# 			saclay[i][j] += (1. + sn2_N_N2 )*one_points[mode, i]*one_points[mode, j]

		y_saclay_simple = np.zeros(lMax)
		y_saclay = np.zeros(lMax)
		y_trento = np.zeros(lMax)
		y_trento_error = np.zeros(lMax)
		for i in range(0, lMax):
			y_trento[i] = trento[i,i]
			y_saclay_simple[i] = saclay_simple[i,i]
			y_saclay[i] = saclay[i,i]
			if ((i == 0) & (mode == 0)):
				y_trento[i] -= trento_one_ml[mode, i] * trento_one_ml[mode, i]
			y_trento_error[i] =abs(trento_error[i,i]) # WHAT ABOUT m=0 ?
		axs[mode][pp].axhline(y=0, color ="black", linestyle="dashed", linewidth=1.0, zorder = 1)
		axs[mode][pp].scatter(l, y_trento, s=100, color = colors[1], label="TrENTo", marker= "x", zorder = 3)
		axs[mode][pp].scatter(l, y_saclay_simple, s=100, color = colors[3], label="magma", marker = "2", zorder = 3)
		axs[mode][pp].scatter(l, y_saclay, s=100, color = colors[2], label="CGC large-$N_c$", marker= "*", zorder = 2)
		axs[mode][pp].set_xticks([1,3,5,7])
		if (mode%2 == 0):
			axs[mode][pp].set_ylim(-0.2*y[-1], 1.2*y[-1])
		else:
			axs[mode][pp].set_ylim(1.2*y[-1], -0.2*y[-1])

filename = "plots/two_point_modules_IPSM"   + ".pdf"
handles, labels = axs[-1][-1].get_legend_handles_labels()
figs.legend(handles, labels, loc='lower center',ncol=4)
plt.savefig(filename, format='pdf', bbox_inches = "tight")










