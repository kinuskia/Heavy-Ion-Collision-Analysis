import numpy as np
import matplotlib.pyplot as plt

modes = [0, 1, 2, 3, 4]
centrality = np.arange(95)
lMax = 5

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']

# Set up figure and image grid
figs, axs = plt.subplots(1, len(modes), figsize=(8.0,1.9),sharex="col", sharey="row", constrained_layout=False)


for mm in range(0, len(modes)):
	mode = modes[mm]
	chi_l = np.zeros((lMax, len(centrality)))
	for p in range(0, len(centrality)):
		###### Read in IPSM Fitting result
		centrality_class =  str(centrality[p]) + '-' + str(centrality[p]+1)
		filename_fit = "../IPSM-Fit/output/" + centrality_class + ".txt"
		One_sw2_N, sn2_N_N2 = np.loadtxt(filename_fit, unpack=True)

		###### Initialize chi profile
		chi = np.zeros((lMax, lMax))

		###### Get Trento data
		source = 'output/two_point_random_' + str(p) + '-' + str(p+1) + '_m' + str(mode) + '_real' +'.txt'
		profile_trento = np.loadtxt(source)
		# import one-point functions
		source_one = 'output/one_point_' + str(p) + '-' + str(p+1)  +'.txt'
		profile_one_ml = np.loadtxt(source_one)
		# subtract disconnected part from m=0 mode
		if (mode == 0):
			profile_trento[0,0] -= profile_one_ml[mode, 0] * profile_one_ml[mode, 0]


		###### Get Large Nc
		source = 'Saclay/output/'+centrality_class+'/Nr41/Nm64/m5.0e-3/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
		profile_largeNc = np.loadtxt(source)
		# # add geometry part
		# source_one_point = 'output/one_point_'+centrality_class+'.txt'
		# one_points = np.loadtxt(source_one_point)
		# for i in range(0, lMax):
		# 	for j in range(0, lMax):
		# 		if ((i==0)&(j==0)&(mode==0)):
		# 			profile_largeNc[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
		# 		else:
		# 			profile_largeNc[i][j] += (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]


		###### Get IPSM
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
		one_points = np.loadtxt(source_one_point)
		profile_IPSM = np.zeros((lMax, lMax))

		clm_old = np.loadtxt("output/clm.txt")
		clm = np.zeros((lMax, modes[-1]+1))
		for l in range(0, lMax):
			for m in range(0, modes[-1]+1):
				if ((m==0)):
					if (l==0):
						clm[l][m] = 1./2
					else:
						clm[l][m] = clm_old[l-1][m]
				else:
					clm[l][m] = clm_old[l][m]

		for i in range(0, lMax):
			for j in range(0, lMax):
				if ((i == 0)&(j==0)&(mode == 0)):
					profile_IPSM[i][j] = 1./np.pi/np.pi*(One_sw2_N+sn2_N_N2)
				elif ((i == j)):
					profile_IPSM[i][j] = (-1.0)**mode/2./np.pi**2/clm[i, mode]*One_sw2_N + (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]
				else:
					profile_IPSM[i][j] = (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]


		###### Get magma
		source = 'Saclay_simplified/output/'+centrality_class+'/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
		profile_magma = np.loadtxt(source)
		# # add geometry part
		# source_one_point = 'output/one_point_'+centrality_class+'.txt'
		# one_points = np.loadtxt(source_one_point)
		# for i in range(0, lMax):
		# 	for j in range(0, lMax):
		# 		if ((i==0)&(j==0)&(mode==0)):
		# 			profile_magma[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
		# 		else:
		# 			profile_magma[i][j] += (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]


		###### Compute chi entries
		for i in range (0, lMax):
			for j in range (0, lMax):
				p1 = profile_trento[i][j]
				p2 = profile_largeNc[i][j]
				p3 = profile_IPSM[i][j]
				p4 = profile_magma[i][j]
				pmean = (p1+p2+p3+p4)/4
				# chi[i][j] = (pow(p1-p2,2) + pow(p1-p3,2) + pow(p1-p4,2) + pow(p2-p3,2) + pow(p2-p4,2) + pow(p3-p4,2))/pmean/pmean
				#chi[i][j] = (pow(p1-pmean,2) + pow(p2-pmean,2) + pow(p3-pmean,2) + pow(p4-pmean,2))/pmean/pmean
				chi[i][j] = (abs(p1-p3) + abs(p2-p3) + abs(p4-p3))/abs(p3)
		#chi = np.sqrt(chi/6)
		#chi = np.sqrt(chi/4)
		chi = chi/3.

		for i in range(0, lMax):
			chi_l[i][p] = chi[i,i]


	####### Plot chi
	#for i in range(lMax-1, -1, -1):
	for i in range(0, lMax):
		if ((i == 0)&(mode==0)):
			continue
		else:
			chi_averaged = chi_l[i]
			# for j in range(1, len(centrality)-1):
			# 	chi_averaged[j] = (chi_l[i][j+1] + chi_l[i][j] + chi_l[i][j-1])/3
			axs[mm].plot(centrality, chi_averaged, label = "$l={0}$".format(i+1), color=colors[i])
	#axs[mm].legend(loc="best")
	axs[mm].set_xlabel("Centrality %")
	axs[mm].set_title("$m={0}$".format(mode))
	axs[mm].set_xticks([0,20,40,60,80])

handles, labels = axs[-1].get_legend_handles_labels()
figs.legend(handles, labels,loc='lower center', bbox_to_anchor=(0.45, -0.04), ncol=5, frameon=False)
figs.subplots_adjust(bottom=0.35)
figs.suptitle("Mean relative deviation of $G_l^{(m)}$ from IPSM", x=0.5, y=1.22, fontsize=14)
filename = "plots/two_point_comparison_mean-diff.pdf"
plt.savefig(filename, format='pdf', bbox_inches = "tight")




