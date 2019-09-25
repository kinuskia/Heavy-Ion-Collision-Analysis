import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

percentiles = np.array([0, 20])
N_values = np.array([223, 436])
s2_N_values = np.zeros(len(percentiles))
s2_w_values = np.array([0.633, 2.18])

for pp in range(0, len(percentiles)):
	p = percentiles[pp]
	modes = [0, 1, 2, 3, 4]
	N = N_values[pp]
	s2_N = s2_N_values[pp]
	s2_w = s2_w_values[pp]

	centrality_class =  str(p) + '-' + str(p+1)

	# Set up figure and image grid
	fig = plt.figure(figsize=(10, 10))
	#plt.rcParams.update({'font.size': 30})
	#plt.rcParams['axes.titlepad'] = 10

	grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
	                 nrows_ncols=(len(modes),1),
	                 axes_pad=0.3,
	                 share_all=True,
	                 cbar_location="bottom",
	                 cbar_mode="edge",
	                 cbar_size="5%",
	                 cbar_pad=0.5,
	                 direction="column"
	                 )

	# Add data to image grid
	# for ax in grid:
	#     im = ax.imshow(np.random.random((6,6)), vmin=0, vmax=0.1)



	from matplotlib.ticker import MaxNLocator
	for ax in grid:
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		ax.yaxis.set_major_locator(MaxNLocator(integer=True))

	counter_fig = 0

	clm = np.loadtxt("output/clm.txt")
	lMax = 10

	trento_gauge = 1

	for mode in modes:
		ax = grid[counter_fig]



		#import one-mode expectation values
		# source_one_point = 'output/one_point_20-21.txt'
		source_one_point = 'output/one_point_' + centrality_class + '.txt'
		one_points = np.loadtxt(source_one_point)
		l = np.zeros(lMax)
		y = np.zeros(lMax)
		y[0] = 1.0/np.pi/np.pi*(s2_w/N+s2_N/N/N)
		for i in range(0, lMax):
			l[i] = i+1
			if mode == 0:
				if i == 0:
					continue
				else:
					y[i] = 1.0/2./np.pi**2/N/clm[i-1, mode]*(1.+s2_w) 
			else:
				y[i] = (-1.0)**mode/2./np.pi**2/N/clm[i, mode]*(1.+s2_w) + (1-1./N+s2_N/N/N)*one_points[mode, i]*one_points[mode, i]

			# # print useful output
			# if mode == 0:
			# 	if i == 0:
			# 		continue
			# 	else:
			# 		print(y[i]*2*np.pi*np.pi*clm[i-1,mode])
					
		ax.scatter(l, y, label="IPSM", s=100, color = "orangered", marker= "+")
		#plt.xlabel("l")
		#plt.ylabel("$G_l^{(" + str(mode)  + ")}$")
		#plt.title("m = " +str(mode))
		# Trento prediction

		filename_trento = "output/two_point_random_" + centrality_class + "_m" + str(mode) + "_real.txt"
		filename_trento_error = "output/two_point_random_" + centrality_class + "_m" + str(mode) + "_real_error.txt"
		filename_trento_one = 'output/one_point_' + centrality_class +'.txt'

		filename_Saclay_simple = "Saclay_simplified/output/"+centrality_class+"/two_point_random_connected_m_" + str(mode) + ".txt"
		filename_Saclay = "Saclay/output/"+centrality_class+"/Nr41/Nm64/m1e-3/two_point_random_connected_m_" + str(mode) + ".txt"

		trento = np.loadtxt(filename_trento)
		trento_error = np.loadtxt(filename_trento_error)
		trento_one_ml = np.loadtxt(filename_trento_one)
		saclay_simple = np.loadtxt(filename_Saclay_simple)
		saclay = np.loadtxt(filename_Saclay)
		# add geometry part
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
		one_points = np.loadtxt(source_one_point)
		for i in range(0, lMax):
			for j in range(0, lMax):
				if ((i==0)&(j==0)&(mode==0)):
					saclay_simple[i][j] += 1./np.pi/np.pi*(s2_N/N/N-1./N)
					saclay[i][j] += 1./np.pi/np.pi*(s2_N/N/N-1./N)
				else:
					saclay_simple[i][j] += (1.-1./N+s2_N/N/N)*one_points[mode, i]*one_points[mode, j]
					saclay[i][j] += (1.-1./N+s2_N/N/N)*one_points[mode, i]*one_points[mode, j]

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
		ax.axhline(y=0, color ="black", linestyle="dashed", linewidth=1.0, zorder = 1)
		ax.scatter(l, y_trento, s=100, color = "green", label="TRENTo", marker= "x", zorder = 2)
		ax.scatter(l, y_saclay_simple, s=100, color = "blue", label="magma", marker= "1", zorder = 2)
		ax.scatter(l, y_saclay, s=100, color = "purple", label="Large-Nc-Glasma", marker= "*", zorder = 2)













		# im = ax.imshow(profile[0:(lMax),0:(lMax)], interpolation=None, cmap=plt.cm.seismic, vmin = -maximal_T, vmax = maximal_T, extent = (-0.5+1, len(profile[0,0:(lMax)])-0.5+1, len(profile[0:(lMax),0])-0.5+1, -0.5+1))
		# #centrality_class =  str(p) + '-' + str(p+1) + '%'
		# if (mode == 0):
		# 	ax.set_title("TRENTo")
		# ax.set_xlabel("$l_2$")
		# ax.set_ylabel("$m={0}$\n$l_1$".format(mode))

		# # colorbar
		# if (mode == 0):
		# 	short = grid[counter_fig].cax
		# 	short.colorbar(im)
		# 	tick_T = float('%.1g' % (maximal_T*0.6))
		# 	short.set_xticks([-tick_T, tick_T])
		# 	short.toggle_label(True)

		counter_fig = counter_fig + 1





	fig.suptitle("$G_{l}^{(m)}$, "+centrality_class+"% class", x=0.5, y=0.94, fontsize=14)
	#plt.tight_layout()   
	#fig.subplots_adjust(left=0.15, top=0.95)
	filename = "plots/two_point_comparison_diagonals_"+centrality_class+".pdf"

	plt.savefig(filename, format='pdf', bbox_inches = "tight")




