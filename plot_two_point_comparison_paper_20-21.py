import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

#percentiles = np.arange(100)
percentiles = np.array([20])
# Read in IPSM fit result
One_sw2_N_values = np.zeros(len(percentiles))
sn2_N_N2_values = np.zeros(len(percentiles))

for k in range(0, len(percentiles)):
	centrality_class =  str(percentiles[k]) + '-' + str(percentiles[k]+1)
	filename_fit = "../IPSM-Fit/output/" + centrality_class + ".txt"
	One_sw2_N, sn2_N_N2 = np.loadtxt(filename_fit, unpack=True)
	One_sw2_N_values[k] = One_sw2_N
	sn2_N_N2_values[k] = sn2_N_N2

for pp in range(0, len(percentiles)):
	p = percentiles[pp]

	# Set up figure and image grid
	fig = plt.figure(figsize=(10, 10))
	#plt.rcParams.update({'font.size': 30})
	#plt.rcParams['axes.titlepad'] = 10

	grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
	                 nrows_ncols=(5,4),
	                 axes_pad=0.3,
	                 share_all=True,
	                 cbar_location="bottom",
	                 cbar_mode="single",
	                 cbar_size="2%",
	                 cbar_pad=None,
	                 direction="column"
	                 )

	# Add data to image grid
	# for ax in grid:
	#     im = ax.imshow(np.random.random((6,6)), vmin=0, vmax=0.1)
	modes = [0, 1, 2, 3, 4]
	One_sw2_N = One_sw2_N_values[pp]
	sn2_N_N2 = sn2_N_N2_values[pp]
	tick_T = 0.005
	tick_N = 0.005
	tick_I = 0.02
	tick_C = 0.02
	# percentiles = [20]
	# N = 112*0+436
	# s2_N = 0*3.1**2+0
	# s2_w = 0*0.85**2+2.18**2
	# tick_T = 0.005
	# tick_N = 0.005
	# tick_I = 0.02
	# tick_C = 0.02
	centrality_class =  str(p) + '-' + str(p+1)
	lMax = 6

	# Find maximal correlator value
	maximal_value_Trento = 0
	maximal_value_CGC = 0
	maximal_value_Nc = 0
	maximal_value_IPSM = 0
	for mode in modes:
		#import two-point functions
		source = 'output/two_point_random_' + str(p) + '-' + str(p+1) + '_m' + str(mode) + '_real' +'.txt'
		profile = np.loadtxt(source)
		# import one-point functions
		source_one = 'output/one_point_' + str(p) + '-' + str(p+1)  +'.txt'
		profile_one_ml = np.loadtxt(source_one)
		# subtract disconnected part from m=0 mode
		if (mode == 0):
			profile[0,0] -= profile_one_ml[mode, 0] * profile_one_ml[mode, 0]

		current_max_Trento = max(np.amax(profile[:lMax,:lMax]), -np.amin(profile[:lMax,:lMax]))
		if (current_max_Trento > maximal_value_Trento):
			maximal_value_Trento = current_max_Trento
		# import CGC simple profiles
		source = 'Saclay_simplified/output/'+centrality_class+'/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
		profile = np.loadtxt(source)
		# add geometry part
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
		one_points = np.loadtxt(source_one_point)
		for i in range(0, lMax):
			for j in range(0, lMax):
				if ((i==0)&(j==0)&(mode==0)):
					profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
				else:
					profile[i][j] += (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]
						
		current_max_CGC = max(np.amax(profile[:lMax,:lMax]), -np.amin(profile[:lMax,:lMax]))
		if (current_max_CGC > maximal_value_CGC):
			maximal_value_CGC = current_max_CGC
		#import Large_Nc profile
		source = 'Saclay/output/'+centrality_class+'/Nr41/Nm64/m1.4e-1/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
		profile = np.loadtxt(source)
		# add geometry part
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
		one_points = np.loadtxt(source_one_point)
		for i in range(0, lMax):
			for j in range(0, lMax):
				if ((i==0)&(j==0)&(mode==0)):
					profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
				else:
					profile[i][j] += (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]

		current_max_Nc = max(np.amax(profile[:lMax,:lMax]), -np.amin(profile[:lMax,:lMax]))
		if (current_max_Nc > maximal_value_Nc):
			maximal_value_Nc = current_max_Nc
		#import IPSM profile
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
		one_points = np.loadtxt(source_one_point)
		profile = np.zeros((lMax, lMax))
		clm_old = np.loadtxt("output/clm.txt")
		clm = np.zeros((lMax, len(modes)))
		for l in range(0, lMax):
			for m in range(0, len(modes)):
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
					profile[i][j] = 1./np.pi/np.pi*(One_sw2_N+sn2_N_N2)
				elif ((i == j)):
					profile[i][j] = (-1.0)**mode/2./np.pi**2/clm[i, mode]*One_sw2_N + (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]
				else:
					profile[i][j] = (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]
		current_max_IPSM = max(np.amax(profile[:lMax,:lMax]), -np.amin(profile[:lMax,:lMax]))
		if (current_max_IPSM > maximal_value_IPSM):
			maximal_value_IPSM = current_max_IPSM
	print("Trento: ", maximal_value_Trento)
	print("CGC: ", maximal_value_CGC)
	print("Nc: ", maximal_value_Nc)
	print("IPSM: ", maximal_value_IPSM)

	maximal_TN = max(maximal_value_Trento, maximal_value_Nc)
	maximal_IC = max(maximal_value_CGC, maximal_value_IPSM)
	maximal_total = max(maximal_TN, maximal_IC)
	maximal_T = maximal_total
	maximal_N = maximal_total
	maximal_I = maximal_total
	maximal_C = maximal_total

	# maximal_T = maximal_value_Trento
	# maximal_N = maximal_value_Nc
	# maximal_I = maximal_value_IPSM
	# maximal_C = maximal_value_CGC

	cmap_Trento = plt.cm.seismic
	cmap_LargeNc = plt.cm.seismic
	cmap_IPSM = plt.cm.seismic
	cmap_magma = plt.cm.seismic


	from matplotlib.ticker import MaxNLocator
	for ax in grid:
		ax.xaxis.set_major_locator(MaxNLocator(integer=True))
		ax.yaxis.set_major_locator(MaxNLocator(integer=True))

	counter_fig = 0
	# plot Trento diagrams
	for mode in modes:
		#import two-point functions
		source = 'output/two_point_random_' + str(p) + '-' + str(p+1) + '_m' + str(mode) + '_real' +'.txt'
		profile = np.loadtxt(source)
		# import one-point functions
		source_one = 'output/one_point_' + str(p) + '-' + str(p+1)  +'.txt'
		profile_one_ml = np.loadtxt(source_one)
		# subtract disconnected part from m=0 mode
		if (mode == 0):
			profile[0,0] -= profile_one_ml[mode, 0] * profile_one_ml[mode, 0]

		ax = grid[counter_fig]
		im = ax.imshow(profile[0:(lMax),0:(lMax)], interpolation=None, cmap=cmap_Trento, vmin = -maximal_T, vmax = maximal_T, extent = (-0.5+1, len(profile[0,0:(lMax)])-0.5+1, len(profile[0:(lMax),0])-0.5+1, -0.5+1))
		#centrality_class =  str(p) + '-' + str(p+1) + '%'
		if (mode == 0):
			ax.set_title("TRENTo")
		ax.set_xlabel("$l_2$")
		ax.set_ylabel("$m={0}$\n$l_1$".format(mode))

		# colorbar
		if (mode == 0):
			short = grid[counter_fig].cax
			short.colorbar(im)
			tick_T = float('%.1g' % (maximal_T*0.6))
			short.set_xticks([-tick_T, tick_T])
			short.toggle_label(True)

		counter_fig = counter_fig + 1

	# plot Saclay diagrams
	for mode in modes:
		#import two-point functions
		source = 'Saclay/output/'+centrality_class+'/Nr41/Nm64/m1.4e-1/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
		profile = np.loadtxt(source)
		# add geometry part
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
		one_points = np.loadtxt(source_one_point)
		for i in range(0, lMax):
			for j in range(0, lMax):
				if ((i==0)&(j==0)&(mode==0)):
					profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
				else:
					profile[i][j] += (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]

		ax = grid[counter_fig]
		im = ax.imshow(profile[0:(lMax),0:(lMax)], interpolation=None, cmap=cmap_LargeNc, vmin = -maximal_N, vmax = maximal_N, extent = (-0.5+1, len(profile[0,0:(lMax)])-0.5+1, len(profile[0:(lMax),0])-0.5+1, -0.5+1))
		#centrality_class =  str(p) + '-' + str(p+1) + '%'
		if (mode == 0):
			ax.set_title("Large $N_c$")
		ax.set_xlabel("$l_2$")
		ax.set_ylabel("$m={0}$\n$l_1$".format(mode))

		# colorbar
		if (mode == modes[-1]):
			short = grid[counter_fig].cax
			short.colorbar(im)
			tick_N = float('%.1g' % (maximal_N*0.6))
			short.set_xticks([-tick_N, tick_N])
			short.toggle_label(True)

		counter_fig = counter_fig + 1


	# plot IPSM diagrams 
	for mode in modes:
		#import one-point functions
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
		one_points = np.loadtxt(source_one_point)
		profile = np.zeros((lMax, lMax))

		clm_old = np.loadtxt("output/clm.txt")
		clm = np.zeros((lMax, len(modes)))
		for l in range(0, lMax):
			for m in range(0, len(modes)):
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
					profile[i][j] = 1./np.pi/np.pi*(One_sw2_N+sn2_N_N2)
				elif ((i == j)):
					profile[i][j] = (-1.0)**mode/2./np.pi**2/clm[i, mode]*One_sw2_N + (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]
				else:
					profile[i][j] = (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]


		ax = grid[counter_fig]
		im = ax.imshow(profile[0:(lMax),0:(lMax)], interpolation=None, cmap=cmap_IPSM, vmin = -maximal_I, vmax = maximal_I, extent = (-0.5+1, len(profile[0,0:(lMax)])-0.5+1, len(profile[0:(lMax),0])-0.5+1, -0.5+1))
		#centrality_class =  str(p) + '-' + str(p+1) + '%'
		if (mode == 0):
			ax.set_title("IPSM")
		ax.set_xlabel("$l_2$")
		ax.set_ylabel("$m={0}$\n$l_1$".format(mode))

		# colorbar
		if (mode == modes[-1]):
			short = grid[counter_fig].cax
			short.colorbar(im)
			tick_I = float('%.1g' % (maximal_I*0.6))
			short.set_xticks([-tick_I, tick_I])
			short.toggle_label(True)

		counter_fig = counter_fig + 1


	# plot Saclay simple diagrams
	for mode in modes:
		#import two-point functions
		source = 'Saclay_simplified/output/'+centrality_class+'/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
		profile = np.loadtxt(source)
		# add geometry part
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
		one_points = np.loadtxt(source_one_point)
		for i in range(0, lMax):
			for j in range(0, lMax):
				if ((i==0)&(j==0)&(mode==0)):
					profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
				else:
					profile[i][j] += (1.+sn2_N_N2)*one_points[mode, i]*one_points[mode, j]

		ax = grid[counter_fig]
		im = ax.imshow(profile[0:(lMax),0:(lMax)], interpolation=None, cmap=cmap_magma, vmin = -maximal_C, vmax = maximal_C, extent = (-0.5+1, len(profile[0,0:(lMax)])-0.5+1, len(profile[0:(lMax),0])-0.5+1, -0.5+1))
		#centrality_class =  str(p) + '-' + str(p+1) + '%'
		if (mode == 0):
			ax.set_title("magma")
		ax.set_xlabel("$l_2$")
		ax.set_ylabel("$m={0}$\n$l_1$".format(mode))

		# colorbar
		if (mode == modes[-1]):
			short = grid[counter_fig].cax
			short.colorbar(im)
			tick_C = float('%.1g' % (maximal_C*0.6))
			short.set_xticks([-tick_C, tick_C])
			short.toggle_label(True)

		counter_fig = counter_fig + 1




	fig.suptitle("$G_{l_1, l_2}^{(m,-m)}$, "+centrality_class+"% class", x=0.5, y=0.94, fontsize=14)
	#plt.tight_layout()   
	#fig.subplots_adjust(left=0.15, top=0.95)
	filename = "plots/two_point_comparison_"+centrality_class+".pdf"

	plt.savefig(filename, format='pdf', bbox_inches = "tight")




