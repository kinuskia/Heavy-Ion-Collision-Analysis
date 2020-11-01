import matplotlib.pyplot as plt
import numpy as np 



#modes = [0, 1, 2, 3, 4, 5]
#modes = [0, 1]
modes = [2, 4]
counter_fig = 0

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']


for i in range(0, len(modes)):
	m = modes[i]
	p = 20
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.figure(figsize=(10,7.3))
	plt.rcParams.update({'font.size': 23})
	plt.rcParams['axes.titlepad'] = 10
	#import TRento
	source_result = 'output/one_point_' + str(p) + '-' + str(p+1)  +'.txt'
	profile = np.loadtxt(source_result)
	source_error = 'output/one_point_' + str(p) + '-' + str(p+1) + '_error' +'.txt'
	profile_error = np.loadtxt(source_error)
	#import CGC
	source_magma = 'Saclay_simplified/output/' + str(p) + '-' + str(p+1) +'/background_coeffs_CGC_' + str(p) + '-' + str(p+1)  +'.txt'
	profile_magma = np.loadtxt(source_magma)
	mMax = 4
	lMax = 10
	l = np.zeros(lMax)
	for i in range(0, lMax):
		l[i] = i+1
	
	plt.errorbar(l, profile[m, 0:lMax], yerr=profile_error[m,0:], linestyle = 'none', elinewidth=2, capsize = 6, capthick = 2, color=colors[1])
	plt.scatter(l, profile[m, 0:lMax], s=30, color=colors[1], label = "TrENTo")
	#plt.plot(l, profile[m, 0:], color=colors[m])
	plt.scatter(l, profile_magma[m, 0:lMax], s=100, color=colors[2], label = "CGC", marker= "*")
	#plt.plot(l, profile_magma[m, 0:], color=colors[1])
	plt.xlabel("$l$")
	plt.ylabel("$\\bar{\\epsilon}_l^{("+str(m)+")}$")
	plt.legend(loc='best')
	centrality_class = str(p) + '-' + str(p+1) + '%'
	plt.title("$\\bar{\\epsilon}_l^{("+str(m)+")}$, "+centrality_class)
	filename = "plots/one_point_comparison_m_"+str(m)+".pdf"
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close(counter_fig)




















# import numpy as np
# import matplotlib.pyplot as plt
# from matplotlib import ticker

# percentiles = [0, 20]

# # Set up figure and image grid
# figs, axs = plt.subplots(2, 2, figsize=(8.0,6.0),sharex="col", sharey="row", constrained_layout=True)




# cmaps = [plt.cm.RdBu,plt.cm.PRGn]
# maximal = [0.02, .2]
# lMax = 10
# mMax = 5

# ims = []

# from matplotlib.ticker import MaxNLocator


# for i in range(0, len(percentiles)):
# 	# plot Trento diagrams
# 	p = percentiles[i]

# 	Trento_one = 'output/one_point_' + str(p) + '-' + str(p+1)  +'.txt'
# 	Trento_one_ml = np.loadtxt(Trento_one)

# 	#im = axs[i][0].imshow(Trento_one_ml[2:mMax:2,0:lMax], interpolation=None, cmap=cmaps[i], vmin = -maximal[i], vmax = maximal[i], extent = (0.5, len(Trento_one_ml[0,0:lMax])-0.5+1, 4, 0))
# 	im = axs[i][0].imshow(Trento_one_ml[2:mMax,0:lMax], interpolation=None, cmap=cmaps[i], vmin = -maximal[i], vmax = maximal[i], extent = (-0.5+1, len(Trento_one_ml[0,0:lMax])+1-0.5, len(Trento_one_ml[2:mMax,0])-0.5, -0.5))
# 	ims.append(im)
# 	# plot CGC diagrams
# 	magma_one = 'Saclay_simplified/output/'+ str(p) + '-' + str(p+1)  +'/background_coeffs_CGC_' + str(p) + '-' + str(p+1)  +'.txt'
# 	magma_one_ml = np.loadtxt(magma_one)

# 	im2 = axs[i][1].imshow(magma_one_ml[2:mMax,0:lMax], interpolation=None, cmap=cmaps[i], vmin = -maximal[i], vmax = maximal[i], extent = (-0.5+1, len(magma_one_ml[0,0:lMax])+1-0.5, len(magma_one_ml[2:mMax,0])-0.5, -0.5))
		
# 	# axs[0][0].set_title("TrENTo")
# 	# axs[1][0].set_xlabel("$l_2$")
# 	# axs[0][0].set_xticks(np.arange(lMax)+1)
# 	#axs[0][mode].set_yticks(np.arange(lMax)+1)






# cb1 = figs.colorbar(ims[0], ax = axs[0,1], shrink = 1, pad=0.0, aspect = 40, location = "right")
# tick_locator = ticker.MaxNLocator(nbins=5)
# cb1.locator = tick_locator
# cb1.update_ticks()
# cb2 = figs.colorbar(ims[1], ax = axs[1,1], shrink = 1, pad=0.0, aspect = 40, location = "right")
# tick_locator = ticker.MaxNLocator(nbins=5)
# cb2.locator = tick_locator
# cb2.update_ticks()
# # axs[0][0].set_ylabel("TrENTo\n$l_1$")
# # axs[1][0].set_ylabel("CGC large-$N_c$\n$l_1$")
# # axs[2][0].set_ylabel("IPSM\n$l_1$")
# # axs[3][0].set_ylabel("magma\n$l_1$")


# figs.suptitle("$\\bar{\\epsilon}_l^{(m)}$", x=0.5, y=1.05, fontsize=14)

# filename = "plots/two_point_comparison_OnePoint.pdf"
# #plt.subplots_adjust(wspace=0.05, hspace=0.05)
# plt.savefig(filename, format='pdf', bbox_inches = "tight")





		










