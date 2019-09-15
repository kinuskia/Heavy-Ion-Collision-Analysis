import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid



modes = [0, 1, 2, 3, 4]
gridpoints = [8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41]
percentiles = [20]
Nm = [128]


# Set up figure and image grid
fig = plt.figure(figsize=(25, 10))
#plt.rcParams.update({'font.size': 30})
#plt.rcParams['axes.titlepad'] = 10

grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(len(modes),len(gridpoints)),
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


tick = 0.005*np.ones(len(gridpoints))


centrality_class =  str(percentiles[0]) + '-' + str(percentiles[0]+1)
lMax = 10

# Find maximal correlator value
maximal_value = np.zeros(len(gridpoints))

for mode in modes:
	for p in percentiles:
		for i in range(0, len(gridpoints)):
			for j in range(0, len(Nm)):
				#import profiles
				source = 'Saclay/output/'+centrality_class+'/Nr'+ str(gridpoints[i]) + '/Nm' + str(Nm[j]) +'/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
				profile = np.loadtxt(source)
				current_max = max(np.amax(profile), -np.amin(profile))
				if (current_max > maximal_value[i]):
					maximal_value[i] = current_max
				print(i, maximal_value[i])
		


from matplotlib.ticker import MaxNLocator
for ax in grid:
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))

counter_fig = 0

# plot diagrams
for i in range(0, len(gridpoints)):
	for j in range(0, len(Nm)):
		for mode in modes:
			for p in percentiles:
				#import two-point functions
				source = 'Saclay/output/'+centrality_class+'/Nr'+ str(gridpoints[i]) + '/Nm' + str(Nm[j]) +'/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
				profile = np.loadtxt(source)

				ax = grid[counter_fig]
				im = ax.imshow(profile[0:(lMax),0:(lMax)], interpolation=None, cmap=plt.cm.seismic, vmin = -maximal_value[i], vmax = maximal_value[i], extent = (-0.5+1, len(profile[0,0:(lMax)])-0.5+1, len(profile[0:(lMax),0])-0.5+1, -0.5+1))
				#centrality_class =  str(p) + '-' + str(p+1) + '%'
				if (mode == 0):
					ax.set_title(str(gridpoints[i]))
				ax.set_xlabel("$l_2$")
				ax.set_ylabel("$m={0}$\n$l_1$".format(mode))

				# colorbar
				if (mode == modes[-1]):
					short = grid[counter_fig].cax
					short.colorbar(im)
					short.set_xticks([-tick[i], tick[i]])
					short.toggle_label(True)

				counter_fig = counter_fig + 1







fig.suptitle("$G_{l_1, l_2}^{(m,-m)}$, "+centrality_class+"% class", x=0.5, y=0.94, fontsize=14)
#plt.tight_layout()   
#fig.subplots_adjust(left=0.15, top=0.95)
filename = "plots/two_point_comparison_"+centrality_class+"_gridpoints.pdf"

plt.savefig(filename, format='pdf', bbox_inches = "tight")




