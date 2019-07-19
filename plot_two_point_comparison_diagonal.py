import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid

# Set up figure and image grid
fig = plt.figure(figsize=(10, 10))
#plt.rcParams.update({'font.size': 30})
#plt.rcParams['axes.titlepad'] = 10

grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(5,2),
                 axes_pad=0.3,
                 share_all=True,
                 cbar_location="bottom",
                 cbar_mode=None,
                 cbar_size="5%",
                 cbar_pad=0.5,
                 direction="column"
                 )

# Add data to image grid
# for ax in grid:
#     im = ax.imshow(np.random.random((6,6)), vmin=0, vmax=0.1)
modes = [0, 1, 2, 3, 4]
percentiles = [0, 20]
N = 172*0+223
s2_N = 0.*1.7**2+0
s2_w = 0.*0.48**2+0.633**2
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
centrality_class =  str(percentiles[0]) + '-' + str(percentiles[0]+1)
clm = np.loadtxt("output/clm.txt")
lMax = 6



from matplotlib.ticker import MaxNLocator
for ax in grid:
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))

counter_fig = 0
# plot Trento diagrams
for p in percentiles:
	for mode in modes:
		centrality_class =  str(p) + '-' + str(p+1)
		ax = grid[counter_fig]
		ax.set_figsize=(5.0,2.0)
		if (p == 0):
			N = 172*0+223
			s2_N = 0.*1.7**2+0
			s2_w = 0.*0.48**2+0.633**2
		elif (p == 20):
			N = 112*0+436
			s2_N = 0*3.1**2+0
			s2_w = 0*0.85**2+2.18**2
		source_one_point = 'output/one_point_'+centrality_class+'.txt'
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


				
		ax.scatter(l, y, label="IPSM", s=100, color = "orangered", marker= "+")
		ax.set_xlabel("l")
		ax.set_ylabel("$G_l^{(" + str(mode)  + ")}$")
		#plt.title("m = " +str(mode))
		# Trento prediction
		# filename_trento = "output/two_point_random_20-21_m" + str(mode) + "_real.txt"
		# filename_trento_error = "output/two_point_random_20-21_m" + str(mode) + "_real_error.txt"
		# filename_trento_one = 'output/one_point_20-21' +'.txt'
		filename_trento = "output/two_point_random_"+centrality_class+"_m" + str(mode) + "_real.txt"
		filename_trento_error = "output/two_point_random_"+centrality_class+"_m" + str(mode) + "_real_error.txt"
		filename_trento_one = 'output/one_point_'+centrality_class+'.txt'

		filename_Saclay_simple = "Saclay_simplified/output/"+centrality_class+"/two_point_random_connected_m_" + str(mode) + ".txt"
		filename_Saclay = "Saclay/output/"+centrality_class+"/two_point_random_connected_m_" + str(mode) + ".txt"

		trento = np.loadtxt(filename_trento)
		trento_error = np.loadtxt(filename_trento_error)
		trento_one_ml = np.loadtxt(filename_trento_one)
		saclay_simple = np.loadtxt(filename_Saclay_simple)
		saclay = np.loadtxt(filename_Saclay)
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
		ax.scatter(l, y_trento, s=100, color = "green", label="TRENTo", marker= "x")
		ax.scatter(l, y_saclay_simple, s=100, color = "blue", label="Generic", marker= "1")
		ax.scatter(l, y_saclay, s=100, color = "purple", label="Large-Nc-Glasma", marker= "*")


		counter_fig = counter_fig + 1

fig.suptitle("$G_{l, l}^{(m,-m)}$", x=0.5, y=0.94, fontsize=14)
#plt.tight_layout()   
#fig.subplots_adjust(left=0.15, top=0.95)
filename = "plots/two_point_comparison_diagonal.pdf"

plt.savefig(filename, format='pdf', bbox_inches = "tight")




