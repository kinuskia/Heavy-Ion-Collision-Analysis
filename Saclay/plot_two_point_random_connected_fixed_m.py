import matplotlib.pyplot as plt
import numpy as np 

modes = [0, 1, 2, 3, 4]
#modes = [0]
counter_fig = 0

lMax = 10
# connected diagrams
for mode in modes:
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.rcParams.update({'font.size': 30})
	plt.rcParams['axes.titlepad'] = 10
	from matplotlib.ticker import MaxNLocator
	ax = plt.figure().gca()
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.figure(figsize=(7,5))
	source = 'output/two_point_random_connected' + '_m_' + str(mode)  +'.txt'
	profile = np.loadtxt(source)

	
	maximal_value = max(np.amax(profile), -np.amin(profile))
	plt.imshow(profile[0:(lMax),0:(lMax)], interpolation=None, cmap=plt.cm.RdYlGn,vmin = -maximal_value, vmax = maximal_value, extent = (-0.5+1, len(profile[0,0:(lMax)])-0.5+1, len(profile[0:(lMax),0])-0.5+1, -0.5+1))
	plt.xlabel("$l_2$")
	plt.ylabel("$l_1$")
	plt.title("$m = $" + str(mode))
	plt.colorbar()
	filename = "plots/two_point_connected_random"  + "_m" + str(mode) + ".pdf"
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close(counter_fig)





