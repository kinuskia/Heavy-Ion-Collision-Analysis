import matplotlib.pyplot as plt
import numpy as np 


percentiles = [0, 5, 10, 20, 30, 40]
counter_fig = 0

# modulus plots
for p in range(1, len(percentiles)):
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.figure(figsize=(10,10))
	plt.rcParams.update({'font.size': 23})
	plt.rcParams['axes.titlepad'] = 20
	source = 'output/two_point_background_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  + '_real' +'.txt'
	profile = np.loadtxt(source)
	maximal_value = max(np.amax(profile), -np.amin(profile))
	plt.imshow(profile, interpolation=None, cmap=plt.cm.RdYlGn,vmin = -maximal_value, vmax = maximal_value, extent = (-0.5+1, len(profile[0,:])-0.5+1, len(profile[:,0])/2, -len(profile[:,0])/2))
	#, extent = (-0.5+1, len(profile[0,:])-0.5+1, len(profile[:,0])-0.5+1, -0.5+1)
	plt.xlabel("$l$")
	plt.ylabel("$m$")
	centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
	plt.title("$\\left\\langle\\epsilon^{(0)}_{0} \\epsilon^{(m)}_{l}\\right\\rangle$, "+centrality_class)
	plt.colorbar()
	filename = "plots/two_point_background_real_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + ".pdf"
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close(counter_fig)

# #phase plots
# for mode in modes:
# 	for p in range(1, len(percentiles)):
# 		counter_fig = counter_fig + 1
# 		plt.figure(counter_fig)
# 		plt.rcParams.update({'font.size': 15})
# 		plt.rcParams['axes.titlepad'] = 10
# 		source = 'output/two_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '_m' + str(mode) + '_phase' +'.txt'
# 		profile = np.loadtxt(source)
# 		plt.imshow(np.abs(profile), interpolation=None, cmap=plt.cm.Blues, extent = (-0.5+1, len(profile[0,:])-0.5+1, len(profile[:,0])-0.5+1, -0.5+1))
# 		plt.xlabel("$l_1$")
# 		plt.ylabel("$l_2$")
# 		centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
# 		plt.title("Phase of $\\left\\langle\\epsilon^{(" + str(mode) + ")}_{l_1} \\epsilon^{(" + str(-mode) + ")}_{l_2}\\right\\rangle$, "+centrality_class)
# 		plt.colorbar()
# 		filename = "plots/two_point_phase_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + "_m" + str(mode) + ".pdf"
# 		plt.savefig(filename, format='pdf', bbox_inches = "tight")
# 		plt.close(counter_fig)





