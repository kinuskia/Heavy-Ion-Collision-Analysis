import matplotlib.pyplot as plt
import numpy as np 

modes = [0, 1, 2, 3, 4]
#modes = [0, 1]
percentiles = [0, 5, 10, 20, 30, 40]
counter_fig = 0

# modulus plots
for mode in modes:
	for p in range(1, len(percentiles)):
		counter_fig = counter_fig + 1
		plt.figure(counter_fig)
		plt.figure(figsize=(10,10))
		plt.rcParams.update({'font.size': 23})
		plt.rcParams['axes.titlepad'] = 20
		source = 'output/two_point_random_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '_m' + str(mode) + '_real' +'.txt'
		profile = np.loadtxt(source)
		maximal_value = max(np.amax(profile), -np.amin(profile))
		plt.imshow(profile, interpolation=None, cmap=plt.cm.RdYlGn,vmin = -maximal_value, vmax = maximal_value, extent = (-0.5+1, len(profile[0,:])-0.5+1, len(profile[:,0])-0.5+1, -0.5+1))
		plt.xlabel("$l_2$")
		plt.ylabel("$l_1$")
		centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
		plt.title("$\\left\\langle\\epsilon^{(" + str(mode) + ")}_{l_1} \\epsilon^{(" + str(-mode) + ")}_{l_2}\\right\\rangle$, "+centrality_class)
		plt.colorbar()
		filename = "plots/two_point_real_random_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + "_m" + str(mode) + ".pdf"
		plt.savefig(filename, format='pdf', bbox_inches = "tight")
		plt.close(counter_fig)

#phase plots
for mode in modes:
	for p in range(1, len(percentiles)):
		counter_fig = counter_fig + 1
		plt.figure(counter_fig)
		plt.rcParams.update({'font.size': 15})
		plt.rcParams['axes.titlepad'] = 10
		source = 'output/two_point_random_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '_m' + str(mode) + '_phase' +'.txt'
		profile = np.loadtxt(source)
		plt.imshow(np.abs(profile), interpolation=None, cmap=plt.cm.Blues, extent = (-0.5+1, len(profile[0,:])-0.5+1, len(profile[:,0])-0.5+1, -0.5+1))
		plt.xlabel("$l_1$")
		plt.ylabel("$l_2$")
		centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
		plt.title("Phase of $\\left\\langle\\epsilon^{(" + str(mode) + ")}_{l_1} \\epsilon^{(" + str(-mode) + ")}_{l_2}\\right\\rangle$, "+centrality_class)
		plt.colorbar()
		filename = "plots/two_point_phase_random_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + "_m" + str(mode) + ".pdf"
		plt.savefig(filename, format='pdf', bbox_inches = "tight")
		plt.close(counter_fig)

# connected diagrams
for mode in modes:
	for p in range(1, len(percentiles)):
		counter_fig = counter_fig + 1
		plt.figure(counter_fig)
		plt.figure(figsize=(10,10))
		plt.rcParams.update({'font.size': 23})
		plt.rcParams['axes.titlepad'] = 20
		source = 'output/two_point_random_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '_m' + str(mode) + '_real' +'.txt'
		profile = np.loadtxt(source)
		# import one-point functions
		source_one = 'output/one_point_random_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
		profile_one_ml = np.loadtxt(source_one)
		# subtract unconnected parts from correlation function
		for i in range(0, len(profile[0:, 0])):
			for j in range(0, len(profile[0, 0:])):
				profile[i, j] -= profile_one_ml[mode, i]*profile_one_ml[mode, j] 
		maximal_value = max(np.amax(profile), -np.amin(profile))
		plt.imshow(profile, interpolation=None, cmap=plt.cm.RdYlGn,vmin = -maximal_value, vmax = maximal_value, extent = (-0.5+1, len(profile[0,:])-0.5+1, len(profile[:,0])-0.5+1, -0.5+1))
		plt.xlabel("$l_2$")
		plt.ylabel("$l_1$")
		centrality_class =  str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
		plt.title("$\\left\\langle\\epsilon^{(" + str(mode) + ")}_{l_1} \\epsilon^{(" + str(-mode) + ")}_{l_2}\\right\\rangle-\\left\\langle\\epsilon^{(" + str(mode) + ")}_{l_1} \\right\\rangle\\left\\langle\\epsilon^{(" + str(mode) + ")}_{l_2} \\right\\rangle$, "+centrality_class)
		plt.colorbar()
		filename = "plots/two_point_real_connected_random_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + "_m" + str(mode) + ".pdf"
		plt.savefig(filename, format='pdf', bbox_inches = "tight")
		plt.close(counter_fig)

# plt.figure(1)
# source = 'output/two_point_' + str(per)
# profile = np.loadtxt('output/two_point_0-5_m0.txt')
# plt.imshow(profile, interpolation=None, cmap=plt.cm.Blues)
# plt.xlabel("$l_1-1$")
# plt.ylabel("$l_2-1$")
# plt.title("$\\left\\|\\left\\langle\\epsilon_{0,l_1}\\epsilon_{0,l_2}^{*}\\right\\rangle\\right\\|$, centrality class 0-5%")
# plt.colorbar()
# plt.savefig("plots/two_point_modules_0-5_m0.pdf", format='pdf', bbox_inches = "tight")




