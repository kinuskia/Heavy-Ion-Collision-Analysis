import matplotlib.pyplot as plt
import numpy as np 

modes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
#modes = [0, 1]
percentiles = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
counter_fig = 0

# modulus plots
for mode in modes:
	for p in range(1, len(percentiles)):
		counter_fig = counter_fig + 1
		plt.figure(counter_fig)
		source = 'output/two_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '_m' + str(mode) + '_modulus' +'.txt'
		profile = np.loadtxt(source)
		plt.imshow(profile, interpolation=None, cmap=plt.cm.Blues)
		plt.xlabel("$l_1-1$")
		plt.ylabel("$l_2-1$")
		centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
		plt.title("$\\left\\|\\left\\langle\\epsilon_{" + str(mode) + ",l_1}\\epsilon_{" + str(mode) + ",l_2}^{*}\\right\\rangle\\right\\|$, "+centrality_class)
		plt.colorbar()
		filename = "plots/two_point_modules_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + "_m" + str(mode) + ".pdf"
		plt.savefig(filename, format='pdf', bbox_inches = "tight")
		plt.close()

#phase plots
for mode in modes:
	for p in range(1, len(percentiles)):
		counter_fig = counter_fig + 1
		plt.figure(counter_fig)
		source = 'output/two_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '_m' + str(mode) + '_phase' +'.txt'
		profile = np.loadtxt(source)
		plt.imshow(np.abs(profile), interpolation=None, cmap=plt.cm.Blues)
		plt.xlabel("$l_1-1$")
		plt.ylabel("$l_2-1$")
		centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
		plt.title("Phase of $\\left\\langle\\epsilon_{" + str(mode) + ",l_1}\\epsilon_{" + str(mode) + ",l_2}^{*}\\right\\rangle$, "+centrality_class)
		plt.colorbar()
		filename = "plots/two_point_phase_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + "_m" + str(mode) + ".pdf"
		plt.savefig(filename, format='pdf', bbox_inches = "tight")
		plt.close()

# plt.figure(1)
# source = 'output/two_point_' + str(per)
# profile = np.loadtxt('output/two_point_0-5_m0.txt')
# plt.imshow(profile, interpolation=None, cmap=plt.cm.Blues)
# plt.xlabel("$l_1-1$")
# plt.ylabel("$l_2-1$")
# plt.title("$\\left\\|\\left\\langle\\epsilon_{0,l_1}\\epsilon_{0,l_2}^{*}\\right\\rangle\\right\\|$, centrality class 0-5%")
# plt.colorbar()
# plt.savefig("plots/two_point_modules_0-5_m0.pdf", format='pdf', bbox_inches = "tight")




