import matplotlib.pyplot as plt
import numpy as np 



#modes = [0, 1, 2, 3, 4, 5]
#modes = [0, 1]
percentiles = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
counter_fig = 0

for p in range(1, len(percentiles)):
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.figure(figsize=(10,7.3))
	plt.rcParams.update({'font.size': 23})
	plt.rcParams['axes.titlepad'] = 10
	source_result = 'output/one_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
	profile = np.loadtxt(source_result)
	source_error = 'output/one_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '_error' +'.txt'
	profile_error = np.loadtxt(source_error)
	mMax = 4
	lMax = 10
	l = np.zeros(lMax)
	for i in range(0, lMax):
		l[i] = i+1
	for m in range(0, mMax+1):
		plt.errorbar(l, profile[m, 0:], yerr=profile_error[m,0:], linestyle = 'none', elinewidth=2, capsize = 6, capthick = 2, label = "$m="+str(m)+"$")
		plt.scatter(l, profile[m, 0:], s=10)
	plt.xlabel("$l$")
	plt.ylabel("$\\bar{\\epsilon}_l^{(m)}$")
	plt.legend(loc='best')
	centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
	plt.title("$\\bar{\\epsilon}_l^{(m)}$, "+centrality_class)
	filename = "plots/one_point_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + ".pdf"
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close()



# modulus plots

# for p in range(1, len(percentiles)):
# 	counter_fig = counter_fig + 1
# 	plt.figure(counter_fig)
# 	source = 'output/one_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  + '_modulus' +'.txt'
# 	profile = np.loadtxt(source)
# 	plt.imshow(profile, interpolation=None, cmap=plt.cm.Blues)
# 	plt.xlabel("$l-1$")
# 	plt.ylabel("$m$")
# 	centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
# 	plt.title("$\\left\\|\\left\\langle\\epsilon_{" + str(0) + ",l}\\right\\rangle\\right\\|$, "+centrality_class)
# 	plt.colorbar()
# 	filename = "plots/one_point_modules_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + ".pdf"
# 	plt.savefig(filename, format='pdf', bbox_inches = "tight")
# 	plt.close()

# #phase plots

# for p in range(1, len(percentiles)):
# 	counter_fig = counter_fig + 1
# 	plt.figure(counter_fig)
# 	source = 'output/one_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) +  '_phase' +'.txt'
# 	profile = np.loadtxt(source)
# 	plt.imshow(np.abs(profile), interpolation=None, cmap=plt.cm.Blues)
# 	plt.xlabel("$l-1$")
# 	plt.ylabel("$m$")
# 	centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
# 	plt.title("Phase of $\\left\\langle\\epsilon_{" + str(0) + ",l}\\right\\rangle$, "+centrality_class)
# 	plt.colorbar()
# 	filename = "plots/one_point_phase_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + ".pdf"
# 	plt.savefig(filename, format='pdf', bbox_inches = "tight")
# 	plt.close()


# percentiles = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# counter_fig = 0

# for p in range(1, len(percentiles)):
# 	counter_fig = counter_fig + 1
# 	plt.figure(counter_fig)
# 	source_modulus = 'output/one_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '_modulus' +'.txt'
# 	eps_0l = np.loadtxt(source_modulus)
# 	l = np.arange(len(eps_0l))+1
# 	plt.scatter(l, eps_0l)
# 	plt.xlabel("$l$")
# 	plt.ylabel("$\\left\\|\\left\\langle\\epsilon_{" + str(0) + ",l}\\right\\rangle\\right\\|$")
# 	centrality_class = "centrality class " + str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
# 	plt.title("$\\left\\|\\left\\langle\\epsilon_{" + str(0) + ",l}\\right\\rangle\\right\\|$, "+centrality_class)
# 	filename = "plots/one_point_modules_" + str(percentiles[p-1]) + "-" + str(percentiles[p])+ ".pdf"
# 	plt.savefig(filename, format='pdf', bbox_inches = "tight")
# 	plt.close()

# plt.figure(1)
# source = 'output/two_point_' + str(per)
# profile = np.loadtxt('output/two_point_0-5_m0.txt')
# plt.imshow(profile, interpolation=None, cmap=plt.cm.Blues)
# plt.xlabel("$l_1-1$")
# plt.ylabel("$l_2-1$")
# plt.title("$\\left\\|\\left\\langle\\epsilon_{0,l_1}\\epsilon_{0,l_2}^{*}\\right\\rangle\\right\\|$, centrality class 0-5%")
# plt.colorbar()
# plt.savefig("plots/two_point_modules_0-5_m0.pdf", format='pdf', bbox_inches = "tight")




