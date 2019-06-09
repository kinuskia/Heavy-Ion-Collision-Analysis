import matplotlib.pyplot as plt
import numpy as np 


percentiles = [0, 1, 10, 11, 20, 21]
counter_fig = 0

# modulus plots
plt.figure(1)
plt.figure(figsize=(10,8))
plt.rcParams.update({'font.size': 23})
plt.rcParams['axes.titlepad'] = 20

for p in range(1, len(percentiles)):
	#two-point functions
	source = 'output/two_point_random_' + str(percentiles[p-1]) + '-' + str(percentiles[p-1]+1) + '_m' + str(0) + '_real' +'.txt'
	profile = np.loadtxt(source)
	#one-point-functions
	source_one = 'output/one_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p-1]+1)  +'.txt'
	profile_one_ml = np.loadtxt(source_one)
	#subtract disconnected part for l=l'=1
	profile[0,0] -= profile_one_ml[0, 0]*profile_one_ml[0,0]

	#two-point error
	source_error = 'output/two_point_random_' + str(percentiles[p-1]) + '-' + str(percentiles[p-1]+1) + '_m' + str(0) + '_real_error' +'.txt'
	profile_error = np.loadtxt(source_error)
	#one-point error
	source_one_error = 'output/one_point_' + str(percentiles[p-1]) + '-' + str(percentiles[p-1]+1)  +'_error.txt'
	profile_one_ml_error = np.loadtxt(source_one_error)

	error = profile_error[0,:]
	# what about l=l'=1 ??
	if (p==1):
		for item in error:
			print(item)

	x = np.arange(1, len(profile[0,:])+1)
	centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p-1]+1) + '%'
	plt.scatter(x, profile[0,:], label=centrality_class)
	plt.errorbar(x, profile[0,:], yerr=error, linestyle = 'none', elinewidth=2, capsize = 6, capthick = 2)
	#plt.plot(x, profile[0,:], label=centrality_class)
	#plt.ylim(-0.00009, 0.00011)
	#maximal_value = max(np.amax(profile), -np.amin(profile))
	#plt.imshow(profile, interpolation=None, cmap=plt.cm.RdYlGn,vmin = -maximal_value, vmax = maximal_value, extent = (-0.5+1, len(profile[0,:])-0.5+1, len(profile[:,0])/2, -len(profile[:,0])/2))
	#, extent = (-0.5+1, len(profile[0,:])-0.5+1, len(profile[:,0])-0.5+1, -0.5+1)
plt.xlabel("$l$")
plt.ylabel("connected two-mode correlation")
	
plt.title("$\\left\\langle\\epsilon^{(0)}_{1} \\epsilon^{(0)}_{l}\\right\\rangle_{\\circ, c}$ ")
plt.legend(loc="best")
filename = "plots/two_point_random_connected_background" + ".pdf"
plt.savefig(filename, format='pdf', bbox_inches = "tight")
plt.close(1)







