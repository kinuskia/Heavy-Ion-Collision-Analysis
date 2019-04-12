import matplotlib.pyplot as plt
import numpy as np 



#modes = [0, 1, 2, 3, 4, 5]
#modes = [0, 1]
percentiles = [0, 5, 10, 20, 30, 40]
counter_fig = 0

for p in range(1, len(percentiles)):
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.figure(figsize=(10,7.3))
	plt.rcParams.update({'font.size': 23})
	plt.rcParams['axes.titlepad'] = 10
	plt.xlabel("N_part")
	plt.ylabel("#")
	centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
	plt.title("N_part, "+centrality_class)
	filename = "plots/N_part_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + ".pdf"
	source = 'output/collision_specs_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
	b, npart, mult = np.loadtxt(source, unpack = True)
	mean = np.mean(npart)
	std = np.std(npart)
	print(mean)
	x = np.arange(0, 500, 1)
	y = np.ones(len(x))
	for i in range(0, len(x)):
		#print(x[i])
		y[i] = 1./np.sqrt(2*np.pi)/std*np.exp(-1./2.*(1.0*x[i]-mean)*(1.0*x[i]-mean)/std/std)
	nbins = 32
	plt.plot(x, y, label = "fit")
	plt.hist(npart, nbins, density = True)
	plt.legend(loc='best')
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close()

	print(str(percentiles[p-1]) + '-' + str(percentiles[p]))
	print("mu: ", mean)
	print("sigma: ", std)
	print("(mu/sigma)^2:", std*std/mean/mean)
	k = 0.8
	print("1/mu/k:", 1/mean/k)






