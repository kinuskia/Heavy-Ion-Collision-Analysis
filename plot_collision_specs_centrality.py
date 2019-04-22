import matplotlib.pyplot as plt
import numpy as np 

def gauss(x, mu, sigma):
	return 1./np.sqrt(2*np.pi)/sigma*np.exp(-1./2./sigma/sigma*(x-mu)*(x-mu))



#modes = [0, 1, 2, 3, 4, 5]
#modes = [0, 1]
percentiles = [0, 5, 10, 20, 30, 40]
counter_fig = 0

for p in range(1, len(percentiles)):
	nbins = 16
	filename = "plots/N_part_" + str(percentiles[p-1]) + "-" + str(percentiles[p]) + ".pdf"
	source = 'output/collision_specs_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
	b, npart, mult = np.loadtxt(source, unpack = True)
	plt.figure()
	bins0, edges0, patches0 = plt.hist(npart, nbins)
	plt.close()
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.figure(figsize=(10,7.3))
	plt.rcParams.update({'font.size': 23})
	plt.rcParams['axes.titlepad'] = 10
	plt.xlabel("$N_{part}$")
	plt.ylabel("#")
	centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
	plt.title("N_part, "+centrality_class)
	mean = np.mean(npart)
	std = np.std(npart)
	print(mean)
	bins, edges, patches = plt.hist(npart, nbins, density = True)
	from scipy.optimize import curve_fit
	popt, pcov = curve_fit(gauss, (edges[1:]+edges[0:-1])/2, bins, sigma=np.sqrt(bins0)*bins/bins0, p0=[mean, std])
	x = np.arange(0, 500, 1)
	y = np.ones(len(x))
	for i in range(0, len(x)):
		#print(x[i])
		y[i] = gauss(x[i], *popt)
	plt.plot(x, y, label = "fit")
	plt.errorbar((edges[1:]+edges[0:-1])/2, bins, yerr = np.sqrt(bins0)/bins0*bins, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1)
	
	plt.legend(loc='best')
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close()

	print(str(percentiles[p-1]) + '-' + str(percentiles[p]))
	print("mu: ", popt[0])
	print("sigma: ", popt[1])
	print("(sigma/mu)^2:", (popt[1]/popt[0])**2)
	k = 0.8
	print("1/mu/k:", 1/popt[0]/k)
	print("1/mu/k+(sigma/mu)^2:", 1/popt[0]/k+(popt[1]/popt[0])**2)








