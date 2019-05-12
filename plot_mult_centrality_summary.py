import matplotlib.pyplot as plt
import numpy as np 

def gauss(x, mu, sigma):
	return 1./np.sqrt(2*np.pi)/sigma*np.exp(-1./2./sigma/sigma*(x-mu)*(x-mu))



#modes = [0, 1, 2, 3, 4, 5]
#modes = [0, 1]
percentiles = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21]
counter_fig = 0

color0001 = '#1f77b4'
color0102 = '#ff7f0e'
color0203 = '#2ca02c'
color0304 = '#d62728'
color0405 = '#9467bd'

nbins = 16
p = 1
filename = "plots/N_mult" +  ".pdf"
source1 = 'output/collision_specs_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
b1, npart1, mult1 = np.loadtxt(source1, unpack = True)
plt.figure()
bins01, edges01, patches01 = plt.hist(mult1, nbins)
plt.close()
counter_fig = 1
plt.figure(counter_fig)
plt.figure(figsize=(10,7.3))
plt.rcParams.update({'font.size': 23})
plt.rcParams['axes.titlepad'] = 10
plt.xlabel("Multiplicity [arb. unit]")
plt.ylabel("probability density")
centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
plt.title("Multiplicity fluctuations")

bins1, edges1, patches1 = plt.hist(mult1, nbins, density = True, color = color0001, alpha = 0.6 , label="0-1%")
plt.errorbar((edges1[1:]+edges1[0:-1])/2, bins1, yerr = np.sqrt(bins01)/bins01*bins1, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1, color = color0001)


p = 2
source2 = 'output/collision_specs_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
b2, npart2, mult2 = np.loadtxt(source2, unpack = True)
plt.figure()
bins02, edges02, patches02 = plt.hist(mult2, nbins)
plt.close()
centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'

bins2, edges2, patches2 = plt.hist(mult2, nbins, density = True, color = color0102, alpha = 0.6, label="1-2%")
plt.errorbar((edges2[1:]+edges2[0:-1])/2, bins2, yerr = np.sqrt(bins02)/bins02*bins2, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1, color = color0102)
plt.legend(loc='best')


p = 3
source3 = 'output/collision_specs_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
b3, npart3, mult3 = np.loadtxt(source3, unpack = True)
plt.figure()
bins03, edges03, patches03 = plt.hist(mult3, nbins)
plt.close()
centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
bins3, edges3, patches3 = plt.hist(mult3, nbins, density = True, color = color0203, alpha = 0.6, label="2-3%")
plt.errorbar((edges3[1:]+edges3[0:-1])/2, bins3, yerr = np.sqrt(bins03)/bins03*bins3, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1, color = color0203)
plt.legend(loc='best')


p = 4
source4 = 'output/collision_specs_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
b4, npart4, mult4 = np.loadtxt(source4, unpack = True)
plt.figure()
bins04, edges04, patches04 = plt.hist(mult4, nbins)
plt.close()
centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'

bins4, edges4, patches4 = plt.hist(mult4, nbins, density = True, color = color0304, alpha = 0.6, label="3-4%")
plt.errorbar((edges4[1:]+edges4[0:-1])/2, bins4, yerr = np.sqrt(bins04)/bins04*bins4, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1, color = color0304)
plt.legend(loc='best')


plt.yscale("log")
plt.savefig(filename, format='pdf', bbox_inches = "tight")
plt.close()







