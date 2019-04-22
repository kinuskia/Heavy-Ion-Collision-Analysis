import matplotlib.pyplot as plt
import numpy as np 

def gauss(x, mu, sigma):
	return 1./np.sqrt(2*np.pi)/sigma*np.exp(-1./2./sigma/sigma*(x-mu)*(x-mu))



#modes = [0, 1, 2, 3, 4, 5]
#modes = [0, 1]
percentiles = [0, 5, 10, 20, 30, 40]
counter_fig = 0

color0005 = '#1f77b4'
color3040 = '#9467bd'

nbins = 16
p = 1
filename = "plots/N_part" +  ".pdf"
source1 = 'output/collision_specs_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
b1, npart1, mult1 = np.loadtxt(source1, unpack = True)
plt.figure()
bins01, edges01, patches01 = plt.hist(npart1, nbins)
plt.close()
counter_fig = 1
plt.figure(counter_fig)
plt.figure(figsize=(10,7.3))
plt.rcParams.update({'font.size': 23})
plt.rcParams['axes.titlepad'] = 10
plt.xlabel("$N_{part}$")
plt.ylabel("probability density")
centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
plt.title("$N_{part}$ fluctuations")
mean1 = np.mean(npart1)
std1 = np.std(npart1)
print(mean1)
bins1, edges1, patches1 = plt.hist(npart1, nbins, density = True, color = color0005, alpha = 0.6)
from scipy.optimize import curve_fit
popt1, pcov1 = curve_fit(gauss, (edges1[1:]+edges1[0:-1])/2, bins1, sigma=np.sqrt(bins01)*bins1/bins01, p0=[mean1, std1])
x = np.arange(0, 500, 1)
y = np.ones(len(x))
for i in range(0, len(x)):
	#print(x[i])
	y[i] = gauss(x[i], *popt1)
plt.plot(x, y, label = "0-5%", color = color0005)
plt.errorbar((edges1[1:]+edges1[0:-1])/2, bins1, yerr = np.sqrt(bins01)/bins01*bins1, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1, color = color0005)


#plt.savefig(filename, format='pdf', bbox_inches = "tight")
#plt.close()
print(str(percentiles[p-1]) + '-' + str(percentiles[p]))
print("mu: ", popt1[0])
print("sigma: ", popt1[1])
print("(sigma/mu)^2:", (popt1[1]/popt1[0])**2)
k = 0.8
print("1/mu/k:", 1/popt1[0]/k)
print("1/mu/k+(sigma/mu)^2:", 1/popt1[0]/k+(popt1[1]/popt1[0])**2)

p = 5
source2 = 'output/collision_specs_' + str(percentiles[p-1]) + '-' + str(percentiles[p])  +'.txt'
b2, npart2, mult2 = np.loadtxt(source2, unpack = True)
plt.figure()
bins02, edges02, patches02 = plt.hist(npart2, nbins)
plt.close()
centrality_class = str(percentiles[p-1]) + '-' + str(percentiles[p]) + '%'
mean2 = np.mean(npart2)
std2 = np.std(npart2)
print(mean2)
bins2, edges2, patches2 = plt.hist(npart2, nbins, density = True, color = color3040, alpha = 0.6)

popt2, pcov2 = curve_fit(gauss, (edges2[1:]+edges2[0:-1])/2, bins2, sigma=np.sqrt(bins02)*bins2/bins02, p0=[mean2, std2])
x = np.arange(0, 500, 1)
y2 = np.ones(len(x))
for i in range(0, len(x)):
	#print(x[i])
	y2[i] = gauss(x[i], *popt2)
plt.plot(x, y2, label = "30-40%", color = color3040)
plt.errorbar((edges2[1:]+edges2[0:-1])/2, bins2, yerr = np.sqrt(bins02)/bins02*bins2, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1, color = color3040)
plt.legend(loc='best')
plt.savefig(filename, format='pdf', bbox_inches = "tight")
plt.close()
print(str(percentiles[p-1]) + '-' + str(percentiles[p]))
print("mu: ", popt2[0])
print("sigma: ", popt2[1])
print("(sigma/mu)^2:", (popt2[1]/popt2[0])**2)
k = 0.8
print("1/mu/k:", 1/popt2[0]/k)
print("1/mu/k+(sigma/mu)^2:", 1/popt2[0]/k+(popt2[1]/popt2[0])**2)






