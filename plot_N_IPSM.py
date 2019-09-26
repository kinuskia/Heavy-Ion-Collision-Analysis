import matplotlib.pyplot as plt
import numpy as np 


# Read in IPSM fit result
percentiles = np.arange(100)
N_values = np.zeros(len(percentiles))


for k in range(0, len(percentiles)):
	centrality_class =  str(percentiles[k]) + '-' + str(percentiles[k]+1)
	filename_fit = "../IPSM-Fit_One/output/" + centrality_class + ".txt"
	N_value = np.loadtxt(filename_fit, unpack=True)
	N_values[k] = N_value


plt.figure(2)

plt.figure(figsize=(10,5))
plt.rcParams.update({'font.size': 23})
plt.rcParams['axes.titlepad'] = 10

plt.plot(percentiles, N_values)
plt.xlabel("centrality class $[p, p+1]$ in %")
plt.ylabel("$\\mu_N$")
plt.title("Number $\\mu_N$ of sources")
plt.savefig("plots/N_IPSM.pdf", format="pdf", bbox_inches="tight")



