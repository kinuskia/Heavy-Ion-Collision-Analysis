import matplotlib.pyplot as plt
import numpy as np 

plt.figure(1)
profile = np.loadtxt('output/decomposition0.txt')
plt.imshow(profile, interpolation=None, cmap=plt.cm.Blues)
plt.xlabel("$l-1$")
plt.ylabel("$m$")
plt.title("$|\\epsilon_{m,l}|$")
plt.colorbar()
plt.savefig("plots/FB_coeffs_module.pdf", format='pdf', bbox_inches = "tight")

plt.figure(2)
l = np.arange(1, len(profile[0,0:-1])+1, 1)
plt.scatter(l,profile[0,0:-1])
plt.title("Decrease of $|\\epsilon_{0,l}|$ with $l$")
plt.xlabel("$l$")
plt.ylabel("$|\\epsilon_{0,l}|$")
plt.savefig("plots/coeffs_decrease.pdf", format = "pdf", bbox_inches = "tight")