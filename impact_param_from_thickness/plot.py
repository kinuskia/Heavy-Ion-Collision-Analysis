import numpy as np 
import matplotlib.pyplot as plt 

b, p = np.loadtxt("b_dist.txt", unpack=True)

plt.plot(b, p)
plt.xlabel("impact parameter [fm]")
plt.ylabel("normed multiplicity [1/fm]")
#plt.yscale("log")
plt.savefig("b_dist.pdf", format="pdf", bbox_inches="tight")