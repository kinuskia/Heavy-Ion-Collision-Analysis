import matplotlib.pyplot as plt 
import numpy as np 

r, phi = np.loadtxt("polars.txt", unpack=True)

plt.scatter(r, phi, s=0.1)
plt.xlabel("$r$")
plt.ylabel("$\\phi$")
plt.savefig("phi-r.pdf", format="pdf", bbox_inches="tight")