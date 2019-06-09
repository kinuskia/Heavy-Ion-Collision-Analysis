import numpy as np 
import matplotlib.pyplot as plt 

r1, rho = np.loadtxt("rho.txt", unpack = True)
r2, W = np.loadtxt("W.txt", unpack = True)

plt.figure(1)
plt.plot(r1, rho)
plt.savefig("rho.pdf", format="pdf", bbox_inches = "tight")
plt.close(1)

plt.figure(2)
plt.plot(r2, W)
plt.savefig("W.pdf", format="pdf", bbox_inches = "tight")
plt.close(2)