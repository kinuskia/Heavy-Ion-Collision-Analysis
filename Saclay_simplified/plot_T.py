import numpy as np 
import matplotlib.pyplot as plt 


r, Glauber, W = np.loadtxt("Saclay_simplified/thickness.txt", unpack=True)

plt.plot(r, Glauber, label ="Glauber")
plt.plot(r, W, label = "Trento")
plt.xlabel("$r [fm]$")
plt.ylabel("$Q_s^2(r)/Q_0^2$")
plt.legend(loc= "best")
plt.title("One-nucleus thickness function")
plt.savefig("Saclay_simplified/thickness.pdf", format="pdf", bbox_inches = "tight")