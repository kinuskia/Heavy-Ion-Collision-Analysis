import matplotlib.pyplot as plt 
import numpy as np 

r0, e01, e02, e03 = np.loadtxt("output/bessel1.txt", unpack=True)





plt.plot(r0, e01, linewidth = 0.5, color = "blue", label = "$m=1, l=1$")
plt.plot(r0, e02, linewidth = 0.5, color = "orange", label = "$m=1, l=2$")
plt.plot(r0, e03, linewidth = 0.5, color = "purple", label = "$m=1, l=3$")
plt.title("$W(r) R^2 J_m(z_m^{(l)} r)$")
plt.xlabel("$r$ in fm")
plt.legend(loc="best")



#plt.errorbar(r, e, yerr = de, linestyle = "none", capsize = 2 , color = "blue")


plt.savefig("plots/bessel1.pdf", format="pdf", bbox_inches="tight")