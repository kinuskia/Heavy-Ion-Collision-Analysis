import matplotlib.pyplot as plt 
import numpy as np 

r0, e0, de0 = np.loadtxt("output/e_ave_fourier0.txt", unpack=True)
r1, e1, de1 = np.loadtxt("output/e_ave_fourier1.txt", unpack=True)
r2, e2, de2 = np.loadtxt("output/e_ave_fourier2.txt", unpack=True)
r3, e3, de3 = np.loadtxt("output/e_ave_fourier3.txt", unpack=True)




plt.plot(r0, e0+de0, linewidth = 0.5, color = "blue")
plt.plot(r0, e0-de0, linewidth = 0.5, color = "blue")
plt.fill_between(r0, e0+de0, e0-de0, color = "blue", label="$m = $0")
plt.plot(r1, e1+de1, linewidth = 0.5, color = "orange")
plt.plot(r1, e1-de1, linewidth = 0.5, color = "orange")
plt.fill_between(r1, e1+de1, e1-de1, color = "orange", label="$m = $1")
plt.plot(r2, e2+de2, linewidth = 0.5, color = "green")
plt.plot(r2, e2-de2, linewidth = 0.5, color = "green")
plt.fill_between(r2, e2+de2, e2-de2, color = "green", label="$m = $2")
plt.plot(r3, e3+de3, linewidth = 0.5, color = "red")
plt.plot(r3, e3-de3, linewidth = 0.5, color = "red")
plt.fill_between(r3, e3+de3, e3-de3, color = "red", label="$m = $3")



#plt.errorbar(r, e, yerr = de, linestyle = "none", capsize = 2 , color = "blue")
plt.xlabel("$r$ [fm]")
plt.ylabel("normalized energy density")
plt.legend(loc = "best")
plt.savefig("plots/energy-e_fou.pdf", format="pdf", bbox_inches="tight")