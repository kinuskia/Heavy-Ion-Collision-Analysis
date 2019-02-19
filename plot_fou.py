import matplotlib.pyplot as plt 
import numpy as np 

r0, e0, de0 = np.loadtxt("output/e_m_mean_m0_0-5.txt", unpack=True)
r1, e1, de1 = np.loadtxt("output/e_m_mean_m1_0-5.txt", unpack=True)
r2, e2, de2 = np.loadtxt("output/e_m_mean_m2_0-5.txt", unpack=True)
r3, e3, de3 = np.loadtxt("output/e_m_mean_m3_0-5.txt", unpack=True)
r4, e4, de4 = np.loadtxt("output/e_m_mean_m4_0-5.txt", unpack=True)

# r0, e0, de0 = np.loadtxt("output/e_m_mean_m0_30-40.txt", unpack=True)
# r1, e1, de1 = np.loadtxt("output/e_m_mean_m1_30-40.txt", unpack=True)
# r2, e2, de2 = np.loadtxt("output/e_m_mean_m2_30-40.txt", unpack=True)
# r3, e3, de3 = np.loadtxt("output/e_m_mean_m3_30-40.txt", unpack=True)
# r4, e4, de4 = np.loadtxt("output/e_m_mean_m4_30-40.txt", unpack=True)






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
plt.plot(r4, e4+de4, linewidth = 0.5, color = "purple")
plt.plot(r4, e4-de4, linewidth = 0.5, color = "purple")
plt.fill_between(r4, e4+de4, e4-de4, color = "purple", label="$m = $4")



#plt.errorbar(r, e, yerr = de, linestyle = "none", capsize = 2 , color = "blue")
plt.xlabel("$r$ [fm]")
plt.ylabel("$\\left\\langle e^{(m)}\\right\\rangle$")
plt.legend(loc = "best")
plt.savefig("plots/energy-e_mean_0-5.pdf", format="pdf", bbox_inches="tight")