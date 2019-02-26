import matplotlib.pyplot as plt 
import numpy as np 

r0, e0, de0 = np.loadtxt("output/e_m_mean_m0_0-5.txt", unpack=True)
r1, e1, de1 = np.loadtxt("output/e_m_mean_m0_5-10.txt", unpack=True)
r2, e2, de2 = np.loadtxt("output/e_m_mean_m0_10-20.txt", unpack=True)
r3, e3, de3 = np.loadtxt("output/e_m_mean_m0_20-30.txt", unpack=True)
r4, e4, de4 = np.loadtxt("output/e_m_mean_m0_30-40.txt", unpack=True)

e0 = e0*np.pi
de0 = de0*np.pi
e1 = e1*np.pi
de1 = de1*np.pi
e2 = e2*np.pi
de2 = de2*np.pi
e3 = e3*np.pi
de3 = de3*np.pi
e4 = e4*np.pi
de4 = de4*np.pi
# r0, e0, de0 = np.loadtxt("output/e_m_mean_m0_30-40.txt", unpack=True)
# r1, e1, de1 = np.loadtxt("output/e_m_mean_m1_30-40.txt", unpack=True)
# r2, e2, de2 = np.loadtxt("output/e_m_mean_m2_30-40.txt", unpack=True)
# r3, e3, de3 = np.loadtxt("output/e_m_mean_m3_30-40.txt", unpack=True)
# r4, e4, de4 = np.loadtxt("output/e_m_mean_m4_30-40.txt", unpack=True)




plt.rcParams.update({'font.size': 23})

plt.figure(figsize=(10,5))
plt.plot(r0, e0+de0, linewidth = 0.5, color = "blue")
plt.plot(r0, e0-de0, linewidth = 0.5, color = "blue")
plt.fill_between(r0, e0+de0, e0-de0, color = "blue", label="0-5%")
plt.plot(r1, e1+de1, linewidth = 0.5, color = "orange")
plt.plot(r1, e1-de1, linewidth = 0.5, color = "orange")
plt.fill_between(r1, e1+de1, e1-de1, color = "orange", label="5-10%")
plt.plot(r2, e2+de2, linewidth = 0.5, color = "green")
plt.plot(r2, e2-de2, linewidth = 0.5, color = "green")
plt.fill_between(r2, e2+de2, e2-de2, color = "green", label="10-20%")
plt.plot(r3, e3+de3, linewidth = 0.5, color = "red")
plt.plot(r3, e3-de3, linewidth = 0.5, color = "red")
plt.fill_between(r3, e3+de3, e3-de3, color = "red", label="20-30%")
plt.plot(r4, e4+de4, linewidth = 0.5, color = "purple")
plt.plot(r4, e4-de4, linewidth = 0.5, color = "purple")
plt.fill_between(r4, e4+de4, e4-de4, color = "purple", label="30-40%")
#plt.errorbar(r, e, yerr = de, linestyle = "none", capsize = 2 , color = "blue")
plt.xlabel("$r$ [fm]")
plt.ylabel("$W(r)$")
plt.title("Background field $W(r)$")
plt.legend(loc = "best")
plt.savefig("plots/W_r_centrality.pdf", format="pdf", bbox_inches="tight")
plt.yscale("log")
plt.xlabel("$r$ [fm]")
plt.ylabel("$W(r)$")
plt.legend(loc = "best")
plt.savefig("plots/W_r_centrality_log.pdf", format="pdf", bbox_inches="tight")


# Check if integral equals 1/2
p = [e0,e1,e2,e3,e4]
for liste in p:
	integral = 0
	dr = r0[1]-r0[0]
	for i in range(0, len(liste)):
		integral += r0[i]*liste[i]*dr
	print(integral)

