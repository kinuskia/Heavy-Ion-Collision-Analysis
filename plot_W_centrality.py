import matplotlib.pyplot as plt
import numpy as np 

r, p1, p2, p3, p4, p5 = np.loadtxt("output/weight_functions.txt", unpack = True)
r, dp1, dp2, dp3, dp4, dp5 = np.loadtxt("output/weight_functions_error.txt", unpack = True)

plt.figure(1)

plt.plot(r, r*(p1+dp1), linewidth = 0.5, color = '#1f77b4')
plt.plot(r, r*(p1-dp1), linewidth = 0.5, color = '#1f77b4')
plt.fill_between(r, r*(p1+dp1), r*(p1-dp1), color = '#1f77b4', label="0-2%")

plt.plot(r, r*(p2+dp2), linewidth = 0.5, color = '#ff7f0e')
plt.plot(r, r*(p2-dp2), linewidth = 0.5, color = '#ff7f0e')
plt.fill_between(r, r*(p2+dp2), r*(p2-dp2), color = '#ff7f0e', label="2-4%")

plt.plot(r, r*(p3+dp3), linewidth = 0.5, color = '#2ca02c')
plt.plot(r, r*(p3-dp3), linewidth = 0.5, color = '#2ca02c')
plt.fill_between(r, r*(p3+dp3), r*(p3-dp3), color = '#2ca02c', label="4-6%")

plt.plot(r, r*(p4+dp4), linewidth = 0.5, color = '#d62728')
plt.plot(r, r*(p4-dp4), linewidth = 0.5, color = '#d62728')
plt.fill_between(r, r*(p4+dp4), r*(p4-dp4), color = '#d62728', label="6-8%")

plt.plot(r, r*(p5+dp5), linewidth = 0.5, color = '#9467bd')
plt.plot(r, r*(p5-dp5), linewidth = 0.5, color = '#9467bd')
plt.fill_between(r, r*(p5+dp5), r*(p5-dp5), color = '#9467bd', label="8-10%")


plt.xlabel("$r$ [fm]")
plt.ylabel("weighting function $r W(r)$")
plt.legend(loc = "best")
plt.savefig("plots/W_centrality.pdf", format="pdf", bbox_inches="tight")

plt.figure(2)

plt.plot(r, (p1+dp1), linewidth = 0.5, color = '#1f77b4')
plt.plot(r, (p1-dp1), linewidth = 0.5, color = '#1f77b4')
plt.fill_between(r, (p1+dp1), (p1-dp1), color = '#1f77b4', label="0-2%")

plt.plot(r, (p2+dp2), linewidth = 0.5, color = '#ff7f0e')
plt.plot(r, (p2-dp2), linewidth = 0.5, color = '#ff7f0e')
plt.fill_between(r, (p2+dp2), (p2-dp2), color = '#ff7f0e', label="2-4%")

plt.plot(r, (p3+dp3), linewidth = 0.5, color = '#2ca02c')
plt.plot(r, (p3-dp3), linewidth = 0.5, color = '#2ca02c')
plt.fill_between(r, (p3+dp3), (p3-dp3), color = '#2ca02c', label="4-6%")

plt.plot(r, (p4+dp4), linewidth = 0.5, color = '#d62728')
plt.plot(r, (p4-dp4), linewidth = 0.5, color = '#d62728')
plt.fill_between(r, (p4+dp4), (p4-dp4), color = '#d62728', label="6-8%")

plt.plot(r, (p5+dp5), linewidth = 0.5, color = '#9467bd')
plt.plot(r, (p5-dp5), linewidth = 0.5, color = '#9467bd')
plt.fill_between(r, (p5+dp5), (p5-dp5), color = '#9467bd', label="8-10%")

plt.yscale("log")

plt.xlabel("$r$ [fm]")
plt.ylabel("weighting function $W(r)$")
plt.legend(loc = "best")
plt.savefig("plots/W_centrality_log.pdf", format="pdf", bbox_inches="tight")



# Check if integral equals 1/2
p = [p1, p2, p3, p4, p5]
for liste in p:
	integral = 0
	dr = r[1]-r[0]
	for i in range(0, len(liste)):
		integral += r[i]*liste[i]*dr
	print(integral)


