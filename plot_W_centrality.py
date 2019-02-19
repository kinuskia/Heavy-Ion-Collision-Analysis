import matplotlib.pyplot as plt
import numpy as np 

r, p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11 = np.loadtxt("output/weight_functions.txt", unpack = True)

plt.figure(1)
plt.plot(r, r*p1, label = "0-5%")
plt.plot(r, r*p2, label = "5-10%")
plt.plot(r, r*p3, label = "10-20%")
plt.plot(r, r*p4, label = "20-30%")
plt.plot(r, r*p5, label = "30-40%")
plt.plot(r, r*p6, label = "40-50%")
plt.plot(r, r*p7, label = "50-60%")
plt.plot(r, r*p8, label = "60-70%")
plt.plot(r, r*p9, label = "70-80%")
plt.plot(r, r*p10, label = "80-90%")
plt.plot(r, r*p11, label = "90-100%")


plt.xlabel("$r$ [fm]")
plt.ylabel("weighting function $r W(r)$")
plt.legend(loc = "best")
plt.savefig("plots/W_centrality.pdf", format="pdf", bbox_inches="tight")

plt.figure(2)

plt.plot(r, p1, label = "0-5%")
plt.plot(r, p2, label = "5-10%")
plt.plot(r, p3, label = "10-20%")
plt.plot(r, p4, label = "20-30%")
plt.plot(r, p5, label = "30-40%")
plt.plot(r, p6, label = "40-50%")
plt.plot(r, p7, label = "50-60%")
plt.plot(r, p8, label = "60-70%")
plt.plot(r, p9, label = "70-80%")
plt.plot(r, p10, label = "80-90%")
plt.plot(r, p11, label = "90-100%")

plt.yscale("log")

plt.xlabel("$r$ [fm]")
plt.ylabel("weighting function $W(r)$")
plt.legend(loc = "best")
plt.savefig("plots/W_centrality_log.pdf", format="pdf", bbox_inches="tight")



# Check if integral equals 1/2
p = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11]
for liste in p:
	integral = 0
	dr = r[1]-r[0]
	for i in range(0, len(liste)):
		integral += r[i]*liste[i]*dr
	print(integral)


