import matplotlib.pyplot as plt
import numpy as np 

r, p1, p2, p3, p4, p5 = np.loadtxt("output/rhos.txt", unpack = True)


plt.figure(2)

plt.rcParams.update({'font.size': 23})

plt.figure(figsize=(10,5))

plt.plot(r, p1, label = "0-5%")
plt.plot(r, p2, label = "5-10%")
plt.plot(r, p3, label = "10-20%")
plt.plot(r, p4, label = "20-30%")
plt.plot(r, p5, label = "30-40%")
# plt.plot(r, p6, label = "40-50%")
# plt.plot(r, p7, label = "50-60%")
# plt.plot(r, p8, label = "60-70%")
# plt.plot(r, p9, label = "70-80%")
# plt.plot(r, p10, label = "80-90%")
# plt.plot(r, p11, label = "90-100%")

plt.xlabel("$r$ [fm]")
plt.ylabel("$\\rho(r)$")
plt.title("Mapping function $\\rho(r)$")
plt.legend(loc = "best")
plt.savefig("plots/rhos.pdf", format="pdf", bbox_inches="tight")



