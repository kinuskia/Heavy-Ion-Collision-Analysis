import numpy as np 
import matplotlib.pyplot as plt 

r1, rho = np.loadtxt("rho.txt", unpack = True)
# r2, W = np.loadtxt("W.txt", unpack = True)

plt.figure(1)
plt.plot(r1, rho)
plt.savefig("rho.pdf", format="pdf", bbox_inches = "tight")
plt.close(1)

# plt.figure(2)
# plt.plot(r2, W)
# plt.savefig("W.pdf", format="pdf", bbox_inches = "tight")
# plt.close(2)

#profile = np.loadtxt("OnePoint_Profile_ml_10-10.txt", unpack = True)
profile = np.loadtxt("../output/profiles_averaged_20-21.txt")

plt.figure(3)
plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
plt.xlabel("x")
plt.ylabel("y")
plt.colorbar()
plt.savefig("OnePoint_Profile.pdf", format="pdf", bbox_inches="tight")
plt.close(3)

plt.figure(4)
x, f = np.loadtxt("OnePoint_cut.txt", unpack = True)
plt.plot(x, f, label="Background")
x4_10, f4_10 = np.loadtxt("OnePoint_cut_ml_4-10.txt", unpack = True)
plt.plot(x4_10, f4_10, label = "$m=4, l=10$")
x10_10, f10_10 = np.loadtxt("OnePoint_cut_ml_10-10.txt", unpack = True)
plt.plot(x10_10, f10_10, label = "$m=10, l=10$")
plt.legend(loc="best")
plt.savefig("OnePoint_cut.pdf", format = "pdf", unpack = True)
plt.figure(4)