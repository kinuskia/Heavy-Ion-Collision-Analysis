import matplotlib.pyplot as plt 
import numpy as np 

b00, mult00 = np.loadtxt("mult-b_x00.txt", unpack=True)
b10, mult10 = np.loadtxt("mult-b_x10.txt", unpack=True)
b20, mult20 = np.loadtxt("mult-b_x20.txt", unpack=True)
b30, mult30 = np.loadtxt("mult-b_x30.txt", unpack=True)
b40, mult40 = np.loadtxt("mult-b_x40.txt", unpack=True)
b50, mult50 = np.loadtxt("mult-b_x50.txt", unpack=True)
b90, mult90 = np.loadtxt("mult-b_x90.txt", unpack=True)
b90f, mult90f = np.loadtxt("mult-b_x90_float.txt", unpack=True)



plt.figure(1)
plt.plot(b00, mult00, label="x=0.0")
plt.plot(b10, mult10, label="x=0.1")
plt.plot(b20, mult20, label="x=0.2")
plt.plot(b30, mult30, label="x=0.3")
plt.plot(b40, mult40, label="x=0.4")
plt.plot(b50, mult50, label="x=0.5")
plt.plot(b90, mult90, label="x=0.9")
plt.plot(b90f, mult90f, label="x=0.9f")

plt.legend(loc="best")

plt.xlabel("impact parameter [fm]")
plt.ylabel("multiplicity")
#plt.legend(loc = "best")
plt.savefig("b-mult_optical_Glauber.pdf", format="pdf", bbox_inches="tight")

