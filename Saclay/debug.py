import numpy as np 
import matplotlib.pyplot as plt 

#x, f = np.loadtxt("debug.txt", unpack = True)
r0001, W0001 = np.loadtxt("weight_functions_0-1.txt", unpack = True)
r1011, W1011 = np.loadtxt("weight_functions_10-11.txt", unpack = True)
r2021, W2021 = np.loadtxt("weight_functions_20-21.txt", unpack = True)

plt.figure(1)
plt.plot(r0001, W0001, label="0-1")
plt.plot(r1011, W1011, label ="10-11")
plt.plot(r2021, W2021, label ="20-21")
plt.legend(loc="best")
plt.savefig("W.pdf", format ="pdf", bbox_inches= "tight")
plt.close(1)

