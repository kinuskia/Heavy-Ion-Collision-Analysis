import numpy as np 
import matplotlib.pyplot as plt 

#x, f = np.loadtxt("debug.txt", unpack = True)
x, f, g = np.loadtxt("debug.txt", unpack = True)

plt.figure(1)
plt.plot(x, f)
plt.savefig("debug.pdf", format ="pdf", bbox_inches= "tight")
plt.close(1)

plt.figure(2)
plt.plot(x, g)
plt.savefig("debug_One.pdf", format ="pdf", bbox_inches= "tight")
plt.close(2)


# plt.figure(1)
# plt.plot(x, f)
# plt.savefig("debug_Bessel.pdf", format ="pdf", bbox_inches= "tight")
# plt.close(1)

# plt.figure(2)
# plt.plot(x, g)
# plt.savefig("debug_log.pdf", format ="pdf", bbox_inches= "tight")
# plt.close(2)