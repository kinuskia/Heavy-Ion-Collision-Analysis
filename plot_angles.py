import matplotlib.pyplot as plt
import numpy as np 

angles = np.loadtxt("output/angles.txt", unpack = True)

plt.figure(1)
plt.hist(angles)
plt.xlabel("angle in rad")
plt.ylabel("#")


plt.savefig("plots/reaction_plane_angles.pdf", format="pdf", bbox_inches="tight")



