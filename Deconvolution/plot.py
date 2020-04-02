import numpy as np 
import matplotlib.pyplot as plt 

raw = np.loadtxt("profiles_averaged_20-21.txt")
convolved = np.loadtxt("outfile.txt")

plt.figure(1)
plt.imshow(raw, interpolation='none', cmap=plt.cm.Blues)
plt.savefig("input.pdf", format="pdf", bbox_inches="tight")
plt.close(1)

plt.figure(2)
plt.imshow(convolved, interpolation='none', cmap=plt.cm.Blues)
plt.savefig("output.pdf", format="pdf", bbox_inches="tight")
plt.close(2)