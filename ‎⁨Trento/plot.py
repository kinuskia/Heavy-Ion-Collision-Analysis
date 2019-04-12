import matplotlib.pyplot as plt
import numpy as np 

# profile = np.loadtxt('002.dat')
# plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
# plt.savefig("profile002.pdf", format='pdf', bbox_inches = "tight")

for n in range(1,10):
	filename = "PbPb10000/000" + str(n)
	plt.figure(n)
	profile = np.loadtxt(filename+'.dat')
	plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
	plt.title("Profile " + str(n))
	plt.savefig("profile"+ str(n)+".pdf", format='pdf', bbox_inches = "tight")
	plt.close(n)