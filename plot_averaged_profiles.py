import matplotlib.pyplot as plt
import numpy as np 

# profile = np.loadtxt('002.dat')
# plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
# plt.savefig("profile002.pdf", format='pdf', bbox_inches = "tight")

percentiles = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

for i in range(0, len(percentiles)-1):
	filename = "output/profiles_averaged_" + str(percentiles[i]) + "-" + str(percentiles[i+1])
	plt.figure(i)
	profile = np.loadtxt(filename + ".txt")
	plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
	plt.title("Profile " + str(percentiles[i]) + "-" + str(percentiles[i+1]) + "%")
	plt.colorbar()
	plt.xlabel("$x$ [0.1 fm]")
	plt.ylabel("$y$ [0.1 fm]")
	plt.savefig("plots/profile_average_" + str(percentiles[i]) + "-" + str(percentiles[i+1]) + ".pdf", format='pdf', bbox_inches = "tight")
	plt.close(i)