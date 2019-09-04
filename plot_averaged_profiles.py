import matplotlib.pyplot as plt
import numpy as np 

# profile = np.loadtxt('002.dat')
# plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
# plt.savefig("profile002.pdf", format='pdf', bbox_inches = "tight")

percentiles = [0, 20]

for i in range(0, len(percentiles)):
	filename = "output/profiles_averaged_" + str(percentiles[i]) + "-" + str(percentiles[i]+1)
	plt.figure(i)
	plt.xticks(np.arange(0, 100, step=20))
	plt.figure(figsize=(10,8))
	plt.rcParams.update({'font.size': 23})
	profile = np.loadtxt(filename + ".txt")
	plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
	plt.title("Profile " + str(percentiles[i]) + "-" + str(percentiles[i]+1) + "%")
	plt.colorbar()
	plt.xlabel("$x$ [0.1 fm]")
	plt.ylabel("$y$ [0.1 fm]")
	plt.savefig("plots/profile_average_" + str(percentiles[i]) + "-" + str(percentiles[i]+1) + ".pdf", format='pdf', bbox_inches = "tight")
	plt.close(i)