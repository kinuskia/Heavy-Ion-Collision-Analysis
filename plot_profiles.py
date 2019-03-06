import matplotlib.pyplot as plt
import numpy as np 

# profile = np.loadtxt('002.dat')
# plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
# plt.savefig("profile002.pdf", format='pdf', bbox_inches = "tight")

plt.figure(1)
filename = "003.dat"
plt.xticks(np.arange(0, 100, step=20))
plt.figure(figsize=(10,8))
plt.rcParams.update({'font.size': 23})
profile = np.loadtxt(filename)
plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
plt.title("initial entropy density")
#plt.colorbar()
plt.xlabel("$x$ [0.1 fm]")
plt.ylabel("$y$ [0.1 fm]")
plt.savefig("plots/initial_profile1.pdf", format='pdf', bbox_inches = "tight")
plt.close(1)

plt.figure(2)
filename = '001.dat'
plt.xticks(np.arange(0, 100, step=20))
plt.figure(figsize=(10,8))
plt.rcParams.update({'font.size': 23})
profile = np.loadtxt(filename)
plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
plt.title("initial entropy density")
#plt.colorbar()
plt.xlabel("$x$ [0.1 fm]")
plt.ylabel("$y$ [0.1 fm]")
plt.savefig("plots/initial_profile2.pdf", format='pdf', bbox_inches = "tight")
plt.close(2)

plt.figure(3)
filename = '000.dat'
plt.xticks(np.arange(0, 100, step=20))
plt.figure(figsize=(10,8))
plt.rcParams.update({'font.size': 23})
profile = np.loadtxt(filename)
plt.imshow(profile, interpolation='none', cmap=plt.cm.Blues)
plt.title("initial entropy density")
#plt.colorbar()
plt.xlabel("$x$ [0.1 fm]")
plt.ylabel("$y$ [0.1 fm]")
plt.savefig("plots/initial_profile3.pdf", format='pdf', bbox_inches = "tight")
plt.close(3)