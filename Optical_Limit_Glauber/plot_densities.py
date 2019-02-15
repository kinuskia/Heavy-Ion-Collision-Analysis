import matplotlib.pyplot as plt 
import numpy as np 

rPb, rhoPb = np.loadtxt("Pb_density.txt", unpack=True)
rAu, rhoAu = np.loadtxt("Au_density.txt", unpack=True)
rCu, rhoCu = np.loadtxt("Cu_density.txt", unpack=True)

RPb = 6.62
RAu = 6.38
RCu = 4.2

# spline interpolation
from scipy.interpolate import interp1d
rnew = np.linspace(rPb.min(), rPb.max(), 300)
fPb = interp1d(rPb, rhoPb, kind='cubic')
fAu = interp1d(rAu, rhoAu, kind='cubic')
fCu = interp1d(rCu, rhoCu, kind='cubic')



plt.figure(1)
plt.plot(rnew, 4*np.pi*rnew*rnew*fPb(rnew), label = 'Pb-208')
plt.plot(rnew, 4*np.pi*rnew*rnew*fAu(rnew), label = 'Au-197')
plt.plot(rnew, 4*np.pi*rnew*rnew*fCu(rnew), label = 'Cu-63')
plt.legend(loc = 'best')
plt.title("Nucleon density profiles")

plt.xlabel("r [fm]")
plt.ylabel("$4\\pi r^2\\rho$")
#plt.legend(loc = "best")
plt.savefig("density_profile.pdf", format="pdf", bbox_inches="tight")

