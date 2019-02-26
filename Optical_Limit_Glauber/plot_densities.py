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
plt.figure(figsize=(10,5))
plt.rcParams.update({'font.size': 15})
plt.plot(rnew, fPb(rnew)/fPb(rnew[0]), label = '$^{208}$Pb ($R=6.624\\,$fm, $a=0.549\\,$fm)')
plt.plot(rnew, fAu(rnew)/fAu(rnew[0]), label = '$^{197}$Au\\ ($R=6.38\\,$fm, $a=0.535\\,$fm)')
plt.plot(rnew, fCu(rnew)/fCu(rnew[0]), label = '$^{63}$Cu\\ ($R=4.2\\,$fm, $a=0.596\\,$fm)')
plt.legend(loc = 'best')
plt.title("Nucleus density profiles")

plt.xlabel("r [fm]")
plt.ylabel("$\\rho/\\rho_0$")
#plt.legend(loc = "best")
plt.savefig("density_profile.pdf", format="pdf", bbox_inches="tight")

