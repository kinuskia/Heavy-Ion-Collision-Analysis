import matplotlib.pyplot as plt 
import numpy as np 

bPb, sigmaPb = np.loadtxt("d_sigma_d_b.txt", unpack=True)



# spline interpolation
from scipy.interpolate import interp1d
bnew = np.linspace(bPb.min(), bPb.max(), 300)
fPb = interp1d(bPb, sigmaPb, kind='cubic')




plt.figure(1)
#plt.figure(figsize=(10,5))
#plt.rcParams.update({'font.size': 15})
plt.plot(bnew, fPb(bnew), label = '$^{208}$Pb')
plt.legend(loc = 'best')
plt.title("$\\sigma(b)$")

plt.xlabel("b [fm]")
plt.ylabel("d^2$\\sigma/$d$b^2$")
#plt.legend(loc = "best")
plt.savefig("cross_section.pdf", format="pdf", bbox_inches="tight")
plt.close(1)

plt.figure(2)
#plt.figure(figsize=(10,5))
#plt.rcParams.update({'font.size': 15})
plt.plot(bnew, 2*np.pi*bnew*fPb(bnew), label = '$^{208}$Pb')
plt.legend(loc = 'best')
plt.title("Impact parameter distribution")

plt.xlabel("b [fm]")
plt.ylabel("$2\\pi b\\,$d$^2$$\\sigma/$d$b^2$")
#plt.legend(loc = "best")
plt.savefig("b_dist.pdf", format="pdf", bbox_inches="tight")
plt.close(2)

# Get percentiles
b_int = np.linspace(0., 12., 100)
y_int = 2*np.pi*b_int*fPb(b_int)
h = (b_int[1]-b_int[0])
norm = 0
# trapezoidal rule
for i in range(0, len(b_int)):
	norm += h*y_int[i]
# I added too much for i = 0 and i = len(b_int):
norm -= y_int[0]*h/2 
norm -= y_int[-1]*h/2 

print ("Pb-Pb total cross section [barn]: ", norm/100)
p0_5 = 0
cumulant = y_int[0]*h/2
for i in range(0, len(b_int)):
	cumulant += h*y_int[i]
	if (cumulant/norm > 0.05):
		p0_5 = b_int[i]
		break;

print("b_max for 0-5%:", p0_5)




