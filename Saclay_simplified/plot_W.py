import numpy as np 
import matplotlib.pyplot as plt 

r, W = np.loadtxt("weight_functions.txt", unpack=True)

def gauss(x, s):
	return 1./2./s/s*np.exp(-x*x/2./s/s)

from scipy.optimize import curve_fit

popt, pcov = curve_fit(gauss, r[r<6], W[r<6], p0=[3])

plt.figure(1)
plt.plot(r, W)
x = np.linspace(np.min(r), np.max(r), 300)
plt.plot(x, gauss(x, *popt))
plt.savefig("W.pdf", format = "pdf", bbox_inches="tight")
plt.close(1)

print("s: ", popt[0], "+-", np.sqrt(pcov[0][0]))

plt.figure(2)
plt.plot(r, W**(1./4))
x = np.linspace(np.min(r), np.max(r), 300)
plt.savefig("W_root.pdf", format = "pdf", bbox_inches="tight")
plt.close(2)

def poly(x, a, b, c, d, e):
	return -np.log(2.*popt[0]*popt[0])+a+b*x+(c-1./2./popt[0]/popt[0])*x*x+d*x*x*x+e*x*x*x*x

popt2, pcov2 = curve_fit(poly, r[r<6], np.log(W[r<6]))

plt.figure(3)
plt.plot(r, np.log(W), label = "Data")
x = np.linspace(np.min(r), np.max(r), 300)
plt.plot(x, poly(x, popt2[0], popt2[1], popt2[2]*1, popt2[3]*1, popt2[4]*1), label = "Fit")
plt.savefig("W_log.pdf", format = "pdf", bbox_inches="tight")
plt.legend(loc = "best")
plt.close(3)

print("a: ", popt2[0], "+-", np.sqrt(pcov2[0][0]))
print("b: ", popt2[1], "+-", np.sqrt(pcov2[1][1]))
print("c: ", popt2[2], "+-", np.sqrt(pcov2[2][2]))
print("d: ", popt2[3], "+-", np.sqrt(pcov2[3][3]))
print("e: ", popt2[4], "+-", np.sqrt(pcov2[4][4]))




