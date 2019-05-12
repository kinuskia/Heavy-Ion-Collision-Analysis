import numpy as np 
import matplotlib.pyplot as plt

N=100000
percent = 0.01
mean, err = np.loadtxt("output/background_coeffs.txt", unpack=True)

var = err*err*N*percent

var_err = var*np.sqrt(2/(N*percent-1))

x = np.arange(0, 21)

plt.figure(figsize=(10,5))
plt.rcParams.update({'font.size': 23})
plt.rcParams['axes.titlepad'] = 20
plt.errorbar(x, var, yerr=var_err, linestyle = 'none', elinewidth=2, capsize = 6, capthick = 2)
plt.ylabel("$\\left\\langle\\epsilon_0^{(0)}\\epsilon_0^{(0)} \\right\\rangle$")
plt.xlabel("$n$")
plt.title("Variance of $\\epsilon_0^{(0)}$ for the $n$-$(n+1)$% class")
plt.savefig("plots/background_coeffs.pdf", format="pdf", bbox_inches = "tight")

print("class", "mean", "s2", "s2err")
for i in range(0, len(mean)):
	print(x[i], mean[i], var[i], var_err[i])