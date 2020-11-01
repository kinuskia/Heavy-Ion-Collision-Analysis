import numpy as np 
import matplotlib.pyplot as plt 

mult, n_mult, n_mult_error = np.loadtxt("output/mult_hist.txt", unpack = True)

plt.figure(1)
plt.xlabel("Multiplicity $m$ [arb. unit]")
plt.ylabel("$p(m)$")
plt.plot(mult, n_mult, color = "black")
plt.errorbar(mult, n_mult, yerr= n_mult_error, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1, color = "blue")
plt.yscale('log')
plt.savefig("plots/mult_hist.pdf", format='pdf', bbox_inches = "tight")
plt.close(1)

b, n_b, n_b_error = np.loadtxt("output/b_hist.txt", unpack = True)
plt.figure(2)
plt.xlabel("impact parameter (fm)")
plt.ylabel("$p(b)$ (1/fm)")
plt.plot(b, n_b, color = "black")
plt.errorbar(b, n_b, yerr= n_b_error, linestyle = 'none', elinewidth=1, capsize = 3, capthick = 1, color = "blue")
plt.savefig("plots/b_hist.pdf", format='pdf', bbox_inches = "tight")
plt.close(2)

# check integral
integral = 0
width = mult[1]-mult[0]
for i in range(0, len(mult)):
	integral += width*n_mult[i]
print("integral multiplicity: ", integral)