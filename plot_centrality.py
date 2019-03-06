import matplotlib.pyplot as plt 
import numpy as np 

b, npart, mult = np.loadtxt("collision_specs.txt", unpack=True)


plt.figure(1)
plt.scatter(b, mult, s=1)

plt.xlabel("impact parameter [fm]")
plt.ylabel("multiplicity")
#plt.legend(loc = "best")
plt.savefig("b-mult.pdf", format="pdf", bbox_inches="tight")

plt.figure(2)
bins = 100
plt.hist(b, bins)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.savefig("b-hist.pdf", format="pdf", bbox_inches="tight")

plt.figure(3)
bins = 100
n, bins, patches = plt.hist(mult, bins)
plt.xlabel("multiplicity [arb. unit]")
plt.ylabel("# events")
plt.yscale("log")
plt.savefig("mult-hist.pdf", format="pdf", bbox_inches="tight")

plt.figure(4)
plt.figure(figsize=(10,3))
plt.rcParams.update({'font.size': 15})
x = (bins[1:]+bins[:-1])/2
dx = bins[1]-bins[0]
integral = np.sum(n)*dx
y = n/integral
d_n = np.sqrt(n)
d_y = d_n/integral 
centrality = np.zeros(len(x))
centrality[0] = 1
for i in range(0,len(centrality)-1):
	centrality[i+1] = centrality[i]-y[i]*dx
plt.plot(x, y-d_y, linewidth = 0.5, color = "black")
plt.plot(x, y+d_y, linewidth = 0.5, color = "black")
plt.fill_between(x, y+d_y, y-d_y, color = "black", label="$b = $0 fm")
plt.xlabel("multiplicity [arb. unit]")
plt.ylabel("probability density")
plt.yscale("log")
plt.savefig("mult-prob.pdf", format="pdf", bbox_inches="tight")
plt.figure(5)
plt.plot(centrality, x)
plt.ylabel("multiplicity")
plt.xlabel("centrality")
plt.savefig("centrality.pdf", format="pdf", bbox_inches="tight")

def centrality_class(p):
	return np.percentile(mult, 100-p)
bins2 = 32

plt.figure(6)
plt.hist(b[mult>centrality_class(5)], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 0-5%")
plt.savefig("b_class005.pdf", format="pdf", bbox_inches="tight")

plt.figure(7)
plt.hist(b[(mult>centrality_class(10)) & (mult<centrality_class(5))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 5-10%")
plt.savefig("b_class010.pdf", format="pdf", bbox_inches="tight")

plt.figure(8)
plt.hist(b[(mult>centrality_class(20)) & (mult<centrality_class(10))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 10-20%")
plt.savefig("b_class020.pdf", format="pdf", bbox_inches="tight")

plt.figure(9)
plt.hist(b[(mult>centrality_class(30)) & (mult<centrality_class(20))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 20-30%")
plt.savefig("b_class030.pdf", format="pdf", bbox_inches="tight")

plt.figure(10)
plt.hist(b[(mult>centrality_class(40)) & (mult<centrality_class(30))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 30-40%")
plt.savefig("b_class040.pdf", format="pdf", bbox_inches="tight")

plt.figure(11)
plt.hist(b[(mult>centrality_class(50)) & (mult<centrality_class(40))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 40-50%")
plt.savefig("b_class050.pdf", format="pdf", bbox_inches="tight")

plt.figure(12)
plt.hist(b[(mult>centrality_class(60)) & (mult<centrality_class(50))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 50-60%")
plt.savefig("b_class060.pdf", format="pdf", bbox_inches="tight")

plt.figure(13)
plt.hist(b[(mult>centrality_class(70)) & (mult<centrality_class(60))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 60-70%")
plt.savefig("b_class070.pdf", format="pdf", bbox_inches="tight")

plt.figure(14)
plt.hist(b[(mult>centrality_class(80)) & (mult<centrality_class(70))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 70-80%")
plt.savefig("b_class080.pdf", format="pdf", bbox_inches="tight")

plt.figure(15)
plt.hist(b[(mult>centrality_class(90)) & (mult<centrality_class(80))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 80-90%")
plt.savefig("b_class090.pdf", format="pdf", bbox_inches="tight")

plt.figure(16)
plt.hist(b[(mult>centrality_class(100)) & (mult<centrality_class(90))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 90-100%")
plt.savefig("b_class100.pdf", format="pdf", bbox_inches="tight")

plt.figure(17)
plt.hist(b[(mult>centrality_class(10))], bins2, alpha = 0.5, label="0-10%")
plt.hist(b[(mult>centrality_class(40)) & (mult<centrality_class(30))], bins2, alpha = 0.8, label="30-40%")
plt.hist(b[(mult>centrality_class(70)) & (mult<centrality_class(60))], bins2, alpha = 0.8, label="60-70%")
plt.hist(b[(mult>centrality_class(100)) & (mult<centrality_class(90))], bins2, alpha = 0.8, label="90-100%")
plt.xlabel("impact parameter [fm]")
plt.legend(loc='best')
plt.ylabel("# events")
plt.title("Impact parameter distribution within different centrality classes")
plt.savefig("b_class_comp.pdf", format="pdf", bbox_inches="tight")

def gauss(x, mu, s, C):
	return C*np.exp(-1./2*(x-mu)*(x-mu)/s/s)/np.sqrt(2*np.pi)/s

from scipy.optimize import curve_fit



plt.figure(18)
plt.figure(figsize=(10,3.5))
plt.rcParams.update({'font.size': 15})
x_values = np.linspace(0, 12, 300)

n_1, bins_1, patches_1 = plt.hist(b[(mult>centrality_class(5))], bins2, alpha = 0.5, label="0-5%", color = '#1f77b4')

popt_1, pcov_1 = curve_fit(gauss, (bins_1[1:]+bins_1[:-1])/2, n_1, p0 = [2.5, 1, 3000], sigma = np.sqrt(n_1))
#plt.plot(x_values, gauss(x_values, *popt_1), color = '#1f77b4')

n_2, bins_2, patches_2 = plt.hist(b[(mult>centrality_class(10)) & (mult<centrality_class(5))], bins2, alpha = 0.8, label="5-10%", color = '#ff7f0e')
popt_2, pcov_2 = curve_fit(gauss, (bins_2[1:]+bins_2[:-1])/2, n_2, p0 = [4.2, 0.8, 6000], sigma = np.sqrt(n_2))
#plt.plot(x_values, gauss(x_values, *popt_2), color = '#ff7f0e')


n_3, bins_3, patches_3 = plt.hist(b[(mult>centrality_class(20)) & (mult<centrality_class(10))], bins2, alpha = 0.8, label="10-20%", color = '#2ca02c')
popt_3, pcov_3 = curve_fit(gauss, (bins_3[1:]+bins_3[:-1])/2, n_3, p0 = [6, 1, 9000], sigma = np.sqrt(n_3))
#plt.plot(x_values, gauss(x_values, *popt_3), color = '#2ca02c')

for item in popt_3:
	print(item)

n_4, bins_4, patches_4 = plt.hist(b[(mult>centrality_class(30)) & (mult<centrality_class(20))], bins2, alpha = 0.8, label="20-30%", color = '#d62728')
popt_4, pcov_4 = curve_fit(gauss, (bins_4[1:]+bins_4[:-1])/2, n_4, p0 = [8, 1, 9000], sigma = np.sqrt(n_4))
#plt.plot(x_values, gauss(x_values, *popt_4), color = '#d62728')

n_5, bins_5, patches_5 = plt.hist(b[(mult>centrality_class(40)) & (mult<centrality_class(30))], bins2, alpha = 0.8, label="30-40%", color = '#9467bd')
popt_5, pcov_5 = curve_fit(gauss, (bins_5[1:]+bins_5[:-1])/2, n_5, p0 = [10, 1, 9000], sigma = np.sqrt(n_5))
#plt.plot(x_values, gauss(x_values, *popt_5), color = '#9467bd')

plt.xlabel("impact parameter [fm]")
plt.legend(loc='best')
plt.ylabel("# events")
plt.title("Impact parameter distribution within different centrality classes")
plt.savefig("b_class_comp_adj.pdf", format="pdf", bbox_inches="tight")

# # compute percentiles
# classes = [0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
# perc = (100-68.27)/2.0
# for i in range(0,len(classes)-1):
# 	print(classes[i], "-", classes[i+1])
# 	print(np.percentile(b[(mult>centrality_class(classes[i+1])) & (mult<centrality_class(i))], perc))
# 	print(np.percentile(b[(mult>centrality_class(classes[i+1])) & (mult<centrality_class(i))], 100-perc))

# print mean and standard deviations
print("Means and standard deviations")

print("0-5%")
print(np.mean(b[(mult>centrality_class(5))]))
print(np.std(b[(mult>centrality_class(5))]))

print("5-10%")
print(np.mean(b[(mult>centrality_class(10)) & (mult<centrality_class(5))]))
print(np.std(b[(mult>centrality_class(10)) & (mult<centrality_class(5))]))

print("10-20%")
print(np.mean(b[(mult>centrality_class(20)) & (mult<centrality_class(10))]))
print(np.std(b[(mult>centrality_class(20)) & (mult<centrality_class(10))]))

print("20-30%")
print(np.mean(b[(mult>centrality_class(30)) & (mult<centrality_class(20))]))
print(np.std(b[(mult>centrality_class(30)) & (mult<centrality_class(20))]))

print("30-40%")
print(np.mean(b[(mult>centrality_class(40)) & (mult<centrality_class(30))]))
print(np.std(b[(mult>centrality_class(40)) & (mult<centrality_class(30))]))

