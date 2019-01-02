import matplotlib.pyplot as plt 
import numpy as np 

b, npart, mult = np.loadtxt("output/collision_specs.txt", unpack=True)


plt.figure(1)
plt.scatter(b, mult, s=1)

plt.xlabel("impact parameter [fm]")
plt.ylabel("multiplicity")
#plt.legend(loc = "best")
plt.savefig("plots/b-mult.pdf", format="pdf", bbox_inches="tight")

plt.figure(2)
bins = 100
plt.hist(b, bins)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.savefig("plots/b-hist.pdf", format="pdf", bbox_inches="tight")

plt.figure(3)
bins = 100
n, bins, patches = plt.hist(mult, bins)
plt.xlabel("multiplicity [fm]")
plt.ylabel("# events")
plt.yscale("log")
plt.savefig("plots/mult-hist.pdf", format="pdf", bbox_inches="tight")

plt.figure(4)
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
plt.plot(x, y-d_y, linewidth = 0.5, color = "blue")
plt.plot(x, y+d_y, linewidth = 0.5, color = "blue")
plt.fill_between(x, y+d_y, y-d_y, color = "blue", label="$b = $0 fm")
plt.xlabel("multiplicity [fm]")
plt.ylabel("probability density")
plt.yscale("log")
plt.savefig("plots/mult-prob.pdf", format="pdf", bbox_inches="tight")
plt.figure(5)
plt.plot(centrality, x)
plt.ylabel("multiplicity")
plt.xlabel("centrality")
plt.savefig("plots/centrality.pdf", format="pdf", bbox_inches="tight")

def centrality_class(p):
	return np.percentile(mult, 100-p)
bins2 = 64

plt.figure(6)
plt.hist(b[mult>centrality_class(5)], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 0-5%")
plt.savefig("plots/b_class005.pdf", format="pdf", bbox_inches="tight")

plt.figure(7)
plt.hist(b[(mult>centrality_class(10)) & (mult<centrality_class(5))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 5-10%")
plt.savefig("plots/b_class010.pdf", format="pdf", bbox_inches="tight")

plt.figure(8)
plt.hist(b[(mult>centrality_class(20)) & (mult<centrality_class(10))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 10-20%")
plt.savefig("plots/b_class020.pdf", format="pdf", bbox_inches="tight")

plt.figure(9)
plt.hist(b[(mult>centrality_class(30)) & (mult<centrality_class(20))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 20-30%")
plt.savefig("plots/b_class030.pdf", format="pdf", bbox_inches="tight")

plt.figure(10)
plt.hist(b[(mult>centrality_class(40)) & (mult<centrality_class(30))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 30-40%")
plt.savefig("plots/b_class040.pdf", format="pdf", bbox_inches="tight")

plt.figure(11)
plt.hist(b[(mult>centrality_class(50)) & (mult<centrality_class(40))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 40-50%")
plt.savefig("plots/b_class050.pdf", format="pdf", bbox_inches="tight")

plt.figure(12)
plt.hist(b[(mult>centrality_class(60)) & (mult<centrality_class(50))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 50-60%")
plt.savefig("plots/b_class060.pdf", format="pdf", bbox_inches="tight")

plt.figure(13)
plt.hist(b[(mult>centrality_class(70)) & (mult<centrality_class(60))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 60-70%")
plt.savefig("plots/b_class070.pdf", format="pdf", bbox_inches="tight")

plt.figure(14)
plt.hist(b[(mult>centrality_class(80)) & (mult<centrality_class(70))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 70-80%")
plt.savefig("plots/b_class080.pdf", format="pdf", bbox_inches="tight")

plt.figure(15)
plt.hist(b[(mult>centrality_class(90)) & (mult<centrality_class(80))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 80-90%")
plt.savefig("plots/b_class090.pdf", format="pdf", bbox_inches="tight")

plt.figure(16)
plt.hist(b[(mult>centrality_class(100)) & (mult<centrality_class(90))], bins2)
plt.xlabel("impact parameter [fm]")
plt.ylabel("# events")
plt.title("centrality class 90-100%")
plt.savefig("plots/b_class100.pdf", format="pdf", bbox_inches="tight")