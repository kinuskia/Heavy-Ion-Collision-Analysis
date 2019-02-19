import matplotlib.pyplot as plt 
import numpy as np 

# r0, e01, e02, e03 = np.loadtxt("output/bessel1.txt", unpack=True)





# plt.plot(r0, e01, linewidth = 0.5, color = "blue", label = "$m=1, l=1$")
# plt.plot(r0, e02, linewidth = 0.5, color = "orange", label = "$m=1, l=2$")
# plt.plot(r0, e03, linewidth = 0.5, color = "purple", label = "$m=1, l=3$")
# plt.title("$W(r) R^2 J_m(z_m^{(l)} r)$")
# plt.xlabel("$r$ in fm")
# plt.legend(loc="best")



# #plt.errorbar(r, e, yerr = de, linestyle = "none", capsize = 2 , color = "blue")


# plt.savefig("plots/bessel1.pdf", format="pdf", bbox_inches="tight")

# import Bessel zeros
zero_ml = np.loadtxt("output/bessel_d_0.txt")



from scipy.special import jv
x = np.linspace(0, 1, 200)
# q_lm
y10 = np.zeros(len(x))
y20 = np.zeros(len(x))
y30 = np.zeros(len(x))
y40 = np.zeros(len(x))

y11 = np.zeros(len(x))
y21 = np.zeros(len(x))
y31 = np.zeros(len(x))
y41 = np.zeros(len(x))

y12 = np.zeros(len(x))
y22 = np.zeros(len(x))
y32 = np.zeros(len(x))
y42 = np.zeros(len(x))

for i in range(0, len(x)):
	y10[i] = jv(0,zero_ml[0, 0]*x[i])
	y20[i] = jv(0,zero_ml[0, 1]*x[i])
	y30[i] = jv(0,zero_ml[0, 2]*x[i])
	y40[i] = jv(0,zero_ml[0, 3]*x[i])

	y11[i] = jv(1,zero_ml[1, 0]*x[i])
	y21[i] = jv(1,zero_ml[1, 1]*x[i])
	y31[i] = jv(1,zero_ml[1, 2]*x[i])
	y41[i] = jv(1,zero_ml[1, 3]*x[i])

	y12[i] = jv(2,zero_ml[2, 0]*x[i])
	y22[i] = jv(2,zero_ml[2, 1]*x[i])
	y32[i] = jv(2,zero_ml[2, 2]*x[i])
	y42[i] = jv(2,zero_ml[2, 3]*x[i])

plt.figure(1)
plt.plot(x, y10, label = "$J_0\\left(z_1^{(0)} \\rho\\right)$")
plt.plot(x, y20, label = "$J_0\\left(z_2^{(0)} \\rho\\right)$")
plt.plot(x, y30, label = "$J_0\\left(z_3^{(0)} \\rho\\right)$")
#plt.plot(x, y40, label = "$J_0\\left(z_4^{(0)} \\rho\\right)$")
plt.xlabel("$\\rho$")
plt.ylabel("$J_m\\left(z_l^{(m)} \\rho\\right)$")
plt.legend(loc = "best")
plt.savefig("plots/Bessel0.pdf", format = "pdf", bbox_inches = "tight")
plt.close(1)

plt.figure(2)
plt.plot(x, y11, label = "$J_1\\left(z_1^{(1)} \\rho\\right)$")
plt.plot(x, y21, label = "$J_1\\left(z_2^{(1)} \\rho\\right)$")
plt.plot(x, y31, label = "$J_1\\left(z_3^{(1)} \\rho\\right)$")
#plt.plot(x, y40, label = "$J_0\\left(z_4^{(0)} \\rho\\right)$")
plt.xlabel("$\\rho$")
plt.ylabel("$J_m\\left(z_l^{(m)} \\rho\\right)$")
plt.legend(loc = "best")
plt.savefig("plots/Bessel1.pdf", format = "pdf", bbox_inches = "tight")
plt.close(2)

plt.figure(3)
plt.plot(x, y12, label = "$J_2\\left(z_1^{(2)} \\rho\\right)$")
plt.plot(x, y22, label = "$J_2\\left(z_2^{(2)} \\rho\\right)$")
plt.plot(x, y32, label = "$J_2\\left(z_3^{(2)} \\rho\\right)$")
#plt.plot(x, y40, label = "$J_0\\left(z_4^{(0)} \\rho\\right)$")
plt.xlabel("$\\rho$")
plt.ylabel("$J_m\\left(z_l^{(m)} \\rho\\right)$")
plt.legend(loc = "best")
plt.savefig("plots/Bessel2.pdf", format = "pdf", bbox_inches = "tight")
plt.close(3)



