import numpy as np 
import matplotlib.pyplot as plt 

# r0001, y0001 = np.loadtxt("TwoPoint_0-1.txt", unpack = True)
# r1011, y1011 = np.loadtxt("TwoPoint_10-11.txt", unpack = True)
# r2021, y2021 = np.loadtxt("TwoPoint_20-21.txt", unpack = True)
# # r0001, W0001 = np.loadtxt("weight_functions_0-1.txt", unpack = True)
# # r1011, W1011 = np.loadtxt("weight_functions_10-11.txt", unpack = True)
# # r2021, W2021 = np.loadtxt("weight_functions_20-21.txt", unpack = True)

# # plt.figure(1)
# # plt.plot(r0001, W0001, label="0-1")
# # plt.plot(r1011, W1011, label ="10-11")
# # plt.plot(r2021, W2021, label ="20-21")
# # plt.legend(loc="best")
# # plt.savefig("W.pdf", format ="pdf", bbox_inches= "tight")
# # plt.close(1)

# plt.figure(1)
# #plt.plot(r0001, y0001, label="0-1")
# #plt.plot(r1011, y1011, label ="10-11")
# plt.plot(r2021, y2021, label ="20-21", linewidth=0.1)
# plt.xlim(0, np.amax(r2021))
# plt.legend(loc="best")
# plt.savefig("TwoPoint.pdf", format ="pdf", bbox_inches= "tight")
# plt.close(1)

m = 1
l = 9
rmax = 9.604-0.2

plt.figure(1)
TwoPoint = np.loadtxt("TwoPoint_mode-"+str(m)+".txt")
plt.imshow(TwoPoint, interpolation='none', cmap=plt.cm.Blues,extent=[0,rmax,rmax,0])
plt.xlabel("r2")
plt.ylabel("r1")
plt.colorbar()
plt.savefig("TwoPoint_mode-"+str(m)+".pdf", format="pdf", bbox_inches="tight")
plt.close(1)

plt.figure(2)
FB = np.loadtxt("FB_weight_ml_"+str(m)+"-"+str(l)+".txt")
maximal_value = np.amax(abs(FB))
plt.imshow(FB, interpolation='none', cmap=plt.cm.seismic, vmin = -maximal_value, vmax = maximal_value, extent=[0,rmax,rmax,0])
plt.xlabel("r2")
plt.ylabel("r1")
plt.colorbar()
plt.savefig("FB_weight_ml_"+str(m)+"-"+str(l)+".pdf", format="pdf", bbox_inches="tight")
plt.close(2)

plt.figure(3)
maximal_value = np.amax(abs(FB*TwoPoint))
plt.imshow(FB*TwoPoint, interpolation='none', cmap=plt.cm.seismic, vmin = -maximal_value, vmax = maximal_value, extent=[0,rmax,rmax,0])
plt.xlabel("r2")
plt.ylabel("r1")
#plt.colorbar()
plt.savefig("FB_TwoPoint_ml_"+str(m)+"-"+str(l)+".pdf", format="pdf", bbox_inches="tight")
plt.close(3)

plt.figure(4)
r_0001, integ_0001 = np.loadtxt("integ_r_0-1.txt", unpack=True)
r_2021, integ_2021 = np.loadtxt("integ_r_20-21.txt", unpack=True)
plt.plot(r_0001, integ_0001, label="0-1")
plt.plot(r_2021, integ_2021, label="20-21")
plt.legend(loc="best")
plt.savefig("integ_r.pdf", format="pdf", bbox_inches="tight")
plt.figure(4)

plt.figure(5)
r_0001, integrand_0001 = np.loadtxt("r_0-1.txt", unpack=True)
r_2021, integrand_2021 = np.loadtxt("r_20-21.txt", unpack=True)
plt.plot(r_0001, integrand_0001, label="0-1")
plt.plot(r_2021, integrand_2021, label="20-21")
plt.legend(loc="best")
plt.savefig("integrand_r.pdf", format="pdf", bbox_inches="tight")
plt.figure(5)

clm = np.loadtxt("../output/clm.txt")

print(m, np.sum(FB*TwoPoint)*rmax/(len(TwoPoint)-1)*rmax/(len(TwoPoint)-1)/clm[l-1,m]/clm[l-1,m])

