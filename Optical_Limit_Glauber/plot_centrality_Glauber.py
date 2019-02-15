import matplotlib.pyplot as plt 
import numpy as np 

b00, mult00 = np.loadtxt("mult-b_x000.txt", unpack=True)
b01, mult01 = np.loadtxt("mult-b_x01.txt", unpack=True)
b02, mult02 = np.loadtxt("mult-b_x02.txt", unpack=True)
b03, mult03 = np.loadtxt("mult-b_x03.txt", unpack=True)
b04, mult04 = np.loadtxt("mult-b_x04.txt", unpack=True)
b05, mult05 = np.loadtxt("mult-b_x05.txt", unpack=True)
b06, mult06 = np.loadtxt("mult-b_x06.txt", unpack=True)
b07, mult07 = np.loadtxt("mult-b_x07.txt", unpack=True)
b08, mult08 = np.loadtxt("mult-b_x08.txt", unpack=True)
b09, mult09 = np.loadtxt("mult-b_x09.txt", unpack=True)
b10, mult10 = np.loadtxt("mult-b_x010.txt", unpack=True)
b11, mult11 = np.loadtxt("mult-b_x11.txt", unpack=True)
b12, mult12 = np.loadtxt("mult-b_x12.txt", unpack=True)
b13, mult13 = np.loadtxt("mult-b_x13.txt", unpack=True)
b14, mult14 = np.loadtxt("mult-b_x14.txt", unpack=True)
b15, mult15 = np.loadtxt("mult-b_x15.txt", unpack=True)
b16, mult16 = np.loadtxt("mult-b_x16.txt", unpack=True)
b17, mult17 = np.loadtxt("mult-b_x17.txt", unpack=True)
b18, mult18 = np.loadtxt("mult-b_x18.txt", unpack=True)
b19, mult19 = np.loadtxt("mult-b_x19.txt", unpack=True)
b20, mult20 = np.loadtxt("mult-b_x020.txt", unpack=True)
b30, mult30 = np.loadtxt("mult-b_x030.txt", unpack=True)
b40, mult40 = np.loadtxt("mult-b_x040.txt", unpack=True)
b50, mult50 = np.loadtxt("mult-b_x050.txt", unpack=True)
b60, mult60 = np.loadtxt("mult-b_x060.txt", unpack=True)
b70, mult70 = np.loadtxt("mult-b_x070.txt", unpack=True)
b80, mult80 = np.loadtxt("mult-b_x080.txt", unpack=True)
b90, mult90 = np.loadtxt("mult-b_x090.txt", unpack=True)
b100, mult100 = np.loadtxt("mult-b_x100.txt", unpack=True)
#b90f, mult90f = np.loadtxt("mult-b_x90_float.txt", unpack=True)
b_Trento, npart_Trento, mult_trento = np.loadtxt("../output/collision_specs.txt", unpack=True)

# spline interpolation
from scipy.interpolate import interp1d
b = np.linspace(b00.min(), b00.max(), 300)
f00 = interp1d(b00, mult00, kind='cubic')
f01 = interp1d(b01, mult01, kind='cubic')
f02 = interp1d(b02, mult02, kind='cubic')
f03 = interp1d(b03, mult03, kind='cubic')
f04 = interp1d(b04, mult04, kind='cubic')
f05 = interp1d(b05, mult05, kind='cubic')
f06 = interp1d(b06, mult06, kind='cubic')
f07 = interp1d(b07, mult07, kind='cubic')
f08 = interp1d(b08, mult08, kind='cubic')
f09 = interp1d(b09, mult09, kind='cubic')
f10 = interp1d(b10, mult10, kind='cubic')
f11 = interp1d(b11, mult11, kind='cubic')
f12 = interp1d(b12, mult12, kind='cubic')
f13 = interp1d(b13, mult13, kind='cubic')
f14 = interp1d(b14, mult14, kind='cubic')
f15 = interp1d(b15, mult15, kind='cubic')
f16 = interp1d(b16, mult16, kind='cubic')
f17 = interp1d(b17, mult17, kind='cubic')
f18 = interp1d(b18, mult18, kind='cubic')
f19 = interp1d(b19, mult19, kind='cubic')
f20 = interp1d(b20, mult20, kind='cubic')
f30 = interp1d(b30, mult30, kind='cubic')
f40 = interp1d(b40, mult40, kind='cubic')
f50 = interp1d(b50, mult50, kind='cubic')
f60 = interp1d(b60, mult60, kind='cubic')
f70 = interp1d(b70, mult70, kind='cubic')
f80 = interp1d(b80, mult80, kind='cubic')
f90 = interp1d(b90, mult90, kind='cubic')
f100 = interp1d(b100, mult100, kind='cubic')


plt.figure(1)
plt.scatter(b_Trento, mult_trento/mult_trento[0]/1.45);
plt.plot(b, f00(b)/f00(b[0]), label="x=0.00")
#plt.plot(b, f01(b)/f01(b[0]), label="x=0.01")
plt.plot(b, f02(b)/f02(b[0]), label="x=0.02")
#plt.plot(b, f03(b)/f03(b[0]), label="x=0.03")
plt.plot(b, f04(b)/f04(b[0]), label="x=0.04")
#plt.plot(b, f05(b)/f05(b[0]), label="x=0.05")
plt.plot(b, f06(b)/f06(b[0]), label="x=0.06")
#plt.plot(b, f07(b)/f07(b[0]), label="x=0.07")
plt.plot(b, f08(b)/f08(b[0]), label="x=0.08")
#plt.plot(b, f09(b)/f09(b[0]), label="x=0.09")
plt.plot(b, f10(b)/f10(b[0]), label="x=0.1")
#plt.plot(b, f11(b)/f11(b[0]), label="x=0.11")
#plt.plot(b, f12(b)/f12(b[0]), label="x=0.12")
#plt.plot(b, f13(b)/f13(b[0]), label="x=0.13")
#plt.plot(b, f14(b)/f14(b[0]), label="x=0.14")
#plt.plot(b, f15(b)/f15(b[0]), label="x=0.15")
#plt.plot(b, f16(b)/f16(b[0]), label="x=0.16")
#plt.plot(b, f17(b)/f17(b[0]), label="x=0.17")
#plt.plot(b, f18(b)/f18(b[0]), label="x=0.18")
#plt.plot(b, f19(b)/f19(b[0]), label="x=0.19")
plt.plot(b, f20(b)/f20(b[0]), label="x=0.2")
#plt.plot(b, f30(b)/f30(b[0]), label="x=0.3")
plt.plot(b, f40(b)/f40(b[0]), label="x=0.4")
#plt.plot(b, f50(b)/f50(b[0]), label="x=0.5")
plt.plot(b, f60(b)/f60(b[0]), label="x=0.6")
#plt.plot(b, f70(b)/f70(b[0]), label="x=0.7")
plt.plot(b, f80(b)/f80(b[0]), label="x=0.8")
#plt.plot(b, f90(b)/f90(b[0]), label="x=0.9")
plt.plot(b, f100(b)/f100(b[0]), label="x=1.0")

#plt.plot(b90f, mult90f, label="x=0.9f")

plt.legend(loc="best")

plt.xlabel("impact parameter b [fm]")
plt.ylabel("$N_{ch}/N_{ch, max}$")
#plt.legend(loc = "best")
plt.savefig("b-mult_optical_Glauber.pdf", format="pdf", bbox_inches="tight")





