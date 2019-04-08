import matplotlib.pyplot as plt
import numpy as np 

modes = [0, 1, 2, 3, 4]
#modes = [0, 1]
counter_fig = 0

N = 200
k = 1


clm = np.loadtxt("output/clm.txt")
lMax = 5

trento_gauge = 1
# modulus plots
for mode in modes:
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.figure(figsize=(10,7.3))
	plt.rcParams.update({'font.size': 20})
	from matplotlib.ticker import MaxNLocator
	ax = plt.figure().gca()
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	#import one-mode expectation values
	source_one_point = 'output/one_point_0-5.txt'
	one_points = np.loadtxt(source_one_point)
	l = np.zeros(lMax)
	y = np.zeros(lMax)
	for i in range(0, lMax):
		l[i] = i+1
		y[i] = (-1)**mode/2./np.pi**2/N/clm[i, mode]*(1.+1./k) + (1-1./N)*one_points[mode, i]*one_points[mode, i]
	plt.scatter(l, y, label="IPSM", s=100, color = "orangered", marker= "+")
	plt.xlabel("l")
	plt.ylabel("$G_l^{(" + str(mode)  + ")}$")
	plt.title("m = " +str(mode))
	# Trento prediction
	filename_trento = "output/two_point_random_0-5_m" + str(mode) + "_real.txt"
	filename_trento_error = "output/two_point_random_0-5_m" + str(mode) + "_real_error.txt"
	trento = np.loadtxt(filename_trento)
	trento_error = np.loadtxt(filename_trento_error)
	y_trento = np.zeros(lMax)
	y_trento_error = np.zeros(lMax)
	for i in range(0, lMax):
		y_trento[i] = trento[i,i]
		y_trento_error[i] =abs(trento_error[i,i])
	plt.scatter(l, y_trento, s=100, color = "green", label="TRENTo", marker= "x")
	#plt.errorbar(l, y_trento, yerr = y_trento_error, linestyle = "none", elinewidth=2, capsize = 6, capthick = 2, color="green")
	plt.legend(loc='best')
	filename = "plots/two_point_modules_IPSM"  + "_m" + str(mode) + ".pdf"
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close()







