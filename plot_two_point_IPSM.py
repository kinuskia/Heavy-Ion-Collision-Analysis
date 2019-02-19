import matplotlib.pyplot as plt
import numpy as np 

modes = [0, 1, 2, 3, 4]
#modes = [0, 1]
counter_fig = 0

clm = np.loadtxt("output/clm.txt")
lMax = 5

trento_gauge = 1
# modulus plots
for mode in modes:
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.rcParams.update({'font.size': 15})
	from matplotlib.ticker import MaxNLocator
	ax = plt.figure().gca()
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	l = np.zeros(lMax)
	y = np.zeros(lMax)
	for i in range(0, lMax):
		l[i] = i+1
		y[i] = clm[0, 0]/clm[i, mode]
	plt.scatter(l, y, label="IPSM", s=30)
	plt.xlabel("l")
	plt.ylabel("$\\left\\|G_l^{(" + str(mode)  + ")}/G_1^{("+str(0)+")}\\right\\|$")
	plt.title("m = " +str(mode))
	# Trento prediction
	filename_trento = "output/two_point_0-5_m" + str(mode) + "_real.txt"
	filename_trento_error = "output/two_point_0-5_m" + str(mode) + "_real_error.txt"
	trento = np.loadtxt(filename_trento)
	trento_error = np.loadtxt(filename_trento_error)
	y_trento = np.zeros(lMax)
	y_trento_error = np.zeros(lMax)
	if mode == 0:
		trento_gauge = trento[0,0]
	for i in range(0, lMax):
		y_trento[i] = abs(trento[i,i]/trento_gauge)
		y_trento_error[i] =abs(trento_error[i,i])/trento_gauge
	plt.scatter(l, y_trento, s=10, color = "green")
	plt.errorbar(l, y_trento, yerr = y_trento_error, linestyle = "none", capsize = 6, color="green", label="TRENTo")
	plt.legend(loc='best')
	filename = "plots/two_point_modules_IPSM"  + "_m" + str(mode) + ".pdf"
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close()







