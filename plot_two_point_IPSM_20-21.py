import matplotlib.pyplot as plt
import numpy as np 

modes = [0, 1, 2, 3, 4]
#modes = [0, 1]
counter_fig = 0


#0-1
# N = 172*0+223
# s2_N = 0.*1.7**2+0
# s2_w = 0.*0.48**2+0.633**2

#20-21
N = 112*0+436
s2_N = 0*3.1**2+0
s2_w = 0*0.85**2+2.18**2

clm = np.loadtxt("output/clm.txt")
lMax = 10

trento_gauge = 1
# modulus plots
for mode in modes:
	counter_fig = counter_fig + 1
	plt.figure(counter_fig)
	plt.rcParams.update({'font.size': 20})
	from matplotlib.ticker import MaxNLocator
	ax = plt.figure().gca()
	ax.xaxis.set_major_locator(MaxNLocator(integer=True))
	ax.yaxis.set_major_locator(MaxNLocator(integer=True))
	plt.figure(figsize=(5.0,2.0))
	#import one-mode expectation values
	source_one_point = 'output/one_point_20-21.txt'
	# source_one_point = 'output/one_point_0-1.txt'
	one_points = np.loadtxt(source_one_point)
	l = np.zeros(lMax)
	y = np.zeros(lMax)
	y[0] = 1.0/np.pi/np.pi*(s2_w/N+s2_N/N/N)
	for i in range(0, lMax):
		l[i] = i+1
		if mode == 0:
			if i == 0:
				continue
			else:
				y[i] = 1.0/2./np.pi**2/N/clm[i-1, mode]*(1.+s2_w) 
		else:
			y[i] = (-1.0)**mode/2./np.pi**2/N/clm[i, mode]*(1.+s2_w) + (1-1./N+s2_N/N/N)*one_points[mode, i]*one_points[mode, i]

		# # print useful output
		# if mode == 0:
		# 	if i == 0:
		# 		continue
		# 	else:
		# 		print(y[i]*2*np.pi*np.pi*clm[i-1,mode])
				
	plt.scatter(l, y, label="IPSM", s=100, color = "orangered", marker= "+")
	plt.xlabel("l")
	plt.ylabel("$G_l^{(" + str(mode)  + ")}$")
	#plt.title("m = " +str(mode))
	# Trento prediction
	filename_trento = "output/two_point_random_20-21_m" + str(mode) + "_real.txt"
	filename_trento_error = "output/two_point_random_20-21_m" + str(mode) + "_real_error.txt"
	filename_trento_one = 'output/one_point_20-21' +'.txt'
	# filename_trento = "output/two_point_random_0-1_m" + str(mode) + "_real.txt"
	# filename_trento_error = "output/two_point_random_0-1_m" + str(mode) + "_real_error.txt"
	# filename_trento_one = 'output/one_point_0-1' +'.txt'

	filename_Saclay_simple = "Saclay_simplified/output/20-21/two_point_random_connected_m_" + str(mode) + ".txt"
	filename_Saclay = "Saclay/output/20-21/29/two_point_random_connected_m_" + str(mode) + ".txt"

	trento = np.loadtxt(filename_trento)
	trento_error = np.loadtxt(filename_trento_error)
	trento_one_ml = np.loadtxt(filename_trento_one)
	saclay_simple = np.loadtxt(filename_Saclay_simple)
	saclay = np.loadtxt(filename_Saclay)
	y_saclay_simple = np.zeros(lMax)
	y_saclay = np.zeros(lMax)
	y_trento = np.zeros(lMax)
	y_trento_error = np.zeros(lMax)
	for i in range(0, lMax):
		y_trento[i] = trento[i,i]
		y_saclay_simple[i] = saclay_simple[i,i]
		y_saclay[i] = saclay[i,i]
		if ((i == 0) & (mode == 0)):
			y_trento[i] -= trento_one_ml[mode, i] * trento_one_ml[mode, i]
		y_trento_error[i] =abs(trento_error[i,i]) # WHAT ABOUT m=0 ?
	plt.axhline(y=0, color ="black", linestyle="dashed", linewidth=1.0, zorder = 1)
	plt.scatter(l, y_trento, s=100, color = "green", label="TRENTo", marker= "x", zorder = 2)
	plt.scatter(l, y_saclay_simple, s=100, color = "blue", label="Generic", marker= "1", zorder = 2)
	plt.scatter(l, y_saclay, s=100, color = "purple", label="Large-Nc-Glasma", marker= "*", zorder = 2)
	plt.xticks([1,3,5,7,9])
	#plt.errorbar(l, y_trento, yerr = y_trento_error, linestyle = "none", elinewidth=2, capsize = 6, capthick = 2, color="green")
	#plt.legend(loc='best')
	# filename = "plots/two_point_modules_IPSM"  + "_m" + str(mode) + "_0-1" + ".pdf"
	filename = "plots/two_point_modules_IPSM"  + "_m" + str(mode) + "_20-21" + ".pdf"
	plt.savefig(filename, format='pdf', bbox_inches = "tight")
	plt.close(counter_fig)








