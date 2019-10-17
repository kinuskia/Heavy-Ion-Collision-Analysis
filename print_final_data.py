import numpy as np 

# useful definitions
location = "Final_data/"
centrality_classes = np.arange(90)
models = ["Trento", "IPSM", "CGCLargeNc", "Magma"]
modes = [0, 1, 2, 3, 4]
lMax = 10

# print normalizations
multiplicities = np.loadtxt("output/average_multiplicity.txt", unpack = True)
outfilename = location + "Normalizations/" + "Normalizations" + ".txt" 
header = "Lower bound of centrality class, upper bound of centrality class, mean multiplicity"
np.savetxt(outfilename, np.array([centrality_classes, centrality_classes+1, multiplicities[0:len(centrality_classes)]]).transpose(), fmt= "%1.6g", header=header)

# print W functions
for c in range(0, len(centrality_classes)):
	centrality_string = str(centrality_classes[c]) + "-" + str(centrality_classes[c]+1)
	infilename = "output/weight_functions_" + centrality_string + ".txt"	
	W_r = np.loadtxt(infilename)
	outfilename = location + "WeightFunctions/" + "WeightFunctions_" + centrality_string + ".txt"
	header = "Weight functions W(r). Left column: r in fm, right column: W(r) in fm^(-2)"
	np.savetxt(outfilename, W_r, fmt= "%1.5e", header=header)


# print background coeffs
for c in range(0, len(centrality_classes)):
	centrality_string = str(centrality_classes[c]) + "-" + str(centrality_classes[c]+1)
	infilename = "output/one_point_" + centrality_string + ".txt"
	background_coeffs = np.loadtxt(infilename)
	outfilename = location + "BackgroundCoeffs/" + "BackgroundCoeffs_" + centrality_string + ".txt"
	header = "Mean Fourier-Bessel coefficients for ensembles with a fixed reaction plane angle of zero, \\bar{\epsilon}_l^{(m)}\nRows count azimuthal modes m = 0,1,2,... and columns count radial modes l = 1,2,3,..."
	np.savetxt(outfilename, background_coeffs, fmt= "%1.5e", header=header)

# print Connected two-mode functions
clm_old = np.loadtxt("output/clm.txt")
clm = np.zeros((lMax, len(modes)))
for l in range(0, lMax):
	for m in range(0, len(modes)):
		if ((m==0)):
			if (l==0):
				clm[l][m] = 1./2
			else:
				clm[l][m] = clm_old[l-1][m]
		else:
			clm[l][m] = clm_old[l][m]
for c in range(0, len(centrality_classes)):
	centrality_string = str(centrality_classes[c]) + "-" + str(centrality_classes[c]+1)
	source_one_point = 'output/one_point_'+centrality_string+'.txt'
	one_points = np.loadtxt(source_one_point)
	filename_fit = "../IPSM-Fit/output/" + centrality_string + ".txt"
	One_sw2_N, sn2_N_N2 = np.loadtxt(filename_fit, unpack=True)
	for m in modes:
		for model in models:
			if (model == models[0]): #Trento
				infilename = "output/two_point_random_" + centrality_string + "_m" + str(m) + "_real.txt"
				profile = np.loadtxt(infilename)
				if (m == 0):
					profile[0,0] -= 1./np.pi/np.pi
				outfilename = location + "ConnectedTwoMode/" + "ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in TrENTo, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)
			elif (model == models[1]): #IPSM
				profile = np.zeros((lMax, lMax))
				for i in range(0, lMax):
					for j in range(0, lMax):
						if ((i == 0)&(j==0)&(m == 0)):
							profile[i][j] = 1./np.pi/np.pi*(One_sw2_N+sn2_N_N2)
						elif ((i == j)):
							profile[i][j] = (-1.0)**m/2./np.pi**2/clm[i, m]*One_sw2_N + (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
						else:
							profile[i][j] = (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
				outfilename = location + "ConnectedTwoMode/" + "ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in the IPSM, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)
			elif (model == models[2]): # Large Nc
				infilename = 'Saclay/output/'+centrality_string+'/Nr41/Nm64/m1.4e-1/two_point_random_connected' + '_m_' + str(m)  +'.txt'
				profile = np.loadtxt(infilename)
				for i in range(0, lMax):
					for j in range(0, lMax):
						if ((i==0)&(j==0)&(m==0)):
							profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
						else:
							profile[i][j] += (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
				outfilename = location + "ConnectedTwoMode/" + "ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in the CGC large-Nc model, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)
			elif (model == models[3]): # magma
				infilename = 'Saclay_simplified/output/'+centrality_string+'/two_point_random_connected' + '_m_' + str(m)  +'.txt'
				profile = np.loadtxt(infilename)
				for i in range(0, lMax):
					for j in range(0, lMax):
						if ((i==0)&(j==0)&(m==0)):
							profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
						else:
							profile[i][j] += (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
				outfilename = location + "ConnectedTwoMode/" + "ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in the magma model, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)









