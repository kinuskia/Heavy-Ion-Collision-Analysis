import numpy as np 

# useful definitions
location = "Final_data/"
centrality_classes = np.arange(90)
models = ["Trento", "IPSM", "CGCLargeNc", "Magma"]
modes = [0, 1, 2, 3, 4, 5, 6]
lMax = 20

# print normalizations
multiplicities = np.loadtxt("output/average_multiplicity.txt", unpack = True)
outfilename = location + "Normalizations/" + "Normalizations" + ".txt" 
header = "Lower bound of centrality class, upper bound of centrality class, mean multiplicity"
np.savetxt(outfilename, np.array([centrality_classes, centrality_classes+1, multiplicities[0:len(centrality_classes)]]).transpose(), fmt= "%1.6g", header=header)

# print impact parameter distribution 
impact_params = np.loadtxt("output/percentiles_b.txt")
outfilename = location + "ImpactParameterDistribution/" + "ImpactParameterDistribution" + ".txt" 
header = "Impact parameter distribution from 1e5 TrENTo events, normalized s.t. \\int p(b) db = 1, sliced up in percentiles: \n percentile edge, probability density"
np.savetxt(outfilename, impact_params, fmt= "%1.6g", header=header)

# print Bessel zeros 
Bessel_zero_ml = np.loadtxt("Saclay_simplified/output/0-1/bessel_d_0.txt")
outfilename = location + "BesselZeros/" + "BesselDerivZeros_ml" + ".txt" 
header = "First l non-negative solutions of J'_m(z) = 0. Rows: m, Columns: l."
Bessel_zero_ml_out = Bessel_zero_ml
# Add z=0
Bessel_zero_ml_out[0,1:] = Bessel_zero_ml[0,0:-1]
Bessel_zero_ml_out[0,0] = 0
np.savetxt(outfilename, Bessel_zero_ml_out, fmt= "%1.12g", header=header)

# print W functions, for TRENTO and CGC
for c in range(0, len(centrality_classes)):
	centrality_string = str(centrality_classes[c]) + "-" + str(centrality_classes[c]+1)
	infilename_Trento = "output/weight_functions_" + centrality_string + ".txt"
	infilename_CGC = "Saclay_simplified/output/"+centrality_string+"/weight_functions_magma_"+centrality_string+".txt"	
	W_r_Trento = np.loadtxt(infilename_Trento)
	W_r_CGC = np.loadtxt(infilename_CGC)
	outfilename_Trento = location + "WeightFunctions/TrENTo/" + "WeightFunctions_Trento_" + centrality_string + ".txt"
	outfilename_CGC = location + "WeightFunctions/CGC/" + "WeightFunctions_CGC_" + centrality_string + ".txt"
	header_Trento = "Weight functions W(r), obtained from TrENTo with 1e5 PbPb collisions at 2.76TeV. Left column: r in fm, right column: W(r) in fm^(-2)"
	header_CGC = "Weight functions W(r), obtained from Woods-Saxon 1-Nucleus thickness functions of Pb. Left column: r in fm, right column: W(r) in fm^(-2)"
	np.savetxt(outfilename_Trento, W_r_Trento, fmt= "%1.5e", header=header_Trento)
	np.savetxt(outfilename_CGC, W_r_CGC, fmt= "%1.5e", header=header_CGC)


# print background coeffs, for TRENTO and CGC
for c in range(0, len(centrality_classes)):
	centrality_string = str(centrality_classes[c]) + "-" + str(centrality_classes[c]+1)
	infilename_Trento = "output/one_point_" + centrality_string + ".txt"
	infilename_CGC = "Saclay_simplified/output/"+centrality_string+"/background_coeffs_CGC_" + centrality_string + ".txt"
	background_coeffs_Trento = np.loadtxt(infilename_Trento)
	background_coeffs_CGC = np.loadtxt(infilename_CGC)
	outfilename_Trento = location + "BackgroundCoeffs/TrENTo/" + "BackgroundCoeffs_Trento_" + centrality_string + ".txt"
	outfilename_CGC = location + "BackgroundCoeffs/CGC/" + "BackgroundCoeffs_CGC_" + centrality_string + ".txt"
	header_Trento = "Mean Fourier-Bessel coefficients for ensembles with a fixed reaction plane angle of zero, \\bar{\epsilon}_l^{(m)}, obtained from TrENTo with 1e5 PbPb collisions at 2.76TeV\nRows count azimuthal modes m = 0,1,2,... and columns count radial modes l = 1,2,3,..."
	header_CGC = "Mean Fourier-Bessel coefficients for ensembles with a fixed reaction plane angle of zero, \\bar{\epsilon}_l^{(m)}, obtained from Woods-Saxon 1-Nucleus thickness functions of Pb\nRows count azimuthal modes m = 0,1,2,... and columns count radial modes l = 1,2,3,..."
	np.savetxt(outfilename_Trento, background_coeffs_Trento, fmt= "%1.5e", header=header_Trento)
	np.savetxt(outfilename_CGC, background_coeffs_CGC, fmt= "%1.5e", header=header_CGC)

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
				if (m > 4):
					continue
				infilename = "output/two_point_random_" + centrality_string + "_m" + str(m) + "_real.txt"
				profile = np.loadtxt(infilename)
				if (m == 0):
					profile[0,0] -= 1./np.pi/np.pi
				outfilename = location + "ConnectedTwoMode/" + model + "/ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in TrENTo, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)
			elif (model == models[1]): #IPSM
				if ((m > 4) | (c > 10)):
					continue
				profile = np.zeros((lMax, lMax))
				for i in range(0, lMax):
					for j in range(0, lMax):
						if ((i == 0)&(j==0)&(m == 0)):
							profile[i][j] = 1./np.pi/np.pi*(One_sw2_N+sn2_N_N2)
						elif ((i == j)):
							profile[i][j] = (-1.0)**m/2./np.pi**2/clm[i, m]*One_sw2_N + (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
						else:
							profile[i][j] = (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
				outfilename = location+ "ConnectedTwoMode/" + model   + "/ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in the IPSM, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)
			elif (model == models[2]): # Large Nc
				if ((c != 0) & (c != 20)):
					continue
				infilename = 'Saclay/output/'+centrality_string+'/Nr41/Nm64/m5.0e-3/two_point_random_connected' + '_m_' + str(m)  +'.txt'
				profile = np.loadtxt(infilename)
				# for i in range(0, lMax):
				# 	for j in range(0, lMax):
				# 		if ((i==0)&(j==0)&(m==0)):
				# 			profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
				# 		else:
				# 			profile[i][j] += (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
				outfilename = location + "ConnectedTwoMode/" +model + "/ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in the CGC large-Nc model, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)
			elif (model == models[3]): # magma
				infilename = 'Saclay_simplified/output/'+centrality_string+'/two_point_random_connected' + '_m_' + str(m)  +'.txt'
				profile = np.loadtxt(infilename)
				# for i in range(0, lMax):
				# 	for j in range(0, lMax):
				# 		if ((i==0)&(j==0)&(m==0)):
				# 			profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
				# 		else:
				# 			profile[i][j] += (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
				outfilename = location + "ConnectedTwoMode/" + model + "/ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in the magma model, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)









