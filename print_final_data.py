import numpy as np 
import os

# Create directories if they do not exist yet
location = "Final_data/"

loc_bg = location + "BackgroundCoeffs/"
loc_bg_CGC = loc_bg + "CGC/"
loc_bg_Trento = loc_bg + "TrENTo/"
os.makedirs(loc_bg_Trento, exist_ok=True)
os.makedirs(loc_bg_CGC, exist_ok=True)

loc_bessel0 = location + "BesselZeros/"
os.makedirs(loc_bessel0, exist_ok=True)

loc_2mode = location + "ConnectedTwoMode/"
loc_2mode_LNc = loc_2mode + "CGCLargeNc/"
loc_2mode_IPSM = loc_2mode + "IPSM/"
loc_2mode_Magma = loc_2mode + "Magma/"
loc_2mode_Trento = loc_2mode + "TrENTo/"
os.makedirs(loc_2mode_LNc, exist_ok=True)
os.makedirs(loc_2mode_IPSM, exist_ok=True)
os.makedirs(loc_2mode_Magma, exist_ok=True)
os.makedirs(loc_2mode_Trento, exist_ok=True)

loc_impact = location + "ImpactParameterDistribution/"
os.makedirs(loc_impact, exist_ok=True)

loc_norm = location + "Normalizations/"
os.makedirs(loc_norm, exist_ok=True)

loc_W = location + "WeightFunctions/"
loc_W_CGC = loc_W + "CGC/"
loc_W_Trento = loc_W + "TrENTo/"
os.makedirs(loc_W_Trento, exist_ok=True)
os.makedirs(loc_W_CGC, exist_ok=True)


# useful definitions
centrality_classes = np.arange(90)
models = ["Trento", "IPSM", "CGCLargeNc", "Magma"]
modes = [0, 1, 2, 3, 4]
lMax = 10

# print normalizations
multiplicities = np.loadtxt("output/average_multiplicity.txt", unpack = True)
outfilename = loc_norm + "Normalizations" + ".txt" 
header = "Lower bound of centrality class, upper bound of centrality class, mean multiplicity"
np.savetxt(outfilename, np.array([centrality_classes, centrality_classes+1, multiplicities[0:len(centrality_classes)]]).transpose(), fmt= "%1.6g", header=header)

# print impact parameter distribution 
impact_params = np.loadtxt("output/percentiles_b.txt")
outfilename = loc_impact + "ImpactParameterDistribution" + ".txt" 
header = "Impact parameter distribution from 1e5 TrENTo events, normalized s.t. \\int p(b) db = 1, sliced up in percentiles: \n percentile edge, probability density"
np.savetxt(outfilename, impact_params, fmt= "%1.6g", header=header)

# print Bessel zeros 
Bessel_zero_ml = np.loadtxt("Saclay_simplified/output/0-1/bessel_d_0.txt")
outfilename = loc_bessel0 + "BesselDerivZeros_ml" + ".txt" 
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
	outfilename_Trento = loc_W_Trento + "WeightFunctions_Trento_" + centrality_string + ".txt"
	outfilename_CGC = loc_W_CGC + "WeightFunctions_CGC_" + centrality_string + ".txt"
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
	outfilename_Trento = loc_bg_Trento + "BackgroundCoeffs_Trento_" + centrality_string + ".txt"
	outfilename_CGC = loc_bg_CGC + "BackgroundCoeffs_CGC_" + centrality_string + ".txt"
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
				infilename = "output/two_point_random_" + centrality_string + "_m" + str(m) + "_real.txt"
				profile = np.loadtxt(infilename)
				if (m == 0):
					profile[0,0] -= 1./np.pi/np.pi
				outfilename = loc_2mode_Trento + "/ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
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
				outfilename = loc_2mode_IPSM  + "/ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in the IPSM, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)
			elif (model == models[2]): # Large Nc
				infilename = 'Saclay/output/'+centrality_string+'/Nr41/Nm64/m5.0e-3/two_point_random_connected' + '_m_' + str(m)  +'.txt'
				profile = np.loadtxt(infilename)
				# for i in range(0, lMax):
				# 	for j in range(0, lMax):
				# 		if ((i==0)&(j==0)&(m==0)):
				# 			profile[i][j] += 1./np.pi/np.pi*(sn2_N_N2)
				# 		else:
				# 			profile[i][j] += (1.+sn2_N_N2)*one_points[m, i]*one_points[m, j]
				outfilename = loc_2mode_LNc + "/ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
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
				outfilename = loc_2mode_Magma + "/ConnectedTwoMode_" + "m" + str(m) + "_" + model + "_" + centrality_string + ".txt"
				header = "Connected two-mode correlation functions G_{{l1,l2}}^{{(m,-m)}}, m = {0}, in the magma model, for ensembles with a randomized reaction plane angle.\nRows count l1 = 1,2,.. and columns count l2=1,2,....".format(m)
				np.savetxt(outfilename, profile, fmt= "%1.5e", header=header)









