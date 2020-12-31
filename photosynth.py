from math import *

# Vanilla max speed P one leaf
###############################################################################################

def kT(T, k, T0, H_a, R):
	return k*exp(-H_a/R * (1/T - 1/T0))

def fT(T, T0, H_d, S, R):
	ans=1 + exp( -H_d/R * (1/T0 - 1/ (H_d/S)))
	ans/=1 + exp(-H_d/R  * (1/T - 1/ (H_d/S)))
	return ans

def P_Max_T(kT, fT):
	return kT * fT

###############################################################################################

#	Light through beer
###############################################################################################

def getProportionBeer(K, LAI, m = 0.1):
	return  (K *  exp(- K * LAI))/ (1 - m)

def IT(I0, propBeer):
	return I0 * propBeer

def LT(L0, propBeer):
	return L0 *(1- propBeer)

###############################################################################################

# Max speed all leaves
###############################################################################################

def kTopen(LAI, kT):
	return LAI * kT

def P_Max_Topen(kTopen, fT):
	return P_Max_T(kTopen,fT)

def P_Max_LT(L, P_Max_T, P_MLT, L_0_5):
	return P_MLT *P_Max_T * L / (L + L_0_5)

###############################################################################################

