from math import *

# Global constants (apply to all model)
g = 9.81
alpha = 0.385
J_MAX_25_Leaf = 210
Theta = 0.7
n_CO2_Air_Stom = 0.67
C_Gamma = 1.7
T_25_K = 298.15
H_J_Pot = 220000
S_J_Pot = 710
E_j = 37000
R_gas = 8.314
R_S_min = 82.0

# Local constants (apply to specific greenhouses)
A_Flr = 4000
C_d = 0.5


# Mutable constants
h_C_Buf = 1
PAR_Can = 100

# Unknown constants/variable
CO2_Out = 0		# meteo.csv dataset?

###############################################################################################

###############################################################################################

# Controller (Blower, Pad, Vent,...) (0<=U<=1)
# Initiate here to avoid compile error only, please edit at #!
U_Blow = 0.5		#!
U_Ext_CO2 = 0.5		#!
U_Pad = 0.5			#!
U_Roof = 0.5		#!
U_Side = 0.5		#!
U_ThScr = 0.5		#!
U_Vent_Forced = 0.5	#!

# Capacity (Blower, External air, Pad,...)
# Initiate here to avoid compile error only, please edit at #!
P_Blow = 7			#!
phi_Ext_CO2 = 7		#!
phi_Pad	= 7			#!
K_ThScr	= 1e-3		#!

###############################################################################################

def MC_Blow_Air(data):
	U_Blow, P_Blow, A_Flr, n_Heat_CO2 = data.U_Blow, data.P_Blow, data.A_Flr, data.n_Heat_CO2

	H_Blow_Air = U_Blow * P_Blow / A_Flr
	return n_Heat_CO2 * H_Blow_Air

###############################################################################################

###############################################################################################

def MC_Ext_Air(data):
	U_Ext_CO2, phi_Ext_CO2, A_Flr = data.U_Ext_CO2, data.phi_Ext_CO2, data.A_Flr

	return U_Ext_CO2 * phi_Ext_CO2 / A_Flr

###############################################################################################

###############################################################################################

def MC_Pad_Air(data, CO2_Air):
	U_Pad, phi_Pad, A_Flr = data.U_Pad, data.phi_Pad, data.A_Flr

	f_Pad = U_Pad * phi_Pad / A_Flr
	return f_Pad * (CO2_Out - CO2_Air)

###############################################################################################

###############################################################################################

def MC_Air_Top(data, CO2_Air, CO2_Top):

	U_ThScr, K_ThScr, T_Air, T_Top = data.U_ThScr, data.K_ThScr, data.T_Air, data.T_Top
	p_Mean_Air, p_Air, p_Top = data.p_Mean_Air, data.p_Air, data.p_Top

	f_ThScr = U_ThScr * K_ThScr * abs(T_Air - T_Top)**(0.66) 
	f_ThScr += (1 - U_ThScr) * (g * (1 - U_ThScr) / (2 * p_Mean_Air) * (p_Air - p_Top))**0.5
	return f_ThScr * (CO2_Air - CO2_Top)

###############################################################################################

def MC_Air_Out(data,CO2_Air):
	CO2_Out = data.CO2_Out

	return (f_Vent_Side(data) + f_Vent_Forced(data)) * (CO2_Air - CO2_Out)

def f_Vent_Roof_Side(data, A_Roof):
	C_d, A_Flr = data.C_d,  data.A_Flr
	U_Roof, U_Side, A_Side = data.U_Roof, data.U_Side, data.A_Side
	g, h_Side_Roof, T_Air, T_Out, T_Mean_Air = data.g, data.h_Side_Roof, data.T_Air, data.T_Out, data.T_Mean_Air
	C_w, v_wind = data.C_w, data.v_wind

	m1 = C_d / A_Flr
	m2 = (U_Roof * U_Side * A_Roof * A_Side)**2 / ((U_Roof * A_Roof)**2 + (U_Side * A_Side)**2)
	m2 *= 2 * g * h_Side_Roof * (T_Air - T_Out) / T_Mean_Air
	m3 = ((U_Roof * A_Roof + U_Side * A_Side) / 2)**2
	m3 *= C_w * v_wind**2
	return m1 * (m2 + m3)**0.5

def n_Ins_Scr(data):
	S_Ins_Scr = data.S_Ins_Scr

	return S_Ins_Scr * (2 - S_Ins_Scr)

def f_leakage(data):
	v_wind, c_leakage = data.v_wind, data.c_leakage 

	if v_wind < 0.25:
		v_wind = 0.25
	return v_wind * c_leakage
	
def f_Vent_Side(data):
	n_Side, n_Side_Thr, n_Ins_Scr, f_leakage, U_ThScr = data.n_Side, data.n_Side_Thr, n_Ins_Scr(data), f_leakage(data), data.U_ThScr

	f_2com_Vent_Side = f_Vent_Roof_Side(data,0)

	if n_Side >= n_Side_Thr:
		return n_Ins_Scr * f_2com_Vent_Side + 0.5 * f_leakage
	else:
		return n_Ins_Scr * (U_ThScr * f_2com_Vent_Side + (1 - U_ThScr) * f_Vent_Roof_Side(data,data.A_Roof) * n_Side) + 0.5* f_leakage

def f_Vent_Forced(data):
	n_Ins_Scr, U_Vent_Forced, phi_Vent_Forced, A_Flr = data.n_Ins_Scr, data.U_Vent_Forced, data.phi_Vent_Forced, data.A_Flr

	return (n_Ins_Scr * U_Vent_Forced * phi_Vent_Forced) / A_Flr

###############################################################################################

###############################################################################################

def MC_Top_Out(data, CO2_Top):
	f_Vent_Roof, CO2_Out = data.f_Vent_Roof, data.CO2_Out

	return f_Vent_Roof * (CO2_Top - CO2_Out)

def f_Vent_Roof(data):
	n_Roof, n_Roof_Thr, U_ThScr, f_leakage = data.n_Roof, data.n_Roof_Thr, data.U_ThScr, data.f_leakage

	if n_Roof >= n_Roof_Thr:
		return n_Ins_Scr * f_2com_Vent_Roof(data) + 0.5 * f_leakage
	else:
		return n_Ins_Scr * (U_ThScr * f_2com_Vent_Roof(data) + (1 - U_ThScr) * f_Vent_Roof_Side(data) * n_Side) + 0.5* f_leakage

def f_2com_Vent_Roof(data):
	C_d, U_Roof, A_Roof, A_Flr = data.C_d, data.U_Roof, data.A_Roof, data.A_Flr
	g, h_Roof, T_Air, T_Out, T_Mean_Air = data.g, data.h_Roof, data.T_Air, data.T_Out, data.T_Mean_Air
	C_w, v_wind = data.C_w, data.v_wind

	m1 = C_d * U_Roof * A_Roof / (2 * A_Flr)
	m2 = g * h_Roof * (T_Air - T_Out) / (2 * T_Mean_Air)
	m3 = C_w * v_wind**2
	return m1 * (m2 + m3) ** 0.5

###############################################################################################

def MC_Air_Can(data):
	M_CH2O, P, R = data.M_CH2O, PhotoSynth(data), PhotoRespi(data)

	return M_CH2O * h_C_Buf * (P - R)

###############################################################################################

def PhotoSynth(data):
	Gamma = GComp(data)
	CO2_Stom = CO2Stomata(data)
	P = JTrans(data) * (CO2_Stom - Gamma)
	P /= 4 * (CO2_Stom + 2 * Gamma)
	return P

def PhotoRespi(data):
	return PhotoSynth(data) * GComp(data) / CO2Stomata(data)

def JTrans(data):
	aP = alpha * PAR_Can
	J_Pot = JPotent(data)
	J = J_Pot + aP - ((J_Pot + aP)**2 - 4 * Theta * J_Pot * aP)**0.5

	return J / (2 * Theta)

def JPotent(data):
	T_Can, LAI = data.T_Can, data.LAI
	J_MAX_25_Can = LAI * J_MAX_25_Leaf

	exp1 = E_j * (T_Can - T_25_K) / (R_gas * T_Can * T_25_K)
	exp1 = exp(exp1)

	exp2 = S_J_Pot * T_25_K - H_J_Pot
	exp2 /= R_gas * T_25_K
	exp2 = 1 + exp(exp2)

	exp3 = S_J_Pot * T_Can - H_J_Pot
	exp3 /= R_gas * T_Can
	exp3 = 1 + exp(exp3)
	
	return J_MAX_25_Can * exp1 * exp2 / exp3

def CO2Stomata(data):
	return n_CO2_Air_Stom * data.CO2_Air

def GComp(data):	# Use Eq. (9.23) (more complex) or Eq. (9.22) (simpler)?
	T_Can, LAI = data.T_Can, data.LAI

	return C_Gamma * T_Can / LAI + 20 * C_Gamma * (1 - 1 / LAI)

###############################################################################################

def dxCO2_Air(data, CO2_Air, CO2_Top):

	return (MC_Blow_Air(data) + MC_Ext_Air(data) + MC_Pad_Air(data, CO2_Air) - MC_Air_Can(data) - MC_Air_Top(data, CO2_Air, CO2_Top) - MC_Air_Out(data, CO2_Air))/ capCO2_Air

def dxCO2_Top(data, CO2_Air, CO2_Top):

	return (MC_Air_Top(data, CO2_Air, CO2_Top) - MC_Top_Out(data, CO2_Top))/ capCO2_Top
	