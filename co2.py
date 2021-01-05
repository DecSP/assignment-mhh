from math import *


data = {}

# Data initialization
data["A_Flr"], data["A_Side"] = 
data["C_d"], data["C_w"] = 
data["phi_Pad"], data["phi_Vent_Forced"], data["phi_Ext_CO2"] = 
data["K_ThScr"] = 
data["P_Blow"] = 
data["zeta_Ins_Scr"] = 
data["T_Air"], data["T_Out"], data["T_Can"], data["T_Top"], data["T_Mean_Air"] = 
data["h_Roof"], data["h_C_Buf"], data["h_Side_Roof"] = 
data["U_Blow"], data["U_Ext_CO2"], data["U_Pad"], data["U_Roof"], data["U_Side"], data["U_ThScr"],  data["U_Vent_Forced"] = 
data["c_leakage"], data["f_leakage"] = 
data["v_wind"] = 

data["rho_Air"] =  
data["rho_Top"] =  
data["rho_Mean_Air"] =  


## main code ##

# :(((
###############

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
eta_Heat_CO2=0.057
eta_Roof_Thr=0.9
eta_Side_Thr=0.9


###############################################################################################

def MC_Blow_Air(data):
	global eta_Heat_CO2
	U_Blow, P_Blow, A_Flr = data["U_Blow"], data["P_Blow"], data["A_Flr"]

	H_Blow_Air = U_Blow * P_Blow / A_Flr
	return eta_Heat_CO2 * H_Blow_Air

###############################################################################################

###############################################################################################

def MC_Ext_Air(data):
	U_Ext_CO2, phi_Ext_CO2, A_Flr = data["U_Ext_CO2"], data["phi_Ext_CO2"], data["A_Flr"]

	return U_Ext_CO2 * phi_Ext_CO2 / A_Flr

###############################################################################################

###############################################################################################

def MC_Pad_Air(data, CO2_Air):
	U_Pad, phi_Pad, A_Flr = data["U_Pad"], data["phi_Pad"], data["A_Flr"]

	f_Pad = U_Pad * phi_Pad / A_Flr
	return f_Pad * (CO2_Out - CO2_Air)

###############################################################################################

###############################################################################################

def MC_Air_Top(data, CO2_Air, CO2_Top):

	U_ThScr, K_ThScr, T_Air, T_Top = data["U_ThScr"], data["K_ThScr"], data["T_Air"], data["T_Top"]
	rho_Mean_Air, rho_Air, rho_Top = data["rho_Mean_Air"], data["rho_Air"], data["rho_Top"]

	f_ThScr = U_ThScr * K_ThScr * abs(T_Air - T_Top)**(0.66) 
	f_ThScr += (1 - U_ThScr) * (g * (1 - U_ThScr) / (2 * rho_Mean_Air) * (rho_Air - rho_Top))**0.5
	return f_ThScr * (CO2_Air - CO2_Top)

###############################################################################################

def MC_Air_Out(data,CO2_Air):
	CO2_Out = data["CO2_Out"]

	return (f_Vent_Side(data) + f_Vent_Forced(data)) * (CO2_Air - CO2_Out)

def f_Vent_Roof_Side(data, A_Roof):
	global g
	C_d, A_Flr = data["C_d"],  data["A_Flr"]
	U_Roof, U_Side, A_Side = data["U_Roof"], data["U_Side"], data["A_Side"]
	h_Side_Roof, T_Air, T_Out, T_Mean_Air = data["h_Side_Roof"], data["T_Air"], data["T_Out"], data["T_Mean_Air"]
	C_w, v_wind = data["C_w"], data["v_wind"]

	m1 = C_d / A_Flr
	m2 = (U_Roof * U_Side * A_Roof * A_Side)**2 / ((U_Roof * A_Roof)**2 + (U_Side * A_Side)**2)
	m2 *= 2 * g * h_Side_Roof * (T_Air - T_Out) / T_Mean_Air
	m3 = ((U_Roof * A_Roof + U_Side * A_Side) / 2)**2
	m3 *= C_w * v_wind**2
	return m1 * (m2 + m3)**0.5

def Compute_eta_Ins_Scr(data):
	zeta_Ins_Scr = data["zeta_Ins_Scr"]

	return zeta_Ins_Scr * (2 - zeta_Ins_Scr)

def Compute_f_leakage(data):
	v_wind, c_leakage = data["v_wind"], data["c_leakage"] 

	if v_wind < 0.25:
		v_wind = 0.25
	return v_wind * c_leakage
	
def f_Vent_Side(data):
	global eta_Side_Thr
	eta_Side, eta_Ins_Scr, f_leakage, U_ThScr = data["eta_Side"], Compute_eta_Ins_Scr(data), Compute_f_leakage(data), data["U_ThScr"]

	f_2com_Vent_Side = f_Vent_Roof_Side(data,0)

	if eta_Side >= eta_Side_Thr:
		return eta_Ins_Scr * f_2com_Vent_Side + 0.5 * f_leakage
	else:
		return eta_Ins_Scr * (U_ThScr * f_2com_Vent_Side + (1 - U_ThScr) * f_Vent_Roof_Side(data,data["A_Roof"]) * eta_Side) + 0.5* f_leakage

def f_Vent_Forced(data):
	U_Vent_Forced, phi_Vent_Forced, A_Flr = data["U_Vent_Forced"], data["phi_Vent_Forced"], data["A_Flr"]

	return (Compute_eta_Ins_Scr(data) * U_Vent_Forced * phi_Vent_Forced) / A_Flr

###############################################################################################

###############################################################################################

def MC_Top_Out(data, CO2_Top):
	CO2_Out = data["CO2_Out"]

	return f_Vent_Roof(data) * (CO2_Top - CO2_Out)

def f_Vent_Roof(data):
	global eta_Roof_Thr
	eta_Side, eta_Roof, U_ThScr, f_leakage = data["eta_Side"], data["eta_Roof"], data["U_ThScr"], Compute_f_leakage(data)

	if eta_Roof >= eta_Roof_Thr:
		return eta_Ins_Scr * f_2com_Vent_Roof(data) + 0.5 * f_leakage
	else:
		return eta_Ins_Scr * (U_ThScr * f_2com_Vent_Roof(data) + (1 - U_ThScr) * f_Vent_Roof_Side(data) * eta_Side) + 0.5* f_leakage

def f_2com_Vent_Roof(data):
	C_d, U_Roof, A_Roof, A_Flr = data["C_d"], data["U_Roof"], data["A_Roof"], data["A_Flr"]
	g, h_Roof, T_Air, T_Out, T_Mean_Air = data["g"], data["h_Roof"], data["T_Air"], data["T_Out"], data["T_Mean_Air"]
	C_w, v_wind = data["C_w"], data["v_wind"]

	m1 = C_d * U_Roof * A_Roof / (2 * A_Flr)
	m2 = g * h_Roof * (T_Air - T_Out) / (2 * T_Mean_Air)
	m3 = C_w * v_wind**2
	return m1 * (m2 + m3) ** 0.5

###############################################################################################

def MC_Air_Can(data):
	M_CH2O, P, R = data["M_CH2O"], data["P"], data["R"]

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
	T_Can, LAI = data["T_Can"], data["LAI"]
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
	return n_CO2_Air_Stom * data["CO2_Air"]

def GComp(data):	# Use Eq. (9.23) (more complex) or Eq. (9.22) (simpler)?
	T_Can, LAI = data["T_Can"], data["LAI"]

	return C_Gamma * T_Can / LAI + 20 * C_Gamma * (1 - 1 / LAI)

###############################################################################################

def dxCO2_Air(data, CO2_Air, CO2_Top):

	return (MC_Blow_Air(data) + MC_Ext_Air(data) + MC_Pad_Air(data, CO2_Air) - MC_Air_Can(data) - MC_Air_Top(data, CO2_Air, CO2_Top) - MC_Air_Out(data, CO2_Air))/ capCO2_Air

def dxCO2_Top(data, CO2_Air, CO2_Top):

	return (MC_Air_Top(data, CO2_Air, CO2_Top) - MC_Top_Out(data, CO2_Top))/ capCO2_Top
	

