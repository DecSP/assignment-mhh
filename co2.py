from math import *

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

def MC_Pad_Air(data):
	U_Pad, phi_Pad, A_Flr = data.U_Pad, data.phi_Pad, data.A_Flr

	f_Pad = U_Pad * phi_Pad / A_Flr
	return f_Pad * (CO2_Out - CO2_Air)

###############################################################################################

###############################################################################################

def MC_Air_Top(data):
	U_ThScr, K_ThScr, T_Air, T_Top = data.U_ThScr, data.K_ThScr, data.T_Air, data.T_Top
	p_Mean_Air, p_Air, p_Top = data.p_Mean_Air, data.p_Air, data.p_Top

	f_ThScr = U_ThScr * K_ThScr * abs(T_Air - T_Top) ** (2/3) 
	f_ThScr+= (1 - U_ThScr) * (g * (1 - U_ThScr)/ (2 * p_Mean_Air)*(p_Air - p_Top)) ** (1/2)
	return 0

###############################################################################################

def MC_Air_Out(data):
	f_Vent_Side, f_Vent_Forced, CO2_Air, CO2_Out = data.f_Vent_Side, data.f_Vent_Forced, data.CO2_Air, data.CO2_Out

	return (f_Vent_Side + f_Vent_Forced) * (CO2_Air - CO2_Out)

def f_Vent_Roof_Side(data, A_Roof):
	C_d, A_Flr = data.C_d,  data.A_Flr
	U_Roof, U_Side, A_Side = data.U_Roof, data.U_Side, data.A_Side
	g, h_Side_Roof, T_Air, T_Out, T_Mean_Air = data.g, data.h_Side_Roof, data.T_Air, data.T_Out, data.T_Mean_Air
	C_w, v_wind = data.C_w, data.v_wind

	m1 = C_d / A_Flr
	m2 = (U_Roof * U_Side * A_Roof * A_Side)**2 / ((U_Roof * A_Roof)**2 + (U_Side * A_Side)**2)
	m2 *= 2 * g * h_Side_Roof * (T_Air - T_Out) / T_Mean_Air
	m3 = ((U_Roof * A_Roof + U_Side * A_Side) / 2)**2
	m3 *=C_w * v_wind**2
	return m1 * (m2 + m3)**(1/2)

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

def MC_Top_Out(data):
	f_Vent_Roof, CO2_Top, CO2_Out = data.f_Vent_Roof, data.CO2_Top, data.CO2_Out

	return f_Vent_Roof * (CO2_Top - CO2_Out)

def f_Vent_Roof(data):
	n_Roof, n_Roof_Thr, U_ThScr, f_leakage = data.n_Roof, data.n_Roof_Thr, data.U_ThScr, data.f_leakage

	if n_Roof >= n_Roof_Thr:
		return n_Ins_Scr * f_2com_Vent_Roof(data) + 0.5 * f_leakage
	else:
		return n_Ins_Scr * (U_ThScr * f_2com_Vent_Roof(data) + (1 - U_ThScr) * f_Vent_Roof_Side * n_Side) + 0.5* f_leakage

def f_2com_Vent_Roof(data):
	C_d, U_Roof, A_Roof, A_Flr = data.C_d, data.U_Roof, data.A_Roof, data.A_Flr
	g, h_Roof, T_Air, T_Out, T_Mean_Air = data.g, data.h_Roof, data.T_Air, data.T_Out, data.T_Mean_Air
	C_w, v_wind = data.C_w, data.v_wind

	m1 = C_d * U_Roof * A_Roof / (2 * A_Flr)
	m2 = g * h_Roof * (T_Air - T_Out) / (2 * T_Mean_Air)
	m3 = C_w * v_wind**2
	return m1 * (m2 + m3) ** 0.5

def MC_Air_Can(data):
	M_CH2O, P, R = data.M_CH2O, data.P, data.R

	return M_CH2O * h_C_Buf(data) * (P - R)

def h_C_Buf(data):
	C_Buf, C_Max_Buf = data.C_Buf, data.C_Max_Buf

	return 0 if C_Buf > C_Max_Buf else 0

###############################################################################################

