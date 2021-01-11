import pandas as pd
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
Delta_H = 2.45 * (10 ** 6)
c_p_Air = 1 * (10 ** 3)
gamma = 65.8
R_B = 275
R_S_min = 82.0
s_MV_1_2 = -0.1
M_Water = 18
eta_Heat_CO2 = 0.057
eta_Roof_Thr = 0.9
eta_Side_Thr = 0.9

data = {}

# Read file
climate = pd.read_csv('Greenhouse_climate.csv')
meteo = pd.read_csv("meteo.csv")
# climate = climate.dropna(how="any")
# meteo = meteo.dropna(how="any")

# Data initialization
data["h_Air"], data["h_Gh"] = 3.8, 4.2
data["A_Flr"], data["A_Side"] = 1.4 * (10 ** 4) , 0
data["A_Roof"] = 0.1 * data["A_Flr"]
data["A_Cov"] = 1.8 * 10 ** 4
data["W_Gutter"] = sqrt(data["A_Flr"] * 0.3)
def Compute_C(GH_C,eta_Sh_Scr_C):
	U_Sh = 0.5
	return GH_C * (1 - eta_Sh_Scr_C * U_Sh)

data["C_d"], data["C_w"] = Compute_C(0.75, 0), Compute_C(0.09, 0)

data["K_ThScr"] = 0.25 * 10 ** -3
data["zeta_Ins_Scr"] = 1.0
data["T_Air"], data["T_Out"], data["T_Can"], data["T_Top"], data["T_Mean_Air"] = 18+273.15, 15.8+273.15, 19.9+273.15, 19.9+273.15, (18+15.8)/2+273.15

def Compute_VP(T, Rh):
    P_Sat = 610.78 * exp((T - 273.15) / (T - 34.85) * 17.2694)

    return (Rh / 100 * P_Sat)

data["Rh_Out"] = 81.7

data["x_Pad"] = 0
data["x_Out"] = data["Rh_Out"]
data["T_ThScr"] = data["T_Air"] + 1
data["T_MechCool"] = data["T_Air"]
data["T_CovIn"] = data["T_Out"]
data["VP_Out"] = Compute_VP(data["T_Out"], data["Rh_Out"])
data["VP_ThScr"] = Compute_VP(data["T_ThScr"], data["Rh_Out"])
data["VP_MechCool"] = Compute_VP(data["T_MechCool"], data["Rh_Out"])
data["VP_CovIn"] = Compute_VP(data["T_CovIn"], data["Rh_Out"])

def update_data(time):
    data["T_Air"] = climate.Tair[time] + 273.15
    data["T_Out"] = meteo.Tout[time] + 273.15
    data["T_Can"] = data["T_Air"] - 1
    data["T_Top"] = data["T_Can"]
    data["T_Mean_Air"] = (data["T_Air"] + data["T_Out"])/2
    data["Rh_Out"] = meteo.Rhout[time]
    # data["U_Side"] = climate.VentLee[time]/100
    # data["U_Vent_Forced"] = climate.Ventwind[time]/100
    data["v_wind"] = meteo.Windsp[time]
    data["rho_Air"] = Compute_rho(data["T_Air"]) 
    data["rho_Top"] = Compute_rho(data["T_Top"])
    data["rho_Mean_Air"] =  (data["rho_Air"]+data["rho_Top"])/2

    data["x_Pad"] = 0
    data["x_Out"] = data["Rh_Out"]
    data["T_ThScr"] = data["T_Air"] + 1
    data["T_MechCool"] = data["T_Air"]
    data["T_CovIn"] = data["T_Out"]
    data["VP_Out"] = Compute_VP(data["T_Out"], data["Rh_Out"])
    data["VP_ThScr"] = Compute_VP(data["T_ThScr"], data["Rh_Out"])
    data["VP_MechCool"] = Compute_VP(data["T_MechCool"], data["Rh_Out"])
    data["VP_CovIn"] = Compute_VP(data["T_CovIn"], data["Rh_Out"])

data["h_Vent"], data["h_C_Buf"], data["h_Side_Roof"] = 0.68, 1, 3.8/2
data["U_Blow"], data["U_Ext_CO2"], data["U_Pad"], data["U_Roof"], data["U_Side"], data["U_ThScr"],  data["U_Vent_Forced"] = 0.0, 0.1, 0.0, 0.1, 0.0, 0.9, 0.0
data["U_Fog"], data["U_MechCool"] = 0.0, 0.0
data["c_leakage"] = 1 * 10 ** - 4
data["c_HECin"] = 1.86
data["v_wind"] = 3.2
data["M_CH2O"] = 30*(10**-3)
data["phi_Pad"], data["phi_Vent_Forced"], data["phi_Ext_CO2"], data["phi_Fog"] = 0.1, 0, 7.2 * (10 ** 4), 0.0
data["COP_MechCool"], data["P_MechCool"] = 0, 0
data["P_Blow"] = 0.5*(10**6)
data["LAI"] = 2.0
data["CO2_Out"] = 668
data["eta_Side"], data["eta_Roof"], data["eta_Pad"], data["eta_HeatVap"] = 1, 0.0, 0, 4.43 * (10 ** -8)

def Compute_rho(T):
	p = 101*1000
	R = 287.058
	return p/(R*T)

data["rho_Air"] = Compute_rho(data["T_Air"]) 
data["rho_Top"] = Compute_rho(data["T_Top"])
data["rho_Mean_Air"] =  (data["rho_Air"]+data["rho_Top"])/2
data["PAR_Can"] = 100