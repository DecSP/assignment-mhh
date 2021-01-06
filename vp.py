from math import *
from co2 import *

# Global constants
Delta_H = 2.45 * 10 ** 6
c_p_Air = 1 * 10 ** 3
gamma = 65.8
R_B = 275
s_MV_1_2 = -0.1
M_Water = 18

data["eta_Pad"] = 0.0
data["Rh_Out"] = 81.7 # Get this from csv, dynamic variable

def Compute_VP(T, Rh):
    P_Sat = 610.78 * exp(T / (T + 238.3) * 17.2694)

    return Rh * P_Sat

data["VP_Out"] = Compute_VP(data["T_Out"], data["Rh_Out"])

###############################################################################################

def MV_Air_Object(data, VP_1, VP_2, HEC_1_2):
    m1 = 1.0 / (1 + exp(s_MV_1_2 * (VP_1 - VP_2)))
    m2 = 6.4 * 10 ** -9 * HEC_1_2 * (VP_1 - VP_2)

    return m1 * m2

###############################################################################################

def MV_Air_Top(data, VP_Air, VP_Top):
    global R_gas, M_Water
    U_ThScr, K_ThScr, T_Air, T_Top = data["U_ThScr"], data["K_ThScr"], data["T_Air"], data["T_Top"]
    rho_Mean_Air, rho_Air, rho_Top = data["rho_Mean_Air"], data["rho_Air"], data["rho_Top"]

    f_ThScr = U_ThScr * K_ThScr * abs(T_Air - T_Top)**(0.66) 
    f_ThScr += ((1 - U_ThScr) * (g * (1 - U_ThScr) / (2 * rho_Mean_Air) * abs(rho_Air - rho_Top))**0.5)

    return (M_Water / R_gas) * f_ThScr * (VP_Air / T_Air - VP_Top / T_Top)

###############################################################################################

def MV_Air_Out(data, VP_Air):
    global R_gas, M_Water
    VP_Out, T_Air, T_Out = data["VP_Out"], data["T_Air"], data["T_Out"]

    f_1_2 = f_Vent_Side(data) + f_Vent_Forced(data)

    return (M_Water / R_gas) * f_1_2 * (VP_Air / T_Air - VP_Out / T_Out)

###############################################################################################

def MV_Can_Air(data, VP_Air):
    return VEC_Can_Air(data) * (0)

def VEC_Can_Air(data):
    R_S = R_S_min

    return (2 * data["rho_Air"] * c_p_Air * data["LAI"]) / (Delta_H * gamma * (R_B + R_S))

###############################################################################################

def MV_Pad_Air(data):
    rho_Air, U_Pad, phi_Pad, A_Flr, eta_Pad = data["rho_Air"], data["U_Pad"], data["phi_Pad"], data["A_Flr"], data["eta_Pad"]
    x_Pad, x_Out = data["x_Pad"], data["x_Out"]

    f_Pad = U_Pad * phi_Pad / A_Flr

    return rho_Air * f_Pad * (eta_Pad * (x_Pad - x_Out) + x_Out)

###############################################################################################

def MV_AirOut_Pad(data, VP_Air):
    global R_gas, M_Water
    T_Air = data["T_Air"]

    f_Pad = U_Pad * phi_Pad / A_Flr

    return f_Pad * (M_Water / R_gas) * (VP_Air / T_Air)
