from math import *
from co2 import *
from data import *

###############################################################################################

def MV_Air_Object(data, VP_1, VP_2, HEC_1_2):
    m1 = 1.0 / (1 + exp(s_MV_1_2 * (VP_1 - VP_2)))
    m2 = 6.4 * 10 ** -9 * HEC_1_2 * (VP_1 - VP_2)

    return m1 * m2

###############################################################################################

def MV_Can_Air(data):
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

def MV_Fog_Air(data, VP_Air):
    U_Fog, phi_Fog, A_Flr = data["U_Fog"], data["phi_Fog"], data["A_Flr"]

    return U_Fog * phi_Fog / A_Flr

###############################################################################################

def MV_Blow_Air(data):
    U_Blow, P_Blow, A_Flr, eta_HeatVap = data["U_Blow"], data["P_Blow"], data["A_Flr"], data["eta_HeatVap"]

    H_Blow_Air = U_Blow * P_Blow / A_Flr

    return eta_HeatVap * H_Blow_Air

###############################################################################################

def MV_Air_ThScr(data, VP_Air):
    return MV_Air_Object(data, VP_Air, data["VP_ThScr"], HEC_Air_ThScr(data))

def HEC_Air_ThScr(data):
    return 1.7 * data["U_ThScr"] * abs(data["T_Air"] - data["T_ThScr"]) ** 0.33

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

def MV_AirOut_Pad(data, VP_Air):
    global R_gas, M_Water
    T_Air, U_Pad, phi_Pad, A_Flr = data["T_Air"], data["U_Pad"], data["phi_Pad"], data["A_Flr"]

    f_Pad = U_Pad * phi_Pad / A_Flr

    return f_Pad * (M_Water / R_gas) * (VP_Air / T_Air)

###############################################################################################

def MV_Air_Mech(data, VP_Air):
    return MV_Air_Object(data, VP_Air, data["VP_MechCool"], -HEC_Mech_Air(data, VP_Air))

def HEC_Mech_Air(data, VP_Air):
    global Delta_H

    U_MechCool, COP_MechCool, P_MechCool, A_Flr, T_Air, T_MechCool, VP_MechCool = data["U_MechCool"], data["COP_MechCool"], data["P_MechCool"], data["A_Flr"], data["T_Air"], data["T_MechCool"], data["VP_MechCool"]
    
    m1 = U_MechCool * COP_MechCool * P_MechCool / A_Flr
    m2 = T_Air - T_MechCool + 6.4 * 10 ** -9 * Delta_H * (VP_Air - VP_MechCool)

    return m1 / m2

###############################################################################################

def MV_Top_CovIn(data, VP_Top):
    return MV_Air_Object(data, VP_Top, data["VP_CovIn"], HEC_Top_CovIn(data, VP_Top))

def HEC_Top_CovIn(data, VP_Top):
    c_HECin, T_Top, T_CovIn, A_Cov, A_Flr = data["c_HECin"], data["T_Top"], data["T_CovIn"], data["A_Cov"], data["A_Cov"]
    
    return c_HECin * (T_Top - T_CovIn) ** 0.33 * A_Cov / A_Flr

###############################################################################################

def MV_Top_Out(data, VP_Top):
    global R_gas, M_Water
    VP_Out, T_Top, T_Out = data["VP_Out"], data["T_Top"], data["T_Out"]

    return (M_Water / R_gas) * f_Vent_Roof(data) * (VP_Top / T_Top - VP_Out / T_Out)

###############################################################################################

def dxVP_Air(data, VP_Air, VP_Top):
    capVP_Air = 3.8
	
    return (MV_Can_Air(data) + MV_Pad_Air(data) + MV_Fog_Air(data, VP_Air) + MV_Blow_Air(data) - MV_Air_ThScr(data, VP_Air) - MV_Air_Top(data, VP_Air, VP_Top) - MV_Air_Out(data, VP_Air) - MV_AirOut_Pad(data, VP_Air) - MV_Air_Mech(data, VP_Air)) / capVP_Air

def dxVP_Top(data, VP_Air, VP_Top):
    capVP_Top = 0.4

    return (MV_Air_Top(data, VP_Air, VP_Top) - MV_Top_CovIn(data, VP_Top) - MV_Top_Out(data, VP_Top)) / capVP_Top