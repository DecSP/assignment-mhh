from math import *

# Global constants
Delta_H = 2.45 * 10 ** 6
c_p_Air = 1 * 10 ** 3
gamma = 65.8
R_B = 275
s_MV_1_2 = -0.1

###############################################################################################

def MV_Can_Air(data, VP_Air):
    return VEC_Can_Air(data) * (0)

def VEC_Can_Air(data):
    R_S = R_S_min

    return (2 * data["rho_Air"] * c_p_Air * data["LAI"]) / (Delta_H * gamma * (R_B + R_S))

###############################################################################################

def MV_Pad_Air(data):
    rho_Air