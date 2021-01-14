from co2 import *
from vp import *
import matplotlib.pyplot as plt

def getco2real(n):
	return [0.0409*climate.CO2air[i]*44.01 for i in range(2,2+n)]

def mse(pred, real):
	return (sum([(pred[i] - real[i])**2 for i in range(len(pred)) ]))/len(pred)

#x0, y0 for Co2air, co2top, t0, t, h mean time start, time end, step
def euler(t0, x0, y0, dx0, dy0, h):
    K1x=h*dx0(data,x0, y0)
    K1y=h*dy0(data, x0,y0)
    K2x=h*dx0(data, x0+K1x, y0+K1y)
    K2y=h*dy0(data, x0+K1x, y0+K1y)
    x0=x0+1.0/2.0*(K1x+K2x)
    y0=y0+1.0/2.0*(K1y+K2y)
    return x0,y0

def eulerN(t0, x0, y0, dx0, dy0, h, n, numCycle = 60, startIndex = 2):
    # n=int((t-t0)/h)
    air=[]
    top=[]
    re=[]
    for i in range(n):
        for j in range(numCycle):
            x0,y0 = euler(t0, x0, y0, dx0, dy0, h)
        update_data(i + startIndex)
        air.append(x0)
        top.append(y0)
        print("After "+str(numCycle*(i+1))+" cycles:",mse(air,getco2real(len(air))))
    return air,top

def rk4(t0, x0, y0, dx0, dy0, h):
    K1x=h*dx0(data, x0, y0)
    K1y=h*dy0(data, x0, y0)
    K2x=h*dx0(data, x0+K1x/2.0, y0+K1y/2.0)
    K2y=h*dy0(data, x0+K1x/2.0, y0+K1y/2.0)
    K3x=h*dx0(data, x0+K2x/2.0, y0+K2y/2.0)
    K3y=h*dy0(data, x0+K2x/2.0, y0+K2y/2.0)
    K4x=h*dx0(data, x0+K3x, y0+K3y)
    K4y=h*dy0(data, x0+K3x, y0+K3y)
    x0=x0+1.0/6.0*(K1x+2*K2x+2*K3x+K4x)
    y0=y0+1.0/6.0*(K1y+2*K2y+2*K3y+K4y)
    return x0,y0


def rk4N(t0, x0, y0, dx0, dy0, h, n, numCycle = 60, startIndex = 2):
    # n=int((t-t0)/h)
    air=[]
    top=[]
    for i in range (n):
        for j in range(numCycle):
            x0,y0 = rk4(t0, x0, y0, dx0, dy0, h)
        update_data(i + startIndex)
        air.append(x0)
        top.append(y0)
        print("After "+str(numCycle*(i+1))+" cycles:",mse(air,getco2real(len(air))))
    return air,top

air, top = rk4N(0, 768.6, 768.6, dxCO2_Air, dxCO2_Top, 5, 4000) # 768.6

# Let's compare by plotting everything

co2air_real = getco2real(len(air))
xline = [5 * k for k in range(1,len(air)+1)]
plt.xlabel("time elapsed (mins)")
plt.ylabel("concentration (mg/m^3)")
plt.title("Changing of CO2 concentration")
plt.plot(xline, air, label = "CO2air")
plt.plot(xline, top, label = "CO2top")
plt.legend()
plt.savefig("CO2_AirTop.png")
plt.clf()
plt.xlabel("time elapsed (mins)")
plt.ylabel("concentration (mg/m^3)")
plt.plot(xline, air, label = "CO2air")
plt.plot(xline, top, label = "CO2top")
plt.title("Compare CO2 concentration")
plt.plot(xline, co2air_real, label = "CO2air_real")
plt.legend()
plt.savefig("CO2_Compare.png")

###############################################################################################
plt.clf()

air,top = rk4N(0, 1525.4, 1525.4, dxVP_Air, dxVP_Top, 5, 4000)

# Let's compare by plotting everything

VPair_real = [Compute_VP(climate.Tair[i] + 273.15, climate.RHair[i]) for i in range(2,2+len(air))]
xline = [5 * k for k in range(1,len(air)+1)]
plt.xlabel("time elapsed (mins)")
plt.ylabel("pressure (Pa)")
plt.title("Changing of vapour pressure")
plt.plot(xline, air, label = "VPair")
plt.plot(xline, top, label = "VPtop")
plt.legend()
plt.savefig("VP_AirTop.png")
plt.clf()
plt.xlabel("time elapsed (mins)")
plt.ylabel("Pressure (Pa)")
plt.plot(xline, air, label = "VPair")
plt.plot(xline, top, label = "VPtop")
plt.title("Compare vapour pressure")
plt.plot(xline, VPair_real, label = "VPair_real")
plt.legend()
plt.savefig("VP_Compare.png")