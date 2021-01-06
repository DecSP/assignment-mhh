from co2 import *
import matplotlib.pyplot as plt 
#x0, y0 for Co2air, co2top, t0, t, h mean time start, time end, step
def euler(t0, x0, y0, h):
    K1x=h*dxCO2_Air(data,x0, y0)
    K1y=h*dxCO2_Top(data, x0,y0)
    K2x=h*dxCO2_Air(data, x0+K1x, y0+K1y)
    K2y=h*dxCO2_Top(data, x0+K1x, y0+K1y)
    x0=x0+1.0/2.0*(K1x+K2x)
    y0=y0+1.0/2.0*(K1y+K2y)
    return x0,y0

def eulerN(t0, x0, y0, h, n, numCycle = 60, startIndex = 2):
    # n=int((t-t0)/h)
    re=[]
    for i in range(n):
        for j in range(numCycle):
            x0,y0 = euler(t0, x0, y0, h)
        update_data(i + startIndex)
        re.append((x0,y0))
    return re

def rk4(t0, x0, y0, h):
    K1x=h*dxCO2_Air(data,x0, y0)
    K1y=h*dxCO2_Top(data, x0,y0)
    K2x=h*dxCO2_Air(data, x0+K1x/2.0, y0+K1y/2.0)
    K2y=h*dxCO2_Top(data, x0+K1x/2.0, y0+K1y/2.0)
    K3x=h*dxCO2_Air(data,x0+K2x/2.0, y0+K2y/2.0)
    K3y=h*dxCO2_Top(data, x0+K2x/2.0,y0+K2y/2.0)
    K4x=h*dxCO2_Air(data, x0+K3x, y0+K3y)
    K4y=h*dxCO2_Top(data, x0+K3x, y0+K3y)
    x0=x0+1.0/6.0*(K1x+2*K2x+2*K3x+K4x)
    y0=y0+1.0/6.0*(K1y+2*K2y+2*K3y+K4y)
    return x0,y0


def rk4N(t0, x0, y0, h, n, numCycle = 60, startIndex = 2):
    # n=int((t-t0)/h)
    re = []
    for i in range (n):
        for j in range(numCycle):
            x0,y0 = rk4(t0, x0, y0, h)
        update_data(i + startIndex)
        re.append((x0,y0))
    return re

ans = rk4N(0,768.6,768.6,5,4000)
co2air = []
co2top = []
for state in ans:
    co2air.append(state[0])
    co2top.append(state[1])
    # print(state)

# Let's compare by plotting everything

co2air_real = [0.0409*climate.CO2air[i]*44.01 for i in range(2,2+len(co2air))]
xline = [5 * k for k in range(1,len(co2air)+1)]
plt.plot(xline,co2air, label = "CO2air")
plt.plot(xline,co2top, label = "CO2top")
plt.plot(xline,co2air_real, label = "CO2air_real")
plt.xlabel("time elapsed (mins)")
plt.ylabel("concentration (mg/m^3)")
plt.title("Changing of CO2air's and CO2top's concentration")
plt.legend()
plt.savefig("mygraph.png")