from co2 import *

#x0, y0 for Co2air, co2top, t0, t, h mean time start, time end, step
def euler(t0, x0, y0, h, t):
    n=int((t-t0)/h)
    for i in range (1, n+1):
        K1x=h*dxCO2_Air(data,x0, y0)
        K1y=h*dxCO2_Top(data, x0,y0)
        K2x=h*dxCO2_Air(data, x0+K1x, y0+K1y)
        K2y=h*dxCO2_Top(data, x0+K1x, y0+K1y)
        x0=x0+1.0/2.0*(K1x+K2x)
        y0=y0+1.0/2.0*(K1y+K2y)
    return x0, y0

def rk4(t0, x0, y0, h, t):
    n=int((t-t0)/h)
    for i in range (1, n+1):
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
    return x0, y0

for i in range(1, 20):
    print(rk4(0,443,443,0.5,5*i))