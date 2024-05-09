import numpy as np
import  structure as st
import  physical as ph
import  friction as fr

def AirCoolerSolver(mdsolve3,detT3,T3in):
    T3out = T3in+detT3
    T3avg = (T3in + T3out)/2

    Pb3 = 9.81 * (ph.A_rho(T3in)*(st.L3[2]-st.L3[3]) -  ph.A_rho(T3out)*st.L3[1] - ph.A_rho(T3avg)*st.L3[0])
    a = fr.PfAirCoolerShellside(mdsolve3,T3avg)
    b = fr.Pfaircoolingtower(mdsolve3,T3avg)
    c = SysAcceleration(mdsolve3,st.D3,T3in,T3out)
    Pf3 = fr.PfAirCoolerShellside(mdsolve3,T3avg)+ fr.Pfaircoolingtower(mdsolve3,T3avg)
    return Pb3,Pf3

def SysAcceleration(mdsolve3,D3,T3in,T3out):
    friction = mdsolve3**2/ ((np.pi/4*(st.D3[1]/2)**2)**2) *(1/ph.A_rho(T3out)-1/ph.A_rho(T3in))

    return friction


def  detTSolver(Tm,T1,t1,T2):
    t2_a = t1 + 1
    detT12_aTm = ((T1-t2_a)-(T2-t1))/ np.log((T1-t2_a)/(T2-t1)) - Tm  #正值
    t2_b = T1 - (1e-15)
    detT12_bTm = ((T1-t2_b)-(T2-t1))/ np.log((T1-t2_b)/(T2-t1)) - Tm  #负值
    t2new = t2_b - (detT12_bTm/(detT12_aTm-detT12_bTm))*(t2_a-t2_b)

    while True:
        detT12_Tm = ((T1-t2new)-(T2-t1))/ np.log((T1-t2new)/(T2-t1)) - Tm
        if abs(detT12_Tm) or abs(T1-t2new) <1e-6:
            t2 = t2new
            break
        elif detT12_Tm > 1e-6:
            t2_a = t2new
            if abs(t2_a-t2_b)<1e-6:
                break
            else:
                t2new = t2_b - (detT12_bTm/(detT12_aTm-detT12_bTm)) * (t2_a-t2_b)
        elif detT12_Tm < 1e-6:
            t2_b = t2new
            if abs(t2_b-t2_a) <1e-6:
                break
            else:
                t2new = t2_b - (detT12_bTm/(detT12_aTm-detT12_bTm)) * (t2_a-t2_b)
    return t2new
