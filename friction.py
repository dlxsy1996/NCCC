import numpy as np
import  structure as st
import  physical as ph
import  friction as fr
#堆芯阻力
def PfCore(md,T1avg):

    friction = 6.14*(md**2) /(2*ph.F_rho(T1avg)*(( np.pi* (st.D1[0]/2)**2 )*st.N1[0])**2)

    return friction

#换热器壳侧阻力(shellside)管侧阻力(tubeside)
def PfHeatExchangeShellside(md,T1avg):

    friction= 2.07 * (md**2) /(2* ph.F_rho(T1avg)*(( np.pi* (st.D1[3]/2)**2 )*st.N1[3])**2)

    return friction

def PfHeatExchangeTubeside(md,T2avg):

    friction= 11.6 * (md**2) / (2 * ph.C_rho(T2avg) * (( np.pi* (st.D2[0]/2)**2 )*st.N2[0]/2)**2)

    return friction

#空冷换热器管侧阻力（tubeside）壳侧阻力（shellside）
def PfAirCoolerTubeside(md,T2avg):
    
    friction = 3.82 * (md**2) / (2* ph.C_rho(T2avg) * (( np.pi* (st.D2[3]/2)**2 )*st.N2[3])**2)
               
    return friction

def PfAirCoolerShellside(md,T3avg):
    s=0.435    #窄隙流通面积与迎风面积比
    friction=5.73 * (md**2) / (2 * ph.A_rho(T3avg) * (2 * 1.86 / s)**2)

    return friction


#空冷塔阻力
def Pfaircoolingtower(md,T3avg):
      
    friction=5.7 * (md**2) /(2 * ph.A_rho(T3avg)  * (np.pi/4 * st.D3[1])**2)

    return friction

