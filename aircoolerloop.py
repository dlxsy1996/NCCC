import numpy as np
import  structure as st
import  physical as ph
import  friction as fr
from aircoolersolver import AirCoolerSolver,detTSolver
def AirCooler(Qinit,T3in,mdsolve2,detT2,T2avg,T2in,T2out,i_T3in,Pb3,Pf3):
    i_T3out = 0
    T3out = 330
    while True:
        
        T3avg = (T3out+T3in)/2
        detT3 = T3out-T3in                               #求空冷器冷测温差
        mdsolve3 = Qinit/(ph.A_Cp(T3avg)*detT3)          #求空冷器冷侧流量
        #空冷器结构
        Sf = 0.0025            #翅片间距
        Sh = 0.016             #翅片高度
        Sr = 0.0005            #翅片厚度
        Sd = 0.062             #管心距
        Sn = 6                 #管排数

        Aa1 = st.N2[3]* np.pi * (st.D2[3]/2)**2 * st.L2[3]         #管内换热面积   
        Aa2 = st.N3[0]* np.pi * (st.D3[0]/2)**2 * st.L3[0]         #管外换热面积
        At = ((st.N3[0]*st.L3[0]/Sf)*(np.pi/4 * ((st.D3[0]+Sh)**2-(st.D3[0]**2 )))*2 + np.pi * (st.D3[0]+Sh)*Sr)*st.N3[0]*st.L3[0]  #翅片管表面积

        #求壁温
        T2m = 0.4 * T2in +0.6 *T2out
        T3m = 0.4 * T3out +0.6 * T3in
        Tw23 = T2m * ph.C_kappa() + T3m * ph.A_kappa(T3avg)

        #管程换热系数   
        u1 = mdsolve2 / (ph.C_rho(T2avg) * (st.N2[3]) * np.pi*(st.D2[3]/2)**2)    #管侧流动速度
        Re1 =st.D2[3] * ph.C_rho(T2avg) * u1 / ph.C_mu(T2avg)
        Pr1 = ph.C_Cp() * ph.C_mu(T2avg) / ph.C_kappa()
        h1 = ph.C_kappa()/st.D3[0] * 1.86 * Re1**(1/3) * (st.D2[3]/st.L2[3])**(1/3) * Pr1**(1/3) *(ph.C_mu(T2avg)/ph.C_mu(Tw23))**0.14

        #壳程换热系数 
        u2 = mdsolve3 / (ph.A_rho(T3avg) * st.L3[0] * (st.N3[0]/Sn) * Sd)    #壳侧流动速度
        Re2 = st.D3[0] * ph.A_rho(T3avg) * u2 /ph.A_mu(T3avg)
        Pr2 = ph.A_Cp(T3avg) * ph.A_mu(T3avg) / ph.A_kappa(T3avg)
        #h2 = 0.1378 * (ph.A_kappa(T3avg)/st.D3[0]) * Re2**0.718 * Pr2**0.333 * (Sf/Sh)**0.296 * (At/Aa2)
        h2 = 0.076*ph.A_kappa(T3avg)/st.D3[0]*Re2**0.683*(At/Aa2)


        #其他热阻
        r1 = 0.000172    #管内结垢热阻
        r2 = 0.00018     #管外结垢热阻
        rp = 0.0000945   #翅片热阻
        rf = 0.000012     #翅片间隙热阻
        h = (1 / ((Aa2/Aa1) * (1/h1) + r1 + 1/h2 + r2 + rp + rf))
        

        detT23 = Qinit/(h*At)
        T3new = detTSolver(detT23,T2in,T3in,T2out)


        if abs(T3new-T3out) <1e-6:
            T3out = T3new
            break
        else:
            T3out = (T3new +T3out)/2
            i_T3out = i_T3out + 1

            
    detT3 = T3out - T3in
    Pb3,Pf3 = AirCoolerSolver(mdsolve3,detT3,T3in)
    detP3 = Pb3 - Pf3

    h2_3 = h   
    h22 = h1
    h33 = h2 
    return T3out,detT3,Pb3,Pf3,mdsolve3,h2_3,h22,h33,i_T3in,detP3
