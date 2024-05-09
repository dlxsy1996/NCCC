import numpy as np
import  structure as st
import  physical as ph
import  friction as fr
from secondloopsolver import SecondLoopSolver,detTSolver
def SecondLoop(Qinit,T2in,mdsolve1,detT1,T1avg,T1in,T1out,i_T2in,Pb2,Pf2):
    i_T2out = 0
    T2out = 650
    while True:
        
        T2avg = (T2out+T2in)/2
        detT2 = T2out-T2in                                         #求主换冷测温差
        mdsolve2 = Qinit/(ph.C_Cp()*detT2)                         #求主换冷侧流量
        Ah1 = st.N1[3]* np.pi * (st.D1[3]/2)**2   * st.L1[3]       #管内换热面积   
        Ah2 = st.N2[0]* np.pi * (st.D2[0]/2)**2   * st.L2[0]       #管外换热面积

        #u1 = mdsolve1 / (ph.F_rho(T1avg)*st.N1[3]*np.pi*(st.D1[3]/2)**2)壳侧流动速度
        u1 = mdsolve1 / 0.137
        p12 = 0.019                                                #管道间距
        Dhe = 4*(p12**2-(np.pi/4)*(st.D1[3]**2))/(np.pi*st.D1[3])  #当量直径
        Re1 = Dhe * ph.F_rho(T1avg) * u1 / ph.F_mu(T1avg)  
        Pr1 = ph.F_Cp() * ph.F_mu(T1avg) / ph.F_kappa()
        #求壁温
        T1m = 0.4 * T1in +0.6 *T1out
        T2m = 0.4 * T2out +0.6 * T2in
        Tw12 = T1m * ph.F_kappa() + T2m * ph.C_kappa()
        #壳程换热系数  
        h1 =  ph.F_kappa() / Dhe * 0.378 * Re1**0.554 * Pr1**(1/3) * (ph.F_mu(T1avg)/ph.F_mu(Tw12))**0.14 * 1.18    #管内径690mm对应旁路挡板传热校正系数1.18


        
        u2 = mdsolve2 / (ph.C_rho(T2avg) * (st.N2[0]/2) * np.pi*(st.D2[0]/2)**2)    #管侧流动速度
        Re2 =st.D2[0] * ph.C_rho(T2avg) * u2 / ph.C_mu(T1avg)
        Pr2 = ph.C_Cp() * ph.C_mu(T2avg) / ph.C_kappa()
        #管程换热系数
        h2 = ph.C_kappa()/st.D1[3] * 0.023 * Re2**0.8 * Pr2**(1/3) *(ph.C_mu(T2avg)/ph.C_mu(Tw12))**0.14

        #其他热阻
        r1 = 0.00018     #管内结垢热阻
        r2 = 0.000172    #管外结垢热阻
        r = 0.000036     #管壁热阻


        h = (1/(Ah2/Ah1*(1/h1+r1)+(1/h2+r2)+r))
        Ah = st.N1[3] * np.pi * st.D1[3] * st.L1[3] *1.08   #1.08为换热管弯头校正系数

        #对数温差输出
        detT12 = Qinit/(h*Ah)
        T2new = detTSolver(detT12,T1in,T2in,T1out)


        if abs(T2new-T2out)<1e-6:
            T2out = T2new
            break
        else:
            T2out = (T2new +T2out)/2
            i_T2out = i_T2out + 1

    detT2 = T2out - T2in
    Pb2,Pf2 = SecondLoopSolver(mdsolve2,detT2,T2in)

    detP2 = Pb2 - Pf2
    h1_2 = h
    h11 = h1
    h21 = h2

    return T2out,detT2,Pb2,Pf2,mdsolve2,h1_2,h11,h21,i_T2in,i_T2out,detP2
