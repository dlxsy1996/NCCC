import numpy as np
import  structure as st
import  physical as ph
import  friction as fr

def SecondLoopSolver(mdsolve2,detT2,T2in):
    T2out = T2in+detT2
    T2avg = (T2in + T2out)/2
    Pb2 = (ph.C_rho(T2in)-ph.C_rho(T2out)) * 9.81 * (st.L2[0]+st.L2[1]-(st.L2[0]+st.L2[3])/2) 
    Pf2 = fr.PfAirCoolerTubeside(mdsolve2,T2avg)+ fr.PfHeatExchangeTubeside(mdsolve2,T2avg) +SysResist(mdsolve2,st.L2,st.D2,st.N2,st.K2)

    return Pb2,Pf2


def SysResist(md,L2,D2,N2,K2): # Calculates Pf(md) 
    T_est = np.array([530,650,650,530,530,530]) 
    v_he = md / (ph.C_rho(590) * st.N2[0]/2 * np.pi/4*(st.D2[0])**2)   
    friction = 0
    u2 = np.zeros(6)  
    j=0
    while j < (len(T_est)): # For every mass flow rate in md  
        # Computing the velocities in every section of the core - Mass Continuity 
        u2[0] = v_he
        for i in range(1,6):
            u2[i] = u2[0] * (ph.C_rho(T_est[0])/ph.C_rho(T_est[i])) * (N2[0]/st.N2[i]) * (D2[0]/D2[i])**2  
        ifric = 0
        for i in range(0,6): 
            # Calculate Reynold's number      
            Re = (ph.C_rho(T_est[i])*u2[i]*D2[i]) / (ph.C_mu(T_est[i]))
            # Calculate friction factor (Blasius Correlation for high Re) 
            if Re < 2100:
                f = 64/Re
            elif Re >= 2100 and Re < 30000:
                f = 0.316 * Re**-0.25
            elif Re >= 30000:
                f = 0.184*Re**-0.20
            # Calculate contribution of major and minor losses for each section  
            ifric += f/2*L2[i]/D2[i]*ph.C_rho(T_est[i])*(u2[i])**2 + K2[i]/2*ph.C_rho(T_est[i])*(u2[i])**2
        friction = ifric
        j += 1
    return friction
 


def  detTSolver(Tm,T1,t1,T2):
    #逆流换热
    t2_a = t1 + 1
    detT12_aTm = ((T1-t2_a)-(T2-t1))/ np.log((T1-t2_a)/(T2-t1)) - Tm  #正值
    t2_b = T1 - 1e-15
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

  