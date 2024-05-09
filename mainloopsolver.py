import numpy as np
import  structure as st
import  physical as ph
import  friction as fr

def MainLoopSolver(detT1, T1in, T1out):
    dx = 0.1
    T1avg= (T1in + T1out)/2                                                                   #定性温度
    m=np.pi/4*(st.D1[4]**2)*ph.F_rho(T1in)*dx                                                 #节块质量
    x = np.arange(0,st.L1[4]+st.L1[5]+dx,dx)                                                  #[4]+[5]长度
    mdl=1000                                                                                  #质量流量数组长度
    dist =int((st.L1[0]/dx+st.L1[1]/dx+st.L1[2]/dx+st.L1[3]/dx+st.L1[4]/dx+st.L1[5]/dx) + 1)  #总数组长度
    Vzong=ZongTiJi(st.L1,st.D1,st.N1)                                                         #总体积计算
    Vjisuan = Vzong - np.pi/4*(st.D1[3]**2)*st.L1[3]*st.N1[3]                                 #计算体积
    Vzongreal = Vjisuan + np.pi/4*(st.D1[3]**2)*st.L1[3]*st.N1[3] 
    # Initializing the arrays out of the loop  初始化数组
    md = np.ones((mdl))
    Pb1 = np.ones((mdl))
    Pf1 = np.ones((mdl))
    Qv = np.ones((mdl))                                                                       #功率密度数组   
    Qinit = np.ones((mdl))
    Tdist = np.ones((dist))
    xdist = np.ones((dist))

    #冷端温度、流量分布计算
    Tcoldleg =  np.ones(((len(x)),(mdl)))*T1in                                                #[4]+[5]即冷端总长温度，流量
    Tout =  np.ones((len(x)))*T1in                                                            #换热器出口温度=T1in
    md[:] = np.linspace(0.1,100,mdl)                                                          #质量流量范围0.1-40

    #换热器换热总功率 对应每个质量流量
    Tcoldleg[0][:] = T1in 
    Qin = (T1out-T1in)*ph.F_Cp()*md[:]                  
    #阻力线求解
    #Pf1 =SysResist(md[:],st.L1,st.D1,st.N1,st.K1) 
    Pf1 =(fr.PfCore(md[:],T1avg)+ fr.PfHeatExchangeShellside(md[:],T1avg) +SysResist(md[:],st.L1,st.D1,st.N1,st.K1))
    
    #循环主体
    while True:
        for k in range(0,mdl):                       #遍历质量流量
            for j in range(1,(len(x))):              #遍历冷端每个节块
                Qv = Qin[k]/Vjisuan                  #不同流量(0.1-40)对应的功率密度 [W/m3]  
                v=np.pi/4*(st.D1[4]**2)*dx           #冷端每个节块的体积
                m1= v*ph.F_rho(T1in)                 #冷端每个节块的质量
                a = (m1*j)*ph.F_Cp()*(Tcoldleg[0][k])
                b = j*Qv*v*(m1*j)/md[k]
                c = a + b
                Tcoldleg[j][k] = (c)/(m1*j)/ph.F_Cp()
            Tdist[:],xdist[:],Pb1[k] = LoopDist(Tcoldleg[:,k],st.N1,st.D1,md[k],Qin[k],T1out,T1in)   #求解整段T，求Pb
        mdsolve1, mdysolve1 = MDSolver(Pb1,Pf1,md)  #求解md
        Qfalse = (T1out-T1in)*ph.F_Cp()*mdsolve1    #错误Q输出
        Qinit = Qfalse * (Vzongreal/Vjisuan)        #实际Q输出
        break
    
    return Qinit, mdsolve1,mdysolve1


def ZongTiJi(L,D,N):
    Vzong=0.0
    for i in range(0,6):
        Vzong+=np.pi/4*(D[i]**2)*L[i]*N[i]
    return Vzong


def LoopDist(Tcoldleg,N,D,md,Q,T1out,T1in): 
    L0a = 0
    L0b = st.L1[0]-0.1 #3.2 
    num_L0 = int(st.L1[0]*10) #33 
    L1a = st.L1[0]
    L1b = st.L1[0] + st.L1[1]-0.1
    num_L1 = int(st.L1[1]*10) #67
    L2a = st.L1[0] + st.L1[1] 
    L2b = st.L1[0] + st.L1[1] + st.L1[2] - 0.1
    num_L2 = int(st.L1[2]*10)
    L3a = st.L1[0] + st.L1[1] + st.L1[2]
    L3b = L3a + st.L1[3] - 0.1
    num_L3 = int(st.L1[3]*10)
    L4a = L3a + st.L1[3]
    L4b = L3a + st.L1[3] + st.L1[4] - 0.1
    num_L4 = int(st.L1[4]*10)
    L5a = L3a + st.L1[3] + st.L1[4] 
    L5b = L3a + st.L1[3] + st.L1[4] +st.L1[5]
    num_L5 = int(st.L1[5]*10+1)

    k12a = 0
    k12b = st.L1[1] +st.L1[2]
    num_k12 = int((st.L1[1]+st.L1[2])*10)
    k3a = 0
    k3b = st.L1[3]
    num_k3 = int(st.L1[3]*10)
    k45a = 0
    k45b = st.L1[4]+st.L1[5]
    num_k45 = int((st.L1[4]+st.L1[5])*10)

    x0 = np.linspace(L0a,L0b,num_L0)
    x1 = np.linspace(L1a,L1b,num_L1)
    x2 = np.linspace(L2a,L2b,num_L2)
    x3 = np.linspace(L3a,L3b,num_L3)
    x4 = np.linspace(L4a,L4b,num_L4)
    x5 = np.linspace(L5a,L5b,num_L5)
    x12a = np.linspace(k12a,k12b,num_k12)
    x3a = np.linspace(k3a,k3b,num_k3)
    x45a = np.linspace(k45a,k45b,num_k45)

    L_hl = int(num_k12 - 1)
    L_cl = int(num_k45 - 1)
    La_Pbjia = 0
    Lb_Pbjia = int(L1b*10)
    Lc_Pbjian = int(L2b*10+1)
    Ld_Pbjian =  int(L4b*10)   
    num_Tcooler = int(k3b*10)
    T1avg = (T1in +T1out)/2

    Thotleg = np.ones(len(x12a))*T1out - detT12(x12a,D,md,Q)


    Tcooler = np.ones((len(x3a))) 
    for i in range(0,(len(x3a))):
        Tcooler[i] = Thotleg[L_hl]-i*(Thotleg[L_hl]-Tcoldleg[0])/num_Tcooler
    

    Tcore = np.ones((len(x0))) 
    for i in range(0,(len(x0))):
        Tcore[i] = Tcoldleg[L_cl]+i*(Thotleg[0]-Tcoldleg[L_cl])/num_L0  

    xdist = np.concatenate((x0,x1,x2,x3,x4,x5))
    Tdist = np.concatenate((Tcore,Thotleg,Tcooler,Tcoldleg[:]))


    Pb1 = ph.F_beta(T1avg)*ph.F_rho(T1avg)*9.81*np.trapz(Tdist[La_Pbjia:Lb_Pbjia],xdist[La_Pbjia:Lb_Pbjia])  
    Pb1 += ph.F_beta(T1avg)*ph.F_rho(T1avg)*9.81*np.trapz(Tdist[Lc_Pbjian:Ld_Pbjian],xdist[Lc_Pbjian:Ld_Pbjian])*(-1)

    return Tdist, xdist , Pb1


def detT12(x12a,D,md,Q):
    Vzong=ZongTiJi(st.L1,st.D1,st.N1)
    Vjisuan = Vzong - np.pi/4*(st.D1[3]**2)*st.L1[3]*st.N1[3]  
    Qv = Q/Vjisuan
    detT12 = np.ones(len(x12a))
    for j in range(0,len(x12a),1):
        dx=0.1
        v=np.pi/4*(D[1]**2)*dx   
        detT12[j] = (len(x12a)-j-1)*Qv*v/ph.F_Cp()/md
    return detT12


def MDSolver(Pb,Pf,md): # Linear interpolation/intersection of Pf = Pb  Pf = Pb 的线性交点
    for j in range(2,(len(md)-1)):
        if (Pb[j] < Pf[j]) and (Pb[j-1] > Pf[j-1]):
            sb = (Pb[j]-Pb[j-1])/(md[j]-md[j-1])
            sf = (Pf[j]-Pf[j-1])/(md[j]-md[j-1])
            mdsolve = (Pb[j-1]-Pf[j-1])/(sf-sb)+md[j-1]
            mdysolve = Pf[j-1]+sf*(mdsolve-md[j-1])
    return mdsolve, mdysolve         #横坐标质量流量 纵坐标压降交点


def SysResist(md,L1,D1,N1,K1): # Calculates Pf(md) 
    T_est = np.array([600,700,700,600,600,600]) 
    v_core = md/(ph.F_rho(700)*N1[0]*np.pi/4*(D1[0])**2)   
    friction = np.ones((len(md)))
    u1 = np.zeros(6)  
    j=0
    while j < (len(md)): # For every mass flow rate in md  
        # Computing the velocities in every section of the core - Mass Continuity 
        u1[0] = v_core[j]
        for i in range(1,6):
            u1[i] = u1[0]*(ph.F_rho(T_est[0])/ph.F_rho(T_est[i]))*(N1[0]/N1[i])*(D1[0]/D1[i])**2  

        ifric = 0
        for i in range(0,6):
            Re = 1  
            # Calculate Reynold's number          
            Re=(ph.F_rho(T_est[i])*u1[i]*D1[i])/(ph.F_mu(T_est[i]))
            # Calculate friction factor (Blasius Correlation for high Re) 
            f = 1
            if Re < 2100:
                f = 64/Re
            elif Re >= 2100 and Re < 30000:
                f = 0.316 * Re**-0.25
            elif Re >= 30000:
                f = 0.184*Re**-0.20
            # Calculate contribution of major and minor losses for each section  
            ifric += f/2*L1[i]/D1[i]*ph.F_rho(T_est[i])*(u1[i])**2 + K1[i]/2*ph.F_rho(T_est[i])*(u1[i])**2
        friction[j] = ifric
        j += 1
    return friction
