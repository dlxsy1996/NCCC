import numpy as np
import  structure as st
import  physical as ph
import  friction as fr
from mainloopsolver import MainLoopSolver
from secondloop import SecondLoop
from aircoolerloop import AirCooler

#“0”初步假设 结构、物性、入口温度
T1in=600                     #恒定量
detT1=30                     #设计变量1
T1out=T1in+detT1             #求解一回路热端温度
T1avg= (T1in + T1out)/2
Qinit, mdsolve1,mdysolve1 = MainLoopSolver(detT1, T1in, T1out)      #求解一回路功率、压降和质量流量


#“1”熔盐自然循环余热排出回路计算循环
i_T2in = 0
Pb2 = 0
Pf2 = 0
T2in_a = 530
T2out,detT2,Pb2,Pf2,mdsolve2, h1_2,h11,h21,i_T2in, i_T2out ,detP2= SecondLoop(Qinit,T2in_a,mdsolve1,detT1,T1avg,T1in,T1out,i_T2in,Pb2,Pf2) 
detP2_a = detP2   #负值
T2in_b = 400
T2out,detT2,Pb2,Pf2,mdsolve2, h1_2,h11,h21,i_T2in, i_T2out ,detP2= SecondLoop(Qinit,T2in_b,mdsolve1,detT1,T1avg,T1in,T1out,i_T2in,Pb2,Pf2) 
detP2_b = detP2   #正值
T2innew = T2in_a - (detP2_a/(detP2_b-detP2_a))*(T2in_b-T2in_a)

while True:

    i_T2in = i_T2in + 1
    T2out,detT2,Pb2,Pf2,mdsolve2, h1_2,h11,h21,i_T2in,i_T2out,detP2= SecondLoop(Qinit,T2innew,mdsolve1,detT1,T1avg,T1in,T1out,i_T2in,Pb2,Pf2) 
    
    if abs(detP2) <1e-6:
        T2in = T2innew
        break
    elif detP2 > 1e-6:
        T2in_b = T2innew
        if abs(T2in_b-T2in_a) <1e-6:
            break
        else:
            T2innew = T2in_a - (detP2_a/(detP2_b-detP2_a)) * (T2in_b-T2in_a)
    elif detP2 < 1e-6:
        T2in_a = T2innew
        if abs(T2in_a-T2in_b) <1e-6:
            break
        else:
            T2innew = T2in_a - (detP2_a/(detP2_b-detP2_a)) * (T2in_b-T2in_a)
T2in =T2innew
T2avg= (T2in + T2out)/2

#“2”空气自然循环余热排出系统计算循环
T3in=20                    
i_T3in=0
i_T3out = 0
Pb3 = 0
Pf3 = 0
T3in_a = 300 
T3out,detT3, Pb3, Pf3 ,mdsolve3,h2_3,h22,h33,i_T3in,detP3 = AirCooler(Qinit,T3in_a,mdsolve2,detT2,T2avg,T2in,T2out,i_T3in,Pb3,Pf3) 
detP3_a = detP3   #负值
T3in_b = -200
T3out,detT3,Pb3,Pf3,mdsolve3,h2_3,h22,h33,i_T3in,detP3 = AirCooler(Qinit,T3in_b,mdsolve2,detT2,T2avg,T2in,T2out,i_T3in,Pb3,Pf3) 
detP3_b = detP3   #正值
T3innew = T3in_a - (detP3_a/(detP3_b-detP3_a))*(T3in_b-T3in_a)

#空冷器换热换热面积计算
Saircooler = st.N3[0] * np.pi * st.D3[0] * st.L3[0]
while True:

    i_T3in= i_T3in+1
    T3out,detT3,Pb3,Pf3,mdsolve3,h2_3,h22,h33,i_T3in,detP3 = AirCooler(Qinit,T3innew,mdsolve2,detT2,T2avg,T2in,T2out,i_T3in,Pb3,Pf3)      #求解三回路压降、流量和温度差
    
    if abs(detP3) <1e-6:
        T3in = T3innew
        break
    elif detP3 > 1e-6:
        T3in_b = T3innew
        if abs(T3in_b-T3in_a) <1e-6:
            break
        else:
            T3innew = T3in_a - (detP3_a/(detP3_b-detP3_a)) * (T3in_b-T3in_a)
    elif detP3 < 1e-6:
        T3in_a = T3innew
        if abs(T3in_a-T3in_b) <1e-6:
            break
        else:
            T3innew = T3in_a - (detP3_a/(detP3_b-detP3_a)) * (T3in_b-T3in_a)
T3in =T3innew
T3avg= (T3in + T3out)/2   
        

#数据汇总
print("回路自然循环能力=",Qinit)
print("---------------------------")
print("一回路自然循环能力=",Qinit)
print("一回路冷端温度=",T1in)
print("一回路热端温度=",T1out)
print("一回路流量=", mdsolve1)
print("---------------------------")
print("二回路自然循环能力=",Qinit)
print("二回路冷端温度=",T2in)
print("二回路热端温度=",T2out)
print("二回路流量=", mdsolve2)
print("---------------------------")
print("空气冷却回路自然循环能力=",Qinit)
print("空气冷却回路冷端温度=",T3in)
print("空气冷却回路热端温度=",T3out)
print("空气冷却回路流量=", mdsolve3)
print("空冷塔高度设计=", st.L3[1]+st.L3[0])
print("---------------------------")
print("一回路压降=", mdysolve1)
print("二回路压降=", Pb2)
print("空气冷却回路压降=",Pb3)
print("---------------------------")
print("主换热器换热系数=", h1_2)
print("空气冷却换热器换热系数=", h2_3)
print("空气冷却换热器换热面积=",Saircooler)
print("---------------------------")
print("二回路冷端温度循环初值=", 530)
print("二回路冷端温度循环终值=", T2in)
print("二回路冷端温度循环范围=", 530,"-",T2in)
print("二回路冷端温度循环次数=", i_T2in)
print("---------------------------")
print("二回路热端温度循环初值=", 650)
print("二回路热端温度循环终值=", T2out)
print("二回路热端温度循环范围=", 650,"-",T2out)
print("二回路热端温度循环次数=", i_T2out)
print("---------------------------")
print("三回路冷端温度循环初值=", 20)
print("三回路冷端温度循环终值=", T3in)
print("三回路冷端温度循环范围=", 20,"-",T3in)
print("三回路冷端温度循环次数=", i_T3in)
print("---------------------------")
print("三回路热端温度循环初值=", 330)
print("三回路热端温度循环终值=", T3out)
print("三回路热端温度循环范围=", 330,"-",T3out)
print("三回路热端温度循环次数=", i_T3out)
print("------------end---------------")

