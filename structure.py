import numpy as np
#确定的设备设计  堆芯、主换、二换、空冷
#PHYSICAL SETUP                                                                        
#Main loop
#[0]Heater;  [1]Hot leg;  [2]Cooler;  [3]Cold leg;  [4]Horz Hot leg;  [5]Horz Cold leg 

#Second loop
#[0]Heater;  [1]Hot leg;  [1]Cooler;  [2]Cold leg;  [4]Horz Hot leg;  [5]Horz Cold leg 

#air cooler
#[1]air cooler;  [2]Top Tower;  [3]Atmosphere; [3]down Tower;


#回路结构H、L、D、N 
#Diameters of each flow path of the loop [m]     环路各流路直径 [m]  !水力直径  当量直径额外计算
global D1,D2,D3
D1=np.zeros(6)
D1[0] = 0.06            #堆芯石墨通道直径
D1[1] = 0.3
D1[2] = 0.3
D1[3] = 0.014             #管外径14mm
D1[4] = 0.3
D1[5] = 0.3
D2=np.zeros(6)
D2[0] = 0.0104             #管内径
D2[1] = 0.3
D2[2] = 0.3
D2[3] = 0.021            #等效面积0.2m2
D2[4] = 0.3
D2[5] = 0.3
D3=np.zeros(4)
D3[0] = 0.027            
D3[1] = 3
D3[2] = 0              #无用值
D3[3] = 3               
# Number of flow paths (tubes) in each section 每个部分的流路（管）数量
global N1,N2,N3
N1=np.zeros(6)
N1[0] = 115
N1[1] = 1
N1[2] = 1
N1[3] = 844
N1[4] = 1
N1[5] = 1
N2=np.zeros(6)
N2[0] = 844
N2[1] = 1
N2[2] = 1
N2[3] = 180
N2[4] = 1
N2[5] = 1 
N3=np.zeros(4)
N3[0] = 180
N3[1] = 1
N3[2] = 1
N3[3] = 1
# Coefficient of friction in each leg  每根管局部阻力
global K1,K2,K3
K1=np.zeros(6)
K1[0] = 0
K1[1] = 0.5
K1[2] = 0.5
K1[3] = 0
K1[4] = 0.5
K1[5] = 0.5
K2=np.zeros(6)
K2[0] = 0
K2[1] = 0.5
K2[2] = 0.5
K2[3] = 0
K2[4] = 0.5
K2[5] = 0.5 
K3=np.zeros(4)
K3[0] = 0
K3[1] = 0
K3[2] = 0
K3[3] = 0
# Lengths of each section of the loop [m]   环路各段的长度 [m]
global L1,L2,L3
L1=np.zeros(6)
L1[0] = 3.3
L1[1] = 6.7
L1[2] = 5
L1[3] = 2.2
L1[4] = 7.8
L1[5] = 5
L2=np.zeros(6)
L2[0] = 2.2
L2[1] = 7.8
L2[2] = 5
L2[3] = 3.0
L2[4] = 7.0
L2[5] = 5 
L3=np.zeros(4)
L3[0] = 3
L3[1] = 17
L3[2] = 20       
L3[3] = 0