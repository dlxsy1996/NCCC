import numpy as np

#物理性质
###FUEL PROPERTIES###
F_rho = lambda T: 2639-0.5286*(T)                                  # Density (kg/m^3)
F_mu = lambda T: 0.0036*np.exp(7043.24/(T+273.15))*0.001           # Dynamic Viscosity (Pa-s) 
F_kappa = lambda : 1.4                                             # Thermal Conductivity (W/m-K) 
F_beta = lambda T: 1/(4992.4-T)                                    # Thermal Expansion Coefficient (1/K) 
F_Cp = lambda : 2000                                               # Specific Heat (J/kg-K) 


###COOL PROPERTIES###
C_rho = lambda T: 2270-0.37*(T)                                    # Density (kg/m^3)
C_mu = lambda T: 0.0346*np.exp(5165.0/(T+273.15))*0.001            # Dynamic Viscosity (Pa-s) 
C_kappa = lambda : 0.87                                            # Thermal Conductivity (W/m-K) 
C_beta = lambda T: 1/(6135-T)                                      # Thermal Expansion Coefficient (1/K) 
C_Cp = lambda : 2177                                               # Specific Heat (J/kg-K) 

###GALLIUM PROPERTIES###  #空气物性                                                                                   
A_rho =lambda T: 3.48816*(10**-3)*86000/(T+273.15)                                                                                                                              # Density (kg/m^3) 
A_mu = lambda T: 5.8449*(10**-8)*(T+273.15)                                                                                                                                     # Dynamic Viscosity (Pa-s)
A_kappa = lambda T: 0.02389+0.100857*(10**-3)*(T)-0.28571*(10**-6)*((T)**2)                                                                                                     # Thermal Conductivity (W/m-K) 
A_Cp = lambda T:(1.00306+0.02413*(10**-3)*(T)+0.4283*(10**-6)*(T**2)+0.03868*(10**-9)*(T**3)-0.95024*(10**-12)*(T**4)+0.89676*(10**-15)*(T**5)-0.25726*(10**-18)*(T**6))*1000   # Specific Heat (J/kg-K) 
A_beta = 0.003               


