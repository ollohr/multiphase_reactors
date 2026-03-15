# -*- coding: utf-8 -*-
"""
"""

#import libraries
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


# function
def funstep1 (z,t):
    c1, c2 = z          # where c1 = [CO_2] and c2 = [OH]

    #Henry_constant for electrolyte solution
    H_L = 3.3e-2
    
    r1 = k1*c1* c2
    co2_ref = H_L * 0.30 # since 30% of co2 in gas phase 
    dc1dt = kla*(co2_ref - c1)-r1
    dc2dt = -r1    

    return [dc1dt,dc2dt]

 
#Temperature
T= 298.15  #K
#diameter column
d = 1.5 #m
#ideal constant
R = 8.314 #J/mol/K;
#gravitational constant
g = 9.81 #m/s²
#considered time interval
t_Max = 10 #10 s

# initial concentration of NaOH
C0_NaOH = 0.04 #mol/l

# flow rate CO2
FCO2 = 70 #l/min

# concentration of CO2 in gas phase --> assumption: ideal gas 
C0_CO2 = (0.3*101325)/(R*T)/1000 #mol/l                     ######################### change this weird as 

print(C0_CO2)
# physical properties
ug = 6.6e-4  #m/s
kla = 0.467*ug**0.82  #1/s

sigma = 0.073 #N/m 
rho_L = 1000 # kg/m^3
k1= 4  #l/(mol*s)


# initial conditions
init= [0, 0.04]                                                           
t_span = np.linspace(0, t_Max) 

# calculation
z = odeint(funstep1, init, t_span)
c1 = z[:, 0]
c2 = z[:, 1]

# Plot
plt.plot(t_span, c1, label='c1(t)')
plt.xlabel("Time (s)", weight = "bold")
plt.ylabel("[CO2] (mol/l)", weight = "bold")
plt.title("Concentration of CO2 over Time")
plt.minorticks_on()
plt.grid(which= "major", linewidth = 1)
plt.grid(which="minor", linewidth = 0.5)
plt.xlim(0,t_Max)
plt.ylim(0, c1[-1])
plt.legend()
plt.show()
plt.plot(t_span, c2, label='c2(t)')
plt.xlabel("Time (s)", weight = "bold")
plt.ylabel("[OH-] (mol/l)", weight = "bold")
plt.title("Concentration of OH- over Time")
plt.minorticks_on()
plt.grid(which= "major", linewidth = 1)
plt.grid(which="minor", linewidth = 0.5)
plt.xlim(0,t_Max)
# plt.ylim(0,c2[-1])
plt.legend()
plt.show()
