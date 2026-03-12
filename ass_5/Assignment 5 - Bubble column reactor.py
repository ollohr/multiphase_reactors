# -*- coding: utf-8 -*-
"""
"""

#import libraries
from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt

# function
def funstep1 (z,t):
    c1, c2 = z

    #Henry_constant for electrolyte solution
    H_L = ***
    r1 = ***
    dc1dt = ***
    dc2dt = ***    
    
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
C0_CO2 = *** #mol/l

# physical properties
ug = ***  #m/s
kla = 0.467*ug**0.82  #1/s

sigma = 0.073 #N/m 
rho_L = 1000 # kg/m^3
k1= 4  #l/(mol*s)


# initial conditions
init= [***, ***]
t_span = np.linspace(0, t_Max) 

# calculation
z = odeint(funstep1, init, t_span)
c1 = z[:, 0]
c2 = z[:, 1]

# Plot
plt.plot(t_span, c1, label='c1(t)')
plt.show()
plt.plot(t_span, c2, label='c2(t)')
plt.show()
