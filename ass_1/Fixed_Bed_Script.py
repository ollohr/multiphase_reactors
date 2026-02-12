# -*- coding: utf-8 -*-
"""

"""

#import libraries
from scipy import integrate
import numpy as np
import matplotlib.pyplot as plt

# function
def funstep1 (z,y):
    FA = y[0]
    FB = y[1]
    FC = y[2]
    FT_v = ****
    CAb  = ****
    CAs  = ****
    
    dy1    = ****
    dy2    = ****
    dy3    = ****
    return np.array([dy1,dy2,dy3])

# constant parameters
R      = 8.31441    # Universal gas constant, J/(mol K)
W      = 20000        # Catalyst weight, kg
rhocat = 915        # Density of the catalyst particle, kgcat/m3
dp     = 910e-6       # Particle diameter, m
K1     = 9.6e-4     # rate constant of reaction, m3/(kgcat s)
kf     = 2.45e-4       # mass transfer coefficient, m/s
P      = 1.01e5     # Total pressure, Pa
Tin    = 920        # Temperatur, K


# Stoichiometric coefficients of the reaction
nuA1 = ****
nuB1 = ****
nuC1 = ****

zspan = np.linspace(0,1,100)

FAin  = 40
FBin  = 0
FCin  = 0

# inital values
y0 = [FAin, FBin, FCin]



# solve
y = np.zeros((len(zspan), len(y0))) #array for solution
y [0,:] = y0
r = integrate.ode(funstep1).set_integrator("dopri5") # choice of method
r.set_initial_value(y0,0) #inital values
for i in range(1, zspan.size):
    y[i,:]=r.integrate(zspan[i])
    if not r.successful():
        raise RuntimeError("Could not integrate")

#Calculate values
FA     = y[:,0]
FB     = y[:,1]
FC     = y[:,2]

FT_v   = ****
CAb    = ****
CBb    = ****
CCb    = ****

CAs    = ****
CBs    = ****
CCs    = ****


# Plot
# Figure 1
plt.figure(1)
plt.plot(zspan*W, y[:,0],'*-', color='blue')
plt.plot(zspan*W, y[:,1], '-+', color='red')
plt.plot(zspan*W, y[:,2],'-^', color='orange')
plt.ylabel('$F_i$ [mol s$^{-1}$]')
plt.xlabel('W [kg$_{cat}$]')
plt.legend(['F$_{HN_3}$','F$_{H_2}$','F$_{N_2}$'],loc='upper right')
plt.show()
# Figure 2
plt.figure(2)
plt.plot(zspan*W, CAs,'*-', color='blue')
plt.plot(zspan*W, CBs, '-+', color='red')
plt.plot(zspan*W, CCs,'-^', color='orange')
plt.ylabel('C$_{surface}$ [mol m$^{-3}$]')
plt.xlabel('W [kg$_{cat}$]')
plt.legend(['C$_{NH_{3},s}$','C$_{H_{2},s}$','C$_{N_{2},s}$'],loc='upper right')
plt.show()
#Figure 3
plt.figure(3)
plt.plot(zspan*W, CAb,'*-', color='blue')
plt.plot(zspan*W, CBb, '-+', color='red')
plt.plot(zspan*W, CCb,'-^', color='orange')
plt.ylabel('C$_{bulk}$ [mol m$^{-3}$]')
plt.xlabel('W [kg$_{cat}$]')
plt.legend(['C$_{NH_{3},b}$','C$_{H_{3},b}$','C$_{N_{2},b}$'],loc='upper right')
plt.show()
#Figure 4
plt.figure(4)
plt.plot(zspan*W, conv,'-', color='red')
plt.ylabel('Conversion [-]')
plt.xlabel('W [kg$_{cat}$]')
plt.show()