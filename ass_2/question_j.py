import numpy as np
import matplotlib.pyplot as plt
"""
This is selectivity analysis on: 
From the techno-economic analysis, it follows that for the process to be profitable the reactor needs to be run at a total current density of 2000 A/m2 (thoughput)

j) Implement th equation from i.) with the equation for FE(y) in Python and plot the conversion and faradaic efficiency over the reactor length. 
"""

ae    = 10e3            # [1/m] specific surface area
w     = 1e-3            # [m]
h     = 1e-3            # [m]
i_tot = 2000            # [A/m2]
nu    = 75e-5/60        # [m3/s]
F     = 96485           # [C/mol]
ne    = 2


L     =  5              #[m]    the height of the reactor (y)



def equation_i(x,y, params): 
    """
    ODE for the equation in question i. 
    
    Arg: 
        params:         (lst): list of all the parameters required 
            params = 
        x               (int):
        y               (int):

        returns: 

        """

    # unpack params
    ae, w, h, i_tot, nu, F, ne = params           

    FE = (F * ne * nu * co2)/ w*L*i_tot
    C_bulk = (ae * w * h) * ( FE * i_tot) / (nu * F * ne)
    
    return 