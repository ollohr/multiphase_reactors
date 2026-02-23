import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Given / constants
w     = 1e-3            # [m]
h     = 1e-3            # [m]
i_tot = 2000            # [A/m2]
nu    = 75e-6/60        # [m3/s]  (75 mL/min)
F     = 96485           # [C/mol]
ne    = 2
kL    = 8e-5            # [m/s]
ae    = 1/h             # [1/m]
L     = 5.0             # [m]

c_ref = 38.0            # [mol/m3] saturation concentration c* = cref

c0 = c_ref              # inlet bulk concentration

def FE_from_cbulk(cbulk, kL, i_tot, F, ne, c_ref):
    """
    FE = c_e/c_ref with c_e found by combining:
      i_CO = (c_e/c_ref)*i_tot  and  i_CO = ne*F*kL*(c_bulk - c_e)
    """
    denom = ne * F * kL * c_ref + i_tot
    FE = (ne * F * kL * cbulk) / denom
    # optional: clip physically (0..1)
    return np.clip(FE, 0.0, 1.0)

def dcbulk_dy(y, cbulk):
    c = cbulk[0]
    FE = FE_from_cbulk(c, kL, i_tot, F, ne, c_ref)
    dcdy = -(ae * w * h) * (FE * i_tot) / (nu * F * ne)
    return [dcdy]

# Solve ODE
y_eval = np.linspace(0, L, 400)
sol = solve_ivp(dcbulk_dy, t_span=(0, L), y0=[c0], t_eval=y_eval)

y = sol.t
c = sol.y[0]

X = 1.0 - c / c0                         # conversion
FE = FE_from_cbulk(c, kL, i_tot, F, ne, c_ref)

# Plot
fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)

ax[0].plot(y, X)
ax[0].set_ylabel("Conversion X(y)")
ax[0].grid(True)

ax[1].plot(y, FE)
ax[1].set_xlabel("Reactor length y [m]")
ax[1].set_ylabel("Faradaic efficiency FE(y)")
ax[1].set_ylim(0, 1.05)
ax[1].grid(True)

plt.tight_layout()
plt.show()