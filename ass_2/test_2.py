import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

ae    = 10e3            # [1/m] specific surface area
w     = 1e-3            # [m]
h     = 1e-3            # [m]
i_tot = 2000            # [A/m2]
nu    = 75e-5/60        # [m3/s]
F     = 96485           # [C/mol]
ne    = 2
kl    = 8e-5            # [m/s]
c_ref = 38              # [mol/m3]
L     = 5               # [m] reactor length

params = [ae, w, h, i_tot, nu, F, ne, kl, c_ref]

def equation_i(y, x, params):
    ae, w, h, i_tot, nu, F, ne, kl, c_ref = params

    C_bulk = x[0]

    FE = (ne * F * kl * C_bulk) / (i_tot + ne * F * kl * c_ref)

    # keep FE physical (optional but strongly recommended)
    FE = np.clip(FE, 0.0, 1.0)

    # FIX 1: negative sign (consumption along reactor)
    dCdy = -(ae * w * h) * (FE * i_tot) / (nu * F * ne)

    return [dCdy]

x0 = [c_ref]
y_span = (0, L)
y_eval = np.linspace(0, L, 1000)

sol = integrate.solve_ivp(equation_i, y_span, x0, t_eval=y_eval, args=(params,))

y = sol.t
c_bulk = sol.y[0]

conversion = 1 - c_bulk / c_ref

FE = (ne * F * kl * c_bulk) / (i_tot + ne * F * kl * c_ref)
FE = np.clip(FE, 0.0, 1.0)

# Plot
fig, ax = plt.subplots(2, 1, figsize=(7, 7), sharex=True)

ax[0].plot(y, conversion)
ax[0].set_ylabel("Conversion X(y)")
ax[0].grid(True)

ax[1].plot(y, FE)
ax[1].set_xlabel("Reactor length y [m]")
ax[1].set_ylabel("Faradaic efficiency FE(y)")
ax[1].set_ylim(0, 1.05)
ax[1].grid(True)

plt.tight_layout()
plt.show()