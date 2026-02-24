import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# --- Constants ---
w     = 1e-3
h     = 1e-3
i_tot = 2000
Vdot  = 75e-6/60
F     = 96485
ne    = 2
L     = 5.0
c_ref = 38.0
c0    = c_ref

kL0 = 8e-5
ae0 = 1/h

# Sensitivity multipliers
mults = [0.25, 0.5, 1, 2, 5, 10]

def FE_from_cbulk(cbulk, kL):
    denom = ne * F * kL * c_ref + i_tot
    FE = (ne * F * kL * cbulk) / denom
    return np.clip(FE, 0.0, 1.0)

def solve_profile(kL, ae):
    def ode(y, cbulk):
        c = cbulk[0]
        FE = FE_from_cbulk(c, kL)
        dcdy = -(ae * w * h) * (FE * i_tot) / (Vdot * F * ne)
        return [dcdy]

    y_eval = np.linspace(0, L, 400)
    sol = solve_ivp(ode, (0, L), [c0], t_eval=y_eval)

    y = sol.t
    c_bulk = sol.y[0]
    X = 1 - c_bulk / c0
    FE = FE_from_cbulk(c_bulk, kL)

    return y, X, FE


# ===============================
# FIGURE 1: Conversion vs y (kL variation)
# ===============================
plt.figure(figsize=(7,5))
for m in mults:
    y, X, FE = solve_profile(kL0*m, ae0)
    plt.plot(y, X, label=f'{m}×')

plt.xlabel("Reactor length y [m]")
plt.ylabel("Conversion X(y)")
plt.title("Conversion sensitivity to $k_L$")
plt.legend(title="Multiplier")
plt.grid(True)
plt.tight_layout()
plt.show()


# ===============================
# FIGURE 2: FE vs y (kL variation)
# ===============================
plt.figure(figsize=(7,5))
for m in mults:
    y, X, FE = solve_profile(kL0*m, ae0)
    plt.plot(y, FE, label=f'{m}×')

plt.xlabel("Reactor length y [m]")
plt.ylabel("Faradaic efficiency FE(y)")
plt.title("Faradaic efficiency sensitivity to $k_L$")
plt.legend(title="Multiplier")
plt.grid(True)
plt.tight_layout()
plt.show()


# ===============================
# FIGURE 3: Conversion vs y (ae variation)
# ===============================
plt.figure(figsize=(7,5))
for m in mults:
    y, X, FE = solve_profile(kL0, ae0*m)
    plt.plot(y, X, label=f'{m}×')

plt.xlabel("Reactor length y [m]")
plt.ylabel("Conversion X(y)")
plt.title("Conversion sensitivity to $a_e$")
plt.legend(title="Multiplier")
plt.grid(True)
plt.tight_layout()
plt.show()


# ===============================
# FIGURE 4: FE vs y (ae variation)
# ===============================
plt.figure(figsize=(7,5))
for m in mults:
    y, X, FE = solve_profile(kL0, ae0*m)
    plt.plot(y, FE, label=f'{m}×')

plt.xlabel("Reactor length y [m]")
plt.ylabel("Faradaic efficiency FE(y)")
plt.title("Faradaic efficiency sensitivity to $a_e$")
plt.legend(title="Multiplier")
plt.grid(True)
plt.tight_layout()
plt.show()