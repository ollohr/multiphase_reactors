import numpy as np
import matplotlib.pyplot as plt

# ============================
# QUESTION A – ΔP vs L fitting
# ============================

# Data estimated from the centers of the square markers
L_data  = np.array([0.1, 0.4, 0.5, 0.6, 0.8, 1.0])          # [m]
dP_data = np.array([1200, 5400, 6800, 7500, 9400, 14200])   # [Pa]

# Linear + quadratic fits
coeff_lin  = np.polyfit(L_data, dP_data, 1)
coeff_quad = np.polyfit(L_data, dP_data, 2)

L_plot = np.linspace(0.0, 1.05, 300)
dP_lin  = np.polyval(coeff_lin,  L_plot)
dP_quad = np.polyval(coeff_quad, L_plot)

def r2(y, yhat):
    ss_res = np.sum((y - yhat)**2)
    ss_tot = np.sum((y - np.mean(y))**2)
    return 1 - ss_res/ss_tot

r2_lin  = r2(dP_data, np.polyval(coeff_lin,  L_data))
r2_quad = r2(dP_data, np.polyval(coeff_quad, L_data))

print(f"Linear fit:    ΔP = {coeff_lin[0]:.4e} L + {coeff_lin[1]:.4e}   (R²={r2_lin:.3f})")
print(f"Quadratic fit: ΔP = {coeff_quad[0]:.4e} L² + {coeff_quad[1]:.4e} L + {coeff_quad[2]:.4e}   (R²={r2_quad:.3f})")

# Plot
plt.figure(figsize=(6.2, 4.4))
plt.plot(L_data, dP_data, 's', label='Estimated data (marker centers)')
plt.plot(L_plot, dP_quad, label=f'Quadratic fit (R²={r2_quad:.3f})')
plt.plot(L_plot, dP_lin, '--', label=f'Linear fit (R²={r2_lin:.3f})')
plt.xlabel('Column length (m)')
plt.ylabel('ΔP (Pa)')
plt.xlim(0, 1.1)
plt.ylim(0, 16000)
plt.grid(True, which='both', alpha=0.4)
plt.legend()
plt.tight_layout()
plt.show()