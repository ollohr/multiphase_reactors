import numpy as np
import matplotlib.pyplot as plt

## Question D

#constants
Ks = 0.3 #g/L
u_max = 0.635/3600 #/s
Y = 0.51 #g bio/g glucose
w = 88E-12 #g / cell
k = 1E-14 #g/s

#initial cond
x_0 = 2 #cells
cs_0 = 1.5 #g/L
v_0 = 2E-9 #L

#fed cond
v_n = 0.8E-9 #L
cs_n = 4 #g/L

#code cond
duur = 24*3600 #s
dt = 1 #s


def ODE(var): #var : X, Cs,V
    dx_dt = u_max*(var[1]/(var[1]+Ks))*var[0]
    dCs_dt = 1/var[2]*(-dx_dt/Y-k*var[0]/w)
    return dx_dt,dCs_dt

x = [x_0*w]
cs = [cs_0]
v = [v_0]

r_growth = [0]
r_product = [0]

count = 0

for i in range(int(duur/dt)):
    # x = x[i]
    # cs = cs[i]
    var = [x[i],cs[i], v[i]]
    dvar = ODE(var)

    r_growth.append(dvar[0] / Y)      # g glucose / s used for biomass growth
    r_product.append(k * x[i] / w)  # g glucose / s used for product formation

    x.append(x[i] + dt*dvar[0]) #assuming feed has no cells only nutrients

    cs_b = cs[i] + dt*dvar[1]

    if cs_b < 0.01:
        v.append(v[i]+v_n)
        cs.append((cs_b*v[i]+cs_n*v_n)/v[i+1])

    else:
        cs.append(cs_b)
        v.append(v[i])

    # if cs[i+1] < 0.0000001: #cs scale does not reach less than
    #     print(i*dt)
    # checking percent <<ks
    if cs[i+1] < Ks:
        count +=1

print(count/int(duur/dt))

#Question e

t = np.arange(0, duur + dt, dt)   # seconds
t_h = t / 3600                    # hours

# Plot 1: X(t)
plt.figure(figsize=(8, 5))
plt.plot(t_h, x)
plt.xlabel('Time [h]')
plt.ylabel('Cell mass, X [g]')
plt.title('Evolution of cell mass X(t)')
plt.xlim(0, 24)
plt.grid(True)
plt.tight_layout()
plt.show()

# Plot 2: Cs(t) with Ks
plt.figure(figsize=(8, 5))
plt.plot(t_h, cs, label=r'$C_s(t)$')
plt.axhline(y=Ks, linestyle='--', label=r'$K_s$')
plt.axhline(y=0.01, linestyle=':', label=r'Feed threshold = 0.01 g/L')
plt.xlabel('Time [h]')
plt.ylabel(r'Glucose concentration, $C_s$ [g/L]')
plt.title(r'Evolution of glucose concentration $C_s(t)$')
plt.xlim(0, 24)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# Plot 3: glucose consumption rates
plt.figure(figsize=(8, 5))
plt.plot(t_h, r_growth, label='Growth-related glucose consumption')
plt.plot(t_h, r_product, label='Product-related glucose consumption')
plt.xlabel('Time [h]')
plt.ylabel('Rate [g/s]')
plt.title('Evolution of glucose consumption rates')
plt.xlim(0, 24)
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()