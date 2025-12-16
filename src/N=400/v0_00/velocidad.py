import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt

# -------------------------
# CONFIG
# -------------------------
hx = 0.10
dt = 5.0e-6
N  = 400

block = 20      # tamaño de ventana para polyfit
step  = 1      # salto entre ventanas (20 => bloques no solapados; 1 => deslizante)

base = Path(__file__).resolve().parent
fname = base / "xcm.dat"
print(fname)
if not fname.exists():
    cand = list(base.rglob("xcm.dat"))
    if cand:
        fname = cand[0]
    else:
        raise FileNotFoundError(
            f"No encuentro xcm_x.dat.\nBuscando en: {base}\n"
            f"Tip: ejecuta `find . -name xcm_x.dat`"
        )

print(f"Leyendo: {fname}")

# -------------------------
# LECTURA
# columnas: t, x_cm, theta_cm, umax, flag
# -------------------------
data  = np.loadtxt(fname)
t     = data[:, 0]
theta = data[:, 2]
flag  = data[:, 4].astype(int)

mask = (flag == 1) & (theta >= 0.0)
t = t[mask]
theta = theta[mask]

if t.size < block + 2:
    raise RuntimeError(f"Muy pocos puntos válidos (flag==1). Necesitas al menos {block+2}.")

# -------------------------
# UNWRAP + x continuo
# -------------------------
theta_u = np.unwrap(theta)
Lx = N * hx
x_u = (theta_u / (2.0*np.pi)) * Lx

# Velocidad instantánea (derivada)
v_inst = np.gradient(x_u, t)
v_inst_star = v_inst * (dt / hx)

# -------------------------
# POLYFIT POR BLOQUES
# -------------------------
t_mid = []
v_fit = []
for k in range(0, len(t) - block + 1, step):
    ts = t[k:k+block]
    xs = x_u[k:k+block]

    # Ajuste lineal en la ventana
    v_phys, b = np.polyfit(ts, xs, 1)

    t_mid.append(0.5 * (ts[0] + ts[-1]))
    v_fit.append(v_phys)

t_mid = np.array(t_mid)
v_fit = np.array(v_fit)
v_fit_star = v_fit * (dt / hx)
print(v_fit)

# -------------------------
# GRÁFICAS
# -------------------------
# Posición
plt.figure()
plt.plot(t, x_u, marker='o', linewidth=1)
plt.xlabel("t")
plt.ylabel("x_cm (unwrapped)")
plt.title("Posición del centro (unwrapped)")
plt.grid(True)

# Velocidad por bloques (polyfit)
plt.figure()
plt.plot(t_mid, v_fit, marker='o', linewidth=1)
plt.axhline(17.4, linestyle='--', linewidth=1)
plt.xlabel("t (centro de ventana)")
plt.ylabel("v* (polyfit por bloque)")
plt.title(f"Velocidad adimensional por bloques | block={block}, step={step}")
plt.grid(True)

plt.show()
