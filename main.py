import numpy as np
import matplotlib.pyplot as plt

from interpolation import computePL, compute_pN, compute_pH


def evaluer_polynome(coeffs, x):
    return sum(c * (x ** i) for i, c in enumerate(coeffs))

def f1(x): return np.exp(x)
def f1_d(x): return np.exp(x)

def f2(x): return 1.0 / (1.0 + x**2)
def f2_d(x): return -2.0 * x / ((1.0 + x**2)**2)

x_fine = np.linspace(-4.0, 4.0, 500)
N_values = [3, 5, 9, 13]
couleurs = ['red', 'orange', 'green', 'blue']


fig, axes = plt.subplots(2, 3, figsize=(18, 10))

# 1. Traitement de f1(x) = exp(x) 
for N, color in zip(N_values, couleurs):
    xn = np.linspace(-4.0, 4.0, N)
    pts = [(x, f1(x)) for x in xn]
    pts_h = [(x, f1(x), f1_d(x)) for x in xn]
    
    # Lagrange
    y_lag = [evaluer_polynome(computePL(pts), x) for x in x_fine]
    axes[0, 0].plot(x_fine, y_lag, color=color, label=f"N={N}")
    
    # Newton
    y_newt = [evaluer_polynome(compute_pN(pts), x) for x in x_fine]
    axes[0, 1].plot(x_fine, y_newt, color=color, label=f"N={N}")
    
    # Hermite
    y_herm = [evaluer_polynome(compute_pH(pts_h), x) for x in x_fine]
    axes[0, 2].plot(x_fine, y_herm, color=color, label=f"N={N}")

# 2. Traitement de f2(x) = 1/(1+x^2) 
for N, color in zip(N_values, couleurs):
    xn = np.linspace(-4.0, 4.0, N)
    pts = [(x, f2(x)) for x in xn]
    pts_h = [(x, f2(x), f2_d(x)) for x in xn]
    
    # Lagrange
    y_lag = [evaluer_polynome(computePL(pts), x) for x in x_fine]
    axes[1, 0].plot(x_fine, y_lag, color=color, label=f"N={N}")
    
    # Newton
    y_newt = [evaluer_polynome(compute_pN(pts), x) for x in x_fine]
    axes[1, 1].plot(x_fine, y_newt, color=color, label=f"N={N}")
    
    # Hermite
    y_herm = [evaluer_polynome(compute_pH(pts_h), x) for x in x_fine]
    axes[1, 2].plot(x_fine, y_herm, color=color, label=f"N={N}")


titres = [
    ["f1(x) - Lagrange", "f1(x) - Newton", "f1(x) - Hermite"],
    ["f2(x) - Lagrange", "f2(x) - Newton", "f2(x) - Hermite"]
]

for i in range(2):
    f_orig = [f1(x) for x in x_fine] if i == 0 else [f2(x) for x in x_fine]
    ylim = (-5, 30) if i == 0 else (-2, 3) 
    
    for j in range(3):
        ax = axes[i, j]
        ax.plot(x_fine, f_orig, 'k--', linewidth=2, label="Originale")
        ax.set_title(titres[i][j], fontsize=12, fontweight='bold')
        ax.set_xlim(-4, 4)
        ax.set_ylim(ylim)
        ax.grid(True, linestyle=":", alpha=0.7)
        ax.legend(fontsize=8)

plt.tight_layout()
plt.show()
