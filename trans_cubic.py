from aesthetic import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

def roots_trans_cubic(b0, T, k_vals, Sigma=1):
    Sigma_prime = 4 * Sigma / 3 
    B = b0**2 / (T**4 * (4/3) * (36 / pi**2)) # compute the rescaled magnetic field
    omega0 = 2 * b0 / (5 * T)

    roots_all = []

    for idx, k in enumerate(k_vals):
        # Define coefficients
        a3 = 1
        a2 = -1j * (Sigma + Sigma_prime)
        a1 = - ( Sigma*Sigma_prime + k**2 / 5 + omega0**2 / 4 )
        a0 = 1j * Sigma_prime * k**2 / 5

        coeffs = [a3, a2, a1, a0]
        roots = np.roots(coeffs)

        if idx == 0:
            prev_roots = roots
            roots_all.append(roots)
        else:
            # Match roots to maintain continuity
            cost_matrix = np.abs(roots[:, None] - prev_roots[None, :])
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            sorted_roots = roots[row_ind[np.argsort(col_ind)]]
            roots_all.append(sorted_roots)
            prev_roots = sorted_roots

    roots_all = np.array(roots_all)
    return np.real(roots_all), np.imag(roots_all)

# === USER INPUT ===
B_values = [0.2, 1.0]  # magnetic field in units of GeV^{-2}
temperature = 0.5 # temperature in units of GeV

# === Computation ===
k_vals = np.linspace(0.01, 5, 500)
colors = ['k', 'b']

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

for idx, B in enumerate(B_values):
    roots_real, roots_imag = roots_trans_cubic(B, temperature, k_vals)
    for j in range(3):  # 4 roots
        axes[0].plot(k_vals, roots_real[:, j], color=colors[idx], label='$B_0$ = '+str(B)+' GeV$^{2}$' if j == 0 else None)
        axes[1].plot(k_vals, roots_imag[:, j], color=colors[idx])

# Final plot formatting
axes[0].set_xlabel("k")
axes[0].set_ylabel("Re (ω)")
axes[0].set_xlim(0, k_vals[-1])
axes[0].legend()
axes[0].text(
    0.05, 0.05, f"T = {temperature} GeV",
    transform=axes[0].transAxes,
    fontsize = 25
)

axes[1].set_xlabel("k")
axes[1].set_ylabel("Im (ω)")
axes[1].set_xlim(0, k_vals[-1])
axes[1].set_ylim(bottom=0)

plt.tight_layout()
plt.savefig('trans-modes-cubic-T='+str(temperature)+'-GeV.pdf')
plt.show()