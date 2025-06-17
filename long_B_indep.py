from aesthetic import *
from math import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import linear_sum_assignment

def roots_long_B_indep(k_vals, Sigma=1):

    # Initialize
    roots_all = []

    for idx, k in enumerate(k_vals):
        # Define coefficients
        a3 = 1
        a2 = -1j * Sigma 
        a1 = - 3 * k**2 / 5
        a0 = 1j * Sigma * k**2 / 3
        coeffs = [a3, a2, a1, a0]
        roots = np.roots(coeffs)
        
        # Sort roots to maintain continuity
        if idx == 0:
            # Store first set as-is
            prev_roots = roots
            roots_all.append(roots)
        else:
            # Cost matrix based on distance to previous roots
            cost_matrix = np.abs(roots[:, None] - prev_roots[None, :])
            row_ind, col_ind = linear_sum_assignment(cost_matrix)
            # Rearranged roots
            sorted_roots = roots[row_ind[np.argsort(col_ind)]]
            roots_all.append(sorted_roots)
            prev_roots = sorted_roots

    roots_all = np.array(roots_all)
    return np.real(roots_all), np.imag(roots_all)

if __name__ == '__main__':
    
    # Parameters
    Sigma = 1
    Sigma_prime = 4 * Sigma / 3
    k_vals = np.linspace(0.01, 5, 500)

    # Plotting
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharex=True)

    roots_real, roots_imag = roots_long_B_indep(k_vals)
    for j in range(3):  # 4 roots
        axes[0].plot(k_vals, roots_real[:, j], color='k')
        axes[1].plot(k_vals, roots_imag[:, j], color='k')

    # Final plot formatting
    axes[0].set_xlabel("k")
    axes[0].set_ylabel("Re (ω)")
    axes[0].set_xlim(0, k_vals[-1])

    axes[1].set_xlabel("k")
    axes[1].set_ylabel("Im (ω)")
    axes[1].set_xlim(0, k_vals[-1])
    axes[1].set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig('long-modes-B-indep.pdf')
    plt.show()