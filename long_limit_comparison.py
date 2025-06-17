from long_B_dep import *
from long_B_indep import *
from aesthetic import *
from math import *
import numpy as np
import matplotlib.pyplot as plt

# Plot options
colors = ['indigo', 'blue', 'darkorange', 'red']
fig, ax = plt.subplots(figsize=(8, 6), sharex=True)

# Parameters
Sigma = 1
Sigma_prime = 4 * Sigma / 3
k_vals = np.linspace(0.01, 5, 500)


# Modes independent of the magnetic field
B_indep_re, B_indep_im = roots_long_B_indep(k_vals, Sigma=1)

# Modes dependent on the magnetic field

B_values = [0.2, 0.5, 1.0, 2.0]  # magnetic field in units of GeV^{-2}
temperature = 0.5 # temperature in units of GeV

# Plots
ax.plot(k_vals, B_indep_im[:, 1], color='k', label='Longitudinal limit')

for idx, b0 in enumerate(B_values):
    B_dep_re, B_dep_im = roots_long_B_dep(b0, temperature, k_vals)
    ax.plot(k_vals, B_dep_im[:, 2], color=colors[idx], label='$B_0$ = '+str(b0)+' GeV$^{2}$')

ax.set_xlabel("k")
ax.set_ylabel("Im (ω)")
ax.set_xlim(0, k_vals[-1])
ax.set_ylim(bottom=0)
#ax.set_ylim(0, .8)
ax.legend(loc=2, fontsize=15)

plt.tight_layout()
plt.savefig('im-long-hydro-modes-comparison-T='+str(temperature)+'-GeV.pdf')
plt.show()