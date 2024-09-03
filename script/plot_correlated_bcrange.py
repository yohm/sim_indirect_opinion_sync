# %%
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
cmap = plt.get_cmap("tab10")

# solve
# h = -2h (1-2mu_a) (1 - hG) + (1 - mu_a)
# h( 1 + 2(1-2mu_a)(1-hG) )= 1 - mu_a
# h = (1 - mu_a) / (1 + 2(1-2mu_a)(1-hG) )
hG = np.linspace(0.5, 1, 100)
mu_a = 0.02
h = (1 - mu_a) / (1 + 2 * (1 - 2 * mu_a) * (1 - hG))
ax.set_xlim(0.49, 1.01)
ax.set_ylim(0.49, 1.01)
ax.plot(hG, h, color=cmap(0))
ax.plot(hG, hG, linestyle='--', color='gray', alpha=0.3)
ax.set_xlabel(r'$h_G$', fontsize=18)
ax.set_ylabel(r'$h$', fontsize=18)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.scatter([0.5], [0.5], color=cmap(1), marker='o', s=84, label='solitary')
ax.scatter([1.0], [0.98], color=cmap(2), marker='s', s=64, label='public')
ax.scatter([0.499971], [0.49999], color=cmap(3), marker='^', s=54, label='simultaneous')   # calculated from simulation
# gossiping model: h = 1 - (1-hG)e^{\tau}
# Define the equations
tau_h = []
tau_hG = []
for tau in [0.1, 0.3, 1.0, 3.0]:
  def equations(vars, tau):
    h, hG = vars
    mu_a = 0.02
    eq1 = h - (1 - mu_a) / (1 + 2 * (1 - 2 * mu_a) * (1 - hG))
    eq2 = h - (1 - (1 - hG) * np.exp(tau))
    return [eq1, eq2]
  solution = fsolve(equations, [0.5, 0.5], args=(tau))
  h, hG = solution
  tau_h.append(h)
  tau_hG.append(hG)
print(tau_h, tau_hG)
ax.scatter(tau_hG, tau_h, color=cmap(4), marker='x', s=64, label='gossiping')
ax.legend(fontsize=14)

# %%
fig.savefig('correlated_L6_h_hG.pdf', bbox_inches='tight')

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
cmap = plt.get_cmap("tab10")
# solve
# h = (1 - mu_a) / [ 2(1-mu_a) - (1-2mu_a)hG ]
hG = np.linspace(0.875, 1, 100)
mu_a = 0.02
h = (1 - mu_a) / (2 * (1 - mu_a) - (1 - 2 * mu_a) * hG)
ax.set_xlim(0.873, 1.002)
ax.set_ylim(0.873, 1.002)
ax.set_xticks([0.9, 0.95, 1])
ax.set_xticklabels(['0.90', '0.95', '1.0'])
ax.set_yticks([0.9, 0.95, 1])
ax.set_yticklabels(['0.90', '0.95', '1.0'])
ax.plot(hG, h, color=cmap(0))
ax.plot(hG, hG, linestyle='-.', color='gray', alpha=0.3)
ax.set_xlabel(r'$h_G$', fontsize=18)
ax.set_ylabel(r'$h$', fontsize=18)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.scatter([0.875], [0.875], color=cmap(1), marker='o', s=84, label='solitary')
ax.scatter([1.0], [0.98], color=cmap(2), marker='s', s=64, label='public')
ax.scatter([0.978208], [0.95972], color=cmap(3), marker='^', s=64, label='simultaneous')   # calculated from simulation
# gossiping model: h = 1 - (1-hG)e^{\tau}
# Define the equations
tau_h = []
tau_hG = []
for tau in [0.1, 0.3, 1.0, 3.0]:
  def equations(vars, tau):
    h, hG = vars
    mu_a = 0.02
    eq1 = h - (1 - mu_a) / (2 * (1 - mu_a) - (1 - 2 * mu_a) * hG)
    eq2 = h - (1 - (1 - hG) * np.exp(tau))
    return [eq1, eq2]
  solution = fsolve(equations, [0.95, 0.95], args=(tau))
  h, hG = solution
  tau_h.append(h)
  tau_hG.append(hG)
print(tau_h, tau_hG)
ax.scatter(tau_hG, tau_h, color=cmap(4), marker='x', s=64, label='gossiping')
ax.legend(fontsize=14)

# %%
fig.savefig('correlated_L3_h_hG.pdf', bbox_inches='tight')

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
mu_a = 0.02
hG = np.linspace(0.5, 1.00, 1000)
b_lower = 1.0 / (2 * hG - 1) / (1 - 2 * mu_a)
ax.set_xlim(0.49, 1.01)
ax.set_ylim(1, 5)
ax.set_yticks([1, 2, 3, 4, 5])
ax.set_yticklabels(['1', '2', '3', '4', '5'])
ax.plot(hG, b_lower, color=cmap(0), label=r'$\mu_a = 0.02$')
ax.fill_between(hG, b_lower, 5, color=cmap(0), alpha=0.5)

ax.set_xlabel(r'$h_G$', fontsize=18)
ax.set_ylabel(r'$b/c$', fontsize=18)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

# %%
fig.savefig('correlated_L6_bc_range.pdf', bbox_inches='tight')

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
mu_a = 0.02
hg_min = math.sqrt(1.0-mu_a)* ( math.sqrt(1.0-mu_a) - math.sqrt(mu_a) ) / (1.0 - 2.0 * mu_a)  # 7/8
hG = np.linspace(hg_min, 1.00, 1000)
# alpha_G > 0
# b/c > 1/(1-2mu_a)hG
bc_lower = 1.0 / ((1.0 - 2.0 * mu_a) * hG)
# alpha_B < 0
# b/c < [(1-mu_a) - (1-2mu_a)hG] / (1-2mu_a)(1-mu_a)(1-hG)
bc_upper = ((1.0 - mu_a) - (1.0 - 2.0 * mu_a) * hG) / (1.0 - 2.0 * mu_a) / (1.0 - mu_a) / (1.0 - hG)
ax.set_xlim(0.873, 1.002)
ax.set_ylim(1, 5)
ax.plot(hG, bc_lower, color=cmap(0), label='lower bound')
ax.plot(hG, bc_upper, color=cmap(0), label='upper bound')
ax.fill_between(hG, bc_lower, bc_upper, color=cmap(0), alpha=0.5)

# set xticks [0.875, 0.900, 0.925, 0.95, 0.975, 1]
ax.set_xticks([0.9, 0.95, 1])
ax.set_xticklabels(['0.90', '0.95', '1.0'])
#ax.set_xticks([0.875, 0.900, 0.925, 0.95, 0.975, 1])
#ax.set_xticklabels(['0.875', '0.900', '0.925', '0.950', '0.975', '1'])
ax.set_yticks([1, 2, 3, 4, 5])
ax.set_yticklabels(['1', '2', '3', '4', '5'])

ax.set_xlabel(r'$h_G$', fontsize=18)
ax.set_ylabel(r'$b/c$', fontsize=18)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

# %%
fig.savefig('correlated_L3_bc_range.pdf', bbox_inches='tight')

