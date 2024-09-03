# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
import os
import re
import json

# %%
# get the path of this script
def run_simulation_return_pcs(resident, mutant):
  script_dir = os.path.dirname(os.path.abspath(__file__))
  executable_path = os.path.normpath(os.path.join(script_dir, "../cmake-build-release/inspect_PrivRepGameTernary"))
  t_measure = 1000000
  param = {"mu_assess": 0.02, "mu_impl": 0.0, "q": 0.0, "seed": 123456789, "t_init": t_measure/10, "t_measure": t_measure}
  print(json.dumps(param))
  command = [executable_path, "-j", json.dumps(param), resident, "99", mutant, "1"]
  print(command)
  result = subprocess.run(command, capture_output=True, text=True)
  print(result.stdout)
  print(result.stderr)
  result
  lines = result.stdout.splitlines()
  start_index = next((i for i, line in enumerate(lines) if "NormCooperationLevels:" in line), None)
  a = [float(num) for num in lines[start_index+1].split()]
  p_res_res = a[0]
  p_res_mut = a[1]
  p_mut_res = float(lines[start_index+2].split()[0])
  return (p_res_res, p_res_mut, p_mut_res)

# %%
resident = "GBGGBN 1 0 0"
p_mut_res = 0.5
mutant = f"BBBBBB {p_mut_res:.1f} {p_mut_res:.1f} {p_mut_res:.1f}"
pcs = run_simulation_return_pcs(resident, mutant)
pcs

# %%
# cell GBGGBN 1 0 0
resident = "GBGGBN 1 0 0"
pc_res_mut_list = []
for mut_p in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
  mutant = f"BBBBBB {mut_p:.1f} {mut_p:.1f} {mut_p:.1f}"
  pcs = run_simulation_return_pcs(resident, mutant)
  pc_res_mut_list.append(pcs[1])
pc_res_mut_list

# %%
# cell GBGGNG 1 0 0
resident = "GBGGNG 1 0 0"
pc_res_mut_list = []
for mut_p in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
  mutant = f"BBBBBB {mut_p:.1f} {mut_p:.1f} {mut_p:.1f}"
  pcs = run_simulation_return_pcs(resident, mutant)
  pc_res_mut_list.append(pcs[1])
pc_res_mut_list

# %%
plt.clf()
fig, ax = plt.subplots()
cmap = plt.get_cmap("tab10")

x = np.linspace(0, 1, 100)
# GBGGBN 1 0 0
# h* = 0.0100211,0.379447,0.610532
y = 0.39 * x + 0.23
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.plot(x, y, linestyle='--', color=cmap(0))
xdata = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
ydata = [0.238591, 0.313516, 0.391558, 0.464644, 0.540442, 0.614334]
ax.scatter(xdata, ydata, color=cmap(0), marker='o', s=32)
ax.scatter([0.383], [0.383], color=cmap(0), marker='s', s=64)

# GBGGNG 1 0 0
y = 0.36 * x + 0.39
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.plot(x, y, linestyle='--', color=cmap(1))
ydata = [0.38922, 0.461884, 0.532643, 0.60452, 0.674638, 0.74838]
ax.scatter(xdata, ydata, color=cmap(1), marker='o', s=32)
ax.scatter([0.606], [0.606], color=cmap(1), marker='s', s=64)

ax.text(0.75, 0.41, 'GBGGBN-100', fontsize=14, ha='center', va='center', color=cmap(0))
ax.text(0.7, 0.73, 'GBGGNG-100', fontsize=14, ha='center', va='center', color=cmap(1))

ax.plot([], [], linestyle='--', color='gray', label='theory')
ax.scatter([], [], color='gray', marker='o', s=32, label='random cooperator')
ax.scatter([], [], color='gray', marker='s', s=64, label='resident')

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel(r'$p_{\rm mut \to res}$', fontsize=18)
ax.set_ylabel(r'$p_{\rm res \to mut}$', fontsize=18)
ax.legend(fontsize=13, frameon=False)


# %%
fig.savefig('pc_mut_res_solitary_ternary.pdf', bbox_inches='tight')
# %%
