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
  executable_path = os.path.normpath(os.path.join(script_dir, "../cmake-build-release/inspect_PrivRepGame"))
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
resident = "L3"
p_mut_res = 0.5
mutant = f"RANDOM-{p_mut_res:.1f}"
pcs = run_simulation_return_pcs(resident, mutant)
pcs

# %%
# cell_L3
resident = "L3"
pc_res_mut_list = []
for mut_strategy in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
  mutant = f"RANDOM-{mut_strategy:.1f}"
  pcs = run_simulation_return_pcs(resident, mutant)
  pc_res_mut_list.append(pcs[1])
pc_res_mut_list

# %%
# cell_L6
resident = "L6"
pc_res_mut_list = []
for mut_strategy in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
  mutant = f"RANDOM-{mut_strategy:.1f}"
  pcs = run_simulation_return_pcs(resident, mutant)
  pc_res_mut_list.append(pcs[1])
pc_res_mut_list

# %%
# cell_Shunning
resident = "SECOND-8"
pc_res_mut_list = []
for mut_strategy in [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]:
  mutant = f"RANDOM-{mut_strategy:.1f}"
  pcs = run_simulation_return_pcs(resident, mutant)
  pc_res_mut_list.append(pcs[1])
pc_res_mut_list

# %%
# cell_L3_vs_SECOND
resident = "L3"
pc_res_mut_list = []
pc_mut_res_list = []
for i in range(0, 16):
  mutant = f"SECOND-{i}"
  pcs = run_simulation_return_pcs(resident, mutant)
  print(f"{mutant}: {pcs[2]} {pcs[1]}")
  pc_res_mut_list.append(pcs[1])
  pc_mut_res_list.append(pcs[2])
pc_res_mut_list, pc_mut_res_list

# %%
# cell_L6_vs_SECOND
resident = "L6"
pc_res_mut_list = []
pc_mut_res_list = []
for i in range(0, 16):
  mutant = f"SECOND-{i}"
  pcs = run_simulation_return_pcs(resident, mutant)
  print(f"{mutant}: {pcs[2]} {pcs[1]}")
  pc_res_mut_list.append(pcs[1])
  pc_mut_res_list.append(pcs[2])
pc_res_mut_list, pc_mut_res_list

# %%
# cell_Shunning_vs_SECOND
resident = "SECOND-8"
pc_res_mut_list = []
pc_mut_res_list = []
for i in range(0, 16):
  mutant = f"SECOND-{i}"
  pcs = run_simulation_return_pcs(resident, mutant)
  print(f"{mutant}: {pcs[2]} {pcs[1]}")
  pc_res_mut_list.append(pcs[1])
  pc_mut_res_list.append(pcs[2])
pc_res_mut_list, pc_mut_res_list


# %%
plt.clf()
fig, ax = plt.subplots()
cmap = plt.get_cmap("tab10")

# L3 :y = 0.84 * x + 0.14
# run cell_L3
x = np.linspace(0, 1, 100)
y = 0.84 * x + 0.14
ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.plot(x, y, linestyle='--', color=cmap(0))
xdata = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
ydata = [0.136821, 0.310654, 0.479187, 0.645467, 0.811017, 0.979835]  # run cell_L3
ax.scatter(xdata, ydata, color=cmap(0), marker='o', s=32)
# run cell_L3_vs_SECOND
xdata = [0.020044, 0.128124, 0.465074, 0.49918, 0.022814, 0.156265, 0.501498, 0.537638, 0.109109, 0.501784, 0.859255, 0.875604, 0.476257, 0.88573, 0.977447, 0.980277]
ydata = [0.154431, 0.247238, 0.533671, 0.562056, 0.158271, 0.273083, 0.562598, 0.59323, 0.230728, 0.563644, 0.861888, 0.875885, 0.541252, 0.883423, 0.960659, 0.963467]
ax.scatter(xdata, ydata, color=cmap(0), marker='v', s=64)
ax.scatter([0.876], [0.876], color=cmap(0), marker='s', s=96)

# L6: y = 0.5
# run cell_L6
x = np.linspace(0, 1, 100)
y = 0.5 * np.ones_like(x)
ax.plot(x, y, linestyle='--', color=cmap(1))
xdata = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
# ydata = [0.498006, 0.496782, 0.498529, 0.501644, 0.501116, 0.504194] for N=200
ydata = [0.49059, 0.494134, 0.497007, 0.502589, 0.505425, 0.510967]
ax.scatter(xdata, ydata, color=cmap(1), marker='o', s=32)
# run cell_L6_vs_SECOND
ydata = [0.489249, 0.49587, 0.498162, 0.501615, 0.489865, 0.500312, 0.500906, 0.502054, 0.492199, 0.499734, 0.499935, 0.50302, 0.49877, 0.509358, 0.509581, 0.507662]
xdata = [0.019654, 0.336337, 0.33962, 0.499461, 0.03742, 0.499029, 0.498658, 0.664357, 0.038616, 0.499591, 0.499769, 0.659579, 0.511966, 0.962772, 0.961977, 0.979967]
ax.scatter(xdata, ydata, color=cmap(1), marker='v', s=64)
ax.scatter([0.5], [0.5], color=cmap(1), marker='s', s=96)

# Shunning: y = 0.0192*x + 0.02
x = np.linspace(0, 1, 100)
y = 0.0192 * x + 0.02
ax.plot(x, y, linestyle='--', color=cmap(2))
xdata = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
ydata = [0.020119, 0.023514, 0.027885, 0.031906, 0.035189, 0.039831]
ax.scatter(xdata, ydata, color=cmap(2), marker='o', s=32)
# run cell_Shunning_vs_SECOND
ydata = [0.020146, 0.029271, 0.020389, 0.029587, 0.024411, 0.038709, 0.027901, 0.039108, 0.020244, 0.029751, 0.020388, 0.029988, 0.028992, 0.038743, 0.032824, 0.039296]
xdata = [0.020101, 0.49438, 0.038774, 0.499703, 0.248908, 0.960527, 0.409751, 0.961194, 0.020701, 0.494665, 0.039314, 0.49981, 0.48695, 0.979816, 0.660204, 0.980014]
ax.scatter(xdata, ydata, color=cmap(2), marker='v', s=64)
ax.scatter([0.02], [0.02], color=cmap(2), marker='s', s=96)

ax.plot([], [], linestyle='--', color='gray', label='theory')
ax.scatter([], [], color='gray', marker='o', s=32, label='random cooperator')
ax.scatter([], [], color='gray', marker='v', s=64, label='second-order norm')
ax.scatter([], [], color='gray', marker='s', s=96, label='resident')

ax.text(0.75, 0.95, 'Simple Standing', fontsize=14, ha='center', va='center', color=cmap(0))
ax.text(0.8, 0.55, 'Stern Judging', fontsize=14, ha='center', va='center', color=cmap(1))
ax.text(0.8, 0.08, 'Shunning', fontsize=14, ha='center', va='center', color=cmap(2))

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)
ax.set_xlabel(r'$p_{\rm mut \to res}$', fontsize=18)
ax.set_ylabel(r'$p_{\rm res \to mut}$', fontsize=18)
ax.legend(fontsize=13, frameon=False)


# %%
fig.savefig('pc_mut_res.pdf', bbox_inches='tight')
# %%
