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
  executable_path = os.path.normpath(os.path.join(script_dir, "../cmake-build-release/inspect_PrivRepGame3rd"))
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
resident = "L5"
mutant = "L1"
pcs = run_simulation_return_pcs(resident, mutant)
print(pcs)

resident = "L5"
mutant = "RANDOM-0.695"
pcs = run_simulation_return_pcs(resident, mutant)
print(pcs)

# %%
# GKGB
# resident = '0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 0.0 1.0 0.0 1.0'
resident = '0.0 1.0 0.0 1.0 0.0 1.0 1.0 1.0 0.0 1.0 0.0 1.0'
mutant = 'AllC'
pcs = run_simulation_return_pcs(resident, mutant)
print(pcs)

# %%
# cell L5
resident = "L5"
pc_res_mut_list = []
pc_mut_res_list = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
pc_mut_res_list = [i/20 for i in range(21)]
pc_self = 0.0
for mut_p in pc_mut_res_list:
  mutant = f"RANDOM-{mut_p:.2f}"
  pcs = run_simulation_return_pcs(resident, mutant)
  pc_res_mut_list.append(pcs[1])
  pc_self = pcs[0]
pc_res_mut_list
print(pc_res_mut_list, pc_self)

plt.clf()
fig, ax = plt.subplots()
cmap = plt.get_cmap("tab10")

ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
ax.scatter(pc_mut_res_list, pc_res_mut_list, color=cmap(0), marker='o', s=32)
ax.scatter([pc_self], [pc_self], color=cmap(0), marker='s', s=64)

# %%
plt.clf()
fig, ax = plt.subplots()
cmap = plt.get_cmap("tab10")

ax.set_xlim(-0.02, 1.02)
ax.set_ylim(-0.02, 1.02)
# x = np.linspace(0, 1, 100)
# y = 0.39 * x + 0.23
# ax.plot(x, y, linestyle='--', color=cmap(0))
xdata = pc_mut_res_list
ydata = [0.387162, 0.411706, 0.432194, 0.452167, 0.471569, 0.491161, 0.510533, 0.528716, 0.545746, 0.560382, 0.57867, 0.594616, 0.611571, 0.623595, 0.638308, 0.651096, 0.665018, 0.677757, 0.69331, 0.704311, 0.716838] # L5
ax.scatter(xdata, ydata, color=cmap(0), marker='o', s=32)
ax.scatter([0.613262], [0.613262], color=cmap(0), marker='s', s=64)

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
fig.savefig('pc_mut_res_solitary_L5.pdf', bbox_inches='tight')
# %%
