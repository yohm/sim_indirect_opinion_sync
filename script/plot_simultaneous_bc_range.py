# %%
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import subprocess
import os
import re
import json
import pickle

# %%
# get the path of this script
def run_simulation_return_pcs(r_bc):
  script_dir = os.path.dirname(os.path.abspath(__file__))
  executable_path = os.path.normpath(os.path.join(script_dir, "../cmake-build-release/main_CorrelatedOpinions"))
  t_measure = 10000
  param = {"N":100, "mu_assess": 0.02, "mu_impl": 0.0, "q": 1.0, "t_init": t_measure/10, "t_measure": t_measure, "seed": 12345}
  print(json.dumps(param))
  num_sample = 100
  command = [executable_path, "-j", json.dumps(param), f"1.0 0.0 {r_bc} 1.0 1.0 0.0", f"{num_sample}"]
  print(command)
  out = subprocess.run(command, capture_output=True, text=True)
  print(out.stdout)
  print(out.stderr)
  results= json.loads(out.stdout)
  return results

# %%
results = run_simulation_return_pcs(0.0)
results

# %%
if os.path.exists("results_all_rbc_dep.pkl"):
  with open("results_all_rbc_dep.pkl", "rb") as f:
    results_all_rbc_dep = pickle.load(f)
else:
  results_all_rbc_dep = []
  for i in range(0, 51):
    r_bc = i / 50.0
    r = run_simulation_return_pcs(r_bc)
    results_all_rbc_dep.append(r)
  with open("results_all_rbc_dep.pkl", "wb") as f:
    pickle.dump(results_all_rbc_dep, f)

results_all_rbc_dep

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
n = len(results_all_rbc_dep)
x = np.array([i/(n-1) for i in range(0, n)], dtype=float)
h = np.array([r["h"] for r in results_all_rbc_dep], dtype=float)
hg = np.array([r["hG"] for r in results_all_rbc_dep], dtype=float)

ax.set_xlim(-0.01, 1.01)
ax.set_ylim(0.48, 1.02)
ax.plot(x, h, marker="o", label=r"$h$")
ax.plot(x, hg, marker="o", label=r"$h_G$")
ax.set_xlabel(r"$R(B,C)$", fontsize=18)
ax.set_ylabel(r"$h, h_G$", fontsize=18)
ax.legend(fontsize=14)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

# %%
fig.savefig("interpolate_h_rbc.pdf", bbox_inches="tight")

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
cmap = mpl.cm.get_cmap("tab10")
# n = len(results_all)
# x = np.array([i/(n-1) for i in range(0, n)], dtype=float)
ax.set_xlim(-0.01, 1.01)
ax.set_ylim(1, 5)
ax.set_xlabel(r"$R(B,C)$", fontsize=18)
ax.set_ylabel(r"$b/c$", fontsize=18)
ax.set_yticks([1, 2, 3, 4, 5])

bc_l_sim = np.array([r["bc_lower_alld"] for r in results_all_rbc_dep], dtype=float)
bc_u_sim = np.array([r["bc_upper_allc"] for r in results_all_rbc_dep], dtype=float)
ax.plot(x, bc_l_sim, marker="o", color=cmap(1))
ax.plot(x, bc_u_sim, marker="o", color=cmap(1))
bc_u_sim_err = np.array([r["bc_upper_allc_err"] for r in results_all_rbc_dep], dtype=float)
ax.errorbar(x, bc_u_sim, yerr=bc_u_sim_err, fmt="none", ecolor=cmap(1))
bc_u_sim_err

bc_l = np.array([r["bc_lower_theory"] for r in results_all_rbc_dep], dtype=float)
bc_u = np.array([r["bc_upper_theory"] for r in results_all_rbc_dep], dtype=float)
ax.plot(x, bc_l, marker=".", color=cmap(0), markersize=2)
ax.plot(x, bc_u, marker=".", color=cmap(0), markersize=2)
s =2
ax.fill_between(np.append(x[s:], 1.05), np.append(bc_u[s:], bc_u[-1]), np.append(bc_l[s:], bc_l[-1]), color=cmap(0), alpha=0.2)

ax.plot([], [], color=cmap(0), label="theory")
ax.plot([], [], marker="o", color=cmap(1), label="simulation")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.legend(fontsize=14)

# %%
fig.savefig("interpolate_bc_rbc.pdf", bbox_inches="tight")
# %%

# get the path of this script
def run_simulation_return_pcs_q_dep(q, r_bc = 0.5):
  script_dir = os.path.dirname(os.path.abspath(__file__))
  executable_path = os.path.normpath(os.path.join(script_dir, "../cmake-build-release/main_CorrelatedOpinions"))
  t_measure = 10000
  param = {"N":100, "mu_assess": 0.02, "mu_impl": 0.0, "q": q, "t_init": t_measure/10, "t_measure": t_measure, "seed": 12345}
  print(json.dumps(param))
  num_sample = 100
  command = [executable_path, "-j", json.dumps(param), f"1.0 0.0 {r_bc} 1.0 1.0 0.0", f"{num_sample}"]
  print(command)
  out = subprocess.run(command, capture_output=True, text=True)
  print(out.stdout)
  print(out.stderr)
  results= json.loads(out.stdout)
  return results


# %%
result = run_simulation_return_pcs_q_dep(0.5, 0.5)
result
# %%
# if file does not exist, run the simulation
# if file exists, load the results
results_all_q_dep = []
if os.path.exists("results_all_q_dep.pkl"):
  with open("results_all_q_dep.pkl", "rb") as f:
    results_all_q_dep = pickle.load(f)
else:
  for i in range(0, 21):
    q = i / 20.0
    print(q)
    r = run_simulation_return_pcs_q_dep(q, 0.5)
    results_all_q_dep.append(r)
  with open("results_all_q_dep.pkl", "wb") as f:
    pickle.dump(results_all_q_dep, f)
results_all_q_dep
# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
n = len(results_all_q_dep)
x = np.array([i/(n-1) for i in range(0, n)], dtype=float)
h = np.array([r["h"] for r in results_all_q_dep], dtype=float)
hg = np.array([r["hG"] for r in results_all_q_dep], dtype=float)

ax.set_xlim(-0.01, 1.01)
ax.set_ylim(0.48, 1.02)
ax.plot(x, h, marker="o", label=r"$h$")
ax.plot(x, hg, marker="o", label=r"$h_G$")
ax.set_xlabel(r"$q$", fontsize=18)
ax.set_ylabel(r"$h, h_G$", fontsize=18)
ax.legend(fontsize=14)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

# %%
fig.savefig("interpolate_h_q.pdf", bbox_inches="tight")

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
cmap = mpl.cm.get_cmap("tab10")
# n = len(results_all)
# x = np.array([i/(n-1) for i in range(0, n)], dtype=float)
ax.set_xlim(-0.01, 1.01)
ax.set_ylim(1, 5)
ax.set_xlabel(r"$q$", fontsize=18)
ax.set_ylabel(r"$b/c$", fontsize=18)
ax.set_yticks([1, 2, 3, 4, 5])

bc_l_sim = np.array([r["bc_lower_alld"] for r in results_all_q_dep], dtype=float)
bc_u_sim = np.array([r["bc_upper_allc"] for r in results_all_q_dep], dtype=float)
ax.plot(x, bc_l_sim, marker="o", color=cmap(1))
ax.plot(x, bc_u_sim, marker="o", color=cmap(1))
bc_u_sim_err = np.array([r["bc_upper_allc_err"] for r in results_all_q_dep], dtype=float)
ax.errorbar(x, bc_u_sim, yerr=bc_u_sim_err, fmt="none", ecolor=cmap(1))
bc_u_sim_err

bc_l = np.array([r["bc_lower_theory"] for r in results_all_q_dep], dtype=float)
bc_u = np.array([r["bc_upper_theory"] for r in results_all_q_dep], dtype=float)
ax.plot(x, bc_l, marker=".", color=cmap(0), markersize=2)
ax.plot(x, bc_u, marker=".", color=cmap(0), markersize=2)
s =2
ax.fill_between(np.append(x[s:], 1.05), np.append(bc_u[s:], bc_u[-1]), np.append(bc_l[s:], bc_l[-1]), color=cmap(0), alpha=0.2)

ax.plot([], [], color=cmap(0), label="theory")
ax.plot([], [], marker="o", color=cmap(1), label="simulation")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

ax.legend(fontsize=14)

# %%
fig.savefig("interpolate_bc_q.pdf", bbox_inches="tight")

# %%
