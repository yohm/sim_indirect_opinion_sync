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
def run_simulation_return_pcs(resident_strategy="L3", mutant_strategy="RANDOM-0.875", q=0.0, t_measure=100000, num_samples=100):
  # run simulation like the following
  # ./cmake-build-release/main_PayoffComparison -j '{"N": 100, "mu_assess":0.02,"mu_impl":0,"q":0.0,"t_init":100000,"t_measure":1000000,"tau":0}' L3 RANDOM-0.875 100
  script_dir = os.path.dirname(os.path.abspath(__file__))
  executable_path = os.path.normpath(os.path.join(script_dir, "../cmake-build-release/main_PayoffComparison"))
  param = {"N": 100, "mu_assess": 0.02, "mu_impl": 0.0, "q": q, "t_init": t_measure//10, "t_measure": t_measure, "seed": 12345}
  print(json.dumps(param))
  command = [executable_path, "-j", json.dumps(param), resident_strategy, mutant_strategy, f"{num_samples}"]
  print(command)
  out = subprocess.run(command, capture_output=True, text=True, env={"OMP_NUM_THREADS": "20"})
  print(out.stdout)
  print(out.stderr)
  # load from out.stdout as numpy array
  dat = np.fromstring(out.stdout, sep=" ")
  return dat.reshape([-1, 9])

# %%
cache_file_name = "payoff_comparison_cache.pkl"
if os.path.exists(cache_file_name):
  dat = pickle.load(open(cache_file_name, "rb"))
  out_random = dat["out_random"]
  out_allc = dat["out_allc"]
  out_q1_allc = dat["out_q1_allc"]
  out_q1_random = dat["out_q1_random"]
else:
  out_random = run_simulation_return_pcs(resident_strategy="L3", mutant_strategy="RANDOM-0.875", q=0.0)
  out_allc = run_simulation_return_pcs(resident_strategy="L3", mutant_strategy="AllC", q=0.0)
  out_q1_allc = run_simulation_return_pcs(resident_strategy="L3", mutant_strategy="AllC", q=1.0, t_measure=10000)
  out_q1_random = run_simulation_return_pcs(resident_strategy="L3", mutant_strategy="RANDOM-0.96", q=1.0, t_measure=10000)
  dat = {"out_random": out_random, "out_allc": out_allc, "out_q1_allc": out_q1_allc, "out_q1_random": out_q1_random}
  pickle.dump(dat, open(cache_file_name, "wb"))

out_random, out_allc

# %%
def calc_payoffs(out, benefit=5.0):
  fmut = out[:, 0]/100.0
  fres = 1.0 - fmut
  num_mut_list = out[:, 0]
  pi_res_list = benefit * (fres * out[:, 1] + fmut * out[:, 5]) - (fres * out[:, 1] + fmut * out[:, 3])
  pi_mut_list = benefit * (fres * out[:, 3] + fmut * out[:, 7]) - (fres * out[:, 5] + fmut * out[:, 7])
  return num_mut_list, pi_res_list, pi_mut_list


# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))

num_mut_list, pi_res_list, pi_mut_list = calc_payoffs(out_random)

ax.set_xlim(0, 100)
ax.set_ylim(3.15, 3.85)
ax.set_yticks(np.arange(3.2, 3.85, 0.2))
ax.plot(num_mut_list, pi_res_list, marker="o", markersize=6, label=r"$\pi_{\rm res}$")
ax.plot(num_mut_list, pi_mut_list, marker="o", markersize=3, label=r"$\pi_{\rm mut}$")
ax.set_xlabel("# of mutants", fontsize=14)
ax.set_ylabel("payoffs", fontsize=14)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.legend(fontsize=14, loc="lower right")

# %%
fig.savefig("L3_RANDOM-0.875_solitary.pdf", bbox_inches="tight")

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))

num_mut_list, pi_res_list_allc, pi_mut_list_allc = calc_payoffs(out_allc)

ax.set_xlim(0, 100)
ax.set_ylim(3.4, 4.1)
ax.set_yticks(np.arange(3.4, 4.1, 0.2))
ax.plot(num_mut_list, pi_res_list_allc, marker="o", markersize=6, label=r"$\pi_{\rm res}$")
ax.plot(num_mut_list, pi_mut_list_allc, marker="o", markersize=3, label=r"$\pi_{\rm mut}$")
ax.set_xlabel("# of mutants", fontsize=14)
ax.set_ylabel("payoffs", fontsize=14)
#ax.legend(fontsize=14, loc="lower right")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

# %%
fig.savefig("L3_ALLC_solitary.pdf", bbox_inches="tight")

# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))
cmap = mpl.cm.get_cmap("tab10")

num_mut_list, pi_res_list_allc, pi_mut_list_allc = calc_payoffs(out_q1_allc)

ax.set_xlim(0, 100)
ax.set_ylim(3.7, 4.1)
ax.set_yticks(np.arange(3.7, 4.11, 0.1))
ax.plot(num_mut_list, pi_res_list_allc, marker="o", markersize=6, label=r"$\pi_{\rm res}$")
ax.plot(num_mut_list, pi_mut_list_allc, marker="o", markersize=3, label=r"$\pi_{\rm mut}$")
ax.set_xlabel("# of mutants", fontsize=14)
ax.set_ylabel("payoffs", fontsize=14)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

#ax.text(10, 3.82, r"L3", fontsize=14, color=cmap(0))
#ax.text(10, 3.94, r"ALLC", fontsize=14, color=cmap(1))
#ax.legend(fontsize=14, loc="lower right")
# %%
fig.savefig("L3_ALLC_q1_solitary.pdf", bbox_inches="tight")
# %%
plt.clf()
fig, ax = plt.subplots(figsize=(6, 4))

num_mut_list, pi_res_list_allc, pi_mut_list_allc = calc_payoffs(out_q1_random)

ax.set_xlim(0, 100)
ax.set_ylim(3.7, 4.1)
ax.set_yticks(np.arange(3.7, 4.11, 0.1))
ax.plot(num_mut_list, pi_res_list_allc, marker="o", markersize=6, label=r"$\pi_{\rm res}$")
ax.plot(num_mut_list, pi_mut_list_allc, marker="o", markersize=3, label=r"$\pi_{\rm mut}$")
ax.set_xlabel("# of mutants", fontsize=14)
ax.set_ylabel("payoffs", fontsize=14)
#ax.legend(fontsize=14, loc="lower right")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='both', which='major', labelsize=14)

# %%
fig.savefig("L3_RANDOM-0.96_q1_solitary.pdf", bbox_inches="tight")

# %%
