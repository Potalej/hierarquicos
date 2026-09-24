import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import f90nml
from colour import Color
import os

color = "#bb8cd1"
out_folder = "img/tree_construction_time"
os.makedirs(out_folder, exist_ok=True)
plt.rcParams["font.family"] = "monospace"
folder = "out/bigger_N"

files = {
    "Monopole":   f"{folder}/test_monopole.txt",
    "Quadrupole": f"{folder}/test_quadrupole.txt",
    "Octupole":   f"{folder}/test_octupole.txt",
}
f90in = f"{folder}/preset.nml"
nml = f90nml.read(f90in)['preset']
num_tests = nml['number_of_tests']

fig, axs = plt.subplots(1, 1, sharex=True, figsize=(10,4))

for i, title in enumerate(files):
    
    file = files[title]
    data = np.loadtxt(file)

    Ns = np.array([N for N in range(nml['nmin'], nml['nmax']+1, nml['nstep'])])
    y_axis_list = []
    
    grouped = {}
    for N in Ns:
        mask = data[:,0] == N
        grouped[N] = data[mask][:,:]
        
        mask = grouped[N][:,1] != -1
        grouped[N] = grouped[N][mask]
    
    for i, N in enumerate(Ns):
        if N < 1e4: continue
        y_axis = np.mean(grouped[N][:,2])
        y_axis_list.append([N, y_axis])
        # for t in grouped[N][:,2]: y_axis_list.append([N, t])
            
    y_axis_list = np.array(y_axis_list)
    xs, ys = y_axis_list[:,0], y_axis_list[:,1]
    # axs.scatter(xs, ys / (xs * np.log(xs)), label=title)
    axs.scatter(xs, ys, label=title)
    
    f_curve_fit = lambda x, a, b: a * x * np.log(x) + b
    # f_curve_fit = lambda x, a, b, c: a * x **2 + b * x + c
    params, cov = curve_fit(f_curve_fit, xs, ys)
    axs.plot(xs, f_curve_fit(xs, *params), linestyle='--')

axs.set_xscale('log')
axs.set_yscale('log')
axs.set_ylabel(r"time in seconds")
axs.set_xlabel(r"$N$")
axs.grid(True)
axs.legend()
axs.set_ylim(2e-3, 1e-1)

plt.tight_layout(rect=[0, 0, 1, 0.93])

fig.suptitle(f"Tree construction time", fontsize=16)
plt.savefig(f"{out_folder}/sequential.png", dpi=200)
plt.close()