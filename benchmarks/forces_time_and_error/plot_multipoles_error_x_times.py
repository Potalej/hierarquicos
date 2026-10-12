import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import f90nml
import os

out_folder = "img/multipoles_error_x_time"
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
num_threads = nml['number_of_threads']

fig, axs = plt.subplots(3, 1, figsize=(6,8), sharex=True, sharey=True)

for i, title in enumerate(files):
    
    file = files[title]
    data = np.loadtxt(file)

    thetas = data[:,1]
    # thetas_uniques = np.unique(thetas)
    thetas_uniques = [0.25, 0.5, 0.75, 1.0, -1]
    grouped = {}
    
    for theta in thetas_uniques:
        mask = data[:,1] == theta
        grouped[theta] = data[mask][:,:]
        
    Ns_uniques = np.unique(data[:,0])
    markers = ["o", "^", "*", "d"]

    for i_theta, theta in enumerate(thetas_uniques):
        if theta < 1e-8: continue
        Ns = grouped[theta][:,0]
        times = grouped[theta][:,3]
        error = grouped[theta][:,4]
        
        x_axis, y_axis = [], []
        times_axis = []
        for N in Ns_uniques:
            y_axis.append(np.mean(
                error[grouped[theta][:,0] == N]
            ))
            times_axis.append(np.mean(
                times[grouped[theta][:,0] == N]
            ))
        
        # cmap = "Greys"
        cmap = "viridis"
        label = rf"$\theta={theta}$"
        # sc = axs[i].scatter(y_axis, times_axis, c=Ns_uniques, cmap=cmap, s=30, label=label, marker=markers[i_theta], edgecolors="black")
        sc = axs[i].scatter(y_axis, times_axis, c=Ns_uniques, cmap=cmap, s=30, label=label, marker=markers[i_theta])

    axs[i].set_xscale('log')
    axs[i].set_yscale('log')
    axs[i].set_ylabel(r"time in seconds")
    axs[i].grid(True)
    axs[i].set_ylim(1e-4, 1e1)
    axs[i].set_title(title)

    cbar = plt.colorbar(sc)
    cbar.set_label(r"$N$", rotation=0)

axs[0].legend(ncol=1)
axs[1].legend(ncol=2)
axs[2].legend(ncol=2)

plt.ylabel("CPU time in seconds")
plt.xlabel(r"$\delta a/a$ (mean)")

plt.tight_layout(rect=[0, 0, 1, 0.93])
fig.suptitle(f"BH with Multipoles error x time\n({num_threads} threads)", fontsize=16)
plt.savefig(f"{out_folder}/multipoles.png", dpi=200)
plt.close()