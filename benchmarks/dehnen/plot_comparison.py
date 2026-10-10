import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import f90nml
from colour import Color
import os

color = "#bb8cd1"
out_folder = "img/direct_vs_bh_vs_dehnen"
os.makedirs(out_folder, exist_ok=True)
plt.rcParams["font.family"] = "monospace"
folder = "out/direct_vs_bh_vs_dehnen"

files = {
  "BHM": f"{folder}/bh_monopole.txt",
  "BHQ": f"{folder}/bh_quadrupole.txt",
  "DM": f"{folder}/dehnen_monopole.txt",
  "DQ": f"{folder}/dehnen_quadrupole.txt",
}
f90in = f"{folder}/preset.nml"
nml = f90nml.read(f90in)['preset']
num_tests = nml['number_of_tests']

for i, title in enumerate(files):
    fig, axs = plt.subplots(1, 2, sharex=True, figsize=(10,4))
    
    file = files[title]
    data = np.loadtxt(file)
    # data= [N, theta, tree_time, time_total, error, time_forces]

    thetas = data[:,1]
    # thetas_uniques = np.unique(thetas)
    thetas_uniques = [0.2, 0.4, 0.6, 0.8, 1.0, -1]
    grouped = {}
    Ns_uniques = np.unique(data[:,0])
    
    for theta in thetas_uniques:
        mask = data[:,1] == theta
        grouped[theta] = data[mask][:,:]
        
    red = Color("violet")
    colors = list(red.range_to(Color("red"),len(thetas_uniques)-1))
    colors = [c.hex for c in colors]
    
    for i_theta, theta in enumerate(thetas_uniques):
        # if theta < 1e-8: continue
        Ns = grouped[theta][:,0]
        times = grouped[theta][:,5]
        error = grouped[theta][:,4]
        
        x_axis, y_axis = [], []
        times_axis = []
        for i_N, N in enumerate(Ns_uniques):
            # if N < 2000: continue
            x_axis.append(N)
            y_axis.append(np.mean(
                error[i_N*num_tests:(i_N+1)*num_tests]
            ))
            times_axis.append(np.mean(
                times[i_N*num_tests:(i_N+1)*num_tests]
            ))
            
        if theta == -1.0:
            # pass
            axs[0].scatter(Ns, times, c='black', label="Dir.", zorder=100, marker='+', s=40)
        else:
            label = rf"$\theta={theta}$"
            # axs[0].scatter(Ns, times, s=30, label=label, c=colors[i_theta])
            axs[0].scatter(x_axis, times_axis, s=30, label=label, c=colors[i_theta])
            axs[1].scatter(x_axis, y_axis, s=30, label=label, c=colors[i_theta])

            # f_curve_fit = lambda x, a, b: a * x * np.log(x) + b
            # f_curve_fit = lambda x, a, b, c: a * x **2 + b * x + c
            x_axis = np.array(x_axis)
            
            f_curve_fit = lambda x, a, b: a * x + b
            params, cov = curve_fit(f_curve_fit, x_axis, times_axis)
            axs[0].plot(x_axis, f_curve_fit(x_axis, *params), linestyle='--', c=colors[i_theta])

            f_curve_fit = lambda x, a, b: a * x * np.log(x) + b
            params, cov = curve_fit(f_curve_fit, x_axis, times_axis)
            axs[0].plot(x_axis, f_curve_fit(x_axis, *params), linestyle=':', c=colors[i_theta])
        
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"time in seconds")
    axs[0].set_xlabel(r"$N$")
    axs[0].grid(True)
    axs[0].legend(ncol=2)
    axs[0].set_ylim(1e-4, 2e1)

    axs[1].set_yscale('log')
    axs[1].set_ylabel(r"$\delta a/a$ (mean)")
    axs[1].set_xlabel(r"$N$")
    axs[1].legend(ncol=2)
    axs[1].grid(True)
    axs[1].set_ylim(1e-7, 5e-1)
    
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.suptitle(f"{title}", fontsize=16)
    plt.savefig(f"{out_folder}/{title}.png", dpi=200)
    plt.close()