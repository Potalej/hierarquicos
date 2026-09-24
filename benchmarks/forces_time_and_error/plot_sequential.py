import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import f90nml
from colour import Color
import os

out_folder = "img/sequential"
os.makedirs(out_folder, exist_ok=True)
plt.rcParams["font.family"] = "monospace"
folder = "out/sequential"

files = {
    "Monopole":   f"{folder}/test_monopole.txt",
    "Quadrupole": f"{folder}/test_quadrupole.txt",
    "Octupole":   f"{folder}/test_octupole.txt",
}
f90in = f"{folder}/preset.nml"
nml = f90nml.read(f90in)['preset']
num_tests = nml['number_of_tests']

for i, title in enumerate(files):
    fig, axs = plt.subplots(1, 2, sharex=True, figsize=(10,4))
    
    file = files[title]
    data = np.loadtxt(file)

    thetas = data[:,1]
    # thetas_uniques = np.unique(thetas)
    thetas_uniques = [0.2, 0.4, 0.6, 0.8, 1.0, -1]
    grouped = {}
    
    for theta in thetas_uniques:
        mask = data[:,1] == theta
        grouped[theta] = data[mask][:,:]
        
    red = Color("violet")
    colors = list(red.range_to(Color("red"),len(thetas_uniques)-1))
    colors = [c.hex for c in colors]
    
    for i_theta, theta in enumerate(thetas_uniques):
        # if theta < 1e-8: continue
        Ns = grouped[theta][:,0]
        times = grouped[theta][:,3]
        error = grouped[theta][:,4]
        
        x_axis, y_axis = [], []
        for test in range(num_tests):
            x_axis.append(Ns[test*num_tests])
            y_axis.append(np.mean(
                    error[test*num_tests:(test+1)*num_tests]
            ))
            
        if theta == -1.0:
            axs[0].scatter(Ns, times, c='black', label="Dir.", zorder=100, marker='+', s=40)
        else:
            label = rf"$\theta={theta}$"
            axs[0].scatter(Ns, times, s=30, label=label, c=colors[i_theta])
            axs[1].scatter(x_axis, y_axis, s=30, label=label, c=colors[i_theta])
        
    axs[0].set_yscale('log')
    axs[0].set_ylabel(r"time in seconds")
    axs[0].set_xlabel(r"$N$")
    axs[0].grid(True)
    axs[0].legend(ncol=2)
    axs[0].set_ylim(1e-4, 1e1)

    axs[1].set_yscale('log')
    axs[1].set_ylabel(r"$\delta a/a$ (mean)")
    axs[1].set_xlabel(r"$N$")
    axs[1].legend(ncol=2)
    axs[1].grid(True)
    axs[1].set_ylim(1e-7, 1e-1)
    
    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.suptitle(f"{title}", fontsize=16)
    plt.savefig(f"{out_folder}/times_{title}_sequential.png", dpi=200)
    plt.close()