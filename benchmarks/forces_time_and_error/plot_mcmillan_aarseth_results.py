import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

out_folder = "img/aarseth"
os.makedirs(out_folder, exist_ok=True)

plt.rcParams["font.family"] = "monospace"

files = {
    "Monopole":   "out/aarseth/test_monopole.txt",
    "Quadrupole": "out/aarseth/test_quadrupole.txt",
    "Octupole":   "out/aarseth/test_octupole.txt",
}

fig, axs = plt.subplots(1, 1, sharex=True, figsize=(5,4))
colors = ["red", "green", "blue"]
symbols = ['^', '+']

points1024 = []
points4096 = []

for i, title in enumerate(files):
    file = files[title]
    data = np.loadtxt(file)

    thetas = data[:,1]
    thetas_uniques = np.unique(thetas)
    grouped = {}
    
    is_set_label1024, is_set_label4096 = False, False
    
    for theta in thetas_uniques:
        mask = data[:,1] == theta
        grouped[theta] = data[mask][:,:]

    for theta in thetas_uniques:
        if theta < 1e-8: continue
        Ns = grouped[theta][:,0]
        times = grouped[theta][:,3]
        error = grouped[theta][:,4]
        
        half_length = len(Ns) // 2
        median1024 = np.median(error[:half_length])
        median4096 = np.median(error[half_length:])
        
        if i == 0:
            points1024.append([theta, median1024])
            points4096.append([theta, median4096])
        
        # label1024 = rf"{title} / $N={1024}$" if not is_set_label1024 else None
        # is_set_label1024 = True
        # axs.scatter(theta, median1024, s=30, c=colors[i], marker=symbols[0], label=label1024)
        axs.scatter(theta, median1024, c=colors[i], marker=symbols[0])
        
        # label4096 = rf"{title} / $N={4096}$"  if not is_set_label4096 else None
        # is_set_label4096 = True
        # axs.scatter(theta, median4096, s=30, c=colors[i], marker=symbols[1], label=label4096)
        axs.scatter(theta, median4096, c=colors[i], marker=symbols[1])

plt.scatter(10, 10, c='black', marker=symbols[0], label=r"$N=1024$")
plt.scatter(10, 10, c='black', marker=symbols[1], label=r"$N=4096$")
plt.text(0.125, 4.7e-7, "octupole")
plt.text(0.112, 2.2e-6, "quadrupole")
plt.text(0.12, 3.9e-5, "monopole")

plt.ylim(2e-7, 5e-2)
plt.xlim(0.1, 1.1)
plt.ylabel(r"$\delta a/a$ (median)")
plt.xlabel(r"$\theta$")
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0], 
          [0.1, "", 0.3, "", "", 0.6, "", "", "", 1.0])
plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig(f"{out_folder}/aarseth_result.png", dpi=200)