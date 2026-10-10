"""
To visualize the file fortran/out/forces_time.txt, with a time comparison
of the generation of quadtrees and evaluate of forces for different 
numbers of bodies and values of theta.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os
from colour import Color
import f90nml

out_folder = "img/forces_timeR3"
os.makedirs(out_folder, exist_ok=True)
plt.rcParams["font.family"] = "monospace"
folder = "out/forces_timeR3"

files = {
    "Monopole": "test_multipoles_thetas_monopole",
    "Quadrupole": "test_multipoles_thetas_quadrupole",
    "Octupole": "test_multipoles_thetas_octupole"
}

f90in = f"{folder}/preset.nml"
nml = f90nml.read(f90in)['preset']
num_tests = nml['number_of_tests']
num_threads = nml['number_of_threads']

for title in files:
    inf = f"{folder}/{files[title]}.txt"
    outf = f"{out_folder}/{files[title]}.png"

    fig, axs = plt.subplots(1, 2, sharex=True, figsize=(10,3))

    data = np.loadtxt(inf)
    thetas = data[:,1]

    # thetas_uniques = np.unique(thetas)
    thetas_uniques = [0.25, 0.5, 0.75, 1.0, -1.0]
    grouped = {}
    for theta in thetas_uniques:
        mask = data[:,1] == theta
        grouped[theta] = data[mask][:,:]

    red = Color("violet")
    colors = list(red.range_to(Color("red"),len(thetas_uniques)-1))
    colors = [c.hex for c in colors]

    for i_theta, theta in enumerate(thetas_uniques):
        Ns = grouped[theta][:,0]
        times = grouped[theta][:,3]
        error = grouped[theta][:,4]

        x_axis, y_axis = [], []
        times_axis = []
        for test in range(num_tests):
            x_axis.append(Ns[test*num_tests])
            y_axis.append(np.mean(
                error[test*num_tests:(test+1)*num_tests]
            ))
            times_axis.append(np.mean(
                times[test*num_tests:(test+1)*num_tests]
            ))

        if theta == -1.0:
            pass
            # axs[0].scatter(Ns, times, c='black', label="Dir.", zorder=100, marker='+')
        else:
            axs[0].scatter(x_axis, times_axis, label=theta, s=30, c=colors[i_theta])
            axs[1].scatter(x_axis, y_axis, label=theta, s=30, c=colors[i_theta])

        # approximate as a*N^2 + b*NlogN
        f = lambda t, a, b: a * t * t + b * t * np.log(t)
        popt, pcov = curve_fit(f, Ns, times)
        Ns = np.unique(Ns)
        coefs = popt / max(np.abs(popt))
        print("{:.4e} {:.4e}".format(coefs[0], coefs[1]))
        # axs[0].plot(Ns, f(Ns, *popt), c='black', linestyle='--')

    axs[0].set_title("Time")
    axs[0].set_ylabel("Time (s)")
    axs[0].set_xlabel(r"$N$")
    axs[0].set_yscale('log')
    axs[0].set_ylim(1e-5, 1.0)
    axs[0].grid(True)
    axs[0].legend(ncol=2)

    axs[1].set_title("Error")
    axs[1].set_ylabel(r"$\delta a/a$ (mean)")
    axs[1].set_xlabel(r"$N$")
    axs[1].set_yscale('log')
    axs[1].set_ylim(1e-7, 1e-1)
    axs[1].grid(True)
    axs[1].legend(ncol=2)

    plt.tight_layout(rect=[0, 0, 1, 0.93])
    fig.suptitle(f"{title} (P={num_threads})", fontsize=16)
    plt.savefig(outf, dpi=200)
    plt.close()