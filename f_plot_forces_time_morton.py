"""
To visualize the file fortran/out/forces_time.txt, with a time comparison
of the generation of quadtrees and evaluate of forces for different 
numbers of bodies and values of theta.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

files = {
    "With Morton (Par. $P=20$)": "forces_time_morton_parallel_20C_N10K",
    "Without Morton (Par. $P=20$)": "forces_time_without_morton_parallel_20C_N10K"
}

fig, axs_2d2 = plt.subplots(2, 2, sharex=True, figsize=(10,6))

for i, title in enumerate(files):
    file = files[title]
    axs = axs_2d2[i,:]
    data = np.loadtxt(f'fortran/out/{file}.txt')

    thetas = data[:,1]
    thetas_uniques = np.unique(thetas)
    grouped = {}
    for theta in thetas_uniques:
        mask = data[:,1] == theta
        grouped[theta] = data[mask][:, [0,2,3]]


    for theta in thetas_uniques:
        Ns = grouped[theta][:,0]
        times = grouped[theta][:,1]
        error = grouped[theta][:,2]

        if theta == -1.0:
            axs[0].scatter(Ns, times, c='black', label="Dir.", zorder=100, marker='+')
        else:
            axs[0].scatter(Ns, times, label=theta, s=5)
            axs[1].scatter(Ns, error, label=theta, s=5)
            
            max_t = np.max(times)

            coef1 = 24.0 / (theta**2 * np.log(1.0 + theta))
            termo1 = coef1 * np.log(theta)
            termo2 = coef1 * np.log(3.0 * Ns / (4.0 * np.pi)) / 3.0
            termo3 = 4.0 * np.pi / (3.0 * theta**3)

            custo_estimado = termo1 + termo2 + termo3
            print(coef1)

        # approximate as a*N^2 + b*NlogN
        f = lambda t, a, b: a * t * t + b * t * np.log(t)
        # f = lambda t, a: a * t * np.log(t)
        popt, pcov = curve_fit(f, Ns, times)
        Ns = np.unique(Ns)
        coefs = popt
        # coefs = popt / max(np.abs(popt))
        print("{:.4e} {:.4e}".format(coefs[0], coefs[1]))
        
        print(f"{1.0/(theta**2)}")
        print()
        axs[0].plot(Ns, f(Ns, *popt), c='black', linestyle='--')

    # axs[0].set_ylim(0,0.0020)
    axs[0].set_title("Time")
    axs[0].set_ylabel("Time (s)")
    axs[0].set_xlabel(r"$N$")
    axs[0].set_yscale('log')
    axs[0].axhline(1e-1, linestyle='--', c='black')
    axs[0].axhline(1e-3, linestyle='--', c='black')
    axs[0].set_ylim(0,5e-1)
    axs[0].legend()

    axs[1].set_title("Error")
    axs[1].set_ylabel(r"$||f_\theta - f_d||/||f_d||$")
    axs[1].set_xlabel(r"$N$")
    axs[1].set_yscale('log')
    axs[1].legend()

fig.text(0.5, 0.95, list(files.keys())[0], ha='center', fontsize=14)
fig.text(0.5, 0.48, list(files.keys())[1], ha='center', fontsize=14)

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.savefig("img/morton_vs_intuitive.png")