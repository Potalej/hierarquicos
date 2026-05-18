"""
To visualize the file fortran/out/forces_time.txt, with a time comparison
of the generation of quadtrees and evaluate of forces for different 
numbers of bodies and values of theta.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.loadtxt('fortran/out/forces_time.txt')
thetas = data[:,1]
thetas_uniques = np.unique(thetas)
grouped = {}
for theta in thetas_uniques:
    mask = data[:,1] == theta
    grouped[theta] = data[mask][:, [0,2,3]]

fig, axs = plt.subplots(1, 2, sharex=True, figsize=(10,3))

for theta in thetas_uniques:
    Ns = grouped[theta][:,0]
    times = grouped[theta][:,1]
    error = grouped[theta][:,2]

    if theta == -1.0:
        axs[0].scatter(Ns, times, c='black', label="Dir.", zorder=100, marker='+')
    else:
        axs[0].scatter(Ns, times, label=theta, s=5)
        axs[1].scatter(Ns, error, label=theta, s=5)

    # approximate as a*N^2 + b*NlogN
    f = lambda t, a, b: a * t * t + b * t * np.log(t)
    popt, pcov = curve_fit(f, Ns, times)
    Ns = np.unique(Ns)
    coefs = popt / max(np.abs(popt))
    print("{:.4e} {:.4e}".format(coefs[0], coefs[1]))
    axs[0].plot(Ns, f(Ns, *popt), c='black', linestyle='--')

axs[0].set_title("Time")
axs[0].set_ylabel("Time (s)")
axs[0].set_xlabel(r"$N$")
axs[0].legend()

axs[1].set_title("Error")
axs[1].set_ylabel(r"$||f_\theta - f_d||_2$")
axs[1].set_xlabel(r"$N$")
axs[1].set_yscale('log')
axs[1].legend()

plt.tight_layout()
plt.savefig("img/f_forces_time_and_error.png")