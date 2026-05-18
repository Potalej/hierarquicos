"""
To estimate the cousts of a quadtree generation.
"""
import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from time import time
from fortran import api

times = []
times_mean = []
N_min, N_max, N_step = 10, 10000, 10
num_tests_N = 50

# testing for each N
for N in range(N_min, N_max + 1, N_step):
  # for each N, tests a quantity of times
  times_N = []
  for _ in range(num_tests_N):
    # random values
    masses = np.random.random(N)
    qs = 5 * (2 * np.random.random((N,2)) - 1)
    x, y = np.array(list(zip(*qs)))
    
    # timing the tree generation
    time_start = time()
    api.api_mod.generate_tree(masses, x, y)
    time_total = time() - time_start
    
    # saving
    times.append([N, time_total])
    times_N.append(time_total)
  
  # save the mean
  times_mean.append([N, np.mean(times_N)])

# ploting
fig, ax = plt.subplots(figsize=(8,3))
Ns, ts = list(zip(*times_mean))
ax.scatter(Ns, ts, c='black')

# curve fit of a N log N
f = lambda t, a: a * t * np.log(t)
Ns, ts = list(zip(*times))
popt, pcov = curve_fit(f, Ns, ts)
Ns = np.linspace(N_min, N_max, 100)
ax.plot(Ns, f(Ns, *popt), c='red', linestyle='--')
ax.text(Ns[len(Ns) // 2] + 10, 0.8*f(Ns[len(Ns) // 2], *popt), rf"$m={popt[0]:.2e}$")

ax.grid(True)
ax.set_ylabel("Time (s)")
ax.set_xlabel(r"$N$")
plt.title(r"Time to create a quadtree by $N$.")
plt.tight_layout()
plt.savefig("img/f2py_tree_times.png")