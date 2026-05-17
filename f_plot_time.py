"""
To visualize the file fortran/out/tree_time.txt, with a time comparison
of the generation of quadtrees for different numbers of bodies.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

data = np.loadtxt('fortran/out/tree_time.txt')
Ns = data[:,0]
times = data[:,1]

tempo_por_N = []

i = 0
for N in np.unique(Ns):
    t_N = []
    while i < len(times) - 1:
        i += 1
        if data[i,0] == N: t_N.append(data[i,1])
        else: break
    tempo_por_N.append([N, np.mean(t_N)])

Ns, times = np.array(list(zip(*tempo_por_N)))
plt.scatter(Ns, times, c='red', marker='+')

f = lambda t, a: a * t * np.log(t)
popt, pcov = curve_fit(f, Ns, times)
print(pcov)
plt.plot(Ns, f(Ns, *popt), c='black', linestyle='--')

plt.savefig("img/f_tree_times.png")