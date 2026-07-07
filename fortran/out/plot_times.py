import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({'font.size': 14})
markersize = 5.
markers = ['o', 'x', 's', 'd', '^']


data_tree_time = np.loadtxt("octree_time.txt", dtype=float)
data_poles_time = np.loadtxt("octpoles_time.txt", dtype=float)
data_forces_time = np.loadtxt("forces_time_without_morton_parallel_20C_N10K.txt", dtype=float)
Ns = np.unique(data_tree_time[:,0])
##avg_time = np.average(data_tree_time[:,1])
avg_time = 1e-8


## Plot tree consctrution and poles comptuation
for i in range(2):
    fig, ax = plt.subplots()

    if i == 0:
        data = np.copy(data_tree_time)
        figname = 'tree_time'
    else:
        data = np.copy(data_poles_time)
        figname = 'poles_time'

    ax.scatter(data[:,0], data[:,1], s = markersize, color = 'C0', label = 'data')
    ##ax.scatter(data_forces_time[:,0], data_forces_time[:,2], s = 1., color = 'C1')
    ax.plot(Ns, avg_time * Ns, color  = 'C1', linestyle = '--', label = r'$\mathcal{O}(N)$')
    ax.plot(Ns, avg_time * Ns * np.log(Ns), color  = 'C2', linestyle = '--', label = r'$\mathcal{O}(N \log N)$')
    ax.plot(Ns, avg_time * Ns**2, color  = 'C3', linestyle = '--', label = r'$\mathcal{O}(N^2)$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$N$')
    ax.set_ylabel("Time (s)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(figname + '.pdf')

## Plot force: time and errors
fig, ax = plt.subplots()
fig2, ax2 = plt.subplots()
for itheta in range(5):
    data = np.copy(data_forces_time[itheta::5,:])
    assert np.all(data[:,1] == data[0,1])

    if itheta == 0:
        label = "direct"
        color = "C0"
    else:
        label = r'$\theta = {}$'.format(data[0,1])
        color = "C" + str(itheta + 3)

    ax.scatter(data[:,0], data[:,2], s = markersize, marker = markers[itheta], color = color, label = label)
    if itheta > 0:
        ax2.scatter(data[:,0], data[:,3], s = markersize, marker = markers[itheta], color = color, label = label)

ax.plot(Ns, avg_time * Ns, color  = 'C1', linestyle = '--', label = r'$\mathcal{O}(N)$')
ax.plot(Ns, avg_time * Ns * np.log(Ns), color  = 'C2', linestyle = '--', label = r'$\mathcal{O}(N \log N)$')
ax.plot(Ns, avg_time * Ns**2, color  = 'C3', linestyle = '--', label = r'$\mathcal{O}(N^2)$')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$N$')
ax.set_ylabel("Time (s)")
ax.legend()
fig.tight_layout()
fig.savefig('forces_time.pdf')

ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlabel(r'$N$')
ax2.set_ylabel("Error")
ax2.legend()
fig2.tight_layout()
fig2.savefig('forces_error.pdf')


plt.show()


