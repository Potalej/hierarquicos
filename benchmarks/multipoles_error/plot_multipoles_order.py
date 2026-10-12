import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcParams["font.family"] = "monospace"
color = "#bb8cd1"
os.makedirs("img", exist_ok=True)

file = "out/multipoles_error.txt"
data = np.loadtxt(file)

thetas = [0.2, 0.4, 0.6, 0.8, 1.0][-1:]
multipoles = [1,4,8]
markers = ['o', "*", "^"]

fig, axs = plt.subplots(1, len(thetas), sharey=True, sharex=True, figsize=(4,3))
if len(thetas) == 1:
  axs = [axs]

grouped_theta = {}
for theta in thetas:
  mask = data[:,1] == theta
  grouped_theta[theta] = data[mask][:,:]

for i_theta, theta in enumerate(thetas):
  for i_multipole, multipole in enumerate(multipoles):
    # multipole, theta, eps, desc, error
    vs = grouped_theta[theta]
    mask = (vs[:,0] == multipole)
    data_list = vs[mask][:,:]

    dist = data_list[:,3]
    dist_unique = np.unique(dist)[1::2]
    error = data_list[:,4]
    error_mean = np.array([
      np.mean(error[dist == d])
      for d in dist_unique
    ])

    axs[i_theta].scatter(dist_unique, error_mean, marker=markers[i_multipole], label=rf"$p={multipole}$", s=15, c=color)

    Y0 = np.log(error_mean[0])
    Y1 = np.log(error_mean[-1])
    X0 = np.log(dist_unique[0])
    X1 = np.log(dist_unique[-1])
    b = (Y1 - Y0)/(X1 - X0)
    A = Y0 - b * X0
    a = np.exp(A)

    print(theta, multipole, A, b)
    b = - i_multipole - 2
    print(b)

    xs = np.linspace(100, 500, 600)
    axs[i_theta].plot(xs, a*xs**b, c='black', linestyle='--')

  # axs[i_theta].set_title(theta)
  axs[i_theta].legend(markerscale=2, scatterpoints=1)
  axs[i_theta].set_yscale('log')
  axs[i_theta].set_xscale('log')
  axs[i_theta].grid(True)
  axs[i_theta].set_xlabel(r"(Log) Distance $d$ between the particle and the node", fontsize=7)
  axs[i_theta].set_ylabel(r"$\log{|a_{bh} - a_{dir.}|/|a_{dir.}|}$", fontsize=9)

  axs[i_theta].text(280, 5e-6, r"$\alpha_1 d^{-2}$", fontsize=10)
  axs[i_theta].text(280, 5e-8, r"$\alpha_2 d^{-3}$", fontsize=10)
  axs[i_theta].text(280, 1e-9, r"$\alpha_3 d^{-4}$", fontsize=10)


plt.title("Relative error of BH with multipoles\n"+r"($\theta=1.0$)", fontsize=10)

plt.tight_layout()
plt.savefig("img/multipoles_order.png", dpi=200)