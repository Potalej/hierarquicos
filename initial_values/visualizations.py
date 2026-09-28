import api as iv
import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

os.makedirs("img", exist_ok=True)
color = "#bb8cd1"
M = 1.0
G = 1.0

########################################################################################################
# UTILITIES
########################################################################################################

def sample (N:int, parameters:list, profile:str, size_qs:float=-1.0, size_ps:float=-1.0, td:bool=False):
  N = np.array(N)
  m = np.zeros(N, order='F')
  qs = np.zeros((3,N), order='F')
  ps = np.zeros((3,N), order='F')
  iv.initial_values_mod.api(size_qs, size_ps, True, td, profile, parameters, m, qs, ps)
  return m, qs, ps

def plot_points (qs:np.array, title:str):
  plt.figure(figsize=(5,4))
  plt.title(title)
  plt.xlim(-5., 5.)
  plt.ylim(-5., 5.)
  plt.scatter(qs[0,:], qs[1,:], c=qs[2,:], cmap='viridis', s=5)
  plt.axis("equal")
  cbar = plt.colorbar()
  cbar.set_label(r"$z$", rotation=0)
  plt.tight_layout()

def plot_histogram_radius (rs, r_axis, y_axis, title, bins):
  plt.figure(figsize=(6,3))
  plt.title(title)
  plt.grid(True)
  plt.ylabel(r"PDF and ocurrence")
  plt.xlabel(r"radius $r$")
  # theory
  plt.plot(r_axis, y_axis, c='black', label="PDF")
  # sample
  plt.hist(rs, bins=bins, color=color, density=True, label="Sample", edgecolor="#212121")
  plt.legend()
  plt.tight_layout()

def plot_sigmar2 (rmid, sigma_sample, rs_theory, sigma_theory, title):
  plt.figure(figsize=(6,3))
  plt.title(title)
  plt.grid(True)  
  plt.plot(rmid, sigma_sample, c=color, label="Sample")
  plt.plot(rs_theory, sigma_theory, c='black', label="Jeans eq. sol.")
  plt.legend()
  plt.tight_layout()

def get_sigma_r2 (r:float, rho, dphi):
  return (1.0/rho(r)) * quad(lambda s: rho(s) * dphi(s), r, np.inf)[0]

def get_sigma_sample (rs, vrs, rbins):
  rmid = 0.5*(rbins[:-1] + rbins[1:])
  sigma_sample = np.zeros(len(rmid))
  for i in range(len(rmid)):
    mask = ((rs >= rbins[i]) & (rs < rbins[i+1]))
    if np.sum(mask) > 2: 
      sigma_sample[i] = np.var(vrs[mask])
  return rmid, sigma_sample

def get_radius (qs:np.array):
  return np.sqrt(np.array([qs[0,i]**2 + qs[1,i]**2 + qs[2,i]**2 for i in range(len(qs[0,:]))]))

def get_vrs (qs, vs, rs):
  return np.array([
    (qs[0,i]*vs[0,i] + qs[1,i]*vs[1,i] + qs[2,i]*vs[2,i])/rs[i]
    for i in range(len(qs[0,:]))])

########################################################################################################
# EXAMPLES
########################################################################################################

def example_homogeneous_sphere (out="homogeneous_sphere"):
  N = 10_000
  R = 3.0
  m, qs, ps = sample(N, [G, R], "homogeneous")

  # points
  plot_points(qs, rf"Homogeneous Sphere Profile $(R={R}, N=10^4)$")
  plt.savefig(f"img/{out}.png")


def example_plummer_sphere (outs:list):
  N = 10_000
  b = 0.5
  
  # rmax = 10 r_h, as suggested by Aarseth
  rmax = 10. * b * 0.5**(1./3.) /np.sqrt(1. - 0.5**(2.3))

  # sampling
  m, qs, ps = sample(N, [G, b], "plummer", size_qs=rmax)

  # points
  plot_points(qs, rf"Plummer Model ($b={b}$, $N=10^4$, cut-off $= 10 r_h$)")
  plt.savefig(f"img/{outs[0]}.png")
  plt.close()

  # histogram of the radius
  rs = get_radius(qs)
  r_axis = np.linspace(0, rmax, 200)
  y_axis = [iv.iv_plummer_mod.plummer_radius_pdf(r, [G, b]) for r in r_axis]

  plot_histogram_radius(rs, r_axis, y_axis, 
    title = r"Plummer radius sample ($N=10^4$, $b=1/2$, cut-off at $r=10r_h$)",
    bins=50)
  plt.savefig(f"img/{outs[1]}.png")
  plt.close()
  
  
  # velocities 
  # get the sigma_2 from the sample
  vrs = get_vrs(qs, ps / m, rs)
  rbins = np.logspace(-0.5, np.log10(rmax), 200)
  rmid, sigma_sample = get_sigma_sample(rs, vrs, rbins)
  print(sigma_sample)
  
  # get sigma_2 from Jeans equations
  rs_theory = rbins
  rho = lambda r: iv.iv_plummer_mod.plummer_density(r, [G, b])
  dphi = lambda r: iv.iv_plummer_mod.plummer_dphi(r, [G, b])
  sigma_theory = np.array([
    get_sigma_r2(r, rho, dphi) for r in rs_theory
  ])

  plot_sigmar2(rmid, sigma_sample, rbins, sigma_theory, 
    title=rf"Plummer $\sigma_r^2$ ($N=10^4$, $b={b}$, cut-off at $r=10r_h$)")
  plt.ylim(0., 0.5)

  plt.savefig(f"img/{outs[2]}.png")
  plt.close()

def example_hernquist_sphere (outs:list):
  N = 10_000
  a = 1.0
  
  # rmax = 10 r_h, as suggested by Aarseth
  rmax = 10. * a * np.sqrt(0.5 * M)/(1.0 - np.sqrt(M*0.5))

  # sampling
  m, qs, ps = sample(N, [G, a], "hernquist", size_qs=rmax)

  # points
  plot_points(qs, rf"Hernquist Model ($a={a}$, $N=10^4$, cut-off $= 10 r_h$)")
  plt.savefig(f"img/{outs[0]}.png")
  plt.close()

  # histogram of the radius
  rs = get_radius(qs)
  r_axis = np.linspace(0, rmax, 200)
  y_axis = [iv.iv_hernquist_mod.hernquist_radius_pdf(r, [G, a]) for r in r_axis]

  plot_histogram_radius(rs, r_axis, y_axis, 
    title = r"Hernquist radius sample ($N=10^4$, $a=1$, cut-off at $r=10r_h$)",
    bins=50)
  plt.savefig(f"img/{outs[1]}.png")
  plt.close()
  

  # velocities 
  # get the sigma_2 from the sample
  vrs = get_vrs(qs, ps / m, rs)
  rbins = np.logspace(-2., np.log10(rmax), 200)
  rmid, sigma_sample = get_sigma_sample(rs, vrs, rbins)
  
  # get sigma_2 from Jeans equations
  rs_theory = rbins
  rho = lambda r: iv.iv_hernquist_mod.hernquist_density(r, [G, a])
  dphi = lambda r: iv.iv_hernquist_mod.hernquist_dphi(r, [G, a])
  sigma_theory = np.array([
    get_sigma_r2(r, rho, dphi) for r in rs_theory
  ])

  plot_sigmar2(rmid, sigma_sample, rbins, sigma_theory, 
    title=rf"Hernquist $\sigma_r^2$ ($N=10^4$, $a={a}$, cut-off at $r=10r_h$)")

  plt.savefig(f"img/{outs[2]}.png")
  plt.close()
    

if __name__ == "__main__":
  # Homogeneous Sphere
  example_homogeneous_sphere("homogeneous_sphere")

  # Plummer Sphere
  example_plummer_sphere(["plummer_sphere", "plummer_density", "plummer_sigmar2"])

  # Hernquist Sphere
  example_hernquist_sphere(["hernquist_sphere", "hernquist_density", "hernquist_sigmar2"])