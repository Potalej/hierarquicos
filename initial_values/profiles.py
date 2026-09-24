import matplotlib.pyplot as plt
import numpy as np
import api as iv
import os

color = "#bb8cd1"

"""
HOMOGENEOUS SPHERE PROFILE
"""
def homogeneous_sphere_profile (N:int, rho0:float, R:float):
  N = np.array(N)
  m = np.zeros(N, order='F')
  qs = np.zeros((3,N), order='F')
  ps = np.zeros((3,N), order='F')
  iv.initial_values_mod.api(-1, -1, True, False, "homogeneous", [1.0, R, rho0], m, qs, ps)
  
  return m, qs, ps
  
def example_homogeneous_sphere ():
  M = 1.0
  G = 1.0

  N = 10000
  m, qs, ps = homogeneous_sphere_profile(N, rho0=0.1, R=3.0)
  x, y, z = qs[0,:], qs[1,:], qs[2,:]

  plt.figure(figsize=(5,4))
  plt.xlim(-5., 5.)
  plt.ylim(-5., 5.)
  plt.title(r"Homogeneous Sphere Profile $(R=3, N=10^4)$")
  plt.scatter(x, y, c=z, cmap='viridis', s=5)
  plt.axis("equal")
  cbar = plt.colorbar()
  cbar.set_label(r'$z$', rotation=0)
  plt.tight_layout()
  plt.savefig("img/homogeneous_sphere.png")

"""
PLUMMER SPHERE PROFILE
"""
def plummer_sphere_profile (N:int, b:float, rmax:float=0.0):
  N = np.array(N)
  m = np.zeros(N, order='F')
  qs = np.zeros((3,N), order='F')
  ps = np.zeros((3,N), order='F')
  iv.initial_values_mod.api(rmax, -1, True, False, "plummer", [1.0, b], m, qs, ps)
  
  return m, qs, ps

def example_plummer_sphere ():
  M = 1.0
  G = 1.0
  
  b = 0.5
  k = 0.5
  rmax = 10. * (k**(1./3.)/np.sqrt(1.-k**(2./3.))) * b

  N = 10000
  m, qs, ps = plummer_sphere_profile(N, b=0.5, rmax=rmax)
  x, y, z = qs[0,:], qs[1,:], qs[2,:]
  
  plt.figure(figsize=(5,4))
  plt.title(r"Plummer Sphere Profile $(b=0.5, N=10^4)$")
  plt.scatter(x, y, c=z, cmap='viridis', s=5)
  plt.axis("equal")
  plt.xlim(-7., 7.)
  plt.ylim(-7., 7.)
  cbar = plt.colorbar()
  cbar.set_label(r'$z$', rotation=0)
  plt.tight_layout()
  plt.savefig("img/plummer_sphere.png")
  plt.close()

  # histogram
  rs = np.sqrt(np.array([x[i]**2 + y[i]**2 + z[i]**2 for i in range(N)]))
  fdp_plummer = lambda r: 3 * b*b * r*r * (b**2 + r**2)**(-5/2)
  r_axis = np.linspace(0, rmax, 200)
  y_axis = fdp_plummer(r_axis)
  plt.figure(figsize=(6,3))
  plt.plot(r_axis, y_axis, c='black', label="PDF")
  plt.hist(rs, bins=50, color=color, density=True, label="Sample", edgecolor="#212121")
  plt.legend()
  plt.grid(True)
  plt.ylabel(r"PDF and ocurrence")
  plt.xlabel(r"radius $r$")
  plt.title(r"Plummer radius sample ($N=10^4$, $b=1/2$, cut-off at $r=10r_h$)")
  plt.tight_layout()
  plt.savefig("img/plummer_density.png")
  
  # velocities
  g = lambda q: q*q*(1.-q*q)**3.5
  pot = lambda q, b: -G * M / np.sqrt(q**2 + b**2)
  v_esc = np.sqrt(- 2.0 * pot(rs, b))
  vs = np.zeros(N)
  qs = np.zeros(N)
  
  for i in range(N):
    # von neumann
    while True:
      qvn, yvn = np.random.uniform(size=2)
      yvn = 0.1 * yvn
      
      if (yvn < qvn*qvn*(1 - qvn**2)**3.5):
        break
    
    vs[i] = qvn * v_esc[i]
    qs[i] = qvn
  
  qs_axis = np.linspace(0, 1.0, 200)
  gs_axis = g(qs_axis)
  gs_axis = gs_axis * 512. / (7. * np.pi)
    
  plt.figure(figsize=(6,3))
  plt.plot(qs_axis, gs_axis, c='black', label="PDF")
  plt.hist(qs, bins=50, density=True, color=color, label="Sample", edgecolor="#212121")
  plt.legend()
  plt.grid(True)
  plt.xlabel(r"$qs$")
  plt.ylabel(r"PDF and ocurrence")
  plt.title("Isotropic Plummer velocities sample \n" + r"($N=10^4$, $b=1/2$, cut-off at $r=10r_h$)")
  plt.tight_layout()
  plt.savefig("img/plummer_isotropic_velocities.png")
  
if __name__ == "__main__":
  os.makedirs("img/", exist_ok=True)
  example_homogeneous_sphere()
  example_plummer_sphere()