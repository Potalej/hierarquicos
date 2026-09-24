import matplotlib.pyplot as plt
import numpy as np

color = "#bb8cd1"

def generate_xyz_from_radius (r, td:bool=False):
  N = len(r)
  phi = 2. * np.pi * np.random.uniform(size=N)
  if td:
    costheta = 0.
    sintheta = 1.
  else:
    x3 = np.random.uniform(size=N)
    costheta = 2. * x3 - 1.
    sintheta = (1. - costheta**2)**0.5

  x = r * sintheta * np.cos(phi)
  y = r * sintheta * np.sin(phi)
  z = r * costheta
  return x, y, z


"""
HOMOGENEOUS SPHERE PROFILE
"""
def homogeneous_sphere_profile (N:int, rho0:float, R:float):
  X1 = np.random.uniform(size=N)
  rs = np.zeros(N)

  i = 0
  while i < N:
    x1 = np.random.uniform()
    r = R * x1 ** (1./3.)
    rs[i] = r
    i += 1

  return rs

def example_homogeneous_sphere ():
  M = 1.0
  G = 1.0

  N = 10000
  rs = homogeneous_sphere_profile(N, rho0=0.1, R=3.0)
  x, y, z = generate_xyz_from_radius(rs)

  plt.figure(figsize=(5,4))
  plt.xlim(-5., 5.)
  plt.ylim(-5., 5.)
  plt.title(r"Homogeneous Sphere Profile $(R=3, N=10^4)$")
  plt.scatter(x, y, c=z, cmap='viridis', s=5)
  plt.axis("equal")
  cbar = plt.colorbar()
  cbar.set_label(r'$z$', rotation=0)
  plt.tight_layout()
  plt.savefig("homogeneous_sphere.png")

"""
PLUMMER SPHERE PROFILE
"""
def plummer_sphere_profile (N:int, b:float, rmax:float=0.0):
  rs = np.zeros(N)
  i = 0
  while i < N:
    x1 = np.random.uniform()
    r = x1**(1./3.) * b / np.sqrt(1. - x1**(2./3.))
    if rmax == 0: 
      rs[i] = r
      i += 1
    else:
      if r <= rmax:
        rs[i] = r
        i += 1
  return rs

def example_plummer_sphere ():
  M = 1.0
  G = 1.0
  
  b = 0.5
  k = 0.5
  rmax = 10. * (k**(1./3.)/np.sqrt(1.-k**(2./3.))) * b

  N = 10000
  rs = plummer_sphere_profile(N, b=0.5, rmax=rmax)
  x, y, z = generate_xyz_from_radius(rs)
  
  plt.figure(figsize=(5,4))
  plt.title(r"Plummer Sphere Profile $(b=0.5, N=10^4)$")
  plt.scatter(x, y, c=z, cmap='viridis', s=5)
  plt.axis("equal")
  plt.xlim(-7., 7.)
  plt.ylim(-7., 7.)
  cbar = plt.colorbar()
  cbar.set_label(r'$z$', rotation=0)
  plt.tight_layout()
  plt.savefig("plummer_sphere.png")
  plt.close()

  # histogram
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
  plt.savefig("plummer_density.png")
  
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
  plt.savefig("plummer_isotropic_velocities.png")
  

if __name__ == "__main__":
  # example_homogeneous_sphere()
  example_plummer_sphere()