import matplotlib.pyplot as plt
import numpy as np
from galpy import potential, df

def plummer (N, b=1.0):
  pot = potential.PlummerPotential(amp=1.0, b=b, normalize=False)
  dfunc = df.isotropicPlummerdf(pot=pot)
  orbits = dfunc.sample(n=N)
  return orbits

def plummer_truncated (N, trunc, b=1.0):
  pot = potential.PlummerPotential(amp=1.0, b=b, normalize=False)
  dfunc = df.isotropicPlummerdf(pot=pot)
  accepted = []

  while len(accepted) < N:
    orb = dfunc.sample(n=N)

    r = np.sqrt(orb.x()**2 + orb.y()**2 + orb.z()**2)

    accepted.extend(orb[r <= trunc])

  return accepted[:N]

if __name__ == '__main__':
    b = 0.5
    M = 1.0
    fdp_plummer = lambda r: 3 * b*b * r*r * (b**2 + r**2)**(-5/2)

    N = 1000
    ms = np.ones(N)/N
    
    # WITHOUT CUT-OFF
    orbits = plummer(N, b)
    xs = orbits.x()
    ys = orbits.y()
    zs = orbits.z()
    qs = np.column_stack((xs, ys, zs))
    rs = np.linalg.norm(qs, axis=1)
    rs.sort()
    
    plt.figure(figsize=(6,3))
    plt.title(r"Random generated $r$ for $N=10^3$ and $b=1/2$")
    plt.hist(rs, bins=N // 5, density=True, color="#bb8cd1", label="Random generate values (histogram)")
    plt.plot(rs, fdp_plummer(rs), c='black', label="PDF")
    plt.legend()
    plt.grid(True)
    plt.xlabel(r"$r$")
    plt.tight_layout()
    plt.savefig("plummer_radius_example_without_cutoff.png")
    
    # WITH CUT-OFF
    orbits = plummer_truncated(N, 10 * 1.3 * b, b)

    xs, ys, zs = np.zeros(N), np.zeros(N), np.zeros(N)
    for i, orb in enumerate(orbits):
        xs[i] = orb.x()
        ys[i] = orb.y()
        zs[i] = orb.z()

    qs = np.column_stack((xs, ys, zs))
    rs = np.linalg.norm(qs, axis=1)
    rs.sort()
    
    plt.figure(figsize=(6,3))
    plt.title(r"Random generated $r$ for $N=10^3$ and $b=1/2$ (cut-off $= 10 r_h$)")
    plt.hist(rs, bins=N // 5, density=True, color="#bb8cd1", label="Random generate values (histogram)")
    plt.plot(rs, fdp_plummer(rs), c='black', label="PDF")
    plt.legend()
    plt.grid(True)
    plt.xlabel(r"$r$")
    plt.tight_layout()
    plt.savefig("plummer_radius_example_cutoff.png")
    
    