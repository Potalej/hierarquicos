import matplotlib.pyplot as plt
import numpy as np
from galpy import potential, df
import ncorpos_utilidades as nut

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
    
    # VELOCITIES
    vxs = orbits.vx()
    vys = orbits.vy()
    vzs = orbits.vz()
    vs = np.linalg.norm(np.column_stack((vxs, vys, vzs)), axis=1)
    
    # evaluating the escape velocity to get the q's
    rs = np.linalg.norm(qs, axis=1)
    pot = potential.PlummerPotential(amp=1.0, b=b, normalize=False)
    v_esc = np.sqrt(-2.0 * pot(rs,0))
    qs = vs / v_esc
    
    g = lambda q: q*q*(1-q*q)**(3.5)
    axis_v = np.linspace(0, 1, 100)
    gs = g(axis_v)
    gs = gs / np.sum(gs)
    plt.plot(axis_v, gs, c='black')
    
    plt.figure(figsize=(6,3))
    plt.title(r"Random generated $v$ for $N=10^3$ and $b=1/2$")
    plt.hist(qs, bins=100, weights=np.ones(N)/N, color="#bb8cd1", label="Random generate values (histogram)")
    plt.plot(axis_v, gs, c='black', label="PDF")
    plt.legend()
    plt.grid(True)
    plt.xlabel(r"$v$")
    plt.tight_layout()
    plt.savefig("plummer_velocity_example.png")
    
    # lets verify some properties
    qs = np.column_stack((xs, ys, zs))
    vs = np.column_stack((vxs, vys, vzs))
    ps = np.array([ms[i] * vs[i] for i in range(N)])
    
    T = nut.energia_cinetica(ms, ps)
    V = nut.energia_potencial(ms, qs, 1.0, 0.0)
    E = T + V
    Q = - 2.0 * T / V
    
    V_expected = - 3.*np.pi*1.0*M**2 / (32. * b)
    T_expected = - V_expected / 2.
    Q_expected = 1.0
    E_expected = - T_expected
    
    print("Dynamical properties (||real - expected||):")
    print(f"Potential V: {abs(V - V_expected)}")
    print(f"Kinect T:    {abs(T - T_expected)}")
    print(f"Total E:     {abs(E - E_expected)}")
    print(f"Virial Q:    {abs(Q - 1.0)}")
    
    
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
    
    