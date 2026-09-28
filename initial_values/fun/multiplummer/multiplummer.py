import sys
sys.path.append("../../")

import api as iv
import numpy as np
import matplotlib.pyplot as plt

def sample (N:int, parameters:list, profile:str, size_qs:float=-1.0, size_ps:float=-1.0, td:bool=False):
  N = np.array(N)
  m = np.zeros(N, order='F')
  qs = np.zeros((3,N), order='F')
  ps = np.zeros((3,N), order='F')
  iv.initial_values_mod.api(size_qs, size_ps, True, td, profile, parameters, m, qs, ps)
  return m, qs, ps

sets = [
  { "N": 300, "b": 0.5 , "qcm": [-5.0,+0.0,+0.0], "size": 2.0},
  { "N": 300, "b": 0.7 , "qcm": [+3.0,+2.0,+0.0], "size": 2.0},
  { "N": 300,  "b": 0.2 , "qcm": [+0.0,-3.0,+0.0], "size": 2.0},
]

Ntotal = sum(preset["N"] for preset in sets)
Mtot = 1.0
m = np.zeros(Ntotal)
qs = np.zeros((3, Ntotal))
ps = np.zeros((3, Ntotal))

plt.figure(figsize=(6,4))
plt.suptitle(r"Three Plummer models with 300 bodies each ($b=0.5, 0.7, 0.2$)")
plt.xlabel(r"$x$", rotation=0)
plt.ylabel(r"$y$", rotation=0)

idx = 0
for i, preset in enumerate(sets):
  m_local, qs_local, ps_local = sample(preset["N"], [1.0, preset["b"]], "plummer", size_qs=preset["size"])
  
  Mtot_local = sum(m_local)
  qs_local[0,:] = qs_local[0,:] - m_local @ qs_local[0,:] / Mtot_local
  qs_local[1,:] = qs_local[1,:] - m_local @ qs_local[1,:] / Mtot_local
  qs_local[2,:] = qs_local[2,:] - m_local @ qs_local[2,:] / Mtot_local

  qs_local[0,:] = qs_local[0,:] + preset["qcm"][0]
  qs_local[1,:] = qs_local[1,:] + preset["qcm"][1]
  qs_local[2,:] = qs_local[2,:] + preset["qcm"][2]

  next_idx = idx + preset["N"]
  m[idx:next_idx] = m_local
  qs[:,idx:next_idx] = qs_local
  ps[:,idx:next_idx] = ps_local

  idx = next_idx

plt.grid(True)

qs[0,:] = qs[0,:] - m @ qs[0,:] / sum(m)
qs[1,:] = qs[1,:] - m @ qs[1,:] / sum(m)
qs[2,:] = qs[2,:] - m @ qs[2,:] / sum(m)

print("qcm_x =", m @ qs[0,:] / sum(m))
print("qcm_y =", m @ qs[1,:] / sum(m))
print("qcm_z =", m @ qs[2,:] / sum(m))
print("Mtot  =", sum(m))

plt.scatter(qs[0,:], qs[1,:], c=qs[2,:], cmap='viridis', s=5)
plt.text(-4.0,  2.4, r"b=0.5", fontsize=11)
plt.text( 1.5,  4.0, r"b=0.7", fontsize=11)
plt.text( 2.2, -3.0, r"b=0.2", fontsize=11)

plt.axis("equal")
plt.xlim(-6, 6)
plt.ylim(-5, 5)
cbar = plt.colorbar()
cbar.set_label(r"$z$", rotation=0)
plt.tight_layout()
plt.savefig("multiplummer.png")
plt.close()


import json

json_obj = {
  "N": Ntotal,
  "valores_iniciais": {
    "massas": m.tolist(),
    "posicoes": [qs[:,i].tolist() for i in range(Ntotal)],
    "momentos": [ps[:,i].tolist() for i in range(Ntotal)]
  }
}

with open("multiplummer_vi.json", "w") as arq:
  json.dump(json_obj, arq)