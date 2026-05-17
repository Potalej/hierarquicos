"""
Para visualizar o quadtree e o calculo das forcas
sobre um corpo.
"""
from srcpy.arvore import Arvore
import numpy as np, matplotlib.pyplot as plt

N = 4
massas = np.ones(N)
qs = 2. * np.random.random((N, 2)) - 1.
x, y = np.array(list(zip(*qs)))


###### Para visualizar o quadtree
fig, ax = plt.subplots(figsize=(6,6))
ax.scatter(x, y, c='red') # plotando as particulas
arvore = Arvore(massas, x, y, plotar=True)
plt.savefig("img/py_quadtree.png")


###### Para visualizar as forcas
fig = plt.figure(figsize=(6,6))
plt.scatter(x, y, c='black') # plotando as particulas
plt.scatter(x[0], y[0], c='blue') # particula observada
# arvore.calcular_forcas(0, theta=0.5)
fx0, fy0 = arvore.calcular_forcas(0, theta=0.0)
fx1, fy1 = arvore.calcular_forcas(1, theta=0.0)
fx2, fy2 = arvore.calcular_forcas(2, theta=0.0)
fx3, fy3 = arvore.calcular_forcas(3, theta=0.0)

print()
print(fx0, fy0)
print(fx1, fy1)
print(fx2, fy2)
print(fx3, fy3)

print()

forcas = np.zeros((N,2))
for p in range(len(x)):
  fx, fy = 0., 0.
  for pb in range(len(x)):
    if pb == p: continue
    dx = x[pb] - x[p]
    dy = y[pb] - y[p]
    dist2 = dx*dx + dy*dy
    dist = np.sqrt(dist2)
    f = 1.0 * massas[p] * massas[pb] / (dist2 * dist)
    fx += dx * f
    fy += dy * f
  forcas[p][0] = fx
  forcas[p][1] = fy

  print(forcas[p])

plt.savefig("img/py_forcas.png")