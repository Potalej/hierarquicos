"""
Para estimar o custo da montagem da quadtree.
"""
import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from time import time
from srcpy.arvore import Arvore

tempos = []
tempos_media = []
N_min, N_max, qnt_em_qnt = 30, 1000, 20
qntd_testes_N = 40

# testando para cada N
for N in range(N_min, N_max + 1, qnt_em_qnt):
  # para cada N, testa uma quantidade de vezes
  tempos_N = []
  for _ in range(qntd_testes_N):
    # sorteia os valores
    massas = np.random.random(N)
    qs = 5 * (2 * np.random.random((N,2)) - 1)
    x, y = np.array(list(zip(*qs)))
    
    # medindo o tempo de criacao da arvore
    tempo_inicial = time()
    arvore = Arvore(massas, x, y)
    tempo_total = time() - tempo_inicial
    
    # salvando
    tempos.append([N, tempo_total])
    tempos_N.append(tempo_total)
  
  # salva a media
  tempos_media.append([N, np.mean(tempos_N)])

# plotando    
fig, ax = plt.subplots(figsize=(8,3))
Ns, ts = list(zip(*tempos_media))
ax.scatter(Ns, ts, c='black')

# para complementar, tambem traco um N log N
f = lambda t, a: a * t * np.log(t)
Ns, ts = list(zip(*tempos))
popt, pcov = curve_fit(f, Ns, ts)
Ns = np.linspace(N_min, N_max, 100)
ax.plot(Ns, f(Ns, *popt), c='red', linestyle='--')
ax.text(Ns[len(Ns) // 2] + 10, 0.8*f(Ns[len(Ns) // 2], *popt), rf"$m={popt[0]:.2e}$")
ax.grid(True)
ax.set_ylabel("Tempo (s)")
ax.set_xlabel(r"$N$")
plt.title("Tempo de criação da quad-tree por $N$.")
plt.tight_layout()
plt.savefig("img/py_custo_arvore.png")