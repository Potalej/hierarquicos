"""
Para estimar o custo do calculo das forcas via Barnes-Hut
"""
import numpy as np, matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from time import time
from fortran import api

thetas = [0.0, 0.1, 0.3, 0.6, 0.9, 1.0]
N_min, N_max, qnt_em_qnt = 10, 1000, 10
qntd_testes_N = 20
eps = 1e-2
G = 1.0

# listas de informacoes
tempos = [[] for theta in thetas]
erros  = [[] for theta in thetas]
temposN2 = []

# testando para cada N
for N in range(N_min, N_max + 1, qnt_em_qnt):
  temposN2_N = []
  tempos_N = [[] for theta in thetas]
  erros_N  = [[] for theta in thetas]
  print(f"N: {N}")
  
  # para cada N, testa uma quantidade de vezes
  for _ in range(qntd_testes_N):
    massas = np.random.random(N)
    qs = 5. * (2. * np.random.random((N,2)) - 1.)
    x, y = np.array(list(zip(*qs)))
    
    # para o calculo direto
    forcas0 = np.zeros((N,2))
    # iniciando o timer
    tempo_0 = time()
    # calcula
    forcas0 = api.api_mod.compute_forces_direct(massas, x, y, G, eps)
    # encerra o timer
    tempo = time() - tempo_0
    temposN2_N.append(tempo)
    
    # para o calculo via arvore
    forcas = np.zeros((N,2))
    
    for i, theta in enumerate(thetas):
      forcas.fill(0.0)
      
      # iniciando o timer
      tempo_0 = time()

      # calculando as forcas
      forcas = api.api_mod.compute_forces(massas, x, y, theta, G, eps)

      # encerra o timer
      tempo = time() - tempo_0
      tempos_N[i].append(tempo)
      
      # agora calcula o erro
      erro = np.max(np.abs(forcas0 - forcas))
      erros_N[i].append(erro / np.max(np.abs(forcas0)))
  
  # salva as medias
  temposN2.append([N, np.mean(temposN2_N)])
  for i in range(len(thetas)):
    tempos[i].append([N, np.mean(tempos_N[i])])
    erros[i].append([N, np.mean(erros_N[i])])

# Analise do tempo de computacao
fig = plt.figure(figsize=(8,3))

x, y = np.array(list(zip(*temposN2)))
plt.scatter(x, y, label="Direto", c='black')
# esperado: O(N^2)
f = lambda t, a: a * t * t
popt, pcov = curve_fit(f, x, y)
plt.plot(x, f(x, *popt), c='black', linestyle='--')

for i, theta in enumerate(thetas):
  x, y = np.array(list(zip(*tempos[i])))
  if i == 0:
    plt.scatter(x, y, label=rf"$\theta={theta}$", c='gray')
    # esperado: O(N^2)
    f = lambda t, a: a * t * t
    popt, pcov = curve_fit(f, x, y)
    plt.plot(x, f(x, *popt), c='gray', linestyle='--')

  else:
    plt.scatter(x, y, label=rf"$\theta={theta}$")
    # esperado: O(N log N)
    f = lambda t, a: a * t * np.log(t)
    popt, pcov = curve_fit(f, x, y)
    plt.plot(x, f(x, *popt), linestyle='--')

plt.legend()
plt.title(r"Custo das forças com diferentes $\theta$")
plt.grid(True)
plt.savefig("img/f2py_custo_forcas_theta.png")
plt.close()


# Analise do erro relativo nas forcas
fig = plt.figure(figsize=(8,3))
for i, theta in enumerate(thetas):
  x, y = np.array(list(zip(*erros[i])))
  if i == 0: plt.scatter(x, y, label=rf"$\theta={theta}$", c='black')
  else: plt.scatter(x, y, label=rf"$\theta={theta}$")

plt.yscale('log')
plt.legend()
plt.title(r"Erros das forças com diferentes $\theta$")
plt.grid(True)
plt.ylim(1e-5, 1e-0)
plt.savefig("img/f2py_erro_forcas_theta.png")
plt.close()