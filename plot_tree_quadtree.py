from srcpy import arvore
import numpy as np
import matplotlib.pyplot as plt

N = 25

np.random.seed(1)

m = np.ones(N)
x = np.random.uniform(-1, 1, N)
y = np.random.uniform(-1, 1, N)

fig = plt.figure(figsize=(4,4))

arvre = arvore.Arvore(m, x, y, plotar=True)
plt.scatter(x,y,c='black',s=20)
plt.axis('off')
plt.tight_layout()
plt.savefig("img/quadtree.png")
plt.close()

def plotar_no (no, x, h, dx):
    pos = x - 1.5*dx
    
    for i, filho in enumerate(no.filhos):
        plt.scatter(pos, h, c='black')
        plt.plot([x, pos], [h-0.5, h], c='black')
        pos = pos + dx
    
    h = h + 0.5
    
    pos = x - 1.5*dx    
    for filho in no.filhos:
        if not filho.eh_folha:
            if h > 1.5:
                plotar_no(filho, pos, h, dx)
            else:
                plotar_no(filho, pos, h, dx / 4)
            pos = pos + dx

def plotar_arvore (tree):
    raiz = tree.raiz
    plt.figure(figsize=(8,3))
    plt.scatter(0,0,c='black')
    plotar_no(raiz, 0.0, 0.5, 1.0)
    plt.axis('off')
    plt.tight_layout()
    plt.savefig("img/tree.png")
    
plotar_arvore(arvre)