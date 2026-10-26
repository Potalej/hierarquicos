import matplotlib.pyplot as plt
import numpy as np
color = "#bb8cd1"

plt.figure(figsize=(4,4))

qs = np.array([
    [0.85, 0.8],
    [0.23, 0.4],
    [0.63, 0.12],
])
qs = qs * 0.8

qcm = [0,0]

for x, y in qs:
    qcm[0] += x
    qcm[1] += y

M = len(qs)
qcm[0] = qcm[0]/M
qcm[1] = qcm[1]/M

p = 1.1

for i in range(len(qs)):
    qs[i][0] -= qcm[0] - 1
    qs[i][1] -= qcm[1] - p

qcm = [0,0]
for x, y in qs:
    qcm[0] += x
    qcm[1] += y

qcm[0] = qcm[0]/M
qcm[1] = qcm[1]/M

lado = 0.85
centro = [1.05, p+1e-1]
plt.plot([centro[0]-lado/2, centro[0]+lado/2],[centro[1]-lado/2, centro[1]-lado/2],c='black')
plt.plot([centro[0]-lado/2, centro[0]+lado/2],[centro[1]+lado/2, centro[1]+lado/2],c='black')
plt.plot([centro[0]-lado/2, centro[0]-lado/2],[centro[1]-lado/2, centro[1]+lado/2],c='black')
plt.plot([centro[0]+lado/2, centro[0]+lado/2],[centro[1]-lado/2, centro[1]+lado/2],c='black')
plt.text(centro[0] + 1.1*lado/2, 0.96*centro[1], r"$s$", fontsize=20)

for x, y in qs:
    plt.scatter(x,y,s=100,c='black')

plt.scatter(*qcm, s=120, c="red", marker='+', zorder=3)

# for x, y in qs:
#     plt.plot([x, qcm[0]], [y, qcm[1]], c='black', linewidth=0.7)

part = [0.0, 0.0]
plt.scatter(*part, c='black', s=100)

theta = 0.7
phi = np.atan(p)
m = np.tan(phi + theta/2)
L = np.sqrt(1**2 + p**2) / np.cos(theta / 2)
x = L / np.sqrt(1 + m**2)
t = np.linspace(0,x,50)
plt.plot(t, m*t, c='black', linestyle='--')
x1 = x
L1 = m*x

m = np.tan(phi - theta/2)
L = np.sqrt(1**2 + p**2) / np.cos(theta / 2)
x = L / np.sqrt(1 + m**2)
t = np.linspace(0,x,50)
plt.plot(t, m*t, c='black', linestyle='--')
x2 = x
L2 = m*x

plt.plot([x1,x2],[L1,L2], c='black', linestyle='--')

t = np.linspace(phi-theta/2,phi+theta/2,50)
r = 0.2
plt.plot(r*np.cos(t), r*np.sin(t), c='black')
plt.text(0.30, 0.23, r"$\theta$", fontsize=20)

plt.plot([0, 1], [0, p], c='black', linestyle='--')
plt.text(0.47, 0.70, r"$d$", fontsize=20)

plt.xlim(-0.2,1.7)
plt.ylim(-0.2,1.7)
plt.axis('off')
plt.tight_layout()
plt.savefig("bh_cell_interaction.png")