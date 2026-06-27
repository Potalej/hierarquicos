"""
This file just generates a visualization to the Barnes [1] suggestion to use 3L-bit keys
to construct the octrees, and the Warren & Salmon [2] suggestion to use this plus the
Morton ordering to optmize the access to memory.

References
[1] BARNES, Joshua E. An efficient N-body algorithm for a fine-grain parallel computer. 1986
[2] WARREN, Michael S; SALMON, John K. A Parallel Hashed Oct-Tree N-Body Algorithm. 1993
"""
import numpy as np
import matplotlib.pyplot as plt

def square (x:float, y:float, side:float, ax=None):
  h = side/2
  xs = [x-h, x+h, x+h, x-h, x-h]
  ys = [y+h, y+h, y-h, y-h, y+h]
  if ax is not None: ax.plot(xs,ys,c='black')
  else:              plt.plot(xs,ys,c='black')

def plot_squares (max_depth:float, q0:list=[], ax=None):
  if len(q0) == 0: squares = [[0.5,0.5,1.0]]
  else:            squares = [q0]
  
  L_min = 1.0 / 2**max_depth

  while len(squares) > 0:
    x,y,L = squares[0]
    square(x,y,L,ax)
    del(squares[0])

    if L > L_min:
      h = L / 2.0
      squares.append([x-h/2,y-h/2,h])
      squares.append([x-h/2,y+h/2,h])
      squares.append([x+h/2,y-h/2,h])
      squares.append([x+h/2,y+h/2,h])

def morton_key (x:float, y:float, L:float)->int:
    # scale is 2**L
    scale = 1 << L
    
    # [0, 2*L]
    ix = int(min(np.floor(x*scale), scale-1))
    iy = int(min(np.floor(y*scale), scale-1))

    key = 0
    for i in range(L):
        # move bits to right
        b = L - 1 - i
        xb = ix >> b
        yb = iy >> b

        # get the last bit
        xb = xb & 1
        yb = yb & 1

        key = key << 2
        key = key | (xb << 1)
        key = key | yb
    
    return key

##########################################################################
if __name__ == '__main__':
    # max depth    
    L = 2

    fig, axs = plt.subplots(1, 2, figsize=(9,4))
    
    plot_squares(L, [], axs[0])
    plot_squares(L, [], axs[1])

    side = 1.0 / 2**L
    bits = []
    morton_points = []
    barnes_point = [0 for i in range(4**L)]
    for i in range(4**L):
        bit_string = f"{i:0{2*L}b}"

        ixs = bit_string[0::2]
        iys = bit_string[1::2]
        ix = int(ixs,2)
        iy = int(iys,2)

        x = (ix + 0.5) / 2**L
        y = (iy + 0.5) / 2**L

        # barnes num
        barnes_num = ""
        for i in range(L-1,-1,-1):
            barnes_num = barnes_num + iys[i]
            barnes_num = barnes_num + ixs[i]
        barnes_point[int(barnes_num,2)] = [x,y]

        # morton num
        morton_num = morton_key(x,y,L)
        morton_num = f"{morton_num:0{2*L}b}"
        morton_points.append([x,y])

        axs[0].text(x, y, rf"$({barnes_num})_2$", ha='center', va='center', zorder=3, fontsize=10)
        axs[1].text(x, y, rf"$({morton_num})_2$", ha='center', va='center', zorder=3, fontsize=10)
        axs[0].text(x, y-0.25*side, rf"${int(barnes_num,2)}$", ha='center', va='center', zorder=3, fontsize=10)
        axs[1].text(x, y-0.25*side, rf"${int(morton_num,2)}$", ha='center', va='center', zorder=3, fontsize=10)

    axs[0].set_title("Barnes method")
    axs[1].set_title("Morton ordering")
    axs[0].set_aspect('equal')
    axs[1].set_aspect('equal')

    x,y = list(zip(*barnes_point))
    axs[0].plot(x, y, alpha=0.8, zorder=-1, c="#bb8cd1")
    x,y = list(zip(*morton_points))
    axs[1].plot(x, y, alpha=0.8, zorder=-1, c="#bb8cd1")

    plt.savefig("img/barnes_morton_keys.png")
    plt.close()