import numpy as np
import matplotlib.pyplot as plt

def vertex_list (cx, cy, half_side):
  vertex_1 = [cx - half_side, cy + half_side]
  vertex_2 = [cx + half_side, cy + half_side]
  vertex_3 = [cx + half_side, cy - half_side]
  vertex_4 = [cx - half_side, cy - half_side]
  
  return [vertex_1, vertex_2, vertex_3, vertex_4, vertex_1]

def plot_boundaries_node (cx, cy, half_side, ax=False):
  xs, ys = list(zip(*vertex_list(cx, cy, half_side)))
  if not ax: plt.plot(xs, ys, c='black')
  else: ax.plot(xs, ys, c='black')

def fill_square_node (cx, cy, half_side, color:str='gray'):
  xs, ys = list(zip(*vertex_list(cx, cy, half_side)))
  if not ax: plt.fill(xs, ys, c=color, zorder=-1)
  else: ax.fill(xs, ys, c=color, zorder=-1)

######################################################

file = "fortran/out/quadtree_test.txt"

with open(file, 'r') as f:
  # fixed part
  N = int(f.readline())
  m = np.fromstring(f.readline(), sep=" ")
  x = np.fromstring(f.readline(), sep=" ")
  y = np.fromstring(f.readline(), sep=" ")

  # root information
  root_x, root_y, root_hs = np.fromstring(f.readline(), sep=" ")

  # tree information
  lines = []
  for line in f:
    line = line.strip()
    if line: lines.append(list(map(float, line.split())))

logs = np.array(lines)

fig = plt.figure(figsize=(5,5))
plt.scatter(x, y, c='black')

# root
plot_boundaries_node(root_x, root_y, root_hs)

# tree
for log in logs:
  subnode_depth, cx, cy = log
  halfsize_node = root_hs * 2**(-subnode_depth)
  plot_boundaries_node(cx, cy, halfsize_node)

plt.savefig("img/f_quadtree_example.pdf")
plt.title("Example of a quadtree generated with Fortran")
plt.savefig("img/f_quadtree_example.png")
