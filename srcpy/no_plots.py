import matplotlib.pyplot as plt
from srcpy.no import No

def vertices (no: No):
  """
  Obtem os vertices de um no.
  """
  vertice_1 = [no.cx - no.meio_lado, no.cy + no.meio_lado]
  vertice_2 = [no.cx + no.meio_lado, no.cy + no.meio_lado]
  vertice_3 = [no.cx + no.meio_lado, no.cy - no.meio_lado]
  vertice_4 = [no.cx - no.meio_lado, no.cy - no.meio_lado]
  
  vertices = [vertice_1, vertice_2, vertice_3, vertice_4, vertice_1]
  return vertices

def plotar_bordas_no (no: No, ax=False):
  """
  Plota o quadrado relativo ao no.
  """
  xs, ys = list(zip(*vertices(no)))
  if not ax: plt.plot(xs, ys, c='black')
  else: ax.plot(xs, ys, c='black')

def preencher_quadrado_no (no: No, ax=False, cor:str='gray'):
  """
  Plota o interior do quadrado relativo ao no.
  """
  xs, ys = list(zip(*vertices(no)))
  if not ax: plt.fill(xs, ys, c=cor, zorder=-1)
  else: ax.fill(xs, ys, c=cor, zorder=-1)
