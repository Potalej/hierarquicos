"""
Implementacao de um no (node). Contem informacoes
de centro, meio lado, massa, centro de massas, 
filhos e profundidade.
"""
import numpy as np

class No:
  __slots__ = (
    "cx", "cy",  # centro (geometrico)
    "meio_lado", # metade do tamanho do lado
    "massa",
    "qcm_x", "qcm_y", # centro de massas
    "particula", # -1 se for galho, inteiro se for folha
    "profundidade",
    "filho_NO", # noroeste
    "filho_NE", # nordeste
    "filho_SO", # sudoeste
    "filho_SE", # sudeste
    "eh_folha"
  )
  
  def __init__ (self, cx:float, cy:float, lado:float, profundidade:int=0):
    self.cx, self.cy = cx, cy
    self.meio_lado = lado / 2
    
    self.massa = 0.0
    self.qcm_x, self.qcm_y = 0.0, 0.0
    
    self.particula = -1
    self.profundidade = profundidade
    
    self.filho_NO = None
    self.filho_NE = None
    self.filho_SO = None
    self.filho_SE = None
    
    self.eh_folha = True # comeca sempre sendo folha
    
  @property
  def filhos (self):
    """
    Retorna os filhos que tiverem sido alocados.
    """
    if self.filho_NO is not None: yield self.filho_NO
    if self.filho_NE is not None: yield self.filho_NE
    if self.filho_SO is not None: yield self.filho_SO
    if self.filho_SE is not None: yield self.filho_SE
    
  def filho (self, ind:int):
    """
    Retorna o filho do no correspondente ao indice, mesmo
    que nao tenha sido alocado (None).
    """
    if   ind == 0: return self.filho_NO
    elif ind == 1: return self.filho_NE
    elif ind == 2: return self.filho_SO
    elif ind == 3: return self.filho_SE
  
def no_tamanho_centro (x:np.array, y:np.array)->float:
  """
  A partir dos vetores de posicao, determina o centro e o tamanho do
  lado que o no deve ter.
  """
  xmin, xmax = np.min(x), np.max(x)
  ymin, ymax = np.min(y), np.max(y)
  tamanho = max(xmax - xmin, ymax - ymin)
  cx = 0.5 * (xmin + xmax)
  cy = 0.5 * (ymin + ymax)
  return cx, cy, tamanho