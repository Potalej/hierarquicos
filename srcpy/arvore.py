from srcpy import no as nos
from srcpy import no_plots
import numpy as np
import matplotlib.pyplot as plt

class Arvore:
  """
  Implementacao do quadtree e do metodo de Barnes-Hut para calculo de forcas de problemas de N-corpos gravitacionais.
  
  Uma vez instanciado,
  ```
  arvore = Arvore(m, x, p)
  ```
  ja monta a arvore, criando as devidas folhas e galhos, e fica pronto para uso.
  """
  
  max_profundidade = 32
  amplificador_lado = 1.2

  """
  Inicializacao da arvore.
  Salva as variaveis principais e ja comeca a criacao da arvore.
  """
  def __init__ (self, m:np.array, x:np.array, y:np.array, plotar:bool=False)->None:
    # variaveis de estado
    self.m, self.x, self.y = m, x, y
    self.N = len(self.m)
    
    # se quiser visualizar
    self.plotar = plotar
    
    # iniciando a raiz
    cx, cy, tamanho = nos.no_tamanho_centro(self.x, self.y)
    self.raiz = nos.No(cx, cy, self.amplificador_lado * tamanho)
    
    # adicionando as particulas
    for i in range(self.N):
      self.__adicionar(self.raiz, i)
  
  """ 
  #-----------------------------#
  # Estrutura interna da arvore # 
  #-----------------------------#
  """ 
  def __alocar_filho (self, no:nos.No, ind:int)->None:
    """
    Aloca o filho no lugar devido.
    """
    h = no.meio_lado
    h_meio = h / 2
    cx, cy = no.cx, no.cy
    
    # noroeste
    if   ind == 0: no.filho_NO = nos.No(cx - h_meio, cy + h_meio, h, no.profundidade + 1)
    # nordeste
    elif ind == 1: no.filho_NE = nos.No(cx + h_meio, cy + h_meio, h, no.profundidade + 1)
    # sudoeste
    elif ind == 2: no.filho_SO = nos.No(cx - h_meio, cy - h_meio, h, no.profundidade + 1)
    # sudeste
    elif ind == 3: no.filho_SE = nos.No(cx + h_meio, cy - h_meio, h, no.profundidade + 1)
   
  def __indice_filho (self, no: nos.No, p:int)->int:
    """
    Identifica o quadrante em que a particula se localiza na arvore,
    retornando o indice do filho correspondente:
      0 -> noroeste
      1 -> nordeste
      2 -> sudoeste
      3 -> sudeste
    """
    cima    = self.y[p] > no.cy
    direita = self.x[p] > no.cx
    
    if cima: return 1 if direita else 0
    else:    return 3 if direita else 2
  
  def __adicionar_ao_filho (self, no:nos.No, p:int)->None:
    """
    Adiciona uma particula a um no.
    """
    ind = self.__indice_filho(no, p)
    filho = no.filho(ind)
    
    # se o filho nao estiver alocado, aloca
    if filho is None:
      self.__alocar_filho(no, ind)
      # filho = no.filho(ind)
    
    # agora adiciona
    self.__adicionar(no.filho(ind), p)
    
    # se quiser plotar
    if self.plotar: 
      no_plots.plotar_bordas_no(no.filho(ind))
  
  def __adicionar (self, no:nos.No, part_ind:int)->None:
    """
    Adiciona uma particula a um no.
    """
    part_m = self.m[part_ind]
    part_x = self.x[part_ind]
    part_y = self.y[part_ind]
    
    # um no vazio vira particula
    if no.particula == -1 and no.eh_folha:
      no.particula = part_ind
      
      # adiciona diretamente o centro de massas
      no.massa = part_m
      no.qcm_x = part_x
      no.qcm_y = part_y
      
      return
    
    # atualiza o centro de massas se nao for vazio
    massa_antiga = no.massa
    no.massa += part_m
    no.qcm_x = (no.qcm_x * massa_antiga + part_x * part_m) / no.massa
    no.qcm_y = (no.qcm_y * massa_antiga + part_y * part_m) / no.massa
    
    # se nao estiver vazio, subdivide
    if no.eh_folha:
      # nao pode estourar o limite de profundidade
      if no.profundidade >= self.max_profundidade:
        print("PROBLEMA SERIO !!!")
        return
      
      # nesse caso, deixa de ser folha
      no.eh_folha = False
      
      # adiciona a particula antiga como particula em si
      part_antiga = no.particula
      no.particula = -1
      self.__adicionar_ao_filho(no, part_antiga)
      
      
    # agora insere nova particula
    self.__adicionar_ao_filho(no, part_ind)
  
  """ 
  #-----------------------------#
  #      Calculo de forcas      # 
  #-----------------------------#
  """
  def calcular_forcas (self, p:int, eps:float=0.0, theta:float=0.0, G:float=1.0)->tuple:
    """
    Calcula as forcas sobre uma particula `p` utilizando a arvore e os parametros de amortecimento `eps` e de escolha (BH) `theta`
    """
    return self.__calcular_forcas(self.raiz, p, eps, theta, G)
  
  def __calcular_forcas (self, no:nos.No, p:int, eps:float, theta:float, G:float)->tuple:
    """
    Calcula as forcas de fato, conforme o potencial newtoniano.
    """
    # se nao tiver massa, nao calcula nada
    if no.massa == 0: return 0.0, 0.0
    
    # se for a propria particula
    if no.eh_folha and no.particula == p: return 0.0, 0.0
    
    # noutro caso, verifica a condicao de Barnes-Hut
    dx = no.qcm_x - self.x[p]
    dy = no.qcm_y - self.y[p]
    
    dist2 = dx*dx + dy*dy
    dist = np.sqrt(dist2)
    l = 2.0 * no.meio_lado
    
    # se for uma folha ou valer o criterio
    if no.eh_folha or (l / dist < theta):
      rab = np.sqrt(dist2 + eps*eps)
      f = G * self.m[p] * no.massa / (rab ** 3)
      
      # caso queira plotar
      if self.plotar: 
        no_plots.preencher_quadrado_no(no)
        no_plots.plotar_bordas_no(no)
        plt.scatter(no.qcm_x, no.qcm_y, marker='+', c='red')
        plt.plot([self.x[p], no.qcm_x], [self.y[p], no.qcm_y], c='black')
      
      return f * dx, f * dy
    
    # se nao valer o criterio nem for folha, vai nos filhos
    fx, fy = 0.0, 0.0
    for filho in no.filhos:
      cfx, cfy = self.__calcular_forcas(filho, p, eps, theta, G)
      fx += cfx
      fy += cfy
    return fx, fy