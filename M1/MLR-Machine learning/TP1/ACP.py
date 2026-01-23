#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACP
---
A fournir:
----------
  X : Ndarray, shape (n,p)
    Tableau individus/variables
Et optionnellement pour les graphiques, les tableaux:
  labels : ndarray of strings, shape (n,)
    Etiquette de chaque individu, sera realisee par une couleur

  indiv : ndarray of strings, shape (n,)
    Les noms des individus, places sur le graphique.
    P.ex. indiv est la marque de voiture et label sa categorie (luxe,...)

  varbs : ndarray of strings, shape (p,)
    Noms de variables, pour le cercle des correlations
"""
import numpy as np
import matplotlib.pyplot as plt
#from importlib import reload; import ACP;reload(ACP)

# Par defaut
labels=None
indiv=[]
varbs=None

# Lecture des donnees de poterie
pot=np.loadtxt("Poteries.csv",skiprows=1,usecols=range(0,10))
labels=pot[:,9] # Les numeros de four
labels=np.array([str(a)[0] for a in labels])
X=pot[:,np.arange(0,9)]
varbs=np.loadtxt("Poteries.csv",dtype='str',usecols=range(0,9))[0,:]
#
def stdise(X):
  """Routine de standardisation
     On pourrrait utiliser scikitlearn: Exemple:
       from sklearn.preprocessing import StandardScaler
       X=np.arange(12).reshape((4,3))
       print(StandardScaler().fit_transform(X))
  """
  Xs=X.astype(float)
  mk=np.mean(Xs,axis=0)
  # Calcul de l'ecart-type avec max pour eviter une division par 0
  sk=np.maximum(np.std(Xs,axis=0),10*np.finfo(float).eps)
  Xs=(X-mk)/sk
  return Xs, mk, sk
# SVD. Axes Composantes
# Apres standardisation les colonnes sont de norme "nb de ligne" et non 1,
# on corrige cela, pour avoir de meilleures echelles.
Xs, mk,sk=stdise(X)

(U,D,VT) = np.linalg.svd(Xs,full_matrices=False)
V=VT.T
# Premieres composantes principales
C1 = D[0]*U[:,0]
C2 = D[1]*U[:,1]
C3 = D[2]*U[:,2]
# Axes principaux modifies pour le cercle des correlations
A1 = D[0]*V[:,0]/np.sqrt(np.shape(X)[0])
A2 = D[1]*V[:,1]/np.sqrt(np.shape(X)[0])
# Graphiques
plt.close('all')
plt.figure()
plt.title('Representation des individus dans le plan (C1,C2)')
if labels is None:
  plt.scatter(C1,C2)
else:
  vlab=np.unique(labels)
  lv=len(vlab)
  #cols=['C0','C1','C2','C3','k'] # un choix de couleurs
  #cols=plt.cm.nipy_spectral(np.arange(lv)/lv) # un choix de couleurs
  for i,vl in enumerate(vlab):
    l=labels==vl
    plt.scatter(C1[l],C2[l],s=47,label=vl)#,color=cols[i])
  plt.legend(title="Four")
for i,nm in enumerate(indiv): plt.text(C1[i],C2[i],nm)
plt.xlabel('C1')
plt.ylabel('C2')

# Inerties
plt.figure()
plt.bar(np.arange(np.shape(D)[0])+1,100*D**2/sum(D**2))
plt.title('Inerties en %')

# Cercle des correlations
if not varbs is None:
  plt.figure()
  plt.title('Cercle des correlations')
  Z = np.linspace(-np.pi, np.pi, 256,endpoint=True)
  C,S = np.cos(Z), np.sin(Z)
  plt.plot(C,S,c='black',lw=.7)
  plt.axvline(c='black',ls='dashed',lw=1)
  plt.axhline(c='black',ls='dashed',lw=1)
  for i, txt in enumerate(varbs):
    plt.arrow(0,0,A1[i],A2[i], length_includes_head=True,
            head_width=0.025, head_length=.05)
    plt.annotate(txt, (A1[i]+.01,A2[i]+.01),fontsize=12)
  plt.xlabel('C1')
  plt.ylabel('C2')


