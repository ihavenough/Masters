#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 12 08:28:18 2018

@author: valerie
"""

import numpy as np
from sklearn.cluster import  KMeans
from sklearn.metrics import confusion_matrix
import matplotlib.pyplot as plt

# Chargement des données ======================================================
dirdat=""
digits=np.loadtxt(dirdat+"digits_extrait_images.csv", delimiter=',',skiprows=1,usecols=range(1,785))
labels=np.loadtxt(dirdat+"digits_extrait_labels.csv", delimiter=',',skiprows=1,usecols=range(1,2))
labels=labels.astype(int)
#from sklearn.preprocessing import StandardScaler
#digits=StandardScaler().fit_transform(digits)
#
# ACP ICI
#
# Clustering avec la fonction KMeans du module scikit learn
# Le nombre de classes est choisi par défaut comme égal au nombre d'étiquettes
nclus=len(np.unique(labels))
nclus=len(np.unique(labels))+2
k_means = KMeans(init='k-means++', n_clusters=nclus, n_init=10)
k_means.fit(digits)
cl = k_means.labels_ # classes "prédites"
print('Premières classes:    ',cl[range(20)])
print('Premières étiquettes: ',labels[range(20)])

#  Images des centres de classe
plt.close('all')
fig=plt.figure(2)
plt.clf()
plt.suptitle("Centres de classe")
for k in range(k_means.n_clusters):
  plt.subplot(1,k_means.n_clusters,k+1)
  plt.imshow(np.reshape(k_means.cluster_centers_[k,:],[28,28]))
  plt.title("Cl. "+str(k))
  plt.axis('off')
plt.show()
# Calcul de l'étiquette majoritaire de chaque classe et du taux d'erreur
# On fabrique le tableau maj_lab qui a un nurmero de classe (p.ex. predit)
# renvoie l'etiquette correspondante.
maj_lab=np.arange(k_means.n_clusters) # Initialisation de tableau
for k in range(k_means.n_clusters):
  counts=np.unique(labels[cl==k],return_counts=True) # Nb d'occurences de chaque label
  imax=np.argmax(counts[1]) # Recherche du majoritaire dans k
  maj_lab[k]=counts[0][imax] # Son étiquette ya majoritaire dans k
print('\n')
print("Classe".ljust(23,'.')+" ",end='')
print(*range(k_means.n_clusters),end='')
print("\n"+"Etiquette majoritaire".ljust(23,'.')+" ",end='')
print(*(maj_lab))
err=sum(labels!=maj_lab[cl])/len(cl)
print("Taux de mal classés:",err.round(3))
# Matrice de confusion
# Rq: La matrice fournie par confusion_matrix est carrée:
#    On retire les lignes de zeros de conf_mat, dues au fait qu'il peut y avoir plus de classes que d'etiquettes
#    Et de meme avec les colonnes s'il y a moins de classes
conf_mat =  confusion_matrix(labels,maj_lab[cl])
from sklearn.metrics import confusion_matrix,ConfusionMatrixDisplay
plt.rcParams.update({'figure.figsize': (3,3),'font.size': 10})
ConfusionMatrixDisplay(conf_mat).plot(cmap='YlOrBr')
im = plt.gca().images[-1].colorbar.remove()
plt.rcdefaults()

# Représentation en barplot
def BarPlotMat(M):
# Fait un barplot pour chaque colonne de M.
# La couleur correspond à l'indice, la hauteur à la valeur
  I=M.shape[0]
  J=M.shape[1]
  ind = np.arange(J)
  haut = 0*M[0,:]
  for i in range(I):
    plt.bar(ind,M[i,:],bottom=haut,color=plt.cm.inferno(i/(I-1)))
    haut += M[i,:]
fig=plt.figure(3)
plt.clf()
conf_mat =  confusion_matrix(labels,cl)
conf_mat=conf_mat[np.sum(conf_mat,axis=1)>0,:]
conf_mat=conf_mat[:,np.sum(conf_mat,axis=0)>0]
BarPlotMat(conf_mat)
plt.xlabel('Classe')
plt.ylabel('Répartition des étiquettes')
plt.title('Répartition dans chaque classe')
plt.legend(['0','1','2','3','4',])
plt.show()
print("Matrice associée (conf_mat[etiq, classe]):\n")
print(conf_mat)
print("   Une ligne = un digit\n   Une colonne = une classe\n")


