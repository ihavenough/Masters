# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from sklearn.cluster import  KMeans

#### Importation des donnees, noms de fromage et noms de variable
fromages=np.loadtxt("fromages.dat",skiprows=1,usecols=range(1,10))
nomfrom=np.loadtxt("fromages.dat",dtype='str',usecols=0)[1:]
nomvar=np.loadtxt("fromages.dat",dtype='str',max_rows=1,usecols=range(1,10))
print(nomfrom)
print(nomvar)

##### Exemple de programation de CAH #####
print("\n******* Classification Ascendante hiérarchique ******* \n")

# Calcul de l'arbre
M=linkage(fromages,method='ward',metric='euclidean')

# Tracé de l'arbre
seuil=40
print("----> !!! A vous de reflechir au seuil !!!\n")
plt.figure()
plt.title('CAH. Visualisation des classes au seuil de '+str(seuil))
d=dendrogram(M,labels=list(nomfrom),orientation='right',color_threshold=seuil)
print(np.round(M[:,2],2))

##### Récupération des groupes
groupes=fcluster(M,t=seuil,criterion='distance')
for k in range(1,np.max(groupes)+1):
    print('Classe '+str(k).ljust(3,' ')+': ', end='')
    print(*nomfrom[np.where(groupes==k)])

#### Décroissance des variances intraclasse
VI=np.cumsum(M[:,2]**2)/2
plt.figure()
plt.plot(np.arange(len(VI))+1,np.flip(VI,axis=0))
plt.xlabel("Nombre de classes")
plt.ylabel("Variance intraclasse")
if 1==0: # Autre figure possible
  plt.figure()
  plt.plot(np.arange(len(VI/max(VI)))+1,np.flip(VI/max(VI),axis=0))
  plt.xlabel("Nombre de classes")
  plt.ylabel("Variance intraclasse/variance totale")

print("\n******* Kmeans ******* \n")

# Comparaison avec les Kmeans
nclus=np.max(groupes)
k_means = KMeans(init='k-means++', n_clusters=nclus, n_init=10)
k_means.fit(fromages)
for k in range(nclus):
    print('Classe '+str(k+1).ljust(3,' ')+': ', end='')
    print(*nomfrom[np.where(k_means.labels_==k)])

print("\n******* Comparaison des inerties ******* \n")

print("Inertie Kmeans",nclus,"centres: ",k_means.inertia_)
print("Inertie CAH",nclus,"classes: ",VI[-nclus])

