#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 25 09:11:35 2018

@author: valerie
"""


# importation des données -----------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
# importation des données -----------------------------------------------------
nomvar=np.loadtxt("HeartData.txt",dtype='str',max_rows=1,delimiter=",")[1:]
heart=np.loadtxt("HeartData.txt",skiprows=1,delimiter=",")[:,1:]
print(*nomvar)
print(*heart[12,:])
X = heart[:,:-1]
y = heart[:,-1]

# Statistiques descriptives: boxplots -----------------------------------------
quanti = np.where((nomvar!="famhist")&(nomvar!="chd"))[0]
plt.close('all')
plt.figure()
for j,nom in enumerate(quanti):
    v = heart[:,quanti[j]]
    data = [v[y==0],v[y==1]]
    plt.subplot(2,4,j+1)
    plt.subplots_adjust(hspace=1,wspace=1)
    plt.boxplot(data, labels=("No","CHD"),widths=.8)
    plt.grid()
    plt.title(nomvar[quanti[j]])

# Analyse discriminante linéaire ----------------------------------------------
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import confusion_matrix,ConfusionMatrixDisplay
lda = LinearDiscriminantAnalysis()
lda.fit(X,y)
yhat = lda.predict(X)
errl=sum(y!=yhat)/len(y)
print("Taux d'erreur: ",round(errl,3))
plt.rcParams.update({'figure.figsize': (3,3),'font.size': 16})
conf_mat =  confusion_matrix(y,yhat)
ConfusionMatrixDisplay(conf_mat).plot(cmap='YlOrBr',colorbar=False)
plt.rcdefaults()

# Pour tracer soi-même la matrice de confusion en format texte
print("Matrice de confusion (sans validation croisée)\n")
conf=conf_mat.astype(str)
conf=np.c_[['y=0 ','y=1 '],conf]
conf=np.r_[[['   ','yh=0','yh=1']],conf]
conft=conf.T
for i,b in enumerate(conft):
  l=max(len(i) for i in b) # Largeur de la colonne i de conf
  conft[i]=np.char.ljust(b,l+1,' ')
conf=conft.T
for row in conf: print(*row)

# La variable discriminante
C = lda.fit_transform(X, y).reshape((-1,)) # On a 2 classes donc 1 seul axe factoriel
tmp = [C[np.where(y==0)[0]],C[np.where(y==1)[0]]]
plt.figure()
plt.boxplot(tmp, labels=("No CHD","CHD"),widths=.8)
plt.title("Distribution de la variable discriminante dans les deux classes")

# Avec validation croisée
ntest=np.floor(len(y)/2).astype(int)
per=np.random.permutation(len(y))
lt,la=per[:ntest], per[ntest:]
Xa,Xt=X[la,:],X[lt,:]
ya,yt=y[la],y[lt]
#lda.fit(????)
#yhat = lda.predict(????)
#errl=
#print("Taux d'erreur: ",round(errl,3))


# Analyse discriminante quadratique ----------------------------------------------
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis

# A vous!

