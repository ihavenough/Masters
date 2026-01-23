#COMMANDES DE LECTURE POUR LES DIGITS AVEC UN EXEMPLE D"ILLUSTRATION

X=np.loadtxt("digits_extrait_images.csv",
  delimiter=',',skiprows=1,usecols=range(1,785))
labels=np.loadtxt("digits_extrait_labels.csv",
  delimiter=',',dtype='str',skiprows=1,usecols=range(1,2))

col=7
for i in range(3*col):
  plt.subplot(3,col,1+i)
  plt.imshow(np.reshape(X[i,:],[28,28]), cmap=plt.get_cmap('gray'))
  plt.title(int(labels[i]))
  plt.axis('off')
plt.show()
