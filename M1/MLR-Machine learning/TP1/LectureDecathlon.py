#COMMANDES DE LECTURE POUR LE DECATHLON

dec=np.loadtxt("Decathlon.csv",dtype='str',delimiter=',')
indiv=dec[:,0][1:]
dec=dec[:,1:-2] # On retire la premiere et les deux dernieres colonnes
varbs=dec[0,:]
X=dec[1:,:].astype(float)

