def proj(x,V,k): # Projection de x sur l'e.v. engendre par les k premieres col de V orthogonale
  M=V[:,:k]
  z=(x.dot(M)).dot(M.T)
  return z
def recon(x,V,mk,sk,k): # Reconstruction prenant en compte la standardisation
  z=(x-mk)/sk
  z=proj(z,V,k)
  z=np.multiply(z,sk)+mk
  return z
plt.figure()
for i in range(9):
  plt.subplot(331+i)
  plt.axis('off')
  imag=np.reshape(recon(X[i+12,:],V,mk,sk,Ncomp),[28,28])
  plt.imshow(imag,cmap=plt.get_cmap('gray'))

