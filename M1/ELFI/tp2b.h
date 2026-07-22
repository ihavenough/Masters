#ifndef TP2B_H
#define TP2B_H

// Intégration de surface
void intElem(int t, float *coord_ref, float *nodes, int nbneel, float eltdif, float cofvar, float **matelm, float *vecelm);

// Intégration sur les arêtes (Neumann/Fourier)
void intAret(int t, float *nodes, float *coords, int nbneel, int *loc, float eltdif, float cofvar, float **matelm, float *vecelm);

// Condition de Dirichlet non homogène
float UD(float *x);

int typeCondition(int ref);

void impCalEl(int K, int typEl, int nbneel, float **MatElem, float *SMbrElem,
              int *NuDElem, float *uDElem);

#endif