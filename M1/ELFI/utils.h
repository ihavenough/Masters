#ifndef UTILS_H
#define UTILS_H

// Allocation mémoire
float **alloctab(int dim1, int dim2);
void freetab(void *ptr);

// Quadrature 
void ppquad(int t, float **w, float ***xhat, int *q);

// Fonctions de forme 
void calFbase(int t, float *coord, float *phi);
void calDerFbase(int t, float *coord, float *dphi);

// Transformations géométriques  
void transFK(int t, float *coord_ref, float *nodes, float *coord);
void matJacob(int t, float *coord_ref, float *dphi, float *nodes, float *jacobian);
int invertM2x2(float *mat, float *det);

// Utilitaires de maillage  
int nbAretEl(int typeEl);
void numNaret(int t, int edge, int *nodes);
void selectPts(int nb, int *num, float **coorEns, float **coorSel);

#endif