#ifndef DSMDASMO_H
#define DSMDASMO_H

void dSMDaSMO(
    int NbLign,
    int *AdPrCoefLi, int *NumCol, int *AdSuccLi, float *Matrice, float *SecMembre,
    int *NumDLDir, float *ValDLDir,
    int *AdPrCoLi0, int *NumCol0, float *Matrice0, float *SecMemb0
);

int dSMOaLongPR(int NbLign, int* AdPrCoefLiO, int* NumColO, float* MatriceO);

void dSMOaPR(int *AdPrCoLi0, float *Matrice0, int *NumCol0,int longProfilMat, float *MatProf, int *Profil, int NbLign);

float solex(float *coor);

void CalSol(int NbLign, float **coord, float *UEX);


#endif