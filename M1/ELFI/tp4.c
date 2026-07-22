/*
SAN Beining
OCHMAN Karolina 
*/

#include <stdio.h>
#include <stdlib.h>

#include "forfun.h"
#include "tp4.h"
#include "solex.h"

/* TP4 */

void dSMDaSMO(int NbLign,
              int *AdPrCoefLi,
              int *NumCol,
              int *AdSuccLi,
              float *Matrice,
              float *SecMembre,
              int *NumDLDir,
              float *ValDLDir,
              int *AdPrCoLi0,
              int *NumCol0,
              float *Matrice0,
              float *SecMemb0)
{
    FORTRANNAME(cdesse)(
        &NbLign,
        AdPrCoefLi,
        NumCol,
        AdSuccLi,
        Matrice,
        SecMembre,
        NumDLDir,
        ValDLDir,
        AdPrCoLi0,
        NumCol0,
        Matrice0,
        SecMemb0
    );
}


/* TP5 */

int dSMOaLongPR(int NbLign,
                int *AdPrCoefLiO,
                int *NumColO,
                float *MatriceO)
{
    int env_size = 0;

    for (int i = 1; i < NbLign; i++) {

        int beg_line = AdPrCoefLiO[i - 1] - 1;
        int end_line = AdPrCoefLiO[i] - 1;

        int min_col = i + 1;

        for (int k = beg_line; k < end_line; k++) {
            int col = NumColO[k];

            if (col >= 1 && col < i + 1 && col < min_col) {
                min_col = col;
            }
        }

        if (min_col < i + 1) {
            env_size += (i + 1) - min_col;
        }
    }

    return NbLign + env_size;
}


void dSMOaPR(int *AdPrCoLi0,
             float *Matrice0,
             int *NumCol0,
             int longProfilMat,
             float *MatProf,
             int *Profil,
             int NbLign)
{
    for (int i = 0; i < longProfilMat; i++) {
        MatProf[i] = 0.0f;
    }

    for (int i = 0; i < NbLign; i++) {
        MatProf[i] = Matrice0[i];
    }

    Profil[0] = 1;

    for (int i = 1; i < NbLign; i++) {

        int beg_line = AdPrCoLi0[i - 1] - 1;
        int end_line = AdPrCoLi0[i] - 1;

        int min_col = i + 1;
        int size_line = 0;

        for (int k = beg_line; k < end_line; k++) {
            int col = NumCol0[k];

            if (col >= 1 && col < i + 1 && col < min_col) {
                min_col = col;
            }
        }

        if (min_col < i + 1) {
            size_line = (i + 1) - min_col;
        }

        Profil[i] = Profil[i - 1] + size_line;
    }

    Profil[NbLign] = Profil[NbLign - 1];

    for (int i = 1; i < NbLign; i++) {

        int beg_line = AdPrCoLi0[i - 1] - 1;
        int end_line = AdPrCoLi0[i] - 1;

        int min_col = i + 1;

        for (int k = beg_line; k < end_line; k++) {
            int col = NumCol0[k];

            if (col >= 1 && col < i + 1 && col < min_col) {
                min_col = col;
            }
        }

        if (min_col < i + 1) {

            for (int k = beg_line; k < end_line; k++) {

                int col = NumCol0[k];

                if (col >= 1 && col < i + 1) {

                    int position_profil =
                        (Profil[i - 1] - 1) + (col - min_col);

                    MatProf[NbLign + position_profil] =
                        Matrice0[NbLign + k];
                }
            }
        }
    }
}


void CalSol(int NbLign, float **coord, float *UEX)
{
    for (int i = 0; i < NbLign; i++) {
        UEX[i] = solex(coord[i]);
    }
}
