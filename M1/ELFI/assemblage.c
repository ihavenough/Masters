/*
SAN Beining
OCHMAN Karolina
gcc main.c assemblage.c lecture.c utils.c cal1Elem.c tp2b.c assmat.f affsmd.f -o prog -lm -lgfortran
./prog
*/

#include <stdlib.h>
#include "forfun.h"
#include "utils.h"
#include "cal1Elem.h"

extern float **g_coord;
extern int **g_ngnel, **g_nRefAr;
extern int g_t, g_p, g_n, g_m;

void Assemblage(
    int NbLign,
    int NbElem,
    float *Matrice,
    float *SecMembre,
    int *AdPrCoefLi,
    int *NumCol,
    int *AdSuccLi,
    int *NumDLDir,
    float *ValDLDir
)
{
    int i, j, K;
    int NextAd = 1;

    for (i = 0; i < NbLign; i++) {
        Matrice[i] = 0.0f;
        SecMembre[i] = 0.0f;
        AdPrCoefLi[i] = 0;
        NumDLDir[i] = i + 1;
        ValDLDir[i] = 0.0f;
    }

    AdPrCoefLi[NbLign] = 0;

    for (K = 0; K < NbElem; K++) {
        int nbneel = g_p;

        float **MatElem = alloctab(nbneel, nbneel);
        float *SMbrElem = calloc(nbneel, sizeof(float));
        int *NuDElem = malloc(nbneel * sizeof(int));
        float *uDElem = calloc(nbneel, sizeof(float));
        float *nodes = malloc(2 * nbneel * sizeof(float));

        if (MatElem == NULL || SMbrElem == NULL || NuDElem == NULL ||
            uDElem == NULL || nodes == NULL) {
            exit(1);
        }

        for (i = 0; i < nbneel; i++) {
            int id = g_ngnel[K][i] - 1;
            nodes[2 * i] = g_coord[id][0];
            nodes[2 * i + 1] = g_coord[id][1];
        }

        for (i = 0; i < nbneel; i++) {
            for (j = 0; j < nbneel; j++) {
                MatElem[i][j] = 0.0f;
            }
        }

        cal1Elem(g_t, nodes, g_nRefAr[K], nbneel, MatElem, SMbrElem, NuDElem, uDElem);

        for (i = 0; i < nbneel; i++) {
            int I = g_ngnel[K][i];
            SecMembre[I - 1] += SMbrElem[i];
        }

        for (i = 0; i < nbneel; i++) {
            int I = g_ngnel[K][i];

            for (j = 0; j <= i; j++) {
                int J = g_ngnel[K][j];
                float X = MatElem[i][j];

                if (I == J) {
                    Matrice[I - 1] += X;
                } else {
                    int II = I;
                    int JJ = J;

                    if (II < JJ) {
                        int tmp = II;
                        II = JJ;
                        JJ = tmp;
                    }

                    FORTRANNAME(assmat)(
                        &II, &JJ, &X,
                        AdPrCoefLi,
                        NumCol,
                        AdSuccLi,
                        Matrice + NbLign,
                        &NextAd
                    );
                }
            }
        }

        for (i = 0; i < nbneel; i++) {
            int I = g_ngnel[K][i];

            if (NuDElem[i] == -1) {
                NumDLDir[I - 1] = -I;
                ValDLDir[I - 1] = uDElem[i];
            } else if (NuDElem[i] == 0) {
                if (NumDLDir[I - 1] != -I) {
                    NumDLDir[I - 1] = 0;
                    ValDLDir[I - 1] = 0.0f;
                }
            }
        }

        free(nodes);
        free(SMbrElem);
        free(NuDElem);
        free(uDElem);
        freetab(MatElem);
    }

    AdPrCoefLi[NbLign] = NextAd;
}
