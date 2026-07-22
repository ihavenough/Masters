/*
SAN Beining
OCHMAN Karolina
gcc main.c assemblage.c lecture.c utils.c cal1Elem.c tp2b.c impcalel.c tp4.c assmat.f affsmd.f cdesse.f affsmo.f tri.f ltlpr.f rsprl.f rspru.f affsol.f impmpr.f solex.c -o prog -lm -lgfortran 
./prog
*/

#include <stdio.h>
#include <stdlib.h>
#include "utils.h"
#include "cal1Elem.h"
#include "forfun.h"
#include "tp4.h"
#include "tp2b.h"
#include "solex.h"



int nucas = 1; // / cas 1 à 6

float **g_coord;
int **g_ngnel, **g_nRefAr;
int g_t, g_p, g_n, g_m;

int lecfima(char *ficmai,
            int *n, float ***coord,
            int *m, int *t, int *p, int *q,
            int ***ngnel, int ***nRefAr);

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
);

#define MAXCOEF 100000

int main(void)
{
    printf("Numero de l'exemple (entre 1 et 6)?\n");
    if (scanf("%d", &nucas) != 1) { 
        return 1;
    }

    if (nucas < 1 || nucas > 6) {
        printf("Numero de cas invalide.\n");
        return 1;
    }
    char ficmai[256];
    printf("Donner le nom du fichier de maillage : ");
    scanf("%255s", ficmai);

    int n, m, t, p, q;
    float **coord;
    int **ngnel, **nRefAr;

    if (lecfima(ficmai, &n, &coord, &m, &t, &p, &q, &ngnel, &nRefAr)) {
        printf("Erreur lecture fichier.\n");
        return 1;
    }

    g_coord = coord;
    g_ngnel = ngnel;
    g_nRefAr = nRefAr;
    g_t = t;
    g_p = p;
    g_n = n;
    g_m = m;

    // TP2

    for (int K = 0; K < m; K++) {

        int nbneel = p;

        // allocate element data
        float **MatElem = alloctab(nbneel, nbneel);
        for (int i = 0; i < nbneel; i++)
            for (int j = 0; j < nbneel; j++)
                MatElem[i][j] = 0.0f;

        float *SMbrElem = calloc(nbneel, sizeof(float));
        int *NuDElem = malloc(nbneel * sizeof(int));
        float *uDElem = calloc(nbneel, sizeof(float));

        // extract node coordinates
        float *nodes = malloc(2 * nbneel * sizeof(float));

        for (int i = 0; i < nbneel; i++) {
            int id = ngnel[K][i] - 1;
            nodes[2*i]   = coord[id][0];
            nodes[2*i+1] = coord[id][1];
        }

        //cal1Elem
        cal1Elem(t, nodes, nRefAr[K], nbneel,
                 MatElem, SMbrElem, NuDElem, uDElem);

        //printing
        impCalEl(K+1, t, nbneel, MatElem, SMbrElem, NuDElem, uDElem);

        // free
        free(nodes);
        free(SMbrElem);
        free(NuDElem);
        free(uDElem);
        freetab(MatElem);
    }

    // TP3
    {
        int NbLign = n;

        float *Matrice = calloc(NbLign + MAXCOEF, sizeof(float));
        float *SecMembre = calloc(NbLign, sizeof(float));

        int *AdPrCoefLi = calloc(NbLign + 1, sizeof(int));
        int *NumCol = calloc(MAXCOEF, sizeof(int));
        int *AdSuccLi = calloc(MAXCOEF, sizeof(int));

        int *NumDLDir = calloc(NbLign, sizeof(int));
        float *ValDLDir = calloc(NbLign, sizeof(float));

        if (Matrice == NULL || SecMembre == NULL || AdPrCoefLi == NULL ||
            NumCol == NULL || AdSuccLi == NULL ||
            NumDLDir == NULL || ValDLDir == NULL) {
            printf("Erreur allocation mémoire.\n");
            return 1;
        }


        Assemblage(
            NbLign,
            m,
            Matrice,
            SecMembre,
            AdPrCoefLi,
            NumCol,
            AdSuccLi,
            NumDLDir,
            ValDLDir
        );
        printf("\nAFFICHAGE DE LA S.M.D \n");
        FORTRANNAME(affsmd)(
            &NbLign,
            AdPrCoefLi,
            NumCol,
            AdSuccLi,
            Matrice,
            SecMembre,
            NumDLDir,
            ValDLDir
        );

        // TP4 

        printf("\nCONSTRUCTION DE LA S.M.O \n");

        int *AdPrCoLi0 = calloc(NbLign + 1, sizeof(int));
        int *NumCol0 = calloc(MAXCOEF, sizeof(int));
        float *Matrice0 = calloc(NbLign + MAXCOEF, sizeof(float));
        float *SecMemb0 = calloc(NbLign, sizeof(float));

        if (AdPrCoLi0 == NULL || NumCol0 == NULL ||
        Matrice0 == NULL || SecMemb0 == NULL) {
        printf("Erreur allocation S.M.O.\n");
        return 1;
    }

        dSMDaSMO(
            NbLign,
            AdPrCoefLi, NumCol, AdSuccLi, Matrice, SecMembre,
            NumDLDir, ValDLDir,
            AdPrCoLi0, NumCol0, Matrice0, SecMemb0
        );

        printf("\nAFFICHAGE DE LA S.M.O \n");
        FORTRANNAME(affsmo)(
            &NbLign,
            AdPrCoLi0, NumCol0, Matrice0, SecMemb0
        );

        // TP5

        int longProfilMat = dSMOaLongPR(NbLign, AdPrCoLi0, NumCol0, Matrice0);

        int *Profil = malloc((NbLign + 1) * sizeof(int));
        float *MatProf = calloc(longProfilMat, sizeof(float));


        //dSMOaPR(AdPrCoLi0, Matrice0, NumCol0, &MatProf, &Profil,  NbLign);
        dSMOaPR( AdPrCoLi0, Matrice0, NumCol0, longProfilMat,MatProf, Profil, NbLign);


        // Cholesky
        printf("\n RESOLUTION PAR CHOLESKY\n");

        float *D = calloc(NbLign, sizeof(float));
        float *L = calloc(longProfilMat - NbLign, sizeof(float));

        float eps = 1e-10;

        if (D == NULL || L == NULL ) {
            printf("Erreur allocation Cholesky.\n");
            return 1;
        }

// A = L * L^T
        FORTRANNAME(ltlpr)(&NbLign, Profil, MatProf, MatProf + NbLign, &eps, D, L);

//  L * y = b 
        FORTRANNAME(rsprl)(&NbLign, Profil, D, L, SecMemb0, SecMemb0);

// L^T * x = y 
        float *U = SecMemb0; // U contient la solution de EF
        FORTRANNAME(rspru)(&NbLign, Profil, D, L, U, U);

        /*
        printf("\n SOLUTION NUMERIQUE U \n");
        for (int i = 0; i < NbLign; i++) {
            printf("U[%d] = %f\n", i + 1, U[i]);
        }

          printf("\nCALCUL DE LA SOLUTION EXACTE\n");
*/
        float *UEX = calloc(NbLign, sizeof(float));

        if (UEX == NULL) {
            printf("Erreur allocation UEX.\n");
            return 1;
        }

        CalSol(NbLign, coord, UEX);

        float *coordFortran = calloc(2 * NbLign, sizeof(float));

        if (coordFortran == NULL) {
            printf("Erreur allocation coordFortran.\n");
            return 1;
        }

        for (int i = 0; i < NbLign; i++) {
            coordFortran[2 * i]     = coord[i][0];
            coordFortran[2 * i + 1] = coord[i][1];
        }

        int IMPFCH = 0; // 0 or 1
        // 1 saves info in a file

        printf("\n ERREUR ENTRE U ET UEX \n");

        FORTRANNAME(affsol)(
            &NbLign,
            coordFortran,
            U,
            UEX,
            &IMPFCH
        );

        // free the memory
        free(coordFortran);
        free(UEX);
        free(L);
        free(D);
        free(MatProf);
        free(Profil);

        free(AdPrCoLi0);
        free(NumCol0);
        free(Matrice0);
        free(SecMemb0);

        free(Matrice);
        free(SecMembre);
        free(AdPrCoefLi);
        free(NumCol);
        free(AdSuccLi);
        free(NumDLDir);
        free(ValDLDir);
    }

    free(ngnel[0]);
    free(ngnel);
    free(nRefAr[0]);
    free(nRefAr);
    freetab(coord);


    return 0;
}

