/*
SAN Beining
OCHMAN Karolina

Compilation :
gcc mainerr.c assemblage.c lecture.c utils.c cal1Elem.c tp2b.c impcalel.c tp4.c assmat.f affsmd.f cdesse.f affsmo.f tri.f ltlpr.f rsprl.f rspru.f affsol.f impmpr.f solex.c -o prog -lm -lgfortran

Execution :
./prog
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "utils.h"
#include "cal1Elem.h"
#include "forfun.h"
#include "tp4.h"
#include "tp2b.h"
#include "solex.h"

int nucas = 1;

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

/* Fonction définie dans tp4.c */
void CalSol(int NbLign, float **coord, float *UEX);

#define MAXCOEF 100000

int main(void)
{
    int liste_cas[] = {1, 2, 3, 4, 5, 6};

char *maillages_d1[] = {
    "d1t1_2",
    "d1t1_4",
    "d1t1_8",
    "d1t1_16",
    "d1t1_32",
    "d1t1_64",

    "d1q1_2",
    "d1q1_4",
    "d1q1_8",
    "d1q1_16",
    "d1q1_32",
    "d1q1_64"
};

char *maillages_d2[] = {
    "d2t1_2",
    "d2t1_4",
    "d2t1_8",
    "d2t1_16",
    "d2t1_32",
    "d2t1_64",

    "d2q1_2",
    "d2q1_4",
    "d2q1_8",
    "d2q1_16",
    "d2q1_32",
    "d2q1_64"
};

    int nb_cas = sizeof(liste_cas) / sizeof(liste_cas[0]);

    FILE *ferr = fopen("erreurs.txt", "w");

    if (ferr == NULL) {
        printf("Impossible de creer erreurs.txt\n");
        return 1;
    }

    for (int ic = 0; ic < nb_cas; ic++) {

        nucas = liste_cas[ic];

        char **liste_maillages;
        int nb_maillages;

        if (nucas <= 3) {
            liste_maillages = maillages_d1;
            nb_maillages = sizeof(maillages_d1) / sizeof(maillages_d1[0]);
        } else {
            liste_maillages = maillages_d2;
            nb_maillages = sizeof(maillages_d2) / sizeof(maillages_d2[0]);
        }

        for (int im = 0; im < nb_maillages; im++) {

            char ficmai[256];
            sprintf(ficmai, "%s", liste_maillages[im]);

            printf("\n====================================\n");
            printf("CAS = %d   MAILLAGE = %s\n", nucas, ficmai);
            printf("====================================\n");

            fprintf(ferr, "\n===== cas=%d maillage=%s =====\n", nucas, ficmai);

            int n, m, t, p, q;
            float **coord;
            int **ngnel, **nRefAr;

            if (lecfima(ficmai, &n, &coord, &m, &t, &p, &q, &ngnel, &nRefAr)) {
                printf("Erreur lecture fichier %s.\n", ficmai);
                fprintf(ferr, "Erreur lecture fichier %s.\n", ficmai);
                continue;
            }

            g_coord = coord;
            g_ngnel = ngnel;
            g_nRefAr = nRefAr;
            g_t = t;
            g_p = p;
            g_n = n;
            g_m = m;

            /* TP2 : calculs elementaires */

            for (int K = 0; K < m; K++) {

                int nbneel = p;

                float **MatElem = alloctab(nbneel, nbneel);

                if (MatElem == NULL) {
                    printf("Erreur allocation MatElem.\n");
                    fprintf(ferr, "Erreur allocation MatElem.\n");
                    fclose(ferr);
                    return 1;
                }

                for (int i = 0; i < nbneel; i++) {
                    for (int j = 0; j < nbneel; j++) {
                        MatElem[i][j] = 0.0f;
                    }
                }

                float *SMbrElem = calloc(nbneel, sizeof(float));
                int *NuDElem = malloc(nbneel * sizeof(int));
                float *uDElem = calloc(nbneel, sizeof(float));
                float *nodes = malloc(2 * nbneel * sizeof(float));

                if (SMbrElem == NULL || NuDElem == NULL ||
                    uDElem == NULL || nodes == NULL) {
                    printf("Erreur allocation TP2.\n");
                    fprintf(ferr, "Erreur allocation TP2.\n");
                    fclose(ferr);
                    return 1;
                }

                for (int i = 0; i < nbneel; i++) {
                    int id = ngnel[K][i] - 1;
                    nodes[2 * i]     = coord[id][0];
                    nodes[2 * i + 1] = coord[id][1];
                }

                cal1Elem(
                    t,
                    nodes,
                    nRefAr[K],
                    nbneel,
                    MatElem,
                    SMbrElem,
                    NuDElem,
                    uDElem
                );

                free(nodes);
                free(SMbrElem);
                free(NuDElem);
                free(uDElem);
                freetab(MatElem);
            }

            /* TP3 : assemblage S.M.D */

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

                    printf("Erreur allocation memoire TP3.\n");
                    fprintf(ferr, "Erreur allocation memoire TP3.\n");
                    fclose(ferr);
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

                /* TP4 : S.M.D vers S.M.O */

                int *AdPrCoLi0 = calloc(NbLign + 1, sizeof(int));
                int *NumCol0 = calloc(MAXCOEF, sizeof(int));
                float *Matrice0 = calloc(NbLign + MAXCOEF, sizeof(float));
                float *SecMemb0 = calloc(NbLign, sizeof(float));

                if (AdPrCoLi0 == NULL || NumCol0 == NULL ||
                    Matrice0 == NULL || SecMemb0 == NULL) {

                    printf("Erreur allocation S.M.O.\n");
                    fprintf(ferr, "Erreur allocation S.M.O.\n");
                    fclose(ferr);
                    return 1;
                }

                dSMDaSMO(
                    NbLign,
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

                /* TP5 : profil et Cholesky */

                int longProfilMat = dSMOaLongPR(
                    NbLign,
                    AdPrCoLi0,
                    NumCol0,
                    Matrice0
                );

                int *Profil = calloc(NbLign + 1, sizeof(int));
                float *MatProf = calloc(longProfilMat, sizeof(float));

                if (Profil == NULL || MatProf == NULL) {
                    printf("Erreur allocation Profil / MatProf.\n");
                    fprintf(ferr, "Erreur allocation Profil / MatProf.\n");
                    fclose(ferr);
                    return 1;
                }

                dSMOaPR(
                    AdPrCoLi0,
                    Matrice0,
                    NumCol0,
                    longProfilMat,
                    MatProf,
                    Profil,
                    NbLign
                );

                float *D = calloc(NbLign, sizeof(float));
                float *L = calloc(longProfilMat - NbLign, sizeof(float));

                if (D == NULL || L == NULL) {
                    printf("Erreur allocation Cholesky.\n");
                    fprintf(ferr, "Erreur allocation Cholesky.\n");
                    fclose(ferr);
                    return 1;
                }

                float eps = 1e-10;

                FORTRANNAME(ltlpr)(
                    &NbLign,
                    Profil,
                    MatProf,
                    MatProf + NbLign,
                    &eps,
                    D,
                    L
                );

                float *Y = calloc(NbLign, sizeof(float));
                float *U = calloc(NbLign, sizeof(float));

                if (Y == NULL || U == NULL) {
                    printf("Erreur allocation Y / U.\n");
                    fprintf(ferr, "Erreur allocation Y / U.\n");
                    fclose(ferr);
                    return 1;
                }

                FORTRANNAME(rsprl)(
                    &NbLign,
                    Profil,
                    D,
                    L,
                    SecMemb0,
                    Y
                );

                FORTRANNAME(rspru)(
                    &NbLign,
                    Profil,
                    D,
                    L,
                    Y,
                    U
                );

                /* Solution exacte */

                float *UEX = calloc(NbLign, sizeof(float));

                if (UEX == NULL) {
                    printf("Erreur allocation UEX.\n");
                    fprintf(ferr, "Erreur allocation UEX.\n");
                    fclose(ferr);
                    return 1;
                }

                CalSol(NbLign, coord, UEX);

                /* Calcul des erreurs */

                double somme_err2 = 0.0;
                double somme_uex2 = 0.0;

                double max_diff = 0.0;
                double max_uex = 0.0;

                for (int i = 0; i < NbLign; i++) {

                    double diff = (double)U[i] - (double)UEX[i];
                    double abs_diff = fabs(diff);
                    double abs_uex = fabs((double)UEX[i]);

                    /* Erreur quadratique relative :
                       ||U - UEX||_2 / ||UEX||_2 */
                    somme_err2 += diff * diff;
                    somme_uex2 += (double)UEX[i] * (double)UEX[i];

                    /* Pour l'erreur maximum relative :
                       ||U - UEX||_inf / ||UEX||_inf */
                    if (abs_diff > max_diff) {
                        max_diff = abs_diff;
                    }

                    if (abs_uex > max_uex) {
                        max_uex = abs_uex;
                    }
                }

                double err_quad_rel;

                if (somme_uex2 > 1e-14) {
                    err_quad_rel = sqrt(somme_err2 / somme_uex2);
                } else {
                    err_quad_rel = sqrt(somme_err2);
                }

                double err_max_rel;

                if (max_uex > 1e-14) {
                    err_max_rel = max_diff / max_uex;
                } else {
                    err_max_rel = max_diff;
                }

                printf("Nombre de noeuds : %d\n", NbLign);
                printf("Erreur quadratique relative : %.9f\n", err_quad_rel);
                printf("Erreur maximum relative     : %.9f\n", err_max_rel);

                fprintf(ferr, "Nombre de noeuds : %d\n", NbLign);
                fprintf(ferr, "Erreur quadratique relative : %.9f\n", err_quad_rel);
                fprintf(ferr, "Erreur maximum relative     : %.9f\n", err_max_rel);

                free(UEX);
                free(Y);
                free(U);
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
        }
    }

    fclose(ferr);

    printf("\nTous les calculs sont termines.\n");
    printf("Les erreurs sont sauvegardees dans erreurs.txt\n");

    return 0;
}
