/*
SAN Beining
OCHMAN Karolina 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "cal1Elem.h"
#include "tp2b.h"
#include "utils.h"

void cal1Elem(
    int typeEl,
    float *nodes,
    int *nRefArEl,
    int nbneel,
    float **MatEl,
    float *SmbEl,
    int *NuDElem,
    float *uDElem
)
{
    int i, l, e;
    int nbaret = nbAretEl(typeEl);

    float *w = NULL;
    float **xhat = NULL;
    int q = 0;

    /* Integration sur l'element */

    ppquad(typeEl, &w, &xhat, &q);

    for (l = 0; l < q; l++) {

        float coord_ref[2] = {0.0f, 0.0f};
        float dphi_ref[8];
        float J[4];
        float detJ;
        float eltdif;
        float cofvar = 1.0f;

        coord_ref[0] = xhat[l][0];

        if (typeEl != 3) {
            coord_ref[1] = xhat[l][1];
        }

        calDerFbase(typeEl, coord_ref, dphi_ref);
        matJacob(typeEl, coord_ref, dphi_ref, nodes, J);

        detJ = J[0] * J[3] - J[1] * J[2];

        if (fabsf(detJ) < 1e-12f) {
            printf("Jacobian singular in cal1Elem\n");
            exit(1);
        }

        eltdif = w[l] * fabsf(detJ);

        intElem(
            typeEl,
            coord_ref,
            nodes,
            nbneel,
            eltdif,
            cofvar,
            MatEl,
            SmbEl
        );
    }

    free(w);
    freetab(xhat);

    /* Integration sur les aretes */

    ppquad(3, &w, &xhat, &q);

    for (e = 1; e <= nbaret; e++) {

        int ref = nRefArEl[e - 1];

        if (typeCondition(ref) == 3) {

            int loc[2];
            float edge_nodes[4];

            numNaret(typeEl, e, loc);

            edge_nodes[0] = nodes[2 * loc[0]];
            edge_nodes[1] = nodes[2 * loc[0] + 1];
            edge_nodes[2] = nodes[2 * loc[1]];
            edge_nodes[3] = nodes[2 * loc[1] + 1];

            for (l = 0; l < q; l++) {

                float coords[2] = {0.0f, 0.0f};
                float dx = edge_nodes[2] - edge_nodes[0];
                float dy = edge_nodes[3] - edge_nodes[1];
                float len = sqrtf(dx * dx + dy * dy);
                float eltdif = w[l] * len;
                float cofvar = 1.0f;

                coords[0] = xhat[l][0];

                intAret(
                    typeEl,
                    edge_nodes,
                    coords,
                    nbneel,
                    loc,
                    eltdif,
                    cofvar,
                    MatEl,
                    SmbEl
                );
            }
        }
    }

    free(w);
    freetab(xhat);

    /* Conditions de Dirichlet */

    for (i = 0; i < nbneel; i++) {
        NuDElem[i] = 1;
        uDElem[i] = 0.0f;
    }

    for (e = 1; e <= nbaret; e++) {

        int ref = nRefArEl[e - 1];
        int type_cond = typeCondition(ref);

        if (type_cond == 1 || type_cond == 2) {

            int loc[2];
            numNaret(typeEl, e, loc);

            if (type_cond == 1) {

                NuDElem[loc[0]] = 0;
                NuDElem[loc[1]] = 0;

                uDElem[loc[0]] = 0.0f;
                uDElem[loc[1]] = 0.0f;
            }

            else if (type_cond == 2) {

                if (NuDElem[loc[0]] != 0) {
                    NuDElem[loc[0]] = -1;
                    uDElem[loc[0]] = UD(&nodes[2 * loc[0]]);
                }

                if (NuDElem[loc[1]] != 0) {
                    NuDElem[loc[1]] = -1;
                    uDElem[loc[1]] = UD(&nodes[2 * loc[1]]);
                }
            }
        }
    }
}
