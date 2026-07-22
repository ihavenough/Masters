/*
SAN Beining
OCHMAN Karolina 
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "utils.h"
#include "tp2b.h"
#include "solex.h"

extern int nucas;

/* Terme de masse */
float a00(float *x)
{
    switch (nucas) {

        case 1:
        case 2:
        case 4:
        case 5:
            return 0.0f;

        case 3:
        case 6:
            return 1.0f;

       // case 4:
        //case 5:
        //    return 0.0f;

        //case 6:
           // return 1.0f;

        default:
            return 0.0f;
    }
}

/* Second membre */
float fomega(float *coor)
{
    switch (nucas) {

        case 1:
        case 4:
            return 32.0f *
                   (coor[1] * (1.0f - coor[1])
                  + coor[0] * (1.0f - coor[0]));

        case 2:
        case 5:
            return 2.0f * M_PI * M_PI
                   * sinf(M_PI * coor[0])
                   * sinf(M_PI * coor[1]);

        case 3:
        case 6:
            return (2.0f * M_PI * M_PI + 1.0f)
                   * cosf(M_PI * coor[0])
                   * cosf(M_PI * coor[1]);
/*
        case 4:
            return 32.0f *
                   (coor[1] * (1.0f - coor[1])
                  + coor[0] * (1.0f - coor[0]));

        case 5:
            return 2.0f * M_PI * M_PI
                   * sinf(M_PI * coor[0])
                   * sinf(M_PI * coor[1]);

        case 6:
            return (2.0f * M_PI * M_PI + 1.0f)
                   * cosf(M_PI * coor[0])
                   * cosf(M_PI * coor[1]);
*/
        default:
            return 0.0f;
    }
}

/* Coefficients de diffusion */
float a11(float *x)
{
    return 1.0f;
}

float a22(float *x)
{
    return 1.0f;
}

float a12(float *x)
{
    return 0.0f;
}

float a21(float *x)
{
    return a12(x);
}

/* Dirichlet */
float UD(float *x)
{
    switch (nucas) {

        case 1:
        case 2:
           // return 0.0f;
        case 3:
        case 4:
        case 5:
        case 6:
            return solex(x);

        default:
            return 0.0f;
    }
}

/* Terme de Robin éventuel */
float bn(float *x)
{
    return 0.0f;
}

/* Neumann */
float fn(float *x)
{
    const float PI = M_PI;
    const float L = 1.0f / 3.0f;
    const float eps = 1e-6f;

    switch (nucas) {

        case 3:
            return 0.0f;

        case 6:
            /*
               Domaine 2 : [0,1/3] x [0,1/3]
               u = cos(pi x) cos(pi y)

               Neumann sur x = 1/3 :
               du/dn = du/dx

               Neumann sur y = 1/3 :
               du/dn = du/dy
            */

            if (fabsf(x[0] - L) < eps) {
                return -PI * sinf(PI * x[0]) * cosf(PI * x[1]);
            }

            if (fabsf(x[1] - L) < eps) {
                return -PI * cosf(PI * x[0]) * sinf(PI * x[1]);
            }

            return 0.0f;

        default:
            return 0.0f;
    }
}

void W(int nbneel, float *fctbas, float eltdif, float cofvar, float *vecelm)
{
    int j;

    for (j = 0; j < nbneel; j++) {
        vecelm[j] = vecelm[j] + eltdif * cofvar * fctbas[j];
    }
}

void WW(int nbneel, float *fctbas, float eltdif, float cofvar, float **matelm)
{
    int i, j;
    float coeff;

    for (i = 0; i < nbneel; i++) {
        coeff = eltdif * cofvar * fctbas[i];

        for (j = 0; j < nbneel; j++) {
            matelm[i][j] = matelm[i][j] + coeff * fctbas[j];
        }
    }
}

void ADWDW(int t,
           int nbneel,
           float *coord_ref,
           float *nodes,
           float eltdif,
           float **matelm)
{
    int i, j;

    float dphi_ref[8];
    float jac[4];
    float det;
    float grad[8];
    float coord[2];

    calDerFbase(t, coord_ref, dphi_ref);
    matJacob(t, coord_ref, dphi_ref, nodes, jac);

    if (!invertM2x2(jac, &det)) {
        printf("Jacobian singular\n");
        exit(1);
    }

    transFK(t, coord_ref, nodes, coord);

    for (i = 0; i < nbneel; i++) {
        grad[2 * i] =
            jac[0] * dphi_ref[i]
          + jac[2] * dphi_ref[i + nbneel];

        grad[2 * i + 1] =
            jac[1] * dphi_ref[i]
          + jac[3] * dphi_ref[i + nbneel];
    }

    for (i = 0; i < nbneel; i++) {
        for (j = 0; j < nbneel; j++) {

            float val = 0.0f;

            val += a11(coord) * grad[2 * i]     * grad[2 * j];
            val += a12(coord) * grad[2 * i]     * grad[2 * j + 1];
            val += a21(coord) * grad[2 * i + 1] * grad[2 * j];
            val += a22(coord) * grad[2 * i + 1] * grad[2 * j + 1];

            matelm[i][j] += eltdif * val;
        }
    }
}

void intElem(int t,
             float *coord_ref,
             float *nodes,
             int nbneel,
             float eltdif,
             float cofvar,
             float **matelm,
             float *vecelm)
{
    float phi[4];
    float dphi_ref[8];
    float jac[4];
    float det;
    float coord[2];

    calFbase(t, coord_ref, phi);
    calDerFbase(t, coord_ref, dphi_ref);

    matJacob(t, coord_ref, dphi_ref, nodes, jac);

    if (!invertM2x2(jac, &det)) {
        printf("Jacobian singular\n");
        exit(1);
    }

    transFK(t, coord_ref, nodes, coord);

    ADWDW(t, nbneel, coord_ref, nodes, eltdif, matelm);

    WW(nbneel, phi, eltdif, a00(coord), matelm);

    W(nbneel, phi, eltdif, fomega(coord), vecelm);
}

void intAret(int t,
             float *nodes,
             float *coords,
             int nbneel,
             int *loc,
             float eltdif,
             float cofvar,
             float **matelm,
             float *vecelm)
{
    float phi_1d[2];
    float phi_element[4];
    float dphi_ref[2];
    float jacobian[1];
    float coord[2];

    calFbase(3, coords, phi_1d);
    calDerFbase(3, coords, dphi_ref);

    for (int i = 0; i < nbneel; i++) {
        phi_element[i] = 0.0f;
    }

    phi_element[loc[0]] = phi_1d[0];
    phi_element[loc[1]] = phi_1d[1];

    matJacob(3, coords, dphi_ref, nodes, jacobian);

    coord[0] = nodes[0] * coords[0]
             + nodes[2] * (1.0f - coords[0]);

    coord[1] = nodes[1] * coords[0]
             + nodes[3] * (1.0f - coords[0]);

    WW(nbneel, phi_element, eltdif, bn(coord), matelm);
    W(nbneel, phi_element, eltdif, fn(coord), vecelm);
}

/*
Retour :
0 = interne
1 = Dirichlet homogène
2 = Dirichlet non homogène
3 = Neumann
*/
int typeCondition(int ref)
{
    if (ref == 0) {
        return 0;
    }

    switch (nucas) {

        case 1:
        case 2:
            if (ref >= 1 && ref <= 4) {
                return 1;
            }
            break;

        case 3:
            if (ref >= 1 && ref <= 4) {
                return 3;
            }
            break;

        case 4:
        case 5:
            if (ref >= 1 && ref <= 4) {
                return 2;
            }
            break;

        case 6:
            /*
               Domaine 2 :
               ref 1 : y = 0       -> Dirichlet
               ref 2 : x = 1/3     -> Neumann
               ref 3 : y = 1/3     -> Neumann
               ref 4 : x = 0       -> Dirichlet
            */
            if (ref == 1 || ref == 4) {
                return 2;
            }

            if (ref == 2 || ref == 3) {
                return 3;
            }

            break;
    }

    return 0;
}
