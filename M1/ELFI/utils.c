/*
SAN Beining
OCHMAN Karolina 
*/

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "utils.h"

int nbAretEl(int typeEl){
    if (typeEl == 1) return 4;
    if (typeEl == 2) return 3;
    if (typeEl == 3) return 1;
    return 0;
}

void ppquad(int t, float **w, float ***xhat, int *q) {
    if (t == 1) { // Quadrangle, ordre 3
        *q = 9;
        *w = malloc(sizeof(float) * 9);
        *xhat = alloctab(9, 2); // 9 points avec 2 valeurs
        
        // poids quadrangle
        for (int i = 0; i < 4; i++) {
            (*w)[i] = 1.0 / 36.0;
        }

        for (int i = 4; i < 8; i++) {
            (*w)[i] = 1.0 / 9.0;
        }

        (*w)[8] = 4.0 / 9.0;
        
        // Points pour quadrangle
        (*xhat)[0][0] = 1.0; (*xhat)[0][1] = 0.0; // (1, 0)
        (*xhat)[1][0] = 1.0; (*xhat)[1][1] = 1.0; // (1, 1)
        (*xhat)[2][0] = 0.0; (*xhat)[2][1] = 1.0; // (0, 1)
        (*xhat)[3][0] = 0.0; (*xhat)[3][1] = 0.0; // (0, 0)
        (*xhat)[4][0] = 1.0; (*xhat)[4][1] = 0.5; // (1, 0.5)
        (*xhat)[5][0] = 0.5; (*xhat)[5][1] = 1.0; // (0.5, 1)
        (*xhat)[6][0] = 0.0; (*xhat)[6][1] = 0.5; // (0, 0.5)
        (*xhat)[7][0] = 0.5; (*xhat)[7][1] = 0.0; // (0.5, 0)
        (*xhat)[8][0] = 0.5; (*xhat)[8][1] = 0.5; // (0.5, 0.5)
    }
    else if (t == 2) { // Triangle, ordre 2
        *q = 3;
        *w = malloc(sizeof(float) * 3);
        *xhat = alloctab(3, 2); //  3 points avec 2 valeurs
        
        // poids pour triangle
        (*w)[0] = (*w)[1] = (*w)[2] = 1.0 / 6.0;
        
        // Points triangle
        (*xhat)[0][0] = 0.5; (*xhat)[0][1] = 0.5; // (0.5, 0.5)
        (*xhat)[1][0] = 0.0; (*xhat)[1][1] = 0.5; // (0, 0.5)
        (*xhat)[2][0] = 0.5; (*xhat)[2][1] = 0.0; // (0.5, 0)
    }
    else if (t == 3) { // Segment, ordre 3
        *q = 3;  
        *w = malloc(sizeof(float) * 3);
        *xhat = alloctab(3, 1); //3 points avec 1 valeur
        
        // poids segment
        (*w)[0] = (*w)[1] = 1.0 / 6.0;
        (*w)[2] = 2.0 / 3.0;
        
        // points segment
        (*xhat)[0][0] = 1.0; // (1)
        (*xhat)[1][0] = 0.0; // (0)
        (*xhat)[2][0] = 0.5; // (0.5)
    }
}

void calFbase(int t, float *coord, float *phi) {

    if (t == 1) { // Quadrangle
        phi[0] = coord[0] * (1 - coord[1]); // P1
        phi[1] = coord[0] * coord[1]; // P2
        phi[2] = (1 - coord[0]) * coord[1]; // P3
        phi[3] = (1 - coord[0]) * (1 - coord[1]); //P4

    }
    else if (t == 2) { // Triangle
        phi[0] = 1 - coord[0] - coord[1]; //P1
        phi[1] = coord[0]; // P2
        phi[2] = coord[1]; // P3
    }
    else if (t == 3) { // Segment
        phi[0] = coord[0]; //P1
        phi[1] = 1 - coord[0]; //P2
    }
}

void calDerFbase(int t, float *coord, float *dphi) {

    if (t == 1) { // Quadrangle
        dphi[0] = 1 - coord[1]; // dP1/dx
        dphi[1] = coord[1];  //dP2/dx
        dphi[2] = -coord[1];  // dP3/dx
        dphi[3] = -1 + coord[1];     // dP4/dx
        
        dphi[4] = -coord[0];     // dP1/dy
        dphi[5] = coord[0];     // dP2/dy
        dphi[6] = 1 - coord[0];      // dP3/dy
        dphi[7] = -1 + coord[0];  // dP4/dy
    }
    else if (t == 2) { // Triangle
        dphi[0] = -1; // dP1/dx
        dphi[1] = 1;  // dP2/dx
        dphi[2] = 0;  // dP3/dx
        
        dphi[3] = -1; // dP1/dy
        dphi[4] = 0;  // dP2/dy
        dphi[5] = 1;  // dP3/dy
    }
    else if (t == 3) { // Segment
        dphi[0] = 1; // dP1/dx
        dphi[1] = -1;  // dP2/dx
    }
}

// node order: BD, HD, HG, BG

void transFK(int t, float *coord_ref, float *nodes, float *coord) {
    if (t == 1) { //Quadrangle
        coord[0] = nodes[0] * coord_ref[0] * (1 - coord_ref[1]) +nodes[2] * coord_ref[0] * coord_ref[1] +
                   nodes[4] * (1 - coord_ref[0]) * coord_ref[1] +nodes[6] * (1 - coord_ref[0]) * (1 - coord_ref[1]);

        coord[1] = nodes[1] * coord_ref[0] * (1 - coord_ref[1]) +nodes[3] * coord_ref[0] * coord_ref[1] +
                   nodes[5] * (1 - coord_ref[0]) * coord_ref[1] +nodes[7] * (1 - coord_ref[0]) * (1 - coord_ref[1]);
    }

    else if (t == 2) { //Triangle
        coord[0] = nodes[0] * (1- coord_ref[0] -coord_ref[1]) +nodes[2] * coord_ref[0] + nodes[4] * coord_ref[1];

        coord[1] = nodes[1] * (1 - coord_ref[0] - coord_ref[1]) +nodes[3] * coord_ref[0] +nodes[5] * coord_ref[1];
    }

    else if (t == 3) { //Segment
        coord[0] = nodes[0] * coord_ref[0] + nodes[1]*(1 - coord_ref[0]) ;
    }
}

void matJacob(int t,float *coord_ref, float *dphi, float *nodes, float *jacobian) {
    if (t == 1) { //quadrangle
        // Jacobian matrix
        jacobian[0] = dphi[0] * nodes[0] + dphi[1] * nodes[2] + dphi[2] * nodes[4] + dphi[3] * nodes[6]; // dX/dxi
        jacobian[1] = dphi[4] * nodes[0] + dphi[5] * nodes[2] + dphi[6] * nodes[4] + dphi[7] * nodes[6]; // dX/deta
        jacobian[2] = dphi[0] * nodes[1] + dphi[1] * nodes[3] + dphi[2] * nodes[5] + dphi[3] * nodes[7]; // dY/dxi
        jacobian[3] = dphi[4] * nodes[1] + dphi[5] * nodes[3] + dphi[6] * nodes[5] + dphi[7] * nodes[7]; // dY/deta
    }
    else if (t == 2) { //triangle
        // Jacobian for triangle (2D): 2x2 matrix
        jacobian[0] = dphi[0] * nodes[0] + dphi[1] * nodes[2] + dphi[2] * nodes[4]; // dX/dxi
        jacobian[1] = dphi[3] * nodes[0] + dphi[4] * nodes[2] + dphi[5] * nodes[4]; // dX/deta
        jacobian[2] = dphi[0] * nodes[1] + dphi[1] * nodes[3] + dphi[2] * nodes[5]; // dY/dxi
        jacobian[3] = dphi[3] * nodes[1] + dphi[4] * nodes[3] + dphi[5] * nodes[5]; // dY/deta
    }
    else if (t == 3) { //segment
        // Jacobian for segment (1D): 1x1 matrix (length of the segment)
        float dx = nodes[2] - nodes[0];
        float dy = nodes[3] - nodes[1];
        jacobian[0] = sqrtf(dx*dx + dy*dy);
    }
}

int invertM2x2(float *mat, float *det) {
    *det = mat[0] * mat[3] - mat[1] * mat[2];   //determinant
    if (fabs(*det) < 1e-8) {
        return 0; // matrice singuliere, pas de inverse
    }
    
    float invDet = 1.0 / (*det);
    
    //inverse of 2x2 matrix
    float invMat[4];
    invMat[0] = mat[3] * invDet;
    invMat[1] = -mat[1] * invDet;
    invMat[2] = -mat[2] * invDet;
    invMat[3] = mat[0] * invDet;
    
    // back to original matrix
    mat[0] = invMat[0];
    mat[1] = invMat[1];
    mat[2] = invMat[2];
    mat[3] = invMat[3];
    
    return 1; // inverse exists
}

void numNaret(int t, int edge, int *nodes) {
    if (t == 1) { // Quadrangle (4 edges)
        if (edge == 1) {
            nodes[0] = 0; nodes[1] = 1; // Edge 1: between nodes 0 and 1
        }
        else if (edge == 2) {
            nodes[0] = 1; nodes[1] = 2; // Edge 2: between nodes 1 and 2
        }
        else if (edge == 3) {
            nodes[0] = 2; nodes[1] = 3; // Edge 3: between nodes 2 and 3
        }
        else if (edge == 4) {
            nodes[0] = 3; nodes[1] = 0; // Edge 4: between nodes 3 and 0
        }
    }
    else if (t == 2) { // Triangle (3 edges)
        if (edge == 1) {
            nodes[0] = 0; nodes[1] = 1; // Edge 1: between nodes 0 and 1
        }
        else if (edge == 2) {
            nodes[0] = 1; nodes[1] = 2; // Edge 2: between nodes 1 and 2
        }
        else if (edge == 3) {
            nodes[0] = 2; nodes[1] = 0; // Edge 3: between nodes 2 and 0
        }
    }
    else if (t == 3) { // Segment (2 nodes)
        nodes[0] = 0; nodes[1] = 1; // Segment has only one edge: between nodes 0 and 1
    }
}

void selectPts(int nb, int *num, float **coorEns, float **coorSel) {
    for (int i = 0; i < nb; i++) {
        coorSel[i] = coorEns[num[i]-1]; // copy the selected point from coorEns to coorSel
    }
}
