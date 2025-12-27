#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <lapack.h>

// Karolina Ochman 
// 26109551

// for matrix A
typedef struct {
    int n;         // number of rows = cols
    float *data; 
    float **row;  
} Matrix;

// for matrix (vector) b
typedef struct {
    int n;
    float *data;
} Vector;
int sizeb = 1;

// declaring the matrix
Matrix create_matrix(int n) {
    Matrix A;
    A.n = n;

    A.data = malloc(n * n * sizeof(float));
    A.row = malloc(n * sizeof(float *));

    for (int i = 0; i < n; i++)
        A.row[i] = A.data + i*n;  
    return A;
}

// 'cleaning' the memory
void free_matrix(Matrix *A) {
    free(A->data);
    free(A->row);
}

// declaring the vector
Vector create_vector(int n) {
    Vector b;
    b.n = n;
    b.data = malloc(n * sizeof(float));
    return b;
}

// 'cleaning' the memory
void free_vector(Vector *b) {
    free(b->data);
}

// initialising both, matrix A and vector b
void init_matrix_A(Matrix *A) {
    int n = A->n;

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            if (i == j)
                A->row[i][j] = 2.0f;
            else if (abs(i - j) == 1)
                A->row[i][j] = -1.0f;
            else
                A->row[i][j] = 0.0f;
        }
}

void init_vector_b(Vector *b) {
    int n = b->n;
    float h = 1.0f / (n + 1);

    for (int i = 0; i < n; i++)
        b->data[i] = h*h;
}

// the product of matrix and vector; Ax = y
void matvec(Matrix *A, Vector *x, Vector *y) {
    int n = A->n;

    for (int i = 0; i < n; i++) {
        y->data[i] = 0.0f;
        for (int j = 0; j < n; j++)
            y->data[i] += (A->row[i][j]) * (x->data[j]);
    }
}

// the remainder after the substraction; y - b = r
void vec_sub(Vector *y, Vector *b, Vector *r) {
    int n = y->n;
    for (int i = 0; i < n; i++)
        r->data[i] = y->data[i] - b->data[i];
}

// norm1 
float norm1(Vector *v){
    float sum = 0.0f;
    for (int i = 0; i < v->n; i++)
    sum += fabsf(v->data[i]);
    return sum;
}


// This type of Lapack function was used when run locally on a university computer. It is kept here for the reference.
//void LAPACK_sgesv_(int * n, int * nrhs, float * a, int * lda, int * ipiv, float * b, int * ldb, int * info){
//} 


int main() {

    FILE *f = fopen("results.tex", "w");
    if (f == NULL) {
        perror("Error opening file");
        return 1;
    }

    fprintf(f, "%% Name: Karolina Ochman\n");
    fprintf(f, "%% Student number: 26109551\n\n");
    // LaTeX table header 
    fprintf(f, "\\begin{tabular}{c c}\n");
    fprintf(f, "\\hline\n");
    fprintf(f, "$n$ & $\\|A\\tilde{x}-b\\|_1/\\|b\\|_1$ \\\\\n");
    fprintf(f, "\\hline\n");


    int n;
    for (n = 5; n<= 705; n+=50){

    // Create matrices and vectors
    Matrix A = create_matrix(n);
    Matrix A_original = create_matrix(n);
    Vector b = create_vector(n);
    Vector b_original = create_vector(n);
    Vector Ax = create_vector(n);
    Vector r = create_vector(n);

    init_matrix_A(&A);
    init_matrix_A(&A_original);
    init_vector_b(&b);
    init_vector_b(&b_original);


    // lapack
    int info;
    int *ipiv = malloc(n * sizeof(int));
    sgesv_(&n, &sizeb, A.data, &n, ipiv, b.data, &n, &info);

    matvec(&A_original, &b, &Ax);

    vec_sub(&Ax, &b_original, &r);

    //norm1
    norm1(&r);
    //printf("%d %0.8e\n", n, norm1(&r)/norm1(&b_original));
    fprintf(f, "%d & %.8e \\\\\n", n, norm1(&r)/norm1(&b_original));

    // Free memory
    free_matrix(&A);
    free_matrix(&A_original);
    free_vector(&b);
    free_vector(&b_original);
    free_vector(&Ax);
    free_vector(&r);

    }

    fprintf(f, "\\hline\n");
    fprintf(f, "\\end{tabular}\n");
    fclose(f);


    return 0;
}