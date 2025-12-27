#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapack.h>

// Karolina Ochman
// 26109551

// for matrix A
typedef struct {
    int n;         // number of rows = cols
    float *data; 
    float **row;  
} Matrix;

// for vector b
typedef struct {
    int n;
    float *data;
} Vector;

int main() {

    int sizeb = 1;

    FILE *f = fopen("results_mono.tex", "w");

    fprintf(f, "%% Name: Karolina Ochman\n");
    fprintf(f, "%% Student number: 26109551\n\n");

    //LaTeX table header 
    fprintf(f, "\\begin{tabular}{c c}\n");
    fprintf(f, "\\hline\n");
    fprintf(f, "$n$ & $\\|A\\tilde{x}-b\\|_1/\\|b\\|_1$ \\\\\n");
    fprintf(f, "\\hline\n");

    for (int n = 5; n <= 705; n += 50) {

        // create matrices and vectors 
        Matrix A, A_original;
        Vector b, b_original, Ax, r;

        A.n = n;
        A.data = malloc(n * n * sizeof(float));
        A.row = malloc(n * sizeof(float *));
        for (int i = 0; i < n; i++) A.row[i] = A.data + i*n;

        A_original.n = n;
        A_original.data = malloc(n * n * sizeof(float));
        A_original.row = malloc(n * sizeof(float *));
        for (int i = 0; i < n; i++) A_original.row[i] = A_original.data + i*n;

        b.n = b_original.n = Ax.n = r.n = n;
        b.data = malloc(n * sizeof(float));
        b_original.data = malloc(n * sizeof(float));
        Ax.data = malloc(n * sizeof(float));
        r.data = malloc(n * sizeof(float));

        // initialize A and A_original 
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                if (i == j) A.row[i][j] = 2.0f;
                else if (abs(i - j) == 1) A.row[i][j] = -1.0f;
                else A.row[i][j] = 0.0f;

                A_original.row[i][j] = A.row[i][j];
            }

        // initialize b and b_original 
        float h = 1.0f / (n + 1);
        for (int i = 0; i < n; i++) {
            b.data[i] = h * h;
            b_original.data[i] = h * h;
        }

        // LAPACK call 
        int info;
        int *ipiv = malloc(n * sizeof(int));
        sgesv_(&n, &sizeb, A.data, &n, ipiv, b.data, &n, &info);

        // compute Ax = A_original * b 
        for (int i = 0; i < n; i++) {
            Ax.data[i] = 0.0f;
            for (int j = 0; j < n; j++)
                Ax.data[i] += A_original.row[i][j] * b.data[j];
        }

        // calculate r = Ax - b_original 
        for (int i = 0; i < n; i++)
            r.data[i] = Ax.data[i] - b_original.data[i];

        // calculate norm1 
        float norm_r = 0.0f;
        float norm_b = 0.0f;
        for (int i = 0; i < n; i++) {
            norm_r += fabsf(r.data[i]);
            norm_b += fabsf(b_original.data[i]);
        }

        fprintf(f, "%d & %.8e \\\\\n", n, norm_r / norm_b);

        // free memory 
        free(A.data);
        free(A.row);
        free(A_original.data);
        free(A_original.row);
        free(b.data);
        free(b_original.data);
        free(Ax.data);
        free(r.data);
        free(ipiv);
    }

    fprintf(f, "\\hline\n");
    fprintf(f, "\\end{tabular}\n");
    fclose(f);

    return 0;
}
