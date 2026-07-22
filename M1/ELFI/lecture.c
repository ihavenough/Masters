/*
SAN Beining
OCHMAN Karolina 
*/

#include <stdio.h>
#include <stdlib.h>

/* ===== Fonctions fournies ===== */
float **alloctab(int dim1, int dim2)
{
    float **ptr = malloc(dim1 * sizeof(float *));
    if (ptr == NULL) return NULL;

    float *tmp = malloc(dim1 * dim2 * sizeof(float));
    if (tmp == NULL) { free(ptr); return NULL; }

    for (int i = 0; i < dim1; i++) {
        ptr[i] = tmp;
        tmp += dim2;
    }
    return ptr;
}

void freetab(void *ptr)
{
    if (ptr == NULL) return;
    void **p = ptr;
    free(p[0]);
    free(p);
}

/* ===== fonction demandée ===== */
int lecfima(char *ficmai,
            int *n, float ***coord,
            int *m, int *t, int *p, int *q,
            int ***ngnel, int ***nRefAr)
{
    FILE *f = fopen(ficmai, "r");
    if (f == NULL) return 1;

    if (fscanf(f, "%d", n) != 1) { fclose(f); return 1; }

    *coord = alloctab(*n, 2);
    if (*coord == NULL) { fclose(f); return 1; }

    for (int i = 0; i < *n; i++)
        fscanf(f, "%f %f", &(*coord)[i][0], &(*coord)[i][1]);

    fscanf(f, "%d %d %d %d", m, t, p, q);

    /* ===== allocation contiguë : ngnel ===== */
    *ngnel = malloc(*m * sizeof(int *));
    if (*ngnel == NULL) { freetab(*coord); fclose(f); return 1; }

    (*ngnel)[0] = malloc((*m) * (*p) * sizeof(int));
    if ((*ngnel)[0] == NULL) {
        free(*ngnel); freetab(*coord); fclose(f); return 1;
    }

    for (int i = 1; i < *m; i++)
        (*ngnel)[i] = (*ngnel)[0] + i * (*p);

    /* ===== allocation contiguë : nRefAr ===== */
    *nRefAr = malloc(*m * sizeof(int *));
    if (*nRefAr == NULL) {
        free((*ngnel)[0]); free(*ngnel);
        freetab(*coord); fclose(f); return 1;
    }

    (*nRefAr)[0] = malloc((*m) * (*q) * sizeof(int));
    if ((*nRefAr)[0] == NULL) {
        free(*nRefAr);
        free((*ngnel)[0]); free(*ngnel);
        freetab(*coord); fclose(f); return 1;
    }

    for (int i = 1; i < *m; i++)
        (*nRefAr)[i] = (*nRefAr)[0] + i * (*q);

    /* ===== lecture des données ===== */
    for (int i = 0; i < *m; i++) {
        for (int j = 0; j < *p; j++)
            fscanf(f, "%d", &(*ngnel)[i][j]);
        for (int j = 0; j < *q; j++)
            fscanf(f, "%d", &(*nRefAr)[i][j]);
    }

    fclose(f);
    return 0;
}

/*
int main(void)
{
    char ficmai[256];
    int n, m, t, p, q;
    float **coord;
    int **ngnel, **nRefAr;

    printf("Nom du fichier maillage :\n");
    scanf("%255s", ficmai);

    if (lecfima(ficmai, &n, &coord, &m, &t, &p, &q, &ngnel, &nRefAr))
        return 1;

    printf("Lecture OK\n");
    printf("n = %d, m = %d, t = %d\n", n, m, t);
    printf("Premier noeud : %.2f %.2f\n", coord[0][0], coord[0][1]);

    //
    free(ngnel[0]);
    free(ngnel);
    free(nRefAr[0]);
    free(nRefAr);
    freetab(coord);

    return 0;
}
*/
