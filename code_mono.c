//Krsinar Xavier - 26102387


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <time.h>

#define N_MIN  5
#define N_MAX  705
#define N_STEP 50
typedef struct {
    int nrow;
    int ncol;
    float *mat;
} Matrix;

int main() {

    FILE *f = fopen("resultats n residu.txt", "w"); // Création de fichiers texte pour l'écriture du résultat
    FILE *g = fopen("resultats n tps execution.txt", "w");

    for (int n = N_MIN; n <= N_MAX; n = n + N_STEP) {

        Matrix A;
        Matrix B;
        Matrix Aprim;
        Matrix Bprim;

        A.nrow = n;
        A.ncol = n;
        A.mat = (float *)calloc(n * n, sizeof(float));

        B.nrow = n;
        B.ncol = 1;
        B.mat = (float *)calloc(n * 1, sizeof(float));

        Aprim.nrow = n;
        Aprim.ncol = n;
        Aprim.mat = (float *)calloc(n * n, sizeof(float));

        Bprim.nrow = n;
        Bprim.ncol = 1;
        Bprim.mat = (float *)calloc(n * 1, sizeof(float));

        float h = 1.0f / (n + 1);
        for (int j = 0; j < n; j++) {
            A.mat[j * n + j] = 2.0f;
            if (j + 1 < n) A.mat[j * n + (j + 1)] = -1.0f;
            if (j - 1 >= 0) A.mat[j * n + (j - 1)] = -1.0f;
        }
        for (int i = 0; i < n; i++) {
            B.mat[i] = h * h;
        }

        for (int j = 0; j < n; j++) {
            Aprim.mat[j * n + j] = 2.0f;
            if (j + 1 < n) Aprim.mat[j * n + (j + 1)] = -1.0f;
            if (j - 1 >= 0) Aprim.mat[j * n + (j - 1)] = -1.0f;
        }
        for (int i = 0; i < n; i++) {
            Bprim.mat[i] = h * h;
        }

        int *ipiv = (int *)malloc(n * sizeof(int)); // Allocation de la mémoire pour le stockage de pivot nécessaire à la fonction LAPACKE_sgesv
        clock_t debut = clock(); // Mesure le "temps" en nombre de ticks CPU au début du calcul

        int info = LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, 1, A.mat, n, ipiv, B.mat, 1); // Ici, le vecteur B sera écrasé par la solution X. Le fonctionnement de la fonction sera explicité dans le rendu
        if (info != 0){ //info est égal à 0 si le calcul s'est déroulé avec succès.
            printf("Erreur, le calcul a échoué pour n=%d \n", n);
        }
        else {
            clock_t fin = clock(); // Mesure le "temps" en nombre de ticks CPU à la fin du calcul

            double temps = ((double)(fin - debut)) / CLOCKS_PER_SEC; // Mesure le temps d'exécution du calcul avec LAPACKE_sgesv, CLOCKS_PER_SEC représente le nombre de ticks CPU en une seconde

            // Calcul du résidu ||Ax - b||_1 / ||b||_1

            Matrix Prod;
            Prod.nrow = n;
            Prod.ncol = 1;
            Prod.mat = (float *)calloc(n * 1, sizeof(float));
            for (int i = 0; i < n; i++) {
                float s = 0.0f;
                for (int j = 0; j < n; j++) {
                    s += Aprim.mat[i * n + j] * B.mat[j];
                }
                Prod.mat[i] = s;
            } //calcul de Ax, x est mis dans la matrice B par la fonction Lapack

            Matrix Diff;
            Diff.nrow = n;
            Diff.ncol = 1;
            Diff.mat = (float *)calloc(n * 1, sizeof(float));
            for (int i = 0; i < n; i++) {
                Diff.mat[i] = Prod.mat[i] - Bprim.mat[i];
            } //calcul de Ax - b, Bprim contient le vecteur original B.

            float s_diff = 0.0f;
            for (int i = 0; i < n; i++) {
                s_diff += fabsf(Diff.mat[i]);
            }
            float s_bprim = 0.0f;
            for (int i = 0; i < n; i++) {
                s_bprim += fabsf(Bprim.mat[i]);
            }

            float res = s_diff / s_bprim;

            fprintf(f, "\n%d         %e\n", n, res); // Résidu pour chaque itération
            fprintf(g, "\n%d         %e\n", n, temps); // Temps d'exécution de la résolution Ax=b par LAPACKE_sgesv pour chaque itération

            // Libération de la mémoire
            free(Prod.mat);
            free(Diff.mat);
        }

        free(A.mat);
        free(B.mat);
        free(Aprim.mat);
        free(Bprim.mat);
        free(ipiv);
    }
}
