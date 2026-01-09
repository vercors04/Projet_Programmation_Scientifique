//KRSINAR Xavier - 26102387

#include "header.h"

int main() {
  
    FILE *f = fopen("resultats n residu.txt", "w"); // Création de fichiers texte pour l'écriture du résultat
    FILE *g = fopen("resultats n tps execution.txt", "w");

    for (int n = N_MIN; n <= N_MAX; n = n + N_STEP) {

        Matrix A;
        Matrix B;
        Matrix Aprim;
        Matrix Bprim;

        mat_create(&A, n, n);
        mat_create(&B, n, 1);
        mat_create(&Aprim, n, n);
        mat_create(&Bprim, n, 1);

        mat_init(&A, &B);
        mat_init(&Aprim, &Bprim);

        int *ipiv = (int *)malloc(n * sizeof(int)); // Allocation de la mémoire pour le stockage de pivot nécessaire à la fonction LAPACKE_sgesv 
        clock_t debut = clock(); // Mesure le "temps" en nombre de ticks CPU au début du calcul
        
        int info = LAPACKE_sgesv(LAPACK_ROW_MAJOR, n, 1, A.mat, n, ipiv, B.mat, 1); // Le vecteur B sera écrasé par la solution X. 
        if (info != 0){ //info est égal à 0 si le calcul s'est déroulé avec succès.
            printf("Erreur, le calcul a échoué pour n=%d \n", n);
        }
        else {
            clock_t fin = clock(); // Mesure le "temps" en nombre de ticks CPU à la fin du calcul

            double temps = ((double)(fin - debut)) / CLOCKS_PER_SEC; // Mesure le temps d'exécution du calcul avec LAPACKE_sgesv, CLOCKS_PER_SEC représente le nombre de ticks CPU en une seconde

            // Calcul du résidu ||Ax - b||_1 / ||b||_1
            Matrix Prod = mat_mul(Aprim, B); //calcul de Ax, x est mis dans la matrice B par la fonction Lapack
            Matrix Diff = mat_dif(Prod, Bprim); //calcul de Ax - b, Bprim contient le vecteur original B. 
            float res = norme1(Diff) / norme1(Bprim);

            fprintf(f, "\n%d         %e\n", n, res); // Résidu pour chaque itération
            fprintf(g, "\n%d         %e\n", n, temps); // Temps d'exécution en secondes de la résolution Ax=b par LAPACKE_sgesv pour chaque itération

            // Libération de la mémoire
            mat_free(&Prod);  
            mat_free(&Diff);
        }
        
        mat_free(&A); 
        mat_free(&B);
        mat_free(&Aprim);
        mat_free(&Bprim);
        free(ipiv);
    }
}


