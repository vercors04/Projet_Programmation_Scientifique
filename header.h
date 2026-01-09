//KRSINAR Xavier - 26102387

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <lapacke.h>
#include <time.h> //necessaire au calcul du temps d'execution

#define N_MIN  5
#define N_MAX  705
#define N_STEP 50

typedef struct {//Structure pour les matrices et vecteurs. 
    int nrow;
    int ncol;
    float *mat;
} Matrix;


void   mat_create(Matrix *A, int nrow, int ncol);

void   mat_init(Matrix *A, Matrix *B);
float  norme1(Matrix A);
Matrix mat_dif(Matrix A, Matrix B);
Matrix mat_mul(Matrix A, Matrix B);


void   mat_free(Matrix *A);
