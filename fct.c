//KRSINAR Xavier - 26102387

#include "header.h"

void mat_create(Matrix *A, int nrow, int ncol) {
    A->nrow = nrow;
    A->ncol = ncol;
    A->mat  = (float *)calloc(nrow * ncol, sizeof(float));
}

void mat_init(Matrix *A, Matrix *B) {
    int n = A->nrow;
    float h = 1.0f / (n + 1);
    for (int j = 0; j < n; j++) {
        A->mat[j * n + j] = 2.0f;
        if (j + 1 < n) A->mat[j * n + (j + 1)] = -1.0f;
        if (j - 1 >= 0) A->mat[j * n + (j - 1)] = -1.0f;
    }
    for (int i = 0; i < n; i++) {
        B->mat[i] = h * h;
    }
}

float norme1(Matrix A) {
    int n = A.nrow;
    float s = 0.0f;
    for (int i = 0; i < n; i++) {
        s += fabsf(A.mat[i]);
    }
    return s;
}


Matrix mat_mul(Matrix A, Matrix B) {
    int n = A.nrow;
    Matrix C;
    mat_create(&C, n, 1);
    for (int i = 0; i < n; i++) {
        float s = 0.0f;
        for (int j = 0; j < n; j++) {
            s += A.mat[i * n + j] * B.mat[j];
        }
        C.mat[i] = s;
    }
    return C;
}

Matrix mat_dif(Matrix A, Matrix B) {
    int n = A.nrow;
    Matrix C;
    mat_create(&C, n, 1);
    for (int i = 0; i < n; i++) {
        C.mat[i] = A.mat[i] - B.mat[i];
    }
    return C;
}


void mat_free(Matrix *A) {
    free(A->mat);
}

