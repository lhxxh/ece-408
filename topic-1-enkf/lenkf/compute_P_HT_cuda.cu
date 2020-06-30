#include <stdlib.h>
#include <assert.h>
#include <cuda.h>
#include <stdio.h>
#include <time.h>

//#include "ensemble.h"
//#include "arg_bundle.h"
//#include "lenkf.h"

// #define BLOCK_SIZE 16
#define BLOCK_SIZE 32

__global__
void compute_P_kernel(double *e, double *C, double *P,
                      int L, int N) {

    int col = threadIdx.x + blockIdx.x * blockDim.x;
    int row = threadIdx.y + blockIdx.y * blockDim.y;

    double e_elem = 0;
    // for (int l = 0; l < L; l ++) {
    //     e_elem += e[col * L + l] * e[row * L + l];
    // }

    //////////////////////////////////////////////////// loop unrolling
    int l;
    for (l = 0 ; l < L%2 ; l++)
      e_elem +=  e[col * L + l] * e[row * L + l];

    for(;l < L; l += 2)
    {
      e_elem +=  e[col * L + l] * e[row * L + l];
      e_elem +=  e[col * L + l + 1] * e[row * L + l + 1];
    }

    /////////////////////////////////////////////////////
    if (col < N && row < N) {
        P[row * N + col] = C[row * N + col] * e_elem;
    }
}

__global__
void compute_P_HT_kernel(double *P, double *H, double *P_HT,
                         int M, int N, int L) {

    int col = threadIdx.x + blockIdx.x * blockDim.x;
    int row = threadIdx.y + blockIdx.y * blockDim.y;

    if (row < N && col < M) {
        double Value = 0;
        // for (int i = 0; i < N; i ++) {
        //     Value += P[row * N + i] * H[col * N + i];
        // }
    /////////////////////////////////////////////////////////////  loop unrolling
        int i;
        for(i = 0; i < N % 2; i++)
          Value += P[row * N + i] * H[col * N + i];
        for(;i < N; i += 2)
        {
          Value += P[row * N + i] * H[col * N + i];
          Value += P[row * N + i + 1] * H[col * N + i + 1];
        }
    ////////////////////////////////////////////////////////////
        P_HT[row * M + col] = Value / double(L - 1);
    }
}

void load_matrix(const char *fname,
                 double *buffer,
                 int dim1,
                 int dim2) {
    FILE *fid;
    int r;

    printf("Loading %s into CPU memory...\n", fname);
    fid = fopen(fname, "r");
    assert(fid);

    r = fseek(fid, 0L, SEEK_END);
    int matrix_size = ftell(fid) / sizeof(double);
    assert(matrix_size == (dim1 * dim2));
    r = fseek(fid, 0L, SEEK_SET);

    r = fread(buffer, sizeof(double), matrix_size, fid);
    assert(r == matrix_size);
    printf("Size of %s is: %d\n", fname, matrix_size);
    fclose(fid);
}

void save_matrix(const char *fname,
                 double *buffer,
                 int dim1,
                 int dim2) {
    FILE *fid;
    int r;
    int matrix_size = dim1 * dim2;

    printf("Saving %s into file...\n", fname);
    fid = fopen(fname, "w");
    assert(fid);

    r = fwrite(buffer, sizeof(double), matrix_size, fid);
    assert(r == matrix_size);
    printf("Size of %s is: %d\n", fname, matrix_size);
    fclose(fid);
}

void do_compute_P_HT(const char *P_HT_fname,
                     const char *e_fname,
                     const char *H_fname,
                     const char *C_fname,
                     const char *N_c,
                     const char *L_c,
                     const char *M_c) {
    int N = atoi(N_c);
    int L = atoi(L_c);
    int M = atoi(M_c);
    double *e = (double *) malloc(N * L * sizeof(double));
    double *H = (double *) malloc(M * N * sizeof(double));
    double *C = (double *) malloc(N * N * sizeof(double));
    double *P_HT = (double *) malloc(N * M * sizeof(double));

    clock_t start, finish;

    printf("0. Problem size: N=%d, L=%d, M=%d\n\n", N, L, M);

    printf("1. Load data into CPU memory.\n");
    load_matrix(e_fname, e, N, L);
    load_matrix(H_fname, H, M, N);
    load_matrix(C_fname, C, N, N);

    printf("2. Allocate GPU memory.\n");
    double *e_device;
    double *H_device;
    double *C_device;

    double *P_device;
    double *P_HT_device;

    cudaMalloc((void **) &e_device, N * L * sizeof(double));
    cudaMalloc((void **) &H_device, M * N * sizeof(double));
    cudaMalloc((void **) &C_device, N * N * sizeof(double));

    cudaMalloc((void **) &P_device, N * N * sizeof(double));
    cudaMalloc((void **) &P_HT_device, N * M * sizeof(double));

    printf("3. Write data into GPU memory.\n");
    start = clock();
    cudaMemcpy(e_device, e, N * L * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(H_device, H, M * N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(C_device, C, N * N * sizeof(double), cudaMemcpyHostToDevice);
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("4. Call GPU cuda kernel.\n");
    start = clock();
    dim3 DimGrid;
    dim3 DimBlock;

    DimGrid = dim3(ceil(N / float(BLOCK_SIZE)), ceil(N / float(BLOCK_SIZE)), 1);
    DimBlock = dim3(BLOCK_SIZE, BLOCK_SIZE, 1);
    compute_P_kernel<<<DimGrid, DimBlock>>>(e_device, C_device, P_device, L, N);
    cudaDeviceSynchronize();

    DimGrid = dim3(ceil(M / float(BLOCK_SIZE)), ceil(N / float(BLOCK_SIZE)), 1);
    DimBlock = dim3(BLOCK_SIZE, BLOCK_SIZE, 1);
    compute_P_HT_kernel<<<DimGrid, DimBlock>>>(P_device, H_device, P_HT_device, M, N, L);
    cudaDeviceSynchronize();
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("5. Read results from GPU memory.\n");
    start = clock();
    cudaMemcpy(P_HT, P_HT_device, N * M * sizeof(double), cudaMemcpyDeviceToHost);
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    //for (int i = 0; i < 3; i ++) {
    //    for (int j = 0; j < M; j ++) {
    //        printf("%.2f, ", *(P_HT + i * M + j));
    //    }
    //    printf("\n");
    //}

    printf("6. Save results to file.\n");
    save_matrix(P_HT_fname, P_HT, N, M);

    printf("7. De-allocate CPU and GPU memory.\n");
    cudaFree(e_device);
    cudaFree(H_device);
    cudaFree(C_device);
    cudaFree(P_device);
    cudaFree(P_HT_device);

    free(e);
    free(H);
    free(C);
    free(P_HT);
}

/*
void compute_P_HT(arg_bundle *ab, const sparse_rcs *H,
                  const int row_H, const char name_H,
                  const int n_rows);


void do_compute_P_HT(const char *P_HT_fname,
                     const char *e_fname,
                     const char *H_fname,
                     const char *C_fname) {
    arg_bundle *ab;
    sparse_rcs *H;
    sparse_rcs *P_HT_rcs;
    int rank, i;

    ab = malloc(sizeof(arg_bundle));
    assert(ab);

    fprintf(stderr, "loading e from %s\n", e_fname);
    ab->e = ensemble_import(e_fname);
    fprintf(stderr, "N=%d  L=%d\n", ab->e->N, ab->e->L);

    ab->config = malloc(sizeof(lenkf_config));
    assert(ab->config);
    ab->config->N = ab->e->N;

    ab->P_HT = NULL;

    fprintf(stderr, "\n");
    fprintf(stderr, "loading H from %s\n", H_fname);
    H = sparse_rcs_import(H_fname);
    fprintf(stderr, "m=%d  n=%d  N=%d\n", H->m, H->n, H->N);

    fprintf(stderr, "\n");
    fprintf(stderr, "loading C from %s\n", C_fname);
    ab->C = sb_toe_r_import(C_fname);
    ab->C_it = sb_toe_r_nz_it_create(ab->C);
    rank = ab->C->dim->rank;
    fprintf(stderr, "rank=%d\n", rank);
    fprintf(stderr, "n_phy=[ ");
    for (i = 0; i < rank; i++) {
        fprintf(stderr, "%d ", ab->C->dim->n_phy[i]);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "n    =[ ");
    for (i = 0; i < rank; i++) {
        fprintf(stderr, "%d ", ab->C->dim->n[i]);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "N_phy=[ ");
    for (i = 0; i < rank; i++) {
        fprintf(stderr, "%d ", ab->C->dim->N_phy[i]);
    }
    fprintf(stderr, "]\n");
    fprintf(stderr, "N    =[ ");
    for (i = 0; i < rank; i++) {
        fprintf(stderr, "%d ", ab->C->dim->N[i]);
    }
    fprintf(stderr, "]\n");

    compute_P_HT(ab, H, 0, 'H', H->m);

    P_HT_rcs = sparse_lil_2_rcs(ab->P_HT);
    sparse_rcs_export(P_HT_fname, P_HT_rcs);

    sparse_rcs_destroy(&H);
    sparse_rcs_destroy(&P_HT_rcs);
    ensemble_destroy(&ab->e);
    sb_toe_r_destroy(&ab->C);
    sb_toe_r_nz_it_destroy(&ab->C_it);
    sparse_lil_destroy(&ab->P_HT);
    free(ab->config);
    free(ab);
}
*/

int main(int argc, char **argv) {
    if (argc != 8) {
        fprintf(stderr, "Usage %s <P_HT_fname> <e_fname> <H_fname> <C_fname>\n", argv[0]);
        return 1;
    }

    do_compute_P_HT(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);

    return 0;
}
