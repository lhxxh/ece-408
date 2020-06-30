#include <stdlib.h>
#include <assert.h>
#include <cuda.h>
#include <stdio.h>
#include <time.h>

#define TILE_SIZE 32
//#define MAX_L 16

__global__
void C_e_H_optim_compute_P_HT_kernel(double *e, double *C, double *H,
                                 int *H_indices, int *H_indptr,
                                 double *P_HT, int M, int N, int L, int n) {
    __shared__ double P_tile[TILE_SIZE][TILE_SIZE];
    //__shared__ double H_tile[TILE_SIZE][TILE_SIZE];

    int tile_row = threadIdx.y;
    int tile_col = threadIdx.x;

    int P_HT_row = blockIdx.y * blockDim.y + tile_row;
    int P_HT_col = blockIdx.x * blockDim.x + tile_col;

    int H_indptr_start = H_indptr[P_HT_col];
    int H_indptr_end = H_indptr[P_HT_col + 1];

    double P_HT_value = 0;
    for (int n = 0; n < N; n += TILE_SIZE) {
        int P_row = P_HT_row;
        int P_col = n + tile_col;
        //int H_row = P_HT_col;
        //int H_col = n + tile_row;

        int e_row_idx = P_row * L; 
        int e_col_idx = P_col * L;

        if (P_row < N && P_col < N) {
            double e_eT_value = 0;
            for (int l = 0; l < L; l ++) {
                e_eT_value += e[e_row_idx + l] * e[e_col_idx + l];
            }

            int diff_n = P_col / n * n - P_row / n * n;
            int C_idx = diff_n * 2 + n - 1 + P_col % n - P_row % n;
            P_tile[tile_row][tile_col] = C[C_idx] * e_eT_value;
        } else {
            P_tile[tile_row][tile_col] = 0;
        }

        __syncthreads();

        //if (H_row < M && H_col < N) {
        //    H_tile[tile_row][tile_col] = H[H_row * N + H_col];
        //} else {
        //    H_tile[tile_row][tile_col] = 0;
        //}

        //__syncthreads();
        
        for (int t = H_indptr_start; t < H_indptr_end; t ++) {
            int H_tile_col = H_indices[t] - n;
            if (H_tile_col < TILE_SIZE) {
                P_HT_value += P_tile[tile_row][H_tile_col] * H[t];
                H_indptr_start += 1;
            } else {
                break;
            }
        }

        //for (int t = 0; t < TILE_SIZE; t ++) {
        //    P_HT_value += P_tile[tile_row][t] * H_tile[t][tile_col];
        //}
    }

    if (P_HT_row < N && P_HT_col < M) {
        P_HT[P_HT_row * M + P_HT_col] = P_HT_value / (L - 1);
    }
}

/*
__global__
void C_e_optim_compute_P_HT_kernel(double *e, double *C, double *H,
                                 double *P_HT, int M, int N, int L) {
    __shared__ double P_tile[TILE_SIZE][TILE_SIZE];
    __shared__ double H_tile[TILE_SIZE][TILE_SIZE];

    int tile_row = threadIdx.y;
    int tile_col = threadIdx.x;

    int P_HT_row = blockIdx.y * blockDim.y + tile_row;
    int P_HT_col = blockIdx.x * blockDim.x + tile_col;

    double P_HT_value = 0;
    for (int n = 0; n < N; n += TILE_SIZE) {
        int P_row = P_HT_row;
        int P_col = n + tile_col;
        int H_row = P_HT_col;
        int H_col = n + tile_row;

        int e_row_idx = P_row * L; 
        int e_col_idx = P_col * L;

        if (P_row < N && P_col < N) {
            double e_eT_value = 0;
            for (int l = 0; l < L; l ++) {
                e_eT_value += e[e_row_idx + l] * e[e_col_idx + l];
            }

            P_tile[tile_row][tile_col] = C[N - 1 + P_col - P_row] * e_eT_value;
        } else {
            P_tile[tile_row][tile_col] = 0;
        }

        __syncthreads();

        if (H_row < M && H_col < N) {
            H_tile[tile_row][tile_col] = H[H_row * N + H_col];
        } else {
            H_tile[tile_row][tile_col] = 0;
        }

        __syncthreads();

        for (int t = 0; t < TILE_SIZE; t ++) {
            P_HT_value += P_tile[tile_row][t] * H_tile[t][tile_col];
        }
    }

    if (P_HT_row < N && P_HT_col < M) {
        P_HT[P_HT_row * M + P_HT_col] = P_HT_value / (L - 1);
    }
}

__global__
void C_optim_compute_P_HT_kernel(double *e, double *C, double *H,
                                   double *P_HT, int M, int N, int L) {
    __shared__ double P_tile[TILE_SIZE][TILE_SIZE];
    __shared__ double H_tile[TILE_SIZE][TILE_SIZE];
    __shared__ double e_row_tile[TILE_SIZE][MAX_L];
    __shared__ double e_col_tile[TILE_SIZE][MAX_L];

    int tile_row = threadIdx.y;
    int tile_col = threadIdx.x;

    int P_HT_row = blockIdx.y * blockDim.y + tile_row;
    int P_HT_col = blockIdx.x * blockDim.x + tile_col;

    double P_HT_value = 0;
    for (int n = 0; n < N; n += TILE_SIZE) {
        int P_row = P_HT_row;
        int P_col = n + tile_col;
        int H_row = P_HT_col;
        int H_col = n + tile_row;

        int e_row_idx = P_row * L; 
        for (int l = 0; l < L; l += TILE_SIZE) {
            e_row_tile[tile_row][l + tile_col] = e[e_row_idx + l + tile_col];
        }

        int e_col_idx = P_col * L;
        for (int l = 0; l < L; l += TILE_SIZE) {
            e_col_tile[tile_col][l + tile_row] = e[e_col_idx + l + tile_row];
        }

        __syncthreads();

        if (P_row < N && P_col < N) {
            double e_eT_value = 0;
            for (int l = 0; l < L; l ++) {
                e_eT_value += e_row_tile[tile_row][l] * e_col_tile[tile_col][l];
            }

            P_tile[tile_row][tile_col] = C[N - 1 + P_col - P_row] * e_eT_value;
        } else {
            P_tile[tile_row][tile_col] = 0;
        }

        if (H_row < M && H_col < N) {
            H_tile[tile_row][tile_col] = H[H_row * N + H_col];
        } else {
            H_tile[tile_row][tile_col] = 0;
        }

        __syncthreads();

        for (int t = 0; t < TILE_SIZE; t ++) {
            P_HT_value += P_tile[tile_row][t] * H_tile[t][tile_col];
        }
    }

    if (P_HT_row < N && P_HT_col < M) {
        P_HT[P_HT_row * M + P_HT_col] = P_HT_value / (L - 1);
    }
}

__global__
void compute_P_kernel(double *e, double *C, double *P,
                      int L, int N) {

    int col = threadIdx.x + blockIdx.x * blockDim.x;
    int row = threadIdx.y + blockIdx.y * blockDim.y;
    
    double e_elem = 0;
    for (int l = 0; l < L; l ++) {
        e_elem += e[col * L + l] * e[row * L + l];
    }

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
        for (int i = 0; i < N; i ++) {
            Value += P[row * N + i] * H[col * N + i];
        }
        P_HT[row * M + col] = Value / double(L - 1);
    }
}
*/

int load_sparse_matrix_nnz(const char *fname,
                           int dim1,
                           int dim2) {
    FILE *fid;
    int sizeof_elem;
    int m, n, nnz;
    int r;

    printf("Loading %s into CPU memory...\n", fname);
    fid = fopen(fname, "r");
    assert(fid);

    r = fread(&sizeof_elem, sizeof(int), 1, fid);
    assert(sizeof_elem == sizeof(double));
    r = fread(&m, sizeof(int), 1, fid);
    assert(m == dim1);
    r = fread(&n, sizeof(int), 1, fid);
    assert(n == dim2);
    r = fread(&nnz, sizeof(int), 1, fid);
    assert(r == 1);
    fclose(fid);
    
    return nnz;
}

void load_sparse_matrix(const char *fname,
                        double *buffer,
                        int *indices,
                        int *indptr,
                        int dim1,
                        int dim2) {
    FILE *fid;
    int sizeof_elem;
    int m, n, nnz;
    int r;

    printf("Loading %s into CPU memory...\n", fname);
    fid = fopen(fname, "r");
    assert(fid);

    r = fread(&sizeof_elem, sizeof(int), 1, fid);
    assert(sizeof_elem == sizeof(double));
    r = fread(&m, sizeof(int), 1, fid);
    assert(m == dim1);
    r = fread(&n, sizeof(int), 1, fid);
    assert(n == dim2);
    r = fread(&nnz, sizeof(int), 1, fid);
    assert(r == 1);

    r = fread(buffer, sizeof(double), nnz, fid);
    assert(r == nnz);
    r = fread(indices, sizeof(int), nnz, fid);
    assert(r == nnz);
    r = fread(indptr, sizeof(int), m + 1, fid);
    assert(r == m + 1);

    printf("Size of %s is: %d\n", fname, nnz);
    fclose(fid);
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

/*
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

    DimGrid = dim3(ceil(N / float(TILE_SIZE)), ceil(N / float(TILE_SIZE)), 1);
    DimBlock = dim3(TILE_SIZE, TILE_SIZE, 1);
    compute_P_kernel<<<DimGrid, DimBlock>>>(e_device, C_device, P_device, L, N);
    cudaDeviceSynchronize();

    DimGrid = dim3(ceil(M / float(TILE_SIZE)), ceil(N / float(TILE_SIZE)), 1);
    DimBlock = dim3(TILE_SIZE, TILE_SIZE, 1);
    compute_P_HT_kernel<<<DimGrid, DimBlock>>>(P_device, H_device, P_HT_device, M, N, L);
    cudaDeviceSynchronize();
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("5. Read results from GPU memory.\n");
    start = clock();
    cudaMemcpy(P_HT, P_HT_device, N * M * sizeof(double), cudaMemcpyDeviceToHost);
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

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

void do_optim_compute_P_HT(const char *P_HT_fname,
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
    double *C = (double *) malloc(2 * N * sizeof(double));
    double *P_HT = (double *) malloc(N * M * sizeof(double));

    clock_t start, finish;

    printf("0. Problem size: N=%d, L=%d, M=%d\n\n", N, L, M);

    printf("1. Load data into CPU memory.\n");
    load_matrix(e_fname, e, N, L);
    load_matrix(H_fname, H, M, N);
    load_matrix(C_fname, C, 2, N);

    printf("2. Allocate GPU memory.\n");
    double *e_device;
    double *H_device;
    double *C_device;
    double *P_HT_device;

    cudaMalloc((void **) &e_device, N * L * sizeof(double));
    cudaMalloc((void **) &H_device, M * N * sizeof(double));
    cudaMalloc((void **) &C_device, 2 * N * sizeof(double));
    cudaMalloc((void **) &P_HT_device, N * M * sizeof(double));

    printf("3. Write data into GPU memory.\n");
    start = clock();
    cudaMemcpy(e_device, e, N * L * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(H_device, H, M * N * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(C_device, C, 2 * N * sizeof(double), cudaMemcpyHostToDevice);
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("4. Call GPU cuda kernel.\n");
    start = clock();
    dim3 DimGrid;
    dim3 DimBlock;

    DimGrid = dim3(ceil(N / float(TILE_SIZE)), ceil(N / float(TILE_SIZE)), 1);
    DimBlock = dim3(TILE_SIZE, TILE_SIZE, 1);
    C_e_optim_compute_P_HT_kernel<<<DimGrid, DimBlock>>>(e_device, C_device, H_device, P_HT_device, M, N, L);
    cudaDeviceSynchronize();
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("5. Read results from GPU memory.\n");
    start = clock();
    cudaMemcpy(P_HT, P_HT_device, N * M * sizeof(double), cudaMemcpyDeviceToHost);
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("6. Save results to file.\n");
    save_matrix(P_HT_fname, P_HT, N, M);

    printf("7. De-allocate CPU and GPU memory.\n");
    cudaFree(e_device);
    cudaFree(H_device);
    cudaFree(C_device);
    cudaFree(P_HT_device);

    free(e);
    free(H);
    free(C);
    free(P_HT);
}
*/

void do_C_e_H_optim_compute_P_HT(const char *P_HT_fname,
                     const char *e_fname,
                     const char *H_fname,
                     const char *C_fname,
                     const char *N_c,
                     const char *L_c,
                     const char *M_c,
                     const char *n_c) {
    int N = atoi(N_c);
    int L = atoi(L_c);
    int M = atoi(M_c);
    int n = atoi(n_c);
    double *e = (double *) malloc(N * L * sizeof(double));
    double *C = (double *) malloc(2 * N * sizeof(double));
    double *P_HT = (double *) malloc(N * M * sizeof(double));

    clock_t start, finish;

    printf("0. Problem size: N=%d, L=%d, M=%d\n\n", N, L, M);

    printf("1. Load data into CPU memory.\n");
    int H_nnz = load_sparse_matrix_nnz(H_fname, M, N);
    double *H = (double *) malloc(H_nnz * sizeof(double));
    int *H_indices = (int *) malloc(H_nnz * sizeof(int));
    int *H_indptr = (int *) malloc((M + 1) * sizeof(int));

    load_matrix(e_fname, e, N, L);
    load_sparse_matrix(H_fname, H, H_indices, H_indptr, M, N);
    load_matrix(C_fname, C, 2, N);

    printf("2. Allocate GPU memory.\n");
    double *e_device;
    double *H_device;
    int *H_indices_device;
    int *H_indptr_device;
    double *C_device;
    double *P_HT_device;

    cudaMalloc((void **) &e_device, N * L * sizeof(double));
    cudaMalloc((void **) &H_device, H_nnz * sizeof(double));
    cudaMalloc((void **) &H_indices_device, H_nnz * sizeof(int));
    cudaMalloc((void **) &H_indptr_device, (M + 1) * sizeof(int));
    cudaMalloc((void **) &C_device, 2 * N * sizeof(double));
    cudaMalloc((void **) &P_HT_device, N * M * sizeof(double));

    printf("3. Write data into GPU memory.\n");
    start = clock();
    cudaMemcpy(e_device, e, N * L * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(H_device, H, H_nnz * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(H_indices_device, H_indices, H_nnz * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(H_indptr_device, H_indptr, (M + 1) * sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(C_device, C, 2 * N * sizeof(double), cudaMemcpyHostToDevice);
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("4. Call GPU cuda kernel.\n");
    start = clock();
    dim3 DimGrid;
    dim3 DimBlock;

    DimGrid = dim3(ceil(M / float(TILE_SIZE)), ceil(N / float(TILE_SIZE)), 1);
    DimBlock = dim3(TILE_SIZE, TILE_SIZE, 1);
    C_e_H_optim_compute_P_HT_kernel<<<DimGrid, DimBlock>>>(e_device, C_device, H_device, H_indices_device, H_indptr_device, P_HT_device, M, N, L, n);
    cudaDeviceSynchronize();
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("5. Read results from GPU memory.\n");
    start = clock();
    cudaMemcpy(P_HT, P_HT_device, N * M * sizeof(double), cudaMemcpyDeviceToHost);
    finish = clock();
    printf("Latency: %fms\n", (double)(finish - start) * 1000 / CLOCKS_PER_SEC);

    printf("6. Save results to file.\n");
    save_matrix(P_HT_fname, P_HT, N, M);

    printf("7. De-allocate CPU and GPU memory.\n");
    cudaFree(e_device);
    cudaFree(H_device);
    cudaFree(H_indices_device);
    cudaFree(H_indptr_device);
    cudaFree(C_device);
    cudaFree(P_HT_device);

    free(e);
    free(H);
    free(H_indices);
    free(H_indptr);
    free(C);
    free(P_HT);
}

int main(int argc, char **argv) {
    //if (argc != 8) {
    //    fprintf(stderr, "Usage %s <P_HT_fname> <e_fname> <H_fname> <C_fname>\n", argv[0]);
    //    return 1;
    //}

    //do_compute_P_HT(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
    //do_optim_compute_P_HT(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7]);
    do_C_e_H_optim_compute_P_HT(argv[1], argv[2], argv[3], argv[4], argv[5], argv[6], argv[7], argv[8]);
    return 0;
}
