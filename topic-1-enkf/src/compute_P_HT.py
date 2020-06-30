#!/usr/bin/env python3

import sys

import sh
import numpy as np
import scipy.sparse
from scipy.linalg import toeplitz

from C_matrix import C_matrix
from libmdb_matrix import export_full_r, export_rcs, export_sb_toe_r, import_rcs, import_sb_toe_r


if __name__ == '__main__':
    np.random.seed(0)

    #LENKF_PATH = '../lenkf'
    LENKF_PATH = '.'

    # Example with 1 physical dimension
    print('##############################')
    print('# 1 Example with 1 Deimension')
    print('##############################')

    N_list = [100, 250, 500, 1000, 2500, 5000, 10000, 25000]
    H_density_list = [1, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01, 0.005, 0.0025, 0.001]

    #N = 10000
    H_density = 0.05

    L = 10
    M = 32

    P_HT_fname = '/tmp/P_HT'
    e_fname = '/tmp/e'
    H_fname = '/tmp/H'
    C_fname = '/tmp/C'

    #for H_density in H_density_list:
    for N in N_list:
        e = np.random.randn(N, L)
        e.tofile(e_fname)

        #H = np.random.randn(M, N)
        #H.tofile(H_fname)
        H = scipy.sparse.random(M, N, density=H_density, format='csr')
        export_rcs(H_fname, H)

        assert N % 2 == 0
        N2 = int(N/2)
        row0 = np.hstack((np.arange(N2, 0, -1), np.zeros((N2,))))
        blocks = row0

        C = toeplitz(row0)
        #C.tofile(C_fname)
        C_compress = np.concatenate((C[-1], C[0, 1:], np.zeros(1)))
        C_compress.tofile(C_fname)

        # RUN Python CODE
        P_HT = np.multiply(C, e @ e.T) @ H.T / (L - 1)

        # RUN CUDA CODE
        print()
        print('##### 1.2 Running Cuda Code')
        print()
        compute_P_HT_cuda = sh.Command(LENKF_PATH + '/compute_P_HT_cuda')
        print(compute_P_HT_cuda(P_HT_fname,
                            e_fname,
                            H_fname,
                            C_fname,
                            N,
                            L,
                            M,
                            _err_to_out=True,
                            _out=sys.stdout.buffer))
    
        P_HT_cuda = np.fromfile(P_HT_fname).reshape([N, M])
        #P_HT_cuda = import_rcs(P_HT_fname)

        #print(P_HT_cuda[0])
        #print(P_HT[0])

        #print('P_HT_cuda == P_HT? {}'.format(np.allclose(P_HT_cuda.toarray(), P_HT)))
        #for n in range(N):
        #    if not np.allclose(P_HT_cuda[n], P_HT[n]):
        #        print(n, 'P_HT_cuda == P_HT? {}'.format(np.allclose(P_HT_cuda[n], P_HT[n])))
        #        print(P_HT_cuda[n])
        #        print(P_HT[n])
        #        pass
        print('P_HT_cuda == P_HT? {}'.format(np.allclose(P_HT_cuda, P_HT)))
        print()