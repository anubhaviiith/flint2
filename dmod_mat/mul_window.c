#include <gmp.h>
#include <stdlib.h>
#include <float.h>
#include "flint.h"
#include "ulong_extras.h"
#include "dmod_vec.h"
#include "dmod_mat.h"

void _dmod_mat_mul_window(dmod_mat_t C, dmod_mat_t A, dmod_mat_t B, slong Cr1, slong Cc1, slong Ar1, slong Ac1, slong Br1, slong Bc1, slong m, slong n, slong k)
{
    slong i, j;

    slong sA, sB, sC;

    /*Starting indices*/
    sA = MATRIX_IDX(A->ld, Ar1, Ac1);
    sB = MATRIX_IDX(B->ld, Br1, Bc1);
    sC = MATRIX_IDX(C->ld, Cr1, Cc1);

 
    cblas_dgemm(101, 111, 111, m, n, k, 1.0, A->rows + sA, A->ld , B->rows + sB, B->ld, 0.0, C->rows + sC, C->ld);
    
    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            C->rows[sC + MATRIX_IDX(C->ld, i, j)] = 
                dmod_reduce(C->rows[sC + MATRIX_IDX(C->ld, i, j)], A->mod);
        }
    }
    
}
