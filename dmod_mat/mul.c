#include <gmp.h>
#include <stdlib.h>
#include <float.h>
#include "flint.h"
#include "ulong_extras.h"
#include "dmod_vec.h"
#include "dmod_mat.h"

void _dmod_mat_mul(dmod_mat_t C, dmod_mat_t A, dmod_mat_t B)
{
    slong m, n, k, i, j;
    m = A->nrows;
    n = B->ncols;
    k = A->ncols;

    if (A->ncols != B->nrows)
        return;

    cblas_dgemm(101, 111, 111, m, n, k, 1.0, A->rows, k, B->rows, n, 0.0, C->rows, n);
    for (i = 0; i < C->nrows; i++)    
    {
        for (j = 0; j < C->ncols; j++)
            C->rows[MATRIX_IDX(C->ncols, i, j)] = dmod_reduce(C->rows[MATRIX_IDX(C->ncols, i, j)], A->mod);
    }
    
}
