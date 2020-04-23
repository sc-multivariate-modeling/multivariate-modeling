/**
 *
 * @file core_zlascal.c
 *
 * @copyright 2009-2016 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 * @copyright 2016-2018 KAUST. All rights reserved.
 *
 ***
 *
 * @brief Chameleon core_zlascal CPU kernel
 *
 * @version 1.0.0
 * @author Dalal Sukkari
 * @date 2015-11-05
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"
#include "coreblas/cblas.h"
#include <math.h>

/**
 *******************************************************************************
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zlascal scales a two-dimensional matrix A. As opposite to
 *  CORE_zlascl(), no checks is performed to prevent under/overflow. This should
 *  have been done at higher level.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the shape of A:
 *          = MorseUpperLower: A is a general matrix.
 *          = MorseUpper: A is an upper trapezoidal matrix.
 *          = MorseLower: A is a lower trapezoidal matrix.
 *
 * @param[in] m is the number of rows of the matrix A. m >= 0
 *
 * @param[in] n is the number of columns of the matrix A. n >= 0
 *
 * @param[in] alpha
 *            The scalar factor.
 *
 * @param[in,out] A is the matrix to be multiplied by alpha
 *
 * @param[in] lda is the leading dimension of the array A. lda >= max(1,m).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */
int
CORE_zlascal( MORSE_enum uplo, int m, int n,
              MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int lda )
{
    int i;

    if ( (uplo != MorseUpperLower) &&
         (uplo != MorseUpper)      &&
         (uplo != MorseLower))
    {
        coreblas_error(1, "illegal value of uplo");
        return -1;
    }

    if (m < 0) {
        coreblas_error(2, "Illegal value of m");
        return -2;
    }
    if (n < 0) {
        coreblas_error(3, "Illegal value of n");
        return -3;
    }
    if ( (lda < chameleon_max(1,m)) && (m > 0) ) {
        coreblas_error(6, "Illegal value of lda");
        return -6;
    }

    switch ( uplo ) {
    case MorseUpper:
        for(i=0; i<n; i++) {
            cblas_zscal( chameleon_min( i+1, m ), CBLAS_SADDR(alpha), A+i*lda, 1 );
        }
        break;

    case MorseLower:
        for(i=0; i<n; i++) {
            cblas_zscal( chameleon_max( m, m-i ), CBLAS_SADDR(alpha), A+i*lda, 1 );
        }
        break;
    default:
        if (m == lda) {
            cblas_zscal( m*n, CBLAS_SADDR(alpha), A, 1 );
        }
        else {
            for(i=0; i<n; i++) {
                cblas_zscal( m, CBLAS_SADDR(alpha), A+i*lda, 1 );
            }
        }
    }

    return MORSE_SUCCESS;
}
