/**
 *
 * @file core_zaxpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2015-09-15
 * @precisions normal z -> c d s
 *
 */
#include "coreblas.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 *  CORE_zaxpy adds to vectors together.
 *
 *       B <- alpha * A  + B
 *
 *******************************************************************************
 *
 * @param[in] M
 *          Number of rows of the vectors A and B.
 *
 * @param[in] alpha
 *          Scalar factor of A.
 *
 * @param[in] A
 *          Vector of size M.
 *
 * @param[in] incA
 *          Offset. incA > 0
 *
 * @param[in,out] B
 *          Vector of size M.
 *
 * @param[in] incB
 *          Offset. incB > 0
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 */

int CORE_zaxpy(int M, MORSE_Complex64_t alpha,
               const MORSE_Complex64_t *A, int incA,
                     MORSE_Complex64_t *B, int incB)
{
    if (M < 0) {
        coreblas_error(1, "Illegal value of M");
        return -1;
    }
    if (incA < 0) {
        coreblas_error(5, "Illegal value of incA");
        return -4;
    }
    if (incB < 0) {
        coreblas_error(5, "Illegal value of incB");
        return -6;
    }

    cblas_zaxpy(M, CBLAS_SADDR(alpha), A, incA, B, incB);

    return MORSE_SUCCESS;
}



