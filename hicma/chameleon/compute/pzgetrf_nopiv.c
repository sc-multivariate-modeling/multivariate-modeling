/**
 *
 * @file pzgetrf_nopiv.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zgetrf_nopiv parallel algorithm
 *
 * @version 1.0.0
 * @author Omar Zenati
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2013-02-01
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n

/**
 *  Parallel tile LU factorization with no pivoting - dynamic scheduling
 */
void morse_pzgetrf_nopiv(MORSE_desc_t *A,
                                MORSE_sequence_t *sequence,
                                MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n, ib;
    int ldak, ldam;
    int tempkm, tempkn, tempmm, tempnn;

    MORSE_Complex64_t zone  = (MORSE_Complex64_t) 1.0;
    MORSE_Complex64_t mzone = (MORSE_Complex64_t)-1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ib = MORSE_IB;

    for (k = 0; k < chameleon_min(A->mt, A->nt); k++) {
        RUNTIME_iteration_push(morse, k);

        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
        ldak = BLKLDD(A, k);

        options.priority = 2*A->nt - 2*k;
        MORSE_TASK_zgetrf_nopiv(
            &options,
            tempkm, tempkn, ib, A->mb,
            A(k, k), ldak, A->mb*k);

        for (m = k+1; m < A->mt; m++) {
            options.priority = 2*A->nt - 2*k - m;
            tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
            ldam = BLKLDD(A, m);
            MORSE_TASK_ztrsm(
                &options,
                MorseRight, MorseUpper, MorseNoTrans, MorseNonUnit,
                tempmm, tempkn, A->mb,
                zone, A(k, k), ldak,
                      A(m, k), ldam);
        }
        for (n = k+1; n < A->nt; n++) {
            tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
            options.priority = 2*A->nt - 2*k - n;
            MORSE_TASK_ztrsm(
                &options,
                MorseLeft, MorseLower, MorseNoTrans, MorseUnit,
                tempkm, tempnn, A->mb,
                zone, A(k, k), ldak,
                      A(k, n), ldak);

            for (m = k+1; m < A->mt; m++) {
                tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
                options.priority = 2*A->nt - 2*k  - n - m;
                ldam = BLKLDD(A, m);
                MORSE_TASK_zgemm(
                    &options,
                    MorseNoTrans, MorseNoTrans,
                    tempmm, tempnn, A->mb, A->mb,
                    mzone, A(m, k), ldam,
                           A(k, n), ldak,
                    zone,  A(m, n), ldam);
            }
        }

        RUNTIME_iteration_pop(morse);
    }

    RUNTIME_options_finalize(&options, morse);
}
