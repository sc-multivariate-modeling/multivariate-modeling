/**
 *
 * @file pzsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

#define A(m,n) A,  m,  n
#define C(m,n) C,  m,  n
/**
 *  Parallel tile symmetric rank-k update - dynamic scheduling
 */
void morse_pzsyrk(MORSE_enum uplo, MORSE_enum trans,
                         MORSE_Complex64_t alpha, MORSE_desc_t *A,
                         MORSE_Complex64_t beta,  MORSE_desc_t *C,
                         MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n, k;
    int ldak, ldam, ldan, ldcm, ldcn;
    int tempnn, tempmm, tempkn, tempkm;

    MORSE_Complex64_t zbeta;
    MORSE_Complex64_t zone = (MORSE_Complex64_t)1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (n = 0; n < C->nt; n++) {
        tempnn = n == C->nt-1 ? C->n-n*C->nb : C->nb;
        ldan = BLKLDD(A, n);
        ldcn = BLKLDD(C, n);
        /*
         *  MorseNoTrans
         */
        if (trans == MorseNoTrans) {
            for (k = 0; k < A->nt; k++) {
                tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                zbeta = k == 0 ? beta : zone;
                MORSE_TASK_zsyrk(
                    &options,
                    uplo, trans,
                    tempnn, tempkn, A->mb,
                    alpha, A(n, k), ldan, /* ldan * K */
                    zbeta, C(n, n), ldcn); /* ldc  * N */
            }
            /*
             *  MorseNoTrans / MorseLower
             */
            if (uplo == MorseLower) {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    ldam = BLKLDD(A, m);
                    ldcm = BLKLDD(C, m);
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        zbeta = k == 0 ? beta : zone;
                        MORSE_TASK_zgemm(
                            &options,
                            trans, MorseTrans,
                            tempmm, tempnn, tempkn, A->mb,
                            alpha, A(m, k), ldam,  /* ldam * K */
                                   A(n, k), ldan,  /* ldan * K */
                            zbeta, C(m, n), ldcm); /* ldc  * N */
                    }
                }
            }
            /*
             *  MorseNoTrans / MorseUpper
             */
            else {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    ldam = BLKLDD(A, m);
                    for (k = 0; k < A->nt; k++) {
                        tempkn = k == A->nt-1 ? A->n-k*A->nb : A->nb;
                        zbeta = k == 0 ? beta : zone;
                        MORSE_TASK_zgemm(
                            &options,
                            trans, MorseTrans,
                            tempnn, tempmm, tempkn, A->mb,
                            alpha, A(n, k), ldan,  /* ldan * K */
                                   A(m, k), ldam,  /* ldam * M */
                            zbeta, C(n, m), ldcn); /* ldc  * M */
                    }
                }
            }
        }
        /*
         *  MorseTrans
         */
        else {
            for (k = 0; k < A->mt; k++) {
                tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                ldak = BLKLDD(A, k);
                zbeta = k == 0 ? beta : zone;
                MORSE_TASK_zsyrk(
                    &options,
                    uplo, trans,
                    tempnn, tempkm, A->mb,
                    alpha, A(k, n), ldak,  /* lda * N */
                    zbeta, C(n, n), ldcn); /* ldc * N */
            }
            /*
             *  MorseTrans / MorseLower
             */
            if (uplo == MorseLower) {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    ldcm = BLKLDD(C, m);
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        MORSE_TASK_zgemm(
                            &options,
                            trans, MorseNoTrans,
                            tempmm, tempnn, tempkm, A->mb,
                            alpha, A(k, m), ldak,  /* lda * M */
                                   A(k, n), ldak,  /* lda * N */
                            zbeta, C(m, n), ldcm); /* ldc * N */
                    }
                }
            }
            /*
             *  MorseTrans / MorseUpper
             */
            else {
                for (m = n+1; m < C->mt; m++) {
                    tempmm = m == C->mt-1 ? C->m-m*C->mb : C->mb;
                    for (k = 0; k < A->mt; k++) {
                        tempkm = k == A->mt-1 ? A->m-k*A->mb : A->mb;
                        ldak = BLKLDD(A, k);
                        zbeta = k == 0 ? beta : zone;
                        MORSE_TASK_zgemm(
                            &options,
                            trans, MorseNoTrans,
                            tempnn, tempmm, tempkm, A->mb,
                            alpha, A(k, n), ldak,  /* lda * K */
                                   A(k, m), ldak,  /* lda * M */
                            zbeta, C(n, m), ldcn); /* ldc * M */
                    }
                }
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
}
