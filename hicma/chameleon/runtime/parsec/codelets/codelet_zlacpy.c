/**
 *
 * @file codelet_zlacpy.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

/**
 *
 * @ingroup CORE_MORSE_Complex64_t
 *
 */
static inline int
CORE_zlacpyx_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    MORSE_enum uplo;
    int M;
    int N;
    int displA;
    MORSE_Complex64_t *A;
    int LDA;
    int displB;
    MORSE_Complex64_t *B;
    int LDB;

    parsec_dtd_unpack_args(
        this_task, &uplo, &M, &N, &displA, &A, &LDA, &displB, &B, &LDB );

    CORE_zlacpy( uplo, M, N, A + (displA), LDA, B + (displB), LDB );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zlacpyx( const MORSE_option_t *options,
                         MORSE_enum uplo, int m, int n, int nb,
                         int displA, const MORSE_desc_t *A, int Am, int An, int lda,
                         int displB, const MORSE_desc_t *B, int Bm, int Bn, int ldb )
{

    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlacpyx_parsec, options->priority, "lacpy",
        sizeof(MORSE_enum),    &uplo,                      VALUE,
        sizeof(int),           &m,                         VALUE,
        sizeof(int),           &n,                         VALUE,
        sizeof(int),           &displA,                    VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),           &lda,                       VALUE,
        sizeof(int),           &displB,                    VALUE,
        PASSED_BY_REF,         RTBLKADDR( B, MORSE_Complex64_t, Bm, Bn ), morse_parsec_get_arena_index( B ) | OUTPUT | AFFINITY,
        sizeof(int),           &ldb,                       VALUE,
        PARSEC_DTD_ARG_END );
    (void)nb;
}

void MORSE_TASK_zlacpy(const MORSE_option_t *options,
                       MORSE_enum uplo, int m, int n, int nb,
                       const MORSE_desc_t *A, int Am, int An, int lda,
                       const MORSE_desc_t *B, int Bm, int Bn, int ldb)
{
    MORSE_TASK_zlacpyx( options, uplo, m, n, nb,
                        0, A, Am, An, lda,
                        0, B, Bm, Bn, ldb );
}
