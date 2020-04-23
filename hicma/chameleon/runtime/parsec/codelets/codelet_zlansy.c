/**
 *
 * @file codelet_zlansy.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlansy PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zlansy_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    int norm;
    MORSE_enum uplo;
    int N;
    MORSE_Complex64_t *A;
    int LDA;
    double *work;
    double *normA;

    parsec_dtd_unpack_args(
        this_task, &norm, &uplo, &N, &A, &LDA, &work, &normA );

    CORE_zlansy( norm, uplo, N, A, LDA, work, normA );

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zlansy(const MORSE_option_t *options,
                       MORSE_enum norm, MORSE_enum uplo, int N, int NB,
                       const MORSE_desc_t *A, int Am, int An, int LDA,
                       const MORSE_desc_t *B, int Bm, int Bn)
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);

    int szeW = chameleon_max( 1, N );

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zlansy_parsec, options->priority, "lansy",
        sizeof(MORSE_enum),            &norm,          VALUE,
        sizeof(MORSE_enum),            &uplo,          VALUE,
        sizeof(int),                   &N,             VALUE,
        PASSED_BY_REF,         RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | INPUT,
        sizeof(int),                   &LDA,           VALUE,
        sizeof(double)*szeW,           NULL,           SCRATCH,
        PASSED_BY_REF,         RTBLKADDR( B, double, Bm, Bn ),            OUTPUT | AFFINITY,
        PARSEC_DTD_ARG_END );

    (void)NB;
}
