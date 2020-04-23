/**
 *
 * @file codelet_zbuild.c
 *
 * @copyright 2009-2015 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild PaRSEC codelet
 *
 * @version 1.0.0
 * @author Reazul Hoque
 * @author Guillaume Sylvand
 * @date 2016-09-05
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_parsec.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

static inline int
CORE_zbuild_parsec( parsec_execution_stream_t *context,
                    parsec_task_t             *this_task )
{
    MORSE_Complex64_t *A;
    int lda;
    void *user_data;
    void (*user_build_callback)( int row_min, int row_max, int col_min, int col_max,
                                 void *buffer, int ld, void *user_data );
    int row_min, row_max, col_min, col_max;

    parsec_dtd_unpack_args(
        this_task, &row_min, &row_max, &col_min, &col_max, &A, &lda, &user_data, &user_build_callback );

    user_build_callback(row_min, row_max, col_min, col_max, A, lda, user_data);

    (void)context;
    return PARSEC_HOOK_RETURN_DONE;
}

void MORSE_TASK_zbuild( const MORSE_option_t *options,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        void *user_data, void* user_build_callback )
{
    parsec_taskpool_t* PARSEC_dtd_taskpool = (parsec_taskpool_t *)(options->sequence->schedopt);
    int row_min, row_max, col_min, col_max;
    row_min = Am*A->mb ;
    row_max = Am == A->mt-1 ? A->m-1 : row_min+A->mb-1 ;
    col_min = An*A->nb ;
    col_max = An == A->nt-1 ? A->n-1 : col_min+A->nb-1 ;

    parsec_dtd_taskpool_insert_task(
        PARSEC_dtd_taskpool, CORE_zbuild_parsec, options->priority, "zbuild",
        sizeof(int),   &row_min,                          VALUE,
        sizeof(int),   &row_max,                          VALUE,
        sizeof(int),   &col_min,                          VALUE,
        sizeof(int),   &col_max,                          VALUE,
        PASSED_BY_REF,  RTBLKADDR( A, MORSE_Complex64_t, Am, An ), morse_parsec_get_arena_index( A ) | OUTPUT | AFFINITY,
        sizeof(int),   &lda,                              VALUE,
        sizeof(void*), &user_data,                        VALUE,
        sizeof(void*), &user_build_callback,              VALUE,
        PARSEC_DTD_ARG_END );
}
