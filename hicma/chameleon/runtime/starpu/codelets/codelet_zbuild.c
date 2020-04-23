/**
 *
 * @file codelet_zbuild.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zbuild StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Guillaume Sylvand
 * @date 2016-09-05
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

void MORSE_TASK_zbuild( const MORSE_option_t *options,
                        const MORSE_desc_t *A, int Am, int An, int lda,
                        void *user_data, void* user_build_callback )
{

  struct starpu_codelet *codelet = &cl_zbuild;
  void (*callback)(void*) = options->profiling ? cl_zbuild_callback : NULL;
  int row_min, row_max, col_min, col_max;

  MORSE_BEGIN_ACCESS_DECLARATION;
  MORSE_ACCESS_W(A, Am, An);
  MORSE_END_ACCESS_DECLARATION;

  row_min = Am*A->mb ;
  row_max = Am == A->mt-1 ? A->m-1 : row_min+A->mb-1 ;
  col_min = An*A->nb ;
  col_max = An == A->nt-1 ? A->n-1 : col_min+A->nb-1 ;
  starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &row_min,                      sizeof(int),
        STARPU_VALUE,    &row_max,                      sizeof(int),
        STARPU_VALUE,    &col_min,                      sizeof(int),
        STARPU_VALUE,    &col_max,                      sizeof(int),
        STARPU_W,         RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,    &lda,                          sizeof(int),
        STARPU_VALUE,    &user_data,                    sizeof(void*),
        STARPU_VALUE,    &user_build_callback,          sizeof(void*),
        STARPU_PRIORITY,  options->priority,
        STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zbuild",
#endif
        0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_zbuild_cpu_func(void *descr[], void *cl_arg)
{
  MORSE_Complex64_t *A;
  int ld;
  void *user_data;
  void (*user_build_callback)(int row_min, int row_max, int col_min, int col_max, void *buffer, int ld, void *user_data) ;
  int row_min, row_max, col_min, col_max;

  A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
  starpu_codelet_unpack_args(cl_arg, &row_min, &row_max, &col_min, &col_max, &ld, &user_data, &user_build_callback );

  /* The callback 'user_build_callback' is expected to build the block of matrix [row_min, row_max] x [col_min, col_max]
   * (with both min and max values included in the intervals, index start at 0 like in C, NOT 1 like in Fortran)
   * and store it at the address 'buffer' with leading dimension 'ld'
   */
  user_build_callback(row_min, row_max, col_min, col_max, A, ld, user_data);

}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zbuild, 1, cl_zbuild_cpu_func)
