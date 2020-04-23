/**
 *
 * @file codelet_zplghe.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplghe StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c
 *
 */
#include "chameleon_starpu.h"
#include "runtime_codelet_z.h"

/*   MORSE_TASK_zplghe - Generate a tile for random hermitian (positive definite if bump is large enough) matrix. */

void MORSE_TASK_zplghe( const MORSE_option_t *options,
                        double bump, int m, int n, const MORSE_desc_t *A, int Am, int An, int lda,
                        int bigM, int m0, int n0, unsigned long long int seed )
{

    struct starpu_codelet *codelet = &cl_zplghe;
    void (*callback)(void*) = options->profiling ? cl_zplghe_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE, &bump,                   sizeof(double),
        STARPU_VALUE,    &m,                      sizeof(int),
        STARPU_VALUE,    &n,                      sizeof(int),
        STARPU_W,         RTBLKADDR(A, MORSE_Complex64_t, Am, An),
        STARPU_VALUE,  &lda,                      sizeof(int),
        STARPU_VALUE, &bigM,                      sizeof(int),
        STARPU_VALUE,   &m0,                      sizeof(int),
        STARPU_VALUE,   &n0,                      sizeof(int),
        STARPU_VALUE, &seed,   sizeof(unsigned long long int),
        STARPU_PRIORITY,    options->priority,
        STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "zplghe",
#endif
        0);
}

/*   cl_zplghe_cpu_func - Generate a tile for random hermitian (positive definite if bump is large enough) matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zplghe_cpu_func(void *descr[], void *cl_arg)
{
    double bump;
    int m;
    int n;
    MORSE_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &bump, &m, &n, &lda, &bigM, &m0, &n0, &seed );
    CORE_zplghe( bump, m, n, A, lda, bigM, m0, n0, seed );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplghe, 1, cl_zplghe_cpu_func)
