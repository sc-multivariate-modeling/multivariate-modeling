/**
 *
 * @file codelet_zlansy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlansy Quark codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.6.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
#include "chameleon_quark.h"
#include "chameleon/morse_tasks_z.h"
#include "coreblas/coreblas_z.h"

void CORE_zlansy_quark(Quark *quark)
{
    double *normA;
    int norm;
    MORSE_enum uplo;
    int N;
    MORSE_Complex64_t *A;
    int LDA;
    double *work;

    quark_unpack_args_7(quark, norm, uplo, N, A, LDA, work, normA);
    CORE_zlansy( norm, uplo, N, A, LDA, work, normA);
}

void MORSE_TASK_zlansy(const MORSE_option_t *options,
                       MORSE_enum norm, MORSE_enum uplo, int N, int NB,
                       const MORSE_desc_t *A, int Am, int An, int LDA,
                       const MORSE_desc_t *B, int Bm, int Bn)
{
    quark_option_t *opt = (quark_option_t*)(options->schedopt);
    DAG_CORE_LANSY;
    int szeW = chameleon_max( 1, N );
    QUARK_Insert_Task(
        opt->quark, CORE_zlansy_quark, (Quark_Task_Flags*)opt,
        sizeof(MORSE_enum),              &norm,  VALUE,
        sizeof(MORSE_enum),              &uplo,  VALUE,
        sizeof(int),                     &N,     VALUE,
        sizeof(MORSE_Complex64_t)*NB*NB, RTBLKADDR(A, MORSE_Complex64_t, Am, An), INPUT,
        sizeof(int),                     &LDA,   VALUE,
        sizeof(double)*szeW,             NULL,   SCRATCH,
        sizeof(double),                  RTBLKADDR(B, double, Bm, Bn), OUTPUT,
        0);
}
