/**
 *
 * @file zsyrk.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zsyrk wrappers
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */
#include "control/common.h"

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t
 *
 *  MORSE_zsyrk - Performs one of the hermitian rank k operations
 *
 *    \f[ C = \alpha [ op( A ) \times conjg( op( A )' )] + \beta C \f],
 *
 *  where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = conjg( X' )
 *
 *  where alpha and beta are real scalars, C is an n-by-n hermitian
 *  matrix and A is an n-by-k matrix in the first case and a k-by-n
 *  matrix in the second case.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MorseUpper: Upper triangle of C is stored;
 *          = MorseLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or conjugate transposed:
 *          = MorseNoTrans:   A is not transposed;
 *          = MorseTrans  :   A is transposed.
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ).
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = MorseNoTrans,
 *          and is N otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA must be at least
 *          max( 1, N ) if trans == MorseNoTrans, otherwise LDA must
 *          be at least max( 1, K ).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max( 1, N ).
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zsyrk_Tile
 * @sa MORSE_csyrk
 * @sa MORSE_dsyrk
 * @sa MORSE_ssyrk
 *
 */
int MORSE_zsyrk( MORSE_enum uplo, MORSE_enum trans, int N, int K,
                 MORSE_Complex64_t alpha, MORSE_Complex64_t *A, int LDA,
                 MORSE_Complex64_t beta,  MORSE_Complex64_t *C, int LDC )
{
    int NB;
    int Am, An;
    int status;
    MORSE_desc_t descAl, descAt;
    MORSE_desc_t descCl, descCt;
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zsyrk", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }

    /* Check input arguments */
    if ((uplo != MorseUpper) && (uplo != MorseLower)) {
        morse_error("MORSE_zsyrk", "illegal value of uplo");
        return -1;
    }
    if ((trans != MorseNoTrans) && (trans != MorseTrans)) {
        morse_error("MORSE_zsyrk", "illegal value of trans");
        return -2;
    }
    if ( trans == MorseNoTrans ) {
        Am = N; An = K;
    } else {
        Am = K; An = N;
    }
    if (N < 0) {
        morse_error("MORSE_zsyrk", "illegal value of N");
        return -3;
    }
    if (K < 0) {
        morse_error("MORSE_zsyrk", "illegal value of K");
        return -4;
    }
    if (LDA < chameleon_max(1, Am)) {
        morse_error("MORSE_zsyrk", "illegal value of LDA");
        return -7;
    }
    if (LDC < chameleon_max(1, N)) {
        morse_error("MORSE_zsyrk", "illegal value of LDC");
        return -10;
    }

    /* Quick return */
    if (N == 0 ||
        ((alpha == (MORSE_Complex64_t)0.0 || K == 0.0) && beta == (MORSE_Complex64_t)1.0))
        return MORSE_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNBSIZE */
    status = morse_tune(MORSE_FUNC_ZSYRK, N, K, 0);
    if (status != MORSE_SUCCESS) {
        morse_error("MORSE_zsyrk", "morse_tune() failed");
        return status;
    }

    /* Set MT & NT & KT */
    NB = MORSE_NB;

    morse_sequence_create( morse, &sequence );

    /* Submit the matrix conversion */
    morse_zlap2tile( morse, &descAl, &descAt, MorseDescInput, MorseUpperLower,
                     A, NB, NB, LDA, An, Am, An, sequence, &request );
    morse_zlap2tile( morse, &descCl, &descCt, MorseDescInout, uplo,
                     C, NB, NB, LDC, N, N,  N, sequence, &request );

    /* Call the tile interface */
    MORSE_zsyrk_Tile_Async( uplo, trans, alpha, &descAt, beta, &descCt, sequence, &request );

    /* Submit the matrix conversion back */
    morse_ztile2lap( morse, &descAl, &descAt,
                     MorseDescInput, MorseUpperLower, sequence, &request );
    morse_ztile2lap( morse, &descCl, &descCt,
                     MorseDescInout, uplo, sequence, &request );

    morse_sequence_wait( morse, sequence );

    /* Cleanup the temporary data */
    morse_ztile2lap_cleanup( morse, &descAl, &descAt );
    morse_ztile2lap_cleanup( morse, &descCl, &descCt );

    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile
 *
 *  MORSE_zsyrk_Tile - Performs rank k update.
 *  Tile equivalent of MORSE_zsyrk().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MorseUpper: Upper triangle of C is stored;
 *          = MorseLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          Specifies whether the matrix A is transposed or conjugate transposed:
 *          = MorseNoTrans:   A is not transposed;
 *          = MorseTrans:     A is transposed.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when trans = MorseNoTrans,
 *          and is N otherwise.
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa MORSE_zsyrk_Tile
 * @sa MORSE_csyrk
 * @sa MORSE_dsyrk
 * @sa MORSE_ssyrk
 *
 */
int MORSE_zsyrk_Tile( MORSE_enum uplo, MORSE_enum trans,
                      MORSE_Complex64_t alpha, MORSE_desc_t *A,
                      MORSE_Complex64_t beta,  MORSE_desc_t *C )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zsyrk_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create( morse, &sequence );

    MORSE_zsyrk_Tile_Async( uplo, trans, alpha, A, beta, C, sequence, &request );

    MORSE_Desc_Flush( A, sequence );
    MORSE_Desc_Flush( C, sequence );

    morse_sequence_wait( morse, sequence );
    status = sequence->status;
    morse_sequence_destroy( morse, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup MORSE_Complex64_t_Tile_Async
 *
 *  MORSE_zsyrk_Tile_Async - Performs rank-k update.
 *  Non-blocking equivalent of MORSE_zsyrk_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa MORSE_zsyrk
 * @sa MORSE_zsyrk_Tile
 * @sa MORSE_csyrk_Tile_Async
 * @sa MORSE_dsyrk_Tile_Async
 * @sa MORSE_ssyrk_Tile_Async
 *
 */
int MORSE_zsyrk_Tile_Async( MORSE_enum uplo, MORSE_enum trans,
                            MORSE_Complex64_t alpha, MORSE_desc_t *A,
                            MORSE_Complex64_t beta,  MORSE_desc_t *C,
                            MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    int N, K;
    int Am, An, Amb;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zsyrk_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zsyrk_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zsyrk_Tile_Async", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS) {
        request->status = MORSE_SUCCESS;
    }
    else {
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_zsyrk_Tile_Async", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(C) != MORSE_SUCCESS) {
        morse_error("MORSE_zsyrk_Tile_Async", "invalid third descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((uplo != MorseUpper) && (uplo != MorseLower)) {
        morse_error("MORSE_zsyrk", "illegal value of uplo");
        return morse_request_fail(sequence, request, -1);
    }
    if ((trans != MorseNoTrans) && (trans != MorseTrans)) {
        morse_error("MORSE_zsyrk", "illegal value of transA");
        return morse_request_fail(sequence, request, -2);
    }

    if ( trans == MorseNoTrans ) {
        Am  = A->m;
        An  = A->n;
        Amb = A->mb;
    } else {
        Am  = A->n;
        An  = A->m;
        Amb = A->nb;
    }

    if (C->mb != C->nb) {
        morse_error("MORSE_zsyrk_Tile_Async", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (Amb != C->mb) {
        morse_error("MORSE_zsyrk_Tile_Async", "tile sizes have to match");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (C->m != C->n) {
        morse_error("MORSE_zsyrk_Tile_Async", "only square matrix C is supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (Am != C->m) {
        morse_error("MORSE_zsyrk_Tile_Async", "sizes of matrices have to match");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    N = C->m;
    K = An;

    /* Quick return */
    if ( N == 0 ||
         ((alpha == (MORSE_Complex64_t)0.0 || K == 0) && beta == (MORSE_Complex64_t)1.0))
        return MORSE_SUCCESS;

    morse_pzsyrk( uplo, trans, alpha, A, beta, C, sequence, request );

    return MORSE_SUCCESS;
}
