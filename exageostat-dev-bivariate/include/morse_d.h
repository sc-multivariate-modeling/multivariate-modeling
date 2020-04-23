/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file morse_d.h
 *
 *  MORSE header file for double routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Azzam Haidar
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @generated d Mon Oct 29 16:30:48 2018
 *
 **/
#ifndef _MORSE_D_H_
#define _MORSE_D_H_

#undef COMPLEX
#define REAL

BEGIN_C_DECLS

/** ****************************************************************************
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 **/
int MORSE_dgeadd(MORSE_enum trans, int M, int N, double alpha, double *A, int LDA, double beta, double *B, int LDB);
//int MORSE_dgebrd(int M, int N, double *A, int LDA, double *D, double *E, MORSE_desc_t *descT);
//int MORSE_dgecon(MORSE_enum norm, int N, double *A, int LDA, double anorm, double *rcond);
//int MORSE_dpocon(MORSE_enum uplo, int N, double *A, int LDA, double anorm, double *rcond);
int MORSE_dgelqf(int M, int N, double *A, int LDA, MORSE_desc_t *descT);
int MORSE_dgelqs(int M, int N, int NRHS, double *A, int LDA, MORSE_desc_t *descT, double *B, int LDB);
int MORSE_dgels(MORSE_enum trans, int M, int N, int NRHS, double *A, int LDA, MORSE_desc_t *descT, double *B, int LDB);
int MORSE_dgemm(MORSE_enum transA, MORSE_enum transB, int M, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
int MORSE_dgeqrf(int M, int N, double *A, int LDA, MORSE_desc_t *descT);
int MORSE_dgeqrs(int M, int N, int NRHS, double *A, int LDA, MORSE_desc_t *descT, double *B, int LDB);
//int MORSE_dgesv(int N, int NRHS, double *A, int LDA, int *IPIV, double *B, int LDB);
int MORSE_dgesv_incpiv(int N, int NRHS, double *A, int LDA, MORSE_desc_t *descL, int *IPIV, double *B, int LDB);
int MORSE_dgesv_nopiv(int N, int NRHS, double *A, int LDA, double *B, int LDB);
int MORSE_dgesvd(MORSE_enum jobu, MORSE_enum jobvt, int M, int N, double *A, int LDA, double *S, MORSE_desc_t *descT, double *U, int LDU, double *VT, int LDVT);
//int MORSE_dgetrf(int M, int N, double *A, int LDA, int *IPIV);
int MORSE_dgetrf_incpiv(int M, int N, double *A, int LDA, MORSE_desc_t *descL, int *IPIV);
int MORSE_dgetrf_nopiv(int M, int N, double *A, int LDA);
//int MORSE_dgetri(int N, double *A, int LDA, int *IPIV);
//int MORSE_dgetrs(MORSE_enum trans, int N, int NRHS, double *A, int LDA, int *IPIV, double *B, int LDB);
int MORSE_dgetrs_incpiv(MORSE_enum trans, int N, int NRHS, double *A, int LDA, MORSE_desc_t *descL, int *IPIV, double *B, int LDB);
int MORSE_dgetrs_nopiv(MORSE_enum trans, int N, int NRHS, double *A, int LDA, double *B, int LDB);
#ifdef COMPLEX
int MORSE_dsymm(MORSE_enum side, MORSE_enum uplo, int M, int N, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
int MORSE_dsyrk(MORSE_enum uplo, MORSE_enum trans, int N, int K, double alpha, double *A, int LDA, double beta, double *C, int LDC);
int MORSE_dsyr2k(MORSE_enum uplo, MORSE_enum trans, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
#endif
//int MORSE_dsyev(MORSE_enum jobz, MORSE_enum uplo, int N, double *A, int LDA, double *W, MORSE_desc_t *descT, double *Q, int LDQ);
int MORSE_dsyevd(MORSE_enum jobz, MORSE_enum uplo, int N, double *A, int LDA, double *W, MORSE_desc_t *descT);
//int MORSE_dhegv(MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, int N, double *A, int LDA, double *B, int LDB, double *W, MORSE_desc_t *descT, double *Q, int LDQ);
//int MORSE_dsygvd(MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, int N, double *A, int LDA, double *B, int LDB, double *W, MORSE_desc_t *descT, double *Q, int LDQ);
//int MORSE_dsygst(MORSE_enum itype, MORSE_enum uplo, int N, double *A, int LDA, double *B, int LDB);
int MORSE_dsytrd(MORSE_enum jobz, MORSE_enum uplo, int N, double *A, int LDA, double *D, double *E, MORSE_desc_t *descT, double *Q, int LDQ);
int MORSE_dlacpy(MORSE_enum uplo, int M, int N, double *A, int LDA, double *B, int LDB);
double MORSE_dlange(MORSE_enum norm, int M, int N, double *A, int LDA);
#ifdef COMPLEX
double MORSE_dlansy(MORSE_enum norm, MORSE_enum uplo, int N, double *A, int LDA);
#endif
double MORSE_dlansy(MORSE_enum norm, MORSE_enum uplo, int N, double *A, int LDA);
double MORSE_dlantr(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, int M, int N, double *A, int LDA);
int MORSE_dlascal(MORSE_enum uplo, int M, int N, double alpha, double *A, int LDA);
int MORSE_dlaset(MORSE_enum uplo, int M, int N, double alpha, double beta, double *A, int LDA);
//int MORSE_dlaswp(int N, double *A, int LDA, int K1, int K2, int *IPIV, int INCX);
//int MORSE_dlaswpc(int N, double *A, int LDA, int K1, int K2, int *IPIV, int INCX);
int MORSE_dlauum(MORSE_enum uplo, int N, double *A, int LDA);
#ifdef COMPLEX
int MORSE_dplgsy( double bump, MORSE_enum uplo, int N, double *A, int LDA, unsigned long long int seed );
#endif
int MORSE_dplgsy( double bump, MORSE_enum uplo, int N, double *A, int LDA, unsigned long long int seed );
int MORSE_dplrnt( int M, int N, double *A, int LDA, unsigned long long int seed );
int MORSE_dposv(MORSE_enum uplo, int N, int NRHS, double *A, int LDA, double *B, int LDB);
int MORSE_dpotrf(MORSE_enum uplo, int N, double *A, int LDA);
int MORSE_dsytrf(MORSE_enum uplo, int N, double *A, int LDA);
int MORSE_dpotri(MORSE_enum uplo, int N, double *A, int LDA);
int MORSE_dpotrimm(MORSE_enum uplo, int N, double *A, int LDA, double *B, int LDB, double *C, int LDC);
int MORSE_dpotrs(MORSE_enum uplo, int N, int NRHS, double *A, int LDA, double *B, int LDB);
int MORSE_dsymm(MORSE_enum side, MORSE_enum uplo, int M, int N, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
int MORSE_dsyrk(MORSE_enum uplo, MORSE_enum trans, int N, int K, double alpha, double *A, int LDA, double beta, double *C, int LDC);
int MORSE_dsyr2k(MORSE_enum uplo, MORSE_enum trans, int N, int K, double alpha, double *A, int LDA, double *B, int LDB, double beta, double *C, int LDC);
int MORSE_dsysv(MORSE_enum uplo, int N, int NRHS, double *A, int LDA, double *B, int LDB);
int MORSE_dsytrs(MORSE_enum uplo, int N, int NRHS, double *A, int LDA, double *B, int LDB);
int MORSE_dtpgqrt( int M, int N, int K, int L, double *V1, int LDV1, MORSE_desc_t *descT1, double *V2, int LDV2, MORSE_desc_t *descT2, double *Q1, int LDQ1, double *Q2, int LDQ2 );
int MORSE_dtpqrt( int M, int N, int L, double *A, int LDA, double *B, int LDB, MORSE_desc_t *descT );
int MORSE_dtradd(MORSE_enum uplo, MORSE_enum trans, int M, int N, double alpha, double *A, int LDA, double beta, double *B, int LDB);
int MORSE_dtrmm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, int N, int NRHS, double alpha, double *A, int LDA, double *B, int LDB);
int MORSE_dtrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, int N, int NRHS, double alpha, double *A, int LDA, double *B, int LDB);
int MORSE_dtrsmpl(int N, int NRHS, double *A, int LDA, MORSE_desc_t *descL, int *IPIV, double *B, int LDB);
int MORSE_dtrsmrv(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, int N, int NRHS, double alpha, double *A, int LDA, double *B, int LDB);
int MORSE_dtrtri(MORSE_enum uplo, MORSE_enum diag, int N, double *A, int LDA);
int MORSE_dorglq(int M, int N, int K, double *A, int LDA, MORSE_desc_t *descT, double *B, int LDB);
int MORSE_dorgqr(int M, int N, int K, double *A, int LDA, MORSE_desc_t *descT, double *B, int LDB);
int MORSE_dormlq(MORSE_enum side, MORSE_enum trans, int M, int N, int K, double *A, int LDA, MORSE_desc_t *descT, double *B, int LDB);
int MORSE_dormqr(MORSE_enum side, MORSE_enum trans, int M, int N, int K, double *A, int LDA, MORSE_desc_t *descT, double *B, int LDB);

//int MORSE_dgecfi(int m, int n, double *A, MORSE_enum fin, int imb, int inb, MORSE_enum fout, int omb, int onb);
//int MORSE_dgetmi(int m, int n, double *A, MORSE_enum fin, int mb,  int nb);

/** ****************************************************************************
 *  Declarations of math functions (tile layout) - alphabetical order
 **/
int MORSE_dgeadd_Tile(MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *B);
//int MORSE_dgebrd_Tile(MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T);
//int MORSE_dgecon_Tile(MORSE_enum norm, MORSE_desc_t *A, double anorm, double *rcond);
//int MORSE_dpocon_Tile(MORSE_enum uplo, MORSE_desc_t *A, double anorm, double *rcond);
int MORSE_dgelqf_Tile(MORSE_desc_t *A, MORSE_desc_t *T);
int MORSE_dgelqs_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_dgels_Tile(MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_dgemm_Tile(MORSE_enum transA, MORSE_enum transB, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C);
int MORSE_dgeqrf_Tile(MORSE_desc_t *A, MORSE_desc_t *T);
int MORSE_dgeqrs_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
//int MORSE_dgesv_Tile(MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B);
int MORSE_dgesv_incpiv_Tile(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B);
int MORSE_dgesv_nopiv_Tile(MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dgesvd_Tile(MORSE_enum jobu, MORSE_enum jobvt, MORSE_desc_t *A, double *S, MORSE_desc_t *T, double *U, int LDU, double *VT, int LDVT);
//int MORSE_dgetrf_Tile(MORSE_desc_t *A, int *IPIV);
int MORSE_dgetrf_incpiv_Tile(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV);
int MORSE_dgetrf_nopiv_Tile(MORSE_desc_t *A);
//int MORSE_dgetri_Tile(MORSE_desc_t *A, int *IPIV);
//int MORSE_dgetrs_Tile(MORSE_enum trans, MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B);
int MORSE_dgetrs_incpiv_Tile(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B);
int MORSE_dgetrs_nopiv_Tile(MORSE_desc_t *A, MORSE_desc_t *B);
#ifdef COMPLEX
int MORSE_dsymm_Tile(MORSE_enum side, MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C);
int MORSE_dsyrk_Tile(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C);
int MORSE_dsyr2k_Tile(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C);
#endif
//int MORSE_dsyev_Tile(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *W, MORSE_desc_t *T, double *Q, int LDQ);
int MORSE_dsyevd_Tile(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *W, MORSE_desc_t *T);
//int MORSE_dhegv_Tile( MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, double *W, MORSE_desc_t *T, MORSE_desc_t *Q);
//int MORSE_dsygvd_Tile(MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, double *W, MORSE_desc_t *T, MORSE_desc_t *Q);
//int MORSE_dsygst_Tile(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dsytrd_Tile(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, double *Q, int LDQ);
int MORSE_dlacpy_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
double MORSE_dlange_Tile(MORSE_enum norm, MORSE_desc_t *A);
#ifdef COMPLEX
double MORSE_dlansy_Tile(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A);
#endif
double MORSE_dlansy_Tile(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A);
double MORSE_dlantr_Tile(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A);
int MORSE_dlascal_Tile(MORSE_enum uplo, double alpha, MORSE_desc_t *A);
int MORSE_dlaset_Tile(MORSE_enum uplo, double alpha, double beta, MORSE_desc_t *A);
//int MORSE_dlaswp_Tile(MORSE_desc_t *A, int K1, int K2, int *IPIV, int INCX);
//int MORSE_dlaswpc_Tile(MORSE_desc_t *A, int K1, int K2, int *IPIV, int INCX);
int MORSE_dlauum_Tile(MORSE_enum uplo, MORSE_desc_t *A);
#ifdef COMPLEX
int MORSE_dplgsy_Tile(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed );
#endif
int MORSE_dplgsy_Tile(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed );
int MORSE_dplrnt_Tile(MORSE_desc_t *A, unsigned long long int seed );
int MORSE_dposv_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dpotrf_Tile(MORSE_enum uplo, MORSE_desc_t *A);
int MORSE_dsytrf_Tile(MORSE_enum uplo, MORSE_desc_t *A);
int MORSE_dpotri_Tile(MORSE_enum uplo, MORSE_desc_t *A);
int MORSE_dpotrimm_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C);
int MORSE_dpotrs_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dsymm_Tile(MORSE_enum side, MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C);
int MORSE_dsyrk_Tile(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C);
int MORSE_dsyr2k_Tile(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C);
int MORSE_dsysv_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dsytrs_Tile(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dtpgqrt_Tile( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2 );
int MORSE_dtpqrt_Tile( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T );
int MORSE_dtradd_Tile(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *B);
int MORSE_dtrmm_Tile(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dtrsm_Tile(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dtrsmpl_Tile(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B);
int MORSE_dtrsmrv_Tile(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *B);
int MORSE_dtrtri_Tile(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A);
int MORSE_dorglq_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_dorgqr_Tile(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_dormlq_Tile(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);
int MORSE_dormqr_Tile(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B);

/** ****************************************************************************
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 **/
int MORSE_dgeadd_Tile_Async(MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dgebrd_Tile_Async(MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dgecon_Tile_Async(MORSE_enum norm, MORSE_desc_t *A, double anorm, double *rcond, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dpocon_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, double anorm, double *rcond, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgelqf_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgelqs_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgels_Tile_Async(MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgemm_Tile_Async(MORSE_enum transA, MORSE_enum transB, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgeqrf_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgeqrs_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dgesv_Tile_Async(MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgesv_incpiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgesv_nopiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgesvd_Tile_Async(MORSE_enum jobu, MORSE_enum jobvt, MORSE_desc_t *A, double *S, MORSE_desc_t *T, double *U, int LDU, double *VT, int LDVT, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dgetrf_Tile_Async(MORSE_desc_t *A, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgetrf_incpiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgetrf_nopiv_Tile_Async(MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dgetri_Tile_Async(MORSE_desc_t *A, int *IPIV, MORSE_desc_t *W, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dgetrs_Tile_Async(MORSE_enum trans, MORSE_desc_t *A, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgetrs_incpiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgetrs_nopiv_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
int MORSE_dsymm_Tile_Async(MORSE_enum side, MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsyrk_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsyr2k_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
//int MORSE_dsyev_Tile_Async(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *W, MORSE_desc_t *T, double *Q, int LDQ, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsyevd_Tile_Async(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *W, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dhegv_Tile_Async( MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, double *W, MORSE_desc_t *T, MORSE_desc_t *Q, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dsygvd_Tile_Async(MORSE_enum itype, MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, double *W, MORSE_desc_t *T, MORSE_desc_t *Q, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dsygst_Tile_Async(MORSE_enum itype, MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsytrd_Tile_Async(MORSE_enum jobz, MORSE_enum uplo, MORSE_desc_t *A, double *D, double *E, MORSE_desc_t *T, double *Q, int LDQ, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dlacpy_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dlange_Tile_Async(MORSE_enum norm, MORSE_desc_t *A, double *value, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
int MORSE_dlansy_Tile_Async(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *value, MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
int MORSE_dlansy_Tile_Async(MORSE_enum norm, MORSE_enum uplo, MORSE_desc_t *A, double *value, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dlantr_Tile_Async(MORSE_enum norm, MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, double *value, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dlascal_Tile_Async(MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dlaset_Tile_Async(MORSE_enum uplo, double alpha, double beta, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dlaswp_Tile_Async(MORSE_desc_t *A, int K1, int K2, int *IPIV, int INCX, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dlaswpc_Tile_Async(MORSE_desc_t *A, int K1, int K2, int *IPIV, int INCX, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dlauum_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
#ifdef COMPLEX
int MORSE_dplgsy_Tile_Async(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
#endif
int MORSE_dplgsy_Tile_Async(double bump, MORSE_enum uplo, MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_dplrnt_Tile_Async(MORSE_desc_t *A, unsigned long long int seed, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_dposv_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dpotrf_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsytrf_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dpotri_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dpotrimm_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dpotrs_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsysv_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsytrs_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsymm_Tile_Async(MORSE_enum side, MORSE_enum uplo, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsyrk_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dsyr2k_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, double beta, MORSE_desc_t *C, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dtpgqrt_Tile_Async( int L, MORSE_desc_t *V1, MORSE_desc_t *T1, MORSE_desc_t *V2, MORSE_desc_t *T2, MORSE_desc_t *Q1, MORSE_desc_t *Q2, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_dtpqrt_Tile_Async( int L, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_desc_t *T, MORSE_sequence_t *sequence, MORSE_request_t *request );
int MORSE_dtradd_Tile_Async(MORSE_enum uplo, MORSE_enum trans, double alpha, MORSE_desc_t *A, double beta, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dtrmm_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dtrsm_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dtrsmpl_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *L, int *IPIV, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dtrsmrv_Tile_Async(MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag, double alpha, MORSE_desc_t *A, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dtrtri_Tile_Async(MORSE_enum uplo, MORSE_enum diag, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dorglq_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dorgqr_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dormlq_Tile_Async(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dormqr_Tile_Async(MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *T, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);

//int MORSE_dgecfi_Async(int m, int n, double *A, MORSE_enum f_in, int imb, int inb, MORSE_enum f_out, int omb, int onb, MORSE_sequence_t *sequence, MORSE_request_t *request);
//int MORSE_dgetmi_Async(int m, int n, double *A, MORSE_enum f_in, int mb,  int inb, MORSE_sequence_t *sequence, MORSE_request_t *request);

/**
 *  Declarations of libhqr dependent functions.
 */
/** ****************************************************************************
 *  Declarations of math functions (LAPACK layout) - alphabetical order
 **/
int MORSE_dgels_param(const libhqr_tree_t *qrtree, MORSE_enum trans, int M, int N, int NRHS, double *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT, double *B, int LDB);
int MORSE_dgelqf_param(const libhqr_tree_t *qrtree, int M, int N, double *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT);
int MORSE_dgelqs_param(const libhqr_tree_t *qrtree, int M, int N, int NRHS, double *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT, double *B, int LDB);
int MORSE_dgeqrf_param(const libhqr_tree_t *qrtree, int M, int N, double *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT);
int MORSE_dgeqrs_param(const libhqr_tree_t *qrtree, int M, int N, int NRHS, double *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *TT, double *B, int LDB);
int MORSE_dorglq_param(const libhqr_tree_t *qrtree, int M, int N, int K, double *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *TT, double *B, int LDB);
int MORSE_dorgqr_param(const libhqr_tree_t *qrtree, int M, int N, int K, double *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *descTT, double *B, int LDB);
int MORSE_dormlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, int M, int N, int K, double *A, int LDA, MORSE_desc_t *descT, MORSE_desc_t *TT, double *B, int LDB);
int MORSE_dormlq_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, int M, int N, int K, double *A, int LDA, MORSE_desc_t *descT, MORSE_desc_t *TT, double *B, int LDB);
int MORSE_dormqr_param(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, int M, int N, int K, double *A, int LDA, MORSE_desc_t *descTS, MORSE_desc_t *TT, double *B, int LDB);
/** ****************************************************************************
 *  Declarations of math functions (tile layout) - alphabetical order
 **/
int MORSE_dgels_param_Tile(const libhqr_tree_t *qrtree, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_dgelqf_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT);
int MORSE_dgelqs_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_dgeqrf_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT);
int MORSE_dgeqrs_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_dorglq_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_dorgqr_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_dorgqr_param_Tile(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_dormlq_param_Tile(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
int MORSE_dormqr_param_Tile(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B);
/** ****************************************************************************
 *  Declarations of math functions (tile layout, asynchronous execution) - alphabetical order
 **/
int MORSE_dgels_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgelqf_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgelqs_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgeqrf_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dgeqrs_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dorglq_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dorgqr_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dormlq_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dormqr_param_Tile_Async(const libhqr_tree_t *qrtree, MORSE_enum side, MORSE_enum trans, MORSE_desc_t *A, MORSE_desc_t *TS, MORSE_desc_t *TT, MORSE_desc_t *B, MORSE_sequence_t *sequence, MORSE_request_t *request);

/**
 *  Declarations of workspace allocation functions (tile layout) - alphabetical order
 */
int MORSE_Alloc_Workspace_dgesv_incpiv(        int N, MORSE_desc_t **descL, int **IPIV, int p, int q);
int MORSE_Alloc_Workspace_dgetrf_incpiv(int M, int N, MORSE_desc_t **descL, int **IPIV, int p, int q);

int MORSE_Alloc_Workspace_dgebrd(int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dgeev( int N,        MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dgehrd(int N,        MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dgelqf(int M, int N, MORSE_desc_t **T, int p, int q);
int MORSE_Alloc_Workspace_dgels( int M, int N, MORSE_desc_t **T, int p, int q);
int MORSE_Alloc_Workspace_dgeqrf(int M, int N, MORSE_desc_t **T, int p, int q);
int MORSE_Alloc_Workspace_dgesvd(int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dsyev( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dsyevd( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dhegv( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dsygvd( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dsytrd(int M, int N, MORSE_desc_t **descT, int p, int q);

/** ****************************************************************************
 *  Declarations of workspace allocation functions (tile layout, asynchronous execution) - alphabetical order
 **/
int MORSE_Alloc_Workspace_dgetri_Tile_Async(MORSE_desc_t *A, MORSE_desc_t *W);

/* Warning: Those functions are deprecated */
int MORSE_Alloc_Workspace_dgelqf_Tile(int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dgels_Tile( int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dgeqrf_Tile(int M, int N, MORSE_desc_t **descT, int p, int q);
int MORSE_Alloc_Workspace_dgesv_incpiv_Tile (int N, MORSE_desc_t **descL, int **IPIV, int p, int q);
int MORSE_Alloc_Workspace_dgetrf_incpiv_Tile(int N, MORSE_desc_t **descL, int **IPIV, int p, int q);

/** ****************************************************************************
 *  Auxiliary function prototypes
 **/
int MORSE_dLapack_to_Tile(double *Af77, int LDA, MORSE_desc_t *A);
int MORSE_dTile_to_Lapack(MORSE_desc_t *A, double *Af77, int LDA);
int MORSE_dLapack_to_Tile_Async(double *Af77, int LDA, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request);
int MORSE_dTile_to_Lapack_Async(MORSE_desc_t *A, double *Af77, int LDA, MORSE_sequence_t *sequence, MORSE_request_t *request);

/** ****************************************************************************
 *  User Builder function prototypes
 **/
int MORSE_dbuild(MORSE_enum uplo, int M, int N, double *A, int LDA, void *user_data, void* user_build_callback);
int MORSE_dbuild_Tile(MORSE_enum uplo,  MORSE_desc_t *A, void *user_data, void* user_build_callback );
int MORSE_dbuild_Tile_Async(MORSE_enum uplo, MORSE_desc_t *A, void *user_data, void* user_build_callback, MORSE_sequence_t *sequence, MORSE_request_t  *request);

END_C_DECLS

#undef COMPLEX

#endif
