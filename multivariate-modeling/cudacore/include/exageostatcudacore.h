/**
 *
 * Copyright (c) 2017-2019  King Abdullah University of Science and Technology
 * All rights reserved.
 *
 * ExaGeoStat is a software package provided by KAUST
 **/
/**
 *
 * @file exageostatcudacore.h
 *
 * CUDA core functions header file.
 *
 * @version 1.0.0
 *
 * @author Sameh Abdulah
 * @date 2019-01-27
 *
 **/
#ifndef _EXAGEOSTATCUDACORE_H_
#define _EXAGEOSTATCUDACORE_H_
#include "cudablas.h"
#include "lapacke.h"
#include "../../misc/include/MLE_misc.h"
#include "../../misc/include/flat_file.h"
int cuda_dsconv(int m, int n,
		cuDoubleComplex *A, int lda,
		cuFloatComplex *B, int ldb, cublasHandle_t handle);


int cuda_sdconv(int m, int n,
		cuFloatComplex *A, int lda,
		cuDoubleComplex *B, int ldb, cublasHandle_t handle);
#ifdef __cplusplus
extern "C"
{
#endif
	void double2float_array(int nrows, int ncols,
			const double *H, int ldh,
			float *F, int ldf,
			cublasOperation_t transa, cudaStream_t stream);

	void float2double_array(int nrows, int ncols,
			const float *F, int ldf,
			double *H, int ldh ,
			cublasOperation_t transa, cudaStream_t stream);


	void dcmg_array(double *A, int m, int n,
			int m0, int n0,
			location  *l1, location *l2,
			double *localtheta, int distance_metric,
			cudaStream_t stream);

        void scmg_array(double *A, int m, int n,
                        int m0, int n0,
                        location  *l1, location *l2,
                        double *localtheta, int distance_metric,
                        cudaStream_t stream);

        void sdcmg_array(double *A, int m, int n,
                        int m0, int n0,
                        location  *l1, location *l2,
                        double *localtheta, int distance_metric,
                        cudaStream_t stream);

#ifdef __cplusplus
}
#endif
#endif


