/**
 *
 * @file MLE.h
 *
 *
 *  ExaGeoStat is a software package provided by KAUST,
 *  King Abdullah University of Science and Technology - ECRC
 *
 * @version 0.1.0
 * @author Sameh Abdulah
 * @date 2016-11-22
 * @generated d Fri Nov 22 15:11:13 2016
 *
 **/
#ifndef _MLE_LR_S_H_
#define _MLE_LR_S_H_
#include "MLE_misc.h"
#include "flat_file.h"
#if defined(EXAGEOSTAT_USE_NETCDF)
        #include "nc_file.h"
#endif
#include "starpu_exageostat_approx.h"
#include "starsh-spatial.h"
#include "hicma_struct.h"
#include "hicma_z.h"
#include "aux/auxcompute_z.h"
#include "aux/auxdescutil.h"


//Generate Observations Vector (Z) for testing Maximum Likelihood function -- MORSE-sync
int HICMA_MLE_szvg_Tile (MLE_data *data,  double * Nrand, double * initial_theta, int n, int dts, int log, int p_grid, int q_grid);
//Generate Observations Vector (Z) for testing Maximum Likelihood function -- MORSE-async
int HICMA_MLE_szvg_Tile_Async (MLE_data *data,  double * Nrand, double * initial_theta, int n, int dts, int log, int p_grid, int q_grid);
//compute MLE function --syn
double HICMA_smle_Tile(unsigned n, const double *theta, double *grad, void *data);
//compute MLE function --async
double HICMA_smle_Tile_Async(unsigned n, const double *theta, double *grad, void *data);
void HICMA_dmle_Predict_Allocate(MLE_data *MORSE_data, int nZmiss, int nZobs, int lts, int p_grid, int q_grid, int mse_flag);
//Predict missing values base on a set of given values and covariance matrix
double HICMA_smle_Predict_Tile(MLE_data *MORSE_data, double * theta, int nZmiss, int nZobs, double *Zobs, double *Zactual, double *Zmiss, int n, int lts);
//init Hicma decriptors
void HICMA_smle_Call(MLE_data  *data, int ncores, int ts, int N, int nZobs, int nZmiss,  int gpus, int p_grid, int q_grid);
void HICMA_MLE_zcpy( MLE_data *data, double *streamdata);

#endif
