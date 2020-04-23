/*! @copyright (c) 2020 King Abdullah University of Science and
 *                      Technology (KAUST). All rights reserved.
 *
 * STARS-H is a software package, provided by King Abdullah
 *             University of Science and Technology (KAUST)
 *
 * Generated from file /home/abdullsm/develop/stars-h/src/applications/spatial/kernel_bivariate_parsimonious.c with NDIM=1
 * Generate different functions for different dimensions. This hack improves
 * performance in certain cases. Value 'n' stands for general case, whereas all
 * other values correspond to static values of dimensionality.
 * During code generation step, each appearance of 1 (including this one)
 * will be replace by proposed values. If you want to use this file outside
 * STARS-H, simply do substitutions yourself.
 *
 * @file /home/abdullsm/develop/stars-h/build/src/applications/spatial/kernel_bivariate_parsimonious_1d.c
 * @version 0.1.1
 * @author Sameh Abdulah
 * @author Aleksandr Mikhalev
 * @date 2020-02-09
 */

#include "common.h"
#include "starsh.h"
#include "starsh-spatial.h"

// If dimensionality is static
#if (1 != n)
//! Replace variable ndim with static integer value
#define ndim 1
#endif

#ifdef GSL

void starsh_ssdata_block_bivariate_parsimonious_kernel_1d(int nrows, int ncols,
        STARSH_int *irow, STARSH_int *icol, void *row_data, void *col_data,
        void *result, int ld)
//! Mat&eacute;rn kernel for 1-dimensional spatial statistics problem
/*! Fills matrix \f$ A \f$ with values
 * \f[
 *      A_{ij} = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu1)} \left( \frac{r_{ij}}
 *      {\beta} \right)^{\nu} K_{\nu} \left( \frac{r_{ij}}{\beta} \right) +
 *      \mu \delta(r_{ij}),
 * \f]
 * where \f$ \Gamma \f$ is the Gamma function, \f$ K_{\nu} \f$ is the modified
 * Bessel function of the second kind, \f$ \delta \f$ is the delta function
 * \f[
 *      \delta(x) = \left\{ \begin{array}{ll} 0, & x \ne 0\\ 1, & x = 0
 *      \end{array} \right.,
 * \f]
 * \f$ r_{ij} \f$ is a distance between \f$i\f$-th and \f$j\f$-th spatial
 * points and variance \f$ \sigma \f$, correlation length \f$ \beta \f$,
 * smoothing parameter \f$ \nu \f$ and noise \f$ \mu \f$ come from \p
 * row_data (\ref STARSH_ssdata object). No memory is allocated in this
 * function!
 *
 * @param[in] nrows: Number of rows of \f$ A \f$.
 * @param[in] ncols: Number of columns of \f$ A \f$.
 * @param[in] irow: Array of row indexes.
 * @param[in] icol: Array of column indexes.
 * @param[in] row_data: Pointer to physical data (\ref STARSH_ssdata object).
 * @param[in] col_data: Pointer to physical data (\ref STARSH_ssdata object).
 * @param[out] result: Pointer to memory of \f$ A \f$.
 * @param[in] ld: Leading dimension of `result`.
 * @sa starsh_ssdata_block_bivariate_parsimonious_kernel_1d(),
 *      starsh_ssdata_block_bivariate_parsimonious_kernel_2d(),
 *      starsh_ssdata_block_bivariate_parsimonious_kernel_3d(),
 *      starsh_ssdata_block_bivariate_parsimonious_kernel_4d(),
 *      starsh_ssdata_block_bivariate_parsimonious_kernel_nd().
 * @ingroup app-spatial-kernels
 * */
{
    int i, j, k;
    STARSH_ssdata *data1 = row_data;
    STARSH_ssdata *data2 = col_data;
    double tmp, dist;
    // Read parameters
// If dimensionality is not static
#if (1 == n)
    int ndim = data1->particles.ndim;
#endif
    double beta = data1->beta;
    double nu1 = data1->nu;
    double noise = data1->noise;
    double sigma1 = data1->sigma;
    double nu2 = data1->nu2;
    double sigma2 = data1->sigma2;
    double corr = data1->corr;



    // Get coordinates
    STARSH_int count1 = data1->particles.count;
    STARSH_int count2 = data2->particles.count;
    double *x1[ndim], *x2[ndim];
    x1[0] = data1->particles.point;
    x2[0] = data2->particles.point;
    //#pragma omp simd
    for(i = 1; i < ndim; i++)
    {
	    x1[i] = x1[0]+i*count1;
	    x2[i] = x2[0]+i*count2;
    }
    double *x1_cur, *x2_cur;
    double *buffer = result;
    //Prepare 6 paramters for bivariate parsimonious kernel
    double expr = 0.0;
    double con1 = 0.0, con2 = 0.0, con12 = 0.0, rho = 0.0, nu12 = 0.0;

    con1 = pow(2,(nu1-1)) * tgamma(nu1);
    con1 = 1.0/con1;
    con1 = sigma1 * con1;

    con2 = pow(2,(nu2-1)) * tgamma(nu2);
    con2 = 1.0/con2;
    con2 = sigma2 * con2;

    //The average
    nu12 = 0.5 * (nu1 + nu2);

    rho = corr * sqrt( (tgamma(nu1 + 1)*tgamma(nu2 + 1)) /
		    (tgamma(nu1) * tgamma(nu2)) ) *
	    tgamma(nu12) / tgamma(nu12 + 1);



    con12 = pow(2,(nu12-1)) * tgamma(nu12);
    con12 = 1.0/con12;
    con12 = rho * sqrt(sigma1 * sigma2) * con12;
    // Fill column-major matrix
    //#pragma omp simd
    for(j = 0; j < ncols; j+=2)
    {
	    for(i = 0; i < nrows; i+=2)
	    {
		    dist = 0.0;
		    for(k = 0; k < ndim; k++)
		    {
			    tmp = x1[k][irow[i/2]]-x2[k][icol[j/2]];
			    dist += tmp*tmp;
		    }
		    dist = sqrt(dist)/beta;
		    if(dist == 0)
		    {
			    buffer[j*(size_t)ld+i] = sigma1+noise;
			    buffer[(j+1)*(size_t)ld+i] = buffer[j*(size_t)ld+(i+1)] = rho * sqrt(sigma1 * sigma2) ;
			    buffer[(j+1)*(size_t)ld+(i+1)] = sigma2;

		    }
		    else
		    {

			    buffer[j*(size_t)ld+i]  = con1 * pow(expr, nu1)*gsl_sf_bessel_Knu(nu1, expr);
			    buffer[(j+1)*(size_t)ld+i] = buffer[j*(size_t)ld+(i+1)] =  con12 * pow(dist, nu12) * gsl_sf_bessel_Knu(nu12, dist);
			    buffer[(j+1)*(size_t)ld+(i+1)] = con2 * pow(expr, nu2)*gsl_sf_bessel_Knu(nu2, dist);

		    }
	    }
    }
}

void starsh_ssdata_block_bivariate_parsimonious_kernel_1d_simd(int nrows, int ncols,
		STARSH_int *irow, STARSH_int *icol, void *row_data, void *col_data,
		void *result, int ld)
	//! Mat&eacute;rn kernel for 1-dimensional spatial statistics problem
	/*! Fills matrix \f$ A \f$ with values
	 * \f[
	 *      A_{ij} = \sigma^2 \frac{2^{1-\nu}}{\Gamma(\nu1)} \left( \frac{r_{ij}}
	 *      {\beta} \right)^{\nu} K_{\nu} \left( \frac{r_{ij}}{\beta} \right) +
	 *      \mu \delta(r_{ij}),
	 * \f]
	 * where \f$ \Gamma \f$ is the Gamma function, \f$ K_{\nu} \f$ is the modified
	 * Bessel function of the second kind, \f$ \delta \f$ is the delta function
	 * \f[
	 *      \delta(x) = \left\{ \begin{array}{ll} 0, & x \ne 0\\ 1, & x = 0
	 *      \end{array} \right.,
	 * \f]
	 * \f$ r_{ij} \f$ is a distance between \f$i\f$-th and \f$j\f$-th spatial
	 * points and variance \f$ \sigma \f$, correlation length \f$ \beta \f$,
	 * smoothing parameter \f$ \nu \f$ and noise \f$ \mu \f$ come from \p
	 * row_data (\ref STARSH_ssdata object). No memory is allocated in this
	 * function!
	 *
	 * Uses SIMD instructions.
	 *
	 * @param[in] nrows: Number of rows of \f$ A \f$.
	 * @param[in] ncols: Number of columns of \f$ A \f$.
	 * @param[in] irow: Array of row indexes.
	 * @param[in] icol: Array of column indexes.
	 * @param[in] row_data: Pointer to physical data (\ref STARSH_ssdata object).
	 * @param[in] col_data: Pointer to physical data (\ref STARSH_ssdata object).
	 * @param[out] result: Pointer to memory of \f$ A \f$.
	 * @param[in] ld: Leading dimension of `result`.
	 * @sa starsh_ssdata_block_bivariate_parsimonious_kernel_1d_simd(),
	 *      starsh_ssdata_block_bivariate_parsimonious_kernel_2d_simd(),
	 *      starsh_ssdata_block_bivariate_parsimonious_kernel_3d_simd(),
	 *      starsh_ssdata_block_bivariate_parsimonious_kernel_4d_simd(),
	 *      starsh_ssdata_block_bivariate_parsimonious_kernel_nd_simd().
	 * @ingroup app-spatial-kernels
	 * */
{
	int i, j, k;
	STARSH_ssdata *data1 = row_data;
	STARSH_ssdata *data2 = col_data;
	double tmp, dist;
	// Read parameters
	// If dimensionality is not static
#if (1 == n)
	int ndim = data1->particles.ndim;
#endif
	double beta = data1->beta;
	double nu1 = data1->nu;
	double noise = data1->noise;
	double sigma1 = data1->sigma;
	double nu2 = data1->nu2;
	double sigma2 = data1->sigma2;
	double corr = data1->corr;
	// Get coordinates
	STARSH_int count1 = data1->particles.count;
	STARSH_int count2 = data2->particles.count;
	double *x1[ndim], *x2[ndim];
	x1[0] = data1->particles.point;
	x2[0] = data2->particles.point;


//	for(k=0;k<4;k++)
//		printf("(%f, %f)\n",x1[0][k], x2[0][k]); 


//	printf("count1=%d, count2=%d\n", count1, count2);

#pragma omp simd
	for(i = 1; i < ndim; i++)
	{
		x1[i] = x1[0]+i*count1;
		x2[i] = x2[0]+i*count2;
//		for(k=0;k<4;k++)
//			printf("i:(%f, %f)\n",x1[i][k], x2[i][k]);
	}
	double *x1_cur, *x2_cur;
	double *buffer = result;
	//Prepare 6 paramters for bivariate parsimonious kernel
	double expr = 0.0;
	double con1 = 0.0, con2 = 0.0, con12 = 0.0, rho = 0.0, nu12 = 0.0;

	con1 = pow(2,(nu1-1)) * tgamma(nu1);
	con1 = 1.0/con1;
	con1 = sigma1 * con1;

	con2 = pow(2,(nu2-1)) * tgamma(nu2);
	con2 = 1.0/con2;
	con2 = sigma2 * con2;

	//The average
	nu12 = 0.5 * (nu1 + nu2);

	rho = corr * sqrt( (tgamma(nu1 + 1)*tgamma(nu2 + 1)) /
			(tgamma(nu1) * tgamma(nu2)) ) *
		tgamma(nu12) / tgamma(nu12 + 1);



	con12 = pow(2,(nu12-1)) * tgamma(nu12);
	con12 = 1.0/con12;
	con12 = rho * sqrt(sigma1 * sigma2) * con12;
	// Fill column-major matrix
	// Fill column-major matrix
	int index_i=0;
	int index_j=0;
//	printf("ncols=%d, nrows=%d\n", ncols, nrows);
#pragma omp simd
	for(j = 0; j < ncols; j+=2)
	{
		index_i=0;
		for(i = 0; i < nrows; i+=2)
		{
		//	printf("j:%d, i:%d\n", j, i);
			dist = 0.0;
			for(k = 0; k < ndim; k++)
			{
				tmp = x1[k][irow[index_i]]-x2[k][icol[index_j]];
		//		printf("(%d)%f - %f, ", k, x1[k][irow[index_i]], x2[k][icol[index_j]]);
				dist += tmp*tmp;
			}
			dist = sqrt(dist)/beta;

//			printf("\n");	

			if(dist == 0)
			{
				buffer[j*(size_t)ld+i] = sigma1+noise;
				buffer[(j+1)*(size_t)ld+i] = buffer[j*(size_t)ld+(i+1)] = rho * sqrt(sigma1 * sigma2) ;
				buffer[(j+1)*(size_t)ld+(i+1)] = sigma2;

			}
			else
			{

				buffer[j*(size_t)ld+i]  = con1 * pow(expr, nu1)*gsl_sf_bessel_Knu(nu1, dist);
				buffer[(j+1)*(size_t)ld+i] = buffer[j*(size_t)ld+(i+1)] =  con12 * pow(dist, nu12) * gsl_sf_bessel_Knu(nu12, dist);
				buffer[(j+1)*(size_t)ld+(i+1)] = con2 * pow(dist, nu2)*gsl_sf_bessel_Knu(nu2, dist);

			}
			index_i++;
		}
		index_j++;
	}
}
#endif // GSL

