/**
 *
 * @file morse_simulate.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon SimGrid simulation header
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2014-10-02
 *
 */
#ifndef _MORSE_SIMULATE_H_
#define _MORSE_SIMULATE_H_

#include "chameleon/chameleon_config.h"

/* we need this when starpu is compiled with simgrid enabled */
#if defined(CHAMELEON_SCHED_STARPU) && defined(CHAMELEON_SIMULATION)
#include <starpu_simgrid_wrap.h>
#if defined(CHAMELEON_USE_MPI)
#define starpu_main smpi_simulated_main_
#endif
#endif

#endif
