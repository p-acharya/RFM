/*
  Mex file for particle flow simulation.

  Based on rfm_stochastic.m

  Compile: mex rfm_stochastic_c.c particle_flow.c -output <name>

  where: <name> is rfm_stochastic_c if particle_flow() function is used below,
         or rfm_stochastic_c_inf if particle_flow_l0inf() function is used below.

  Note: The current mex option file is in ~/.matlab/R2012a/mexopts.sh
        This file was modifie as follows (all in the maci64 section):
        1. Changed CC=llvm-gcc-4.2 to CC=llvm-gcc.
        2. Changed CC=llvm-g-4.2++ to CC=llvm-g++.
        3. Changed CFLAGS="$CFLAGS  -fexceptions" to CFLAGS="$CFLAGS  -fexceptions -Dchar16_t=UINT16_T".


  Yoram Zarai, 11/2/14
 */

#include "matrix.h"
#include "mex.h"

#include "particle_flow.h"

#define NUM_IN_ARGUMENTS   (4)


/* Matlab interface function */
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
  /* check for proper number of arguments */
  if( nrhs != NUM_IN_ARGUMENTS ) {
    mexErrMsgIdAndTxt( "rfm_stochastic_c:num_in_arg", "%d inputs are required", NUM_IN_ARGUMENTS );
  }

  /* input arguments */
  double*  lambda = mxGetPr( prhs[ 0 ] );
  double   sim_step = mxGetScalar( prhs[ 1 ] );
  double   sim_time = mxGetScalar( prhs[ 2 ] );
  uint32_t num_nodes = mxGetScalar( prhs[ 3 ] );

  /* output arguments */
  plhs[ 0 ] = mxCreateDoubleMatrix( 1, num_nodes + 2, mxREAL );
  double *occu = mxGetPr( plhs[ 0 ] );

  /* maximum number of particles that can traverse the chain of num_nodes nodes */
  //uint32_t max_num_part = (uint32_t)( sim_time / sim_step ); 
  uint32_t max_num_partc = estimate_terminate_partc( sim_time, sim_step, num_nodes, lambda );
  plhs[ 1 ] = mxCreateDoubleMatrix( 1, max_num_partc, mxREAL );
  double *delay = mxGetPr( plhs[ 1 ] );

  //int32_t last_id = particle_flow_l0inf( sim_time, sim_step, num_nodes, lambda, delay, occu );
  int32_t last_id = particle_flow( sim_time, sim_step, num_nodes, lambda, delay, occu );
  if( last_id == (int32_t)EXIT_ERROR_VAL ) {
    fprintf( stderr, "Error executing rfm_stochastic_c !!\n" );
  }
  plhs[ 2 ] = mxCreateDoubleScalar( (double)last_id );
}
