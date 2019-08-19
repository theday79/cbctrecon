#ifndef ERR_CODE_H
#define ERR_CODE_H
/*----------------------------------------------------------------------------
 *
 * Name:     err_code()
 *
 * Purpose:  Function to output descriptions of errors for an input error code
 *           and quit a program on an error with a user message
 *
 *
 * RETURN:   echoes the input error code / echos user message and exits
 *
 * HISTORY:  Written by Tim Mattson, June 2010
 *           This version automatically produced by genErrCode.py
 *           script written by Tom Deakin, August 2013
 *           Modified by Bruce Merry, March 2014
 *           Updated by Tom Deakin, October 2014
 *               Included the checkError function written by
 *               James Price and Simon McIntosh-Smith
 *
 *----------------------------------------------------------------------------
 */

/*
 *
 * This code is released under the "attribution CC BY" creative commons license.
 * In other words, you can use it in any way you see fit, including
 * commercially, but please retain an attribution for the original authors: the
 * High Performance Computing Group at the University of Bristol. Contributors
 * include Simon McIntosh-Smith, James Price, Tom Deakin and Mike O'Connor.
 *
 */

/*
 * Modified by A. Gravgaard
 */
#ifdef __cplusplus
#include <cstdio>
#else
#include <stdio.h>
#include <stdlib.h>
#endif

#if CBCTRECON_OPENCL_VERSION >= 210
#define CL_HPP_MINIMUM_OPENCL_VERSION 200
#define CL_HPP_TARGET_OPENCL_VERSION 210
#else
#define CL_HPP_MINIMUM_OPENCL_VERSION 120
#define CL_HPP_TARGET_OPENCL_VERSION 120
#endif

#include "OpenCL/cl2.hpp"

const char *ocl_err_code(cl_int err_in);

void check_ocl_error(cl_int err, const char *operation, const char *filename,
                     int line);

#define checkError(E, S) check_ocl_error(E, S, __FILE__, __LINE__)

#endif
