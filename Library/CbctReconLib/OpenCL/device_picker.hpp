/*------------------------------------------------------------------------------
 *
 * Name:       device_picker.h
 *
 * Purpose:    Provide a simple CLI to specify an OpenCL device at runtime
 *
 * Note:       Must be included AFTER the relevant OpenCL header
 *             See one of the Matrix Multiply exercises for usage
 *
 * HISTORY:    Method written by James Price, October 2014
 *             Extracted to a common header by Tom Deakin, November 2014
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

#ifndef DEVICE_PICKER_HPP
#define DEVICE_PICKER_HPP

#include <string>
#include <vector>

#include "cbctrecon_config.h"

#include "OpenCL/err_code.hpp"

#ifndef CL_DEVICE_BOARD_NAME_AMD
#define CL_DEVICE_BOARD_NAME_AMD 0x4038
#endif

#define MAX_INFO_STRING 256

CBCTRECON_API std::vector<cl::Device> OpenCL_getDeviceList();

CBCTRECON_API std::string OpenCL_getDeviceName(const cl::Device &device);

CBCTRECON_API void OpenCL_parseArguments(int argc, char *argv[],
                                         cl_uint *deviceIndex);

#endif
