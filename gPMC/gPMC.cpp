/*
 *
 * file  example.cpp
 * brief example program for using goPMC.
 *
 * author Nan Qin
 *
 * last update on 9/21/2016
 *
 */
#include <algorithm>
#include <random>
#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <fstream>

#define NDOSECOUNTERS 1
#include "goPMC.h"
#define N 100000

#include "gPMC_ggo.h"

template < class TArgsInfo, class TCleanupFunction = void(*)(TArgsInfo*) >
class args_info_manager
{
public:
	args_info_manager(TArgsInfo & args_info, TCleanupFunction cf)
	{
		this->args_info_pointer = &args_info;
		this->cleanup_function = cf;
	}
	~args_info_manager()
	{
		this->cleanup_function(this->args_info_pointer);
	}
private:
	TArgsInfo * args_info_pointer;
	TCleanupFunction cleanup_function;
};
//--------------------------------------------------------------------
/** \brief Process gengetopt with config file option (shamelessly stolen from RTK)
*
* \author Simon Rit
*
* \ingroup Macro
*/
#define GGO(ggo_filename, args_info)                                                                       \
	args_info_##ggo_filename args_info;                                                                    \
	cmdline_parser_##ggo_filename##_params args_params;                                                    \
	cmdline_parser_##ggo_filename##_params_init(&args_params);                                             \
	args_params.print_errors = 1;                                                                          \
	args_params.check_required = 0;                                                                        \
	args_params.override = 1;                                                                              \
	args_params.initialize = 1;                                                                            \
if (0 != cmdline_parser_##ggo_filename##_ext(argc, argv, &args_info, &args_params))                        \
	{                                                                                                      \
	std::cerr << "Error in cmdline_parser_" #ggo_filename "_ext" << std::endl;                             \
	exit(1);                                                                                               \
	}                                                                                                      \
	std::string configFile;                                                                                \
if (args_info.config_given)                                                                                \
	configFile = args_info.config_arg;                                                                     \
	cmdline_parser_##ggo_filename##_free(&args_info);                                                      \
if (configFile != "")                                                                                      \
	{                                                                                                      \
if (0 != cmdline_parser_##ggo_filename##_config_file(configFile.c_str(), &args_info, &args_params))        \
	  {                                                                                                    \
	  std::cerr << "Error in cmdline_parser_" #ggo_filename "_config_file" << std::endl;                   \
	  exit(1);                                                                                             \
	  }                                                                                                    \
	  args_params.initialize = 0;                                                                          \
	}                                                                                                      \
	args_params.check_required = 1;                                                                        \
if (0 != cmdline_parser_##ggo_filename##_ext(argc, argv, &args_info, &args_params))                        \
	{                                                                                                      \
	std::cerr << "Error in cmdline_parser_" #ggo_filename "_ext" << std::endl;                             \
	exit(1);                                                                                               \
	}                                                                                                      \
	args_info_manager< args_info_##ggo_filename >                                                          \
	manager_object(args_info, cmdline_parser_##ggo_filename##_free);
//--------------------------------------------------------------------

// A function to initialize source protons. Should be replaced by real beams.
void initSource(cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight){
	unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
	std::minstd_rand0 g1(seed1);  // minstd_rand0 is a standard linear_congruential_engine
	std::fill_n(T, N, 120.0f);
	std::fill_n(weight, N, 1.0f);
	for (int i = 0; i < N; i++){
		pos[i].s[0] = 5 * float(g1()) / g1.max(); // -15; // -15 to -10
		pos[i].s[1] = 120; // -20;   // ^--- between 0 and the largest possible max = 2147483646
		pos[i].s[2] = 5 * float(g1()) / g1.max(); // +25; //  25 to  30
	} //                                ^------- UP TO the largest possible max = 2147483646
	const cl_float3 temp2 = { 0.0f, 1.0f, 0.0f };
	std::fill_n(dir, N, temp2);
}

int main(int argc, char * argv[])
{
	GGO(gPMC, args_info);


	FILE *stream;
	if ((stream = freopen("D:\stdout.txt", "w", stdout)) == NULL)
		exit(-1);
	FILE *stream_err;
	if ((stream_err = freopen("D:\stderr.txt", "w", stderr)) == NULL)
		exit(-1);
	// Get OpenCL platform and device.
	cl::Platform platform;
	cl::Platform::get(&platform);
	std::vector<cl::Device> devs;
	platform.getDevices(CL_DEVICE_TYPE_CPU, &devs);

	cl::Device device;
	try{
		device = devs.at(0); // throws exception in contrary to []
	}
	catch (const std::exception& e) {
		std::cout << "Well, this happened: " << e.what() << std::endl << "That usually means you tried to compile for CPU-device with CUDA" << std::endl;
		std::cout << "OR you compiled for GPU-device didn't and didn't have a GPU" << std::endl;
		std::cin.ignore();
		return -1;
	}
	
	// Initialize simulation engine.
	goPMC::MCEngine mcEngine;
	mcEngine.initializeComputation(platform, device);
	
	// Read and process physics data.
	mcEngine.initializePhysics("input");
	
	// Read and process patient Dicom CT data. ITK has origin at upper left
	mcEngine.initializePhantom("zzzCetphan504"); // "090737"); // "directoryToDicomData");
	
	// Initialize source protons with arrays of energy (T), position (pos), direction (dir) and weight (weight) of each proton.
	// Position and direction should be defined in Dicom CT coordinate.
	cl_float * T = new cl_float[N];     //Energy(MeV?) = [120.0, ..., 120.0]
	cl_float3 * pos = new cl_float3[N]; //Position    = [(5*rand_1-15, -20, 5*rand_1+25), ..., (5*rand_N-15, -20, 5*rand_N+25)]
	cl_float3 * dir = new cl_float3[N]; //Direction   = [(0, 1, 0), ..., (0, 1, 0)] = y?          ^-------- 0 <= rand_X <= 1
	cl_float * weight = new cl_float[N];//Weight      = [1.0, ..., 1.0]
	initSource(T, pos, dir, weight);
	
	// Choose a physics quantity to score for this simulation run.
	// Scoring quantity could be one of {DOSE2MEDIUM, DOSE2WATER, FLUENCE, LETD}.
	// LETD is dose weighted LET, to get dose averaged LET, divide it by DOSE2MEDIUM from another simulation run.
	std::string quantity("DOSE2WATER"); 
	
	// Run simulation.
	mcEngine.simulate(T, pos, dir, weight, N, quantity);
	
	// Get simulation results.
	std::vector<cl_float> doseMean, doseStd;
	mcEngine.getResult(doseMean, doseStd);
	
	// Do something with doseMean and doseStd //
	std::cout << "doseMean size: " << doseMean.size();
	cl_float mean = 0;
	for (std::vector<cl_float>::iterator it = doseMean.begin(); it != doseMean.end(); ++it)
		mean += *it;
	std::cout << " sum: " << mean;
	std::cout << " mean: " << mean / doseMean.size() << std::endl;

	std::cout << "doseStd size: " << doseStd.size();
	mean = 0;
	for (std::vector<cl_float>::iterator it = doseStd.begin(); it != doseStd.end(); ++it)
		mean += *it;
	std::cout << " sum: " << mean;
	std::cout << " mean: " << mean / doseStd.size() << std::endl;
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //

	// Clear the scoring counters in previous simulation runs.
	mcEngine.clearCounter();


	delete[] T;
	delete[] pos;
	delete[] dir;
	delete[] weight;
	stream = freopen("CON", "w", stdout);
	stream_err = freopen("CON", "w", stderr);
	return 0;
}

