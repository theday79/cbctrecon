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

#include "itkImage.h"
#include "itkImageFileWriter.h"

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
#define GGO(ggo_filename, args_info)                                                                     \
	args_info_##ggo_filename args_info;                                                                    \
	cmdline_parser_##ggo_filename##_params args_params;                                                    \
	cmdline_parser_##ggo_filename##_params_init(&args_params);                                             \
	args_params.print_errors = 1;                                                                          \
	args_params.check_required = 0;                                                                        \
	args_params.override = 1;                                                                              \
	args_params.initialize = 1;                                                                            \
if (0 != cmdline_parser_##ggo_filename##_ext(argc, argv, &args_info, &args_params))                      \
	{                                                                                                      \
	std::cerr << "Error in cmdline_parser_" #ggo_filename "_ext" << std::endl;                             \
	exit(1);                                                                                               \
	}                                                                                                      \
	std::string configFile;                                                                                \
if (args_info.config_given)                                                                              \
	configFile = args_info.config_arg;                                                                     \
	cmdline_parser_##ggo_filename##_free(&args_info);                                                      \
if (configFile != "")                                                                                    \
	{                                                                                                      \
if (0 != cmdline_parser_##ggo_filename##_config_file(configFile.c_str(), &args_info, &args_params))      \
	  {                                                                                                    \
	  std::cerr << "Error in cmdline_parser_" #ggo_filename "_config_file" << std::endl;                   \
	  exit(1);                                                                                             \
	  }                                                                                                    \
	  args_params.initialize = 0;                                                                          \
	}                                                                                                      \
	args_params.check_required = 1;                                                                        \
if (0 != cmdline_parser_##ggo_filename##_ext(argc, argv, &args_info, &args_params))                      \
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

	const std::string stdout_file = std::string(args_info.path_arg) + "\\..\\gPMCstdout.txt";
	const std::string stderr_file = std::string(args_info.path_arg) + "\\..\\gPMCstderr.txt";
	FILE *stream;
	if ((stream = freopen(stdout_file.c_str(), "w", stdout)) == NULL)
		exit(-1);
	FILE *stream_err;
	if ((stream_err = freopen(stderr_file.c_str(), "w", stderr)) == NULL)
		exit(-1);
	// Get OpenCL platform and device.
	cl::Platform platform;
	cl::Platform::get(&platform);
	std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;
	std::vector<cl::Device> devs;
	if (!strcmp(args_info.hardware_arg, "gpu"))
	  std::cout << "Getting device GPU returned: " << platform.getDevices(CL_DEVICE_TYPE_GPU, &devs) << std::endl;
	else if (!strcmp(args_info.hardware_arg, "cpu"))
		std::cout << "Getting devices CPU returned: " << platform.getDevices(CL_DEVICE_TYPE_CPU, &devs) << std::endl;
  else if (!strcmp(args_info.hardware_arg, "acc"))
		std::cout << "Getting devices ACCELERATOR returned: " << platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devs) << std::endl;

	cl::Device device;
	try{
		device = devs.at(0); // throws exception in contrary to []
		std::cout << "Using device: " << device.getInfo<CL_DEVICE_NAME>() << std::endl;
	}
	catch (const std::exception& e) {
		std::cout << "Error getting device: " << e.what() << std::endl;
		std::cout << "Ususally means the program wasn't compiled for desired device!" << std::endl;
		return -1;
	}

	// Initialize simulation engine.
	goPMC::MCEngine mcEngine;
	mcEngine.initializeComputation(platform, device);
  std::cout << "Context created! Initialising physics... ";

	// Read and process physics data.
	mcEngine.initializePhysics("../lut"); // install-prefix/lut relative to install-prefix/bin/gPMC.exe
  std::cout << "Physics initialised! Now reading dicom... ";
	// Read and process patient Dicom CT data. ITK has origin at upper left
	mcEngine.initializePhantom(args_info.path_arg); // "090737"); // "directoryToDicomData");
  std::cout << "Dicom read!" << std::endl;

	// Initialize source protons with arrays of energy (T), position (pos), direction (dir) and weight (weight) of each proton.
	// Position and direction should be defined in Dicom CT coordinate.
	cl_float * T = new cl_float[N];     //Energy(MeV?) = [120.0, ..., 120.0]
	cl_float3 * pos = new cl_float3[N]; //Position    = [(5*rand_1-15, -20, 5*rand_1+25), ..., (5*rand_N-15, -20, 5*rand_N+25)]
	cl_float3 * dir = new cl_float3[N]; //Direction   = [(0, 1, 0), ..., (0, 1, 0)] = y?          ^-------- 0 <= rand_X <= 1
	cl_float * weight = new cl_float[N];//Weight      = [1.0, ..., 1.0]
	initSource(T, pos, dir, weight);

	std::cout << "Source initialised, now simulating... ";
	// Choose a physics quantity to score for this simulation run.
	// Scoring quantity could be one of {DOSE2MEDIUM, DOSE2WATER, FLUENCE, LETD}.
	// LETD is dose weighted LET, to get dose averaged LET, divide it by DOSE2MEDIUM from another simulation run.
	std::string quantity("DOSE2WATER");

	// Run simulation.
	mcEngine.simulate(T, pos, dir, weight, N, quantity);
  std::cout << "Simulation complete! Now getting results..." << std::endl;
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
	typedef itk::Image<float, 3> ImageType;
	ImageType::Pointer doseImage = ImageType::New();

	float * im_spc = new float[3];
	im_spc[0] = args_info.spacing_arg[0];
	im_spc[1] = args_info.spacing_arg[1];
	im_spc[2] = args_info.spacing_arg[2];
	doseImage->SetSpacing(im_spc);

	itk::Index<3> im_org;
	im_org[0] = args_info.origin_arg[0]; //manually unrolled to help compiler optimize cache
	im_org[1] = args_info.origin_arg[1];
	im_org[2] = args_info.origin_arg[2];

	// float * im_dim = new float[3];
	itk::Size<3> im_dim;
	im_dim[0] = args_info.dimension_arg[0];
	im_dim[1] = args_info.dimension_arg[1];
	im_dim[2] = args_info.dimension_arg[2];
	ImageType::RegionType region(im_org, im_dim);
	doseImage->SetRegions(region);
	doseImage->Allocate();

	for (size_t i = 0; i < im_dim[0]; i++)
	{
		for (size_t j = 0; j < im_dim[1]; j++)
		{
			for (size_t k = 0; k < im_dim[2]; k++)
			{
				ImageType::IndexType pixelIndex;
				pixelIndex[0] = i;
				pixelIndex[1] = j;
				pixelIndex[2] = k;
				doseImage->SetPixel(pixelIndex, doseMean[i + j * im_dim[0] + k * im_dim[0] * im_dim[1]]);
			}
		}
	}

	std::cout << "Writing output... " << std::endl;
	typedef  itk::ImageFileWriter<ImageType> WriterType;
	WriterType::Pointer outputWriter = WriterType::New();
	outputWriter->SetFileName(args_info.output_arg);
	outputWriter->SetInput(doseImage);
	outputWriter->Update();

	// Clear the scoring counters in previous simulation runs.
	mcEngine.clearCounter();

	delete[] T;
	delete[] pos;
	delete[] dir;
	delete[] weight;
	// "Real" classes has destructors?

	stream = freopen("CON", "w", stdout);
	stream_err = freopen("CON", "w", stderr);
	return 0;
}
