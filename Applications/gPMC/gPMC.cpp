/*
 *
 * file  gPMC.cpp
 * brief program for using goPMC from commandline.
 * specifically made for use with CBCTrecon
 *
 * author Andreas Gravgaard Andersen
 *
 * last update on 6/6/2017
 *
 */
#include <algorithm>
#include <random>
#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <fstream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#include <itkEuler3DTransform.h>
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "gdcmReader.h"
#include "gdcmAttribute.h"

#define NDOSECOUNTERS 1
#include "goPMC.h"
#define N 100000
#define M_PI 3.1415926535897932384626433832795028841971693993751

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

void initSourceFromDicom(const char * dicom_path, cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight){
	// check dicom integrity
	gdcm::Reader reader;
	reader.SetFileName(dicom_path);
	if (!reader.Read())
	{
		std::cout << "Reading dicom plan failed!" << std::endl;
		return;
	}
	gdcm::File &file = reader.GetFile();
	gdcm::DataSet &ds = file.GetDataSet();
	//if (ds.FindDataElement(gdcm::Tag(0x10, 0x20)))
	//	const gdcm::DataElement &pt_id = ds.GetDataElement(gdcm::Tag(0x10, 0x20));
	
	const gdcm::DataElement &beam_seq_tag = ds.GetDataElement(gdcm::Tag(0x300a, 0x3a2));
	const gdcm::SequenceOfItems *beam_seq = beam_seq_tag.GetValueAsSQ();

	typedef itk::Euler3DTransform< double > TransformType;
	TransformType::ParametersType fixedParam(3); //rotation center
	fixedParam.put(0, 0);
	fixedParam.put(1, 0);
	fixedParam.put(2, 0);

	const double halfC = M_PI / 180.0;
	int i = 0;
	int i_beam = 0;
	gdcm::SequenceOfItems::ConstIterator it_beams = beam_seq->Begin();
	do { // LOOPING BEAMS, assumes there are at least one beam
		gdcm::Attribute<0x300a, 0xC6> at_rt_type;
		at_rt_type.SetFromDataElement(it_beams->GetDataElement(at_rt_type.GetTag()));
		const std::string rt_type(at_rt_type.GetValue());

		float one_over_charge = 1.0; // Energy is float
		if (rt_type.find("ION") != std::string::npos)
		{
			gdcm::Attribute<0x300a, 0x306> at_charge;
			at_charge.SetFromDataElement(it_beams->GetDataElement(at_charge.GetTag()));
			one_over_charge = 1.0 / at_charge.GetValue();
			std::cout << "Warning: not proton plan! -> dividing energy by charge for approximation: " << at_charge.GetValue() << std::endl;
		}
		const gdcm::DataElement &cp_seq_tag = it_beams->GetDataElement(gdcm::Tag(0x300a, 0x3a8));
		const gdcm::SequenceOfItems *control_point_seq = cp_seq_tag.GetValueAsSQ();

		gdcm::SequenceOfItems::ConstIterator it_control_points = control_point_seq->Begin();
		// While we're on the first control point we get the data that may only be defined here
		gdcm::Attribute<0x300a, 0x11e> at_gantry;
		at_gantry.SetFromDataElement(it_control_points->GetDataElement(at_gantry.GetTag()));
		const double gantry = at_gantry.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x122> at_couch;
		at_couch.SetFromDataElement(it_control_points->GetDataElement(at_couch.GetTag()));
		const double couch = at_couch.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x12c> at_isocenter;
		at_isocenter.SetFromDataElement(it_control_points->GetDataElement(at_isocenter.GetTag()));
		const double* isocenter = at_isocenter.GetValues();

		const gdcm::DataElement &range_shifter_seq_tag = it_control_points->GetDataElement(gdcm::Tag(0x300a, 0x360));
		const gdcm::SequenceOfItems *range_shifter_seq = range_shifter_seq_tag.GetValueAsSQ();
		gdcm::SequenceOfItems::ConstIterator it_rng_shifter = range_shifter_seq->Begin();
		gdcm::Attribute<0x300a, 0x364> at_sid; // Isocenter to Rangeshifter Distance
		at_sid.SetFromDataElement(it_rng_shifter->GetDataElement(at_sid.GetTag()));
		const float sid = at_sid.GetValue();

		const std::vector<double> direction = { 
			std::sin(gantry)*std::cos(couch),
			std::cos(gantry),
			std::sin(gantry)*std::sin(couch) 
		};

		//            ( x )                              ( x' )       ( x )             [a b c]
		// IF point = ( y ) THEN point after rotation is ( y' ) = A * ( y ) , WHERE A = [d e f]
		//            ( z )                              ( z' )       ( z )             [g h i]
		TransformType::Pointer transform = TransformType::New();
		transform->SetRotation(0, couch, gantry);
		transform->SetFixedParameters(fixedParam); //Center of the Transform
		const itk::Matrix<double, 3U, 3U> A = transform->GetMatrix();

		const vnl_vector_fixed<double, 3> beam_offset(
			-isocenter[0] - direction[0] * sid,
			-isocenter[1] - direction[1] * sid,
			isocenter[2] - direction[2] * sid
			);

		std::cout << "Beam# " << ++i_beam << " Distance to Rng Shifter: " << sid;
		std::cout << " Gantry,Couch: " << gantry << ", " << couch;
		std::cout << " Direction: (" << direction[0] << ", " << direction[1] << ", " << direction[2] << ") " << std::endl;

		do { // LOOPING-CONTROL POINTS, assumes there are at least one control point per beam.
			gdcm::Attribute<0x300a, 0x114> at_nom_beam_energy;
			at_nom_beam_energy.SetFromDataElement(it_control_points->GetDataElement(at_nom_beam_energy.GetTag()));
			const float nom_beam_energy = at_nom_beam_energy.GetValue() * one_over_charge; // gPMC uses float and it's not used for anything else.

			// Below is only necessary to check the sum of Scan Spot Meterset Weights, we don't do that
			// gdcm::Attribute<0x300a, 0x134> at_cumsum_meterset_w; // weights of spots
			// at_cumsum_meterset_w.SetFromDataElement(it_control_points->GetDataElement(at_cumsum_meterset_w.GetTag()));
			// const double* p_cumsum_w = at_cumsum_meterset_w.GetValues(); 

			gdcm::Attribute< 0x300a, 0x392> at_n_spots;
			at_n_spots.SetFromDataElement(it_control_points->GetDataElement(at_n_spots.GetTag()));
			const int n_spots = at_n_spots.GetValue();

			gdcm::Attribute<0x300a, 0x394> at_spot_pos_map; // x, y - positions of spots
			at_spot_pos_map.SetFromDataElement(it_control_points->GetDataElement(at_spot_pos_map.GetTag()));
			const float* p_spot_pos_map = at_spot_pos_map.GetValues(); // gPMC uses float and it's not used for anything else.
			const std::vector<float> spot_pos_map{ p_spot_pos_map, p_spot_pos_map + at_spot_pos_map.GetNumberOfValues() };

			gdcm::Attribute<0x300a, 0x396> at_spot_meterset_w;
			at_spot_meterset_w.SetFromDataElement(it_control_points->GetDataElement(at_spot_meterset_w.GetTag()));
			const float* p_spot_w = at_spot_meterset_w.GetValues(); // gPMC uses float and it's not used for anything else.
			const std::vector<float> spot_meterset_w{ p_spot_w, p_spot_w + at_spot_meterset_w.GetNumberOfValues() };

			// gdcm::Attribute<0x300a, 0x398> at_spot_size;
			// at_spot_size.SetFromDataElement(it_control_points->GetDataElement(at_spot_size.GetTag()));
			// const double* spot_size = at_spot_size.GetValues();
			
#pragma omp parallel for private(i_spot_xy)
			for (int i_spot_xy = 0; i_spot_xy < 2*n_spots; i_spot_xy+=2)
			{
				const unsigned int i_local = i + i_spot_xy / 2;
				// Dir(ection) is the opposite (unit)vector of that pointing from iso->source
				// Theta is the angle from z->phi and phi is angle from x->y
				// theta is gantry if z is towards anterior and phi is couch if y is cranial and x is lateral right
				dir[i_local].s[0] = direction[0]; // x-lateral   = -sin(theta)*cos(phi)
				dir[i_local].s[1] = direction[1]; // y-AP        = -cos(theta)
				dir[i_local].s[2] = direction[2]; // z-CC        = -sin(theta)*sin(phi)

				const vnl_vector_fixed<double, 3> point(
					spot_pos_map[i_spot_xy], 
					0.0, 
					spot_pos_map[i_spot_xy + 1]
					);

				const vnl_vector_fixed<double, 3> point_after_transform = (A * point) + beam_offset;
				pos[i_local].s[0] = point_after_transform[0]; // x-lateral
				pos[i_local].s[1] = point_after_transform[1]; // y-AP
				pos[i_local].s[2] = point_after_transform[2]; // z-CC

				T[i_local] = nom_beam_energy;
				weight[i_local] = spot_meterset_w[i_local - i];
			} // end for and omp parallel for
			i += n_spots; // for enabeling full parallelism

			++it_control_points;
		} while (it_control_points < control_point_seq->End());

		++it_beams;
	} while (it_beams < beam_seq->End());
}

size_t GetNFromDicom(const char * dicom_path){
	// check dicom integrity
	if (!dicom_path) return 0;

	gdcm::Reader reader;
	reader.SetFileName(dicom_path);
	if (!reader.Read())
	{
		std::cout << "Reading dicom plan failed!" << std::endl;
		return 69;
	}
	gdcm::File &file = reader.GetFile();
	gdcm::DataSet &ds = file.GetDataSet();

	size_t total_n_spots = 0;
	const gdcm::DataElement &beam_seq_tag = ds.GetDataElement(gdcm::Tag(0x300a, 0x3a2));
	const gdcm::SequenceOfItems *beam_seq = beam_seq_tag.GetValueAsSQ();

	gdcm::SequenceOfItems::ConstIterator it_beams = beam_seq->Begin();
	while (it_beams < beam_seq->End())
	{
		const gdcm::DataElement &cp_seq_tag = it_beams->GetDataElement(gdcm::Tag(0x300a, 0x3a8));
		const gdcm::SequenceOfItems *control_point_seq = cp_seq_tag.GetValueAsSQ();

		gdcm::SequenceOfItems::ConstIterator it_control_points = control_point_seq->Begin();
		while (it_control_points < control_point_seq->End())
		{
			gdcm::Attribute< 0x300a, 0x392> at_n_spots;
			at_n_spots.SetFromDataElement(it_control_points->GetDataElement(at_n_spots.GetTag()));
			total_n_spots += at_n_spots.GetValue();
			++it_control_points;
		}
		++it_beams;
	}
	return total_n_spots;
}


int main(int argc, char * argv[])
{
	GGO(gPMC, args_info);

	const std::string stdout_file = std::string(args_info.path_arg) + "\\..\\gPMCstdout.txt";
	const std::string stderr_file = std::string(args_info.path_arg) + "\\..\\gPMCstderr.txt";
	//FILE *stream;
	//if ((stream = freopen(stdout_file.c_str(), "w", stdout)) == NULL)
	//	exit(-1);
	//FILE *stream_err;
	//if ((stream_err = freopen(stderr_file.c_str(), "w", stderr)) == NULL)
	//	exit(-1);
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
	// mcEngine.__autoclassinit2(651384); // just for hacks

	mcEngine.initializeComputation(platform, device);
    std::cout << "Context created! Initialising physics... " << std::endl;

	std::cout << "Using LookUpTable path: " << (!args_info.lut_arg ? "../lut" : std::string(args_info.lut_arg)) << std::endl;
	// Read and process physics data.
	mcEngine.initializePhysics(!args_info.lut_arg ? "../lut" : std::string(args_info.lut_arg)); // Look Up Tables path, relative path accepted
	std::cout << "Physics initialised! Now reading dicom... " << std::endl;
	// Read and process patient Dicom CT data. ITK has origin at upper left
	std::string dicom_path = args_info.path_arg; // to help debugging

	std::cout << dicom_path << std::endl;
	try{
		mcEngine.initializePhantom(dicom_path); // "090737"); // "directoryToDicomData");
	}
	catch (const std::exception& e) {
		std::cout << "Error reading dicom: " << e.what() << std::endl;
		std::cout << "Ususally means the implemented writer didn't give a compatible dicom image!" << std::endl;
		return -1;
	}
	
	std::cout << "Dicom images read!";
	std::cout << (!args_info.plan_arg ? "" : "Reading dicom RT plan...") << std::endl;

	// Initialize source protons with arrays of energy (T), position (pos), direction (dir) and weight (weight) of each proton.
	// Position and direction should be defined in Dicom CT coordinate.

	const size_t N_dicom = GetNFromDicom(args_info.plan_arg);
	if (N_dicom == 69) return -1; // yes, 69 is a joke, but 69 is a highly improbable and specific number of spots.
	std::cout << " Dicom plan readable! " << N_dicom << " spots will be simulated." << std::endl;
	cl_float * T;
	cl_float3 * pos;
	cl_float3 * dir;
	cl_float * weight;

	if (!args_info.plan_arg)
	{
		T = new cl_float[N];     //Energy(MeV?)= [120.0, ..., 120.0]
		pos = new cl_float3[N];  //Position    = [(5*rand_1-15, -20, 5*rand_1+25), ..., (5*rand_N-15, -20, 5*rand_N+25)]
		dir = new cl_float3[N];  //Direction   = [(0, 1, 0), ..., (0, 1, 0)] = y?          ^-------- 0 <= rand_X <= 1
		weight = new cl_float[N];//Weight      = [1.0, ..., 1.0]

		initSource(T, pos, dir, weight);
		std::cout << "Source initialised, now simulating... ";
	}
	else 
	{
		T = new cl_float[N_dicom];     //Energy(MeV?)= nominal control point energy
		pos = new cl_float3[N_dicom];  //Position    = Scan Spot Position transformed to gPMC coordinates
		dir = new cl_float3[N_dicom];  //Direction   = vector in unit sphere defined by the gantry and couch angle 
		weight = new cl_float[N_dicom];//Weight      = Scan Spot Meterset Weights of control point

		initSourceFromDicom(args_info.plan_arg, T, pos, dir, weight);
		if (args_info.spacing_arg[2])
		{
#pragma omp parallel for private(i)
			for (size_t i = 0; i < N_dicom; i++)
			{
				pos[i].s[0] /= 10.0; //2.0*args_info.spacing_arg[0];
				pos[i].s[1] /= 10.0; //2.0*args_info.spacing_arg[1];
				pos[i].s[2] /= 10.0; // 2.0*args_info.spacing_arg[2];
			}
		}
	}

	// Choose a physics quantity to score for this simulation run.
	// Scoring quantity could be one of {DOSE2MEDIUM, DOSE2WATER, FLUENCE, LETD}.
	// LETD is dose weighted LET, to get dose averaged LET, divide it by DOSE2MEDIUM from another simulation run.
	std::string quantity("DOSE2WATER");

	// Run simulation.
	mcEngine.simulate(T, pos, dir, weight, (!args_info.plan_arg ? N : N_dicom), quantity);
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
	im_spc[0] = args_info.spacing_arg[0] * 2.0;
	im_spc[1] = args_info.spacing_arg[1] * 2.0;
	im_spc[2] = args_info.spacing_arg[2];
	doseImage->SetSpacing(im_spc);

	itk::Index<3U> im_org;
	im_org[0] = args_info.origin_arg[0];
	im_org[1] = args_info.origin_arg[1];
	im_org[2] = args_info.origin_arg[2];

	itk::Size<3U> im_dim;
	im_dim[0] = args_info.dimension_arg[0] / 2;
	im_dim[1] = args_info.dimension_arg[1] / 2;
	im_dim[2] = args_info.dimension_arg[2];
	ImageType::RegionType region(im_org, im_dim);
	doseImage->SetRegions(region);
	if (args_info.direction_arg[0] && args_info.direction_arg[8]){
		ImageType::DirectionType ImDirection;
		ImDirection[0][0] = args_info.direction_arg[0];
		ImDirection[0][1] = args_info.direction_arg[1];
		ImDirection[0][2] = args_info.direction_arg[2];
		ImDirection[1][0] = args_info.direction_arg[3];
		ImDirection[1][1] = args_info.direction_arg[4];
		ImDirection[1][2] = args_info.direction_arg[5];
		ImDirection[2][0] = args_info.direction_arg[6];
		ImDirection[2][1] = args_info.direction_arg[7];
		ImDirection[2][2] = args_info.direction_arg[8];
		doseImage->SetDirection(ImDirection);
	}
	doseImage->Allocate();

	unsigned int i = 0;
	itk::ImageRegionIterator<ImageType> imIter(doseImage, region);
	while (!imIter.IsAtEnd())
	{
		imIter.Set(doseMean[i]);
		++imIter;
		++i;
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
	//doseImage->Delete();

	//stream = freopen("CON", "w", stdout);
	//stream_err = freopen("CON", "w", stderr);
	return 0;
}
