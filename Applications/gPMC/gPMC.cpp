/*
 *
 * file  gPMC.cpp
 * brief program for using goPMC from commandline.
 * specifically made for use with CBCTrecon
 *
 * author Andreas Gravgaard Andersen
 *
 * last update on 26/7/2017
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
#include <itkImageSource.h>
#include "gdcmReader.h"
#include "gdcmAttribute.h"

#define NDOSECOUNTERS 1
#include "goPMC.h"
#define N 1000000
#define N_per_spot 100

#define M_PI 3.1415926535897932384626433832795028841971693993751

#include "gPMC_ggo.h"
#include "range_modulator_data.hxx"

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

// prototypes, real functions is defined after main.
void write_dose_to_mha(std::vector<cl_float> dose, args_info_gPMC &args_info);
void print_sum_and_mean(std::vector<cl_float> dose);

size_t initFromSpotScanning(gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq,
	cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight, size_t &n_total_spots);

size_t initFromPassiveScatter(gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq,
	cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight, size_t &n_total_spots);




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

size_t initSourceFromDicom(args_info_gPMC args_info, cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight){
	size_t total_spots = 0;

	size_t size = 1;
	if (args_info.Nplans_given)
		size = args_info.Nplans_arg;

	for (size_t i = 0; i < size; i++){
		const char * dicom_path = args_info.plan_arg[i];
		// check dicom integrity
		gdcm::Reader reader;
		reader.SetFileName(dicom_path);
		if (!reader.Read())
		{
			std::cout << "Reading dicom plan failed!" << std::endl;
			return 0;
		}
		gdcm::File &file = reader.GetFile();
		gdcm::DataSet &ds = file.GetDataSet();
		//if (ds.FindDataElement(gdcm::Tag(0x10, 0x20)))
		//	const gdcm::DataElement &pt_id = ds.GetDataElement(gdcm::Tag(0x10, 0x20));

		const gdcm::DataElement &beam_seq_tag = ds.GetDataElement(gdcm::Tag(0x300a, 0x3a2));
		gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq = beam_seq_tag.GetValueAsSQ();


		const gdcm::Item &it_beams = beam_seq->GetItem(1);
		gdcm::Attribute< 0x300a, 0x308> at_scan_mode;
		at_scan_mode.SetFromDataElement(it_beams.GetDataElement(at_scan_mode.GetTag()));
		const char* scan_mode = at_scan_mode.GetValue();

		if (!strcmp(scan_mode, "NONE")){ // passive scatter
			size_t total_spots_cur = initFromPassiveScatter(beam_seq, T, pos, dir, weight, total_spots);
			std::cout << total_spots_cur << std::endl;
		}
		else // spot scanning
		{
			size_t total_spots_cur = initFromSpotScanning(beam_seq, T, pos, dir, weight, total_spots);
			std::cout << total_spots_cur << std::endl;
		}
		if (total_spots == 0)
			std::cout << "Something went wrong!! No points returned from dicom reader!!" << std::endl;
	}
	return total_spots;

}

//                        N points      x array       y array      x point      y point
int point_in_aperture(int nvert, float *vertx, float *verty, float testx, float testy)
{
	int i, j, c = 0;
	for (i = 0, j = nvert - 1; i < nvert; j = i++) {
		if (((verty[i]>testy) != (verty[j]>testy)) &&
			(testx < (vertx[j] - vertx[i]) * (testy - verty[i]) / (verty[j] - verty[i]) + vertx[i]))
			c = !c;
	}
	return c; // 0=outside, 1=inside
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

	gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq = beam_seq_tag.GetValueAsSQ();
	for (unsigned int i = 0; i < beam_seq->GetNumberOfItems(); ++i){
		const gdcm::Item &it_beams = beam_seq->GetItem(i + 1);
		gdcm::Attribute< 0x300a, 0x308> at_scan_mode;
		at_scan_mode.SetFromDataElement(it_beams.GetDataElement(at_scan_mode.GetTag()));
		const char* scan_mode = at_scan_mode.GetValue();
		
		if (!strcmp(scan_mode, "NONE")){ // passive scatter
			// Range compensators:
			const gdcm::DataElement &rng_comp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x2ea));
			gdcm::SmartPointer<gdcm::SequenceOfItems> rng_comp_seq = rng_comp_seq_tag.GetValueAsSQ();
			const gdcm::Item &it_rng_comp = rng_comp_seq->GetItem(1);

			gdcm::Attribute< 0x300a, 0xe7> at_n_rows; // number of rows in rng comp. -> x-direction.
			at_n_rows.SetFromDataElement(it_rng_comp.GetDataElement(at_n_rows.GetTag()));
			const unsigned int n_rows = at_n_rows.GetValue(); // value can't be signed, ignore warning
			gdcm::Attribute< 0x300a, 0xe8> at_n_cols; // number of columns in rng comp. -> y-direction.
			at_n_cols.SetFromDataElement(it_rng_comp.GetDataElement(at_n_cols.GetTag()));
			const unsigned int n_cols = at_n_cols.GetValue(); // value can't be signed, ignore warning

			gdcm::Attribute< 0x300a, 0xe9> at_spacing; // spacing of rows and columns in rng comp.
			at_spacing.SetFromDataElement(it_rng_comp.GetDataElement(at_spacing.GetTag()));
			const double* spacing = at_spacing.GetValues(); // mm
			gdcm::Attribute< 0x300a, 0xea> at_offset; // position of 0,0 in rng comp. -> upper lefthand.
			at_offset.SetFromDataElement(it_rng_comp.GetDataElement(at_offset.GetTag()));
			const double* offset = at_offset.GetValues(); // mm

			// Block:
			const gdcm::DataElement &block_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a6));
			gdcm::SmartPointer<gdcm::SequenceOfItems> block_seq = block_seq_tag.GetValueAsSQ();
			const gdcm::Item &it_block = block_seq->GetItem(1);

			gdcm::Attribute< 0x300a, 0x104> at_block_n_points;
			at_block_n_points.SetFromDataElement(it_block.GetDataElement(at_block_n_points.GetTag()));
			const unsigned int block_n_points = at_block_n_points.GetValue(); // number of x-y pairs defining block

			gdcm::Attribute< 0x300a, 0x106> at_block_points;
			at_block_points.SetFromDataElement(it_block.GetDataElement(at_block_points.GetTag()));
			const double* block_points = at_block_points.GetValues(); // x-y pairs defining block
			float* block_points_x = new float[block_n_points + 1];
			float* block_points_y = new float[block_n_points + 1];

#pragma omp parallel for
			for (int j = 0; j < block_n_points; j++){
				block_points_x[j] = block_points[j * 2] + (n_rows / 2 * spacing[0] - offset[0]);
				block_points_y[j] = block_points[j * 2 + 1] + (n_cols / 2 * spacing[1] - offset[1]);
			}
			block_points_x[block_n_points] = block_points[0] + (n_rows / 2 * spacing[0] - offset[0]);; //close the loop
			block_points_y[block_n_points] = block_points[0 + 1] + (n_cols / 2 * spacing[1] - offset[1]);;

			size_t n_spots_in_aperture = 0;
			for (int j = 0; j < n_cols; j++){
				const cl_float cur_col_pos = j * spacing[1] - offset[1];
				for (int k = 0; k < n_rows; k++){ // k + j*n_rows
					if (1 == point_in_aperture(block_n_points, block_points_x, block_points_y,
						k * spacing[0] - offset[0], cur_col_pos))
						n_spots_in_aperture++;
				}
			}

			delete[] block_points_x;
			delete[] block_points_y;
			// end block

			size_t n_mod_wepl = 0;
			const gdcm::DataElement &cp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a8));
			gdcm::SmartPointer<gdcm::SequenceOfItems> control_point_seq = cp_seq_tag.GetValueAsSQ();

			const gdcm::Item &it_control_points = control_point_seq->GetItem(1);

			// Range modulator settings:
			const gdcm::DataElement &rng_mod_set_seq_tag = it_control_points.GetDataElement(gdcm::Tag(0x300a, 0x380));
			gdcm::SmartPointer<gdcm::SequenceOfItems> rng_mod_set_seq = rng_mod_set_seq_tag.GetValueAsSQ();

			for (unsigned int k = 0; k < rng_mod_set_seq->GetNumberOfItems(); ++k){
				const gdcm::Item &it_rng_mod_set = rng_mod_set_seq->GetItem(k + 1);
				gdcm::Attribute< 0x300a, 0x382> at_mod_start; // range modulator start index
				at_mod_start.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_start.GetTag()));
				const unsigned int mod_start = at_mod_start.GetValue();

				gdcm::Attribute< 0x300a, 0x384> at_mod_end; // range modulator end index
				at_mod_end.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_end.GetTag()));
				const unsigned int mod_end = at_mod_end.GetValue();



				// Range modulator:
				const gdcm::DataElement &rng_mod_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x342));
				gdcm::SmartPointer<gdcm::SequenceOfItems> rng_mod_seq = rng_mod_seq_tag.GetValueAsSQ();

				const gdcm::Item &it_rng_mod = rng_mod_seq->GetItem(1);

				gdcm::Attribute< 0x300a, 0x346> at_mod_ID; // range modulator type
				at_mod_ID.SetFromDataElement(it_rng_mod.GetDataElement(at_mod_ID.GetTag()));
				const char* mod_ID = at_mod_ID.GetValue();
				if (!isdigit(mod_ID[3])){
					std::cout << "\a" << "Non-standard format of range modulator id: " << mod_ID << " source code edit may be necessary." << std::endl;
					return 0;
				}

				const size_t ID_num = mod_ID[3] - '0'; // apparently the best way to convert from char to int

				//         range modulator angles    range modulator heights   material     lead or alu angles        lead or alu heights
				std::tuple<std::vector<const double>, std::vector<const double>, const char*, std::vector<const double>, std::vector<const double>>
					range_modulator = get_range_modulator(ID_num);

				size_t n_lead_alu_only = 0;
				for (size_t j = 0; j < std::get<3>(range_modulator).size(); j++)
					if (std::get<3>(range_modulator)[j] < std::get<0>(range_modulator)[0] && std::get<3>(range_modulator)[j] > mod_start)
						n_lead_alu_only++;

				const size_t n_rng_mod = std::get<0>(range_modulator).size() + n_lead_alu_only;

				size_t n_non_zero_weights = 0;
				// Range modulator energy difference
#pragma omp parallel for
				for (int j = 0; j < n_lead_alu_only; j++){
					if (j < (n_lead_alu_only - 1)){
						if (mod_start <= std::get<3>(range_modulator)[j]){
							if (std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j + 1])
								n_non_zero_weights++;
							else if (std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j])
								n_non_zero_weights++;
						}
						else if ((std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j + 1]) && (mod_start < std::get<3>(range_modulator)[j + 1]))
							n_non_zero_weights++;
					}
					else if ((std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j]) && (mod_start <= std::get<3>(range_modulator)[j]))
							n_non_zero_weights++;
				}

#pragma omp parallel for
				for (int j = n_lead_alu_only; j < n_rng_mod; j++){
					size_t j_RM = j - n_lead_alu_only;
					if (j < (n_rng_mod - 1)){
						if (mod_start <= std::get<0>(range_modulator)[j_RM]){
							if (mod_end >= std::get<0>(range_modulator)[j_RM + 1])
								n_non_zero_weights++;
							else if (mod_end >= std::get<0>(range_modulator)[j_RM])
								n_non_zero_weights++;
						}
						else if ((mod_end >= std::get<0>(range_modulator)[j_RM + 1]) && (mod_start < std::get<0>(range_modulator)[j_RM + 1]))
							n_non_zero_weights++;
					}
					else if ((mod_start <= std::get<0>(range_modulator)[j_RM]) && (mod_end >= std::get<0>(range_modulator)[j_RM]))
							n_non_zero_weights++;
				}

				n_mod_wepl += (n_non_zero_weights * N_per_spot);
			}
			total_n_spots += n_spots_in_aperture * n_mod_wepl;
		}
		else // spot scanning
		{
			const gdcm::DataElement &cp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a8));
			gdcm::SmartPointer<gdcm::SequenceOfItems> control_point_seq = cp_seq_tag.GetValueAsSQ();
			
			for (unsigned int j = 0; j < control_point_seq->GetNumberOfItems(); ++j){
				const gdcm::Item &it_control_points = control_point_seq->GetItem(j + 1);
				gdcm::Attribute< 0x300a, 0x392> at_n_spots;
				at_n_spots.SetFromDataElement(it_control_points.GetDataElement(at_n_spots.GetTag()));
				total_n_spots += at_n_spots.GetValue()  * N_per_spot;
			}
		}
	}
	return total_n_spots;
}

std::tuple<std::vector<cl_float>, std::vector<cl_float>> simulate(goPMC::MCEngine &mcEngine, args_info_gPMC &args_info, const size_t N_dicom)
{
	cl_float * T;
	cl_float3 * pos;
	cl_float3 * dir;
	cl_float * weight;

	if (!args_info.plan_given)
	{
		T = new cl_float[N];     //Energy(MeV?)= [120.0, ..., 120.0]
		pos = new cl_float3[N];  //Position    = [(5*rand_1-15, -20, 5*rand_1+25), ..., (5*rand_N-15, -20, 5*rand_N+25)]
		dir = new cl_float3[N];  //Direction   = [(0, 1, 0), ..., (0, 1, 0)] = y?          ^-------- 0 <= rand_X <= 1
		weight = new cl_float[N];//Weight      = [1.0, ..., 1.0]

		initSource(T, pos, dir, weight);
		if (args_info.verbose_flag)
			std::cout << "Source initialised, now simulating... ";
	}
	else
	{
		T = new cl_float[N_dicom];     //Energy(MeV?)= nominal control point energy
		pos = new cl_float3[N_dicom];  //Position    = Scan Spot Position transformed to gPMC coordinates
		dir = new cl_float3[N_dicom];  //Direction   = vector in unit sphere defined by the gantry and couch angle 
		weight = new cl_float[N_dicom];//Weight      = Scan Spot Meterset Weights of control point
		
		size_t N_result = initSourceFromDicom(args_info, T, pos, dir, weight);
		if (N_result != N_dicom){
			std::cout << "\a" << "ONE OF THE COUNTERS ARE WRONG!!" << std::endl;
			std::cout << "dicom: " << N_dicom << ", result: " << N_result << std::endl;
		}

		if (args_info.spacing_arg[2]) // not necessary as it is now, but more tests are needed, so lets keep it a little longer...
		{
#pragma omp parallel for
			for (int i = 0; i < N_dicom; i++)
			{
				pos[i].s[0] *= 0.1; //2.0*args_info.spacing_arg[0];
				pos[i].s[1] *= 0.1; //2.0*args_info.spacing_arg[1];
				pos[i].s[2] *= 0.1; // 2.0*args_info.spacing_arg[2];
			}
		}
	}

	// Choose a physics quantity to score for this simulation run.
	// Scoring quantity could be one of {DOSE2MEDIUM, DOSE2WATER, FLUENCE, LETD}.
	// LETD is dose weighted LET, to get dose averaged LET, divide it by DOSE2MEDIUM from another simulation run.
	//if (!strcmp(args_info.metric_arg, "dose2water"))
	std::string quantity("DOSE2WATER");
	//	quantity = "DOSE2WATER";
	if (!strcmp(args_info.metric_arg, "dose2medium"))
		quantity = "DOSE2MEDIUM";
	else if (!strcmp(args_info.metric_arg, "let"))
		quantity = "LETD";
	else if (!strcmp(args_info.metric_arg, "fluence"))
		quantity = "FLUENCE";

	// Run simulation.
	mcEngine.simulate(T, pos, dir, weight,
		(!args_info.plan_given ? N : N_dicom),
		quantity);

	if (args_info.verbose_flag)
		std::cout << "Simulation complete! Now getting results..." << std::endl;
	// Get simulation results.
	std::vector<cl_float> doseMean, doseStd;
	mcEngine.getResult(doseMean, doseStd);

	if (!strcmp(args_info.metric_arg, "let")){
		quantity = "DOSE2MEDIUM";
		mcEngine.simulate(T, pos, dir, weight,
			(!args_info.plan_given ? N : N_dicom),
			quantity);
		// Get simulation results.
		std::vector<cl_float> dmedMean, dmedStd;
		mcEngine.getResult(dmedMean, dmedStd);

#pragma omp parallel for
		for (int i = 0; i < dmedMean.size(); i++)
			doseMean[i] /= dmedMean[i];
	}

	delete[] T;
	delete[] pos;
	delete[] dir;
	delete[] weight;

	return std::make_tuple(doseMean, doseStd);
}

int main(int argc, char * argv[])
{
	GGO(gPMC, args_info);

	// Get OpenCL platform and device.
	cl::Platform platform;
	cl::Platform::get(&platform);

	if (args_info.verbose_flag)
		std::cout << "Using platform: " << platform.getInfo<CL_PLATFORM_NAME>() << std::endl;

	std::vector<cl::Device> devs;

	if (!strcmp(args_info.hardware_arg, "gpu"))
		std::cout << "Getting device GPU returned: " << platform.getDevices(CL_DEVICE_TYPE_GPU, &devs) << std::endl;
	else if (!strcmp(args_info.hardware_arg, "cpu"))
		std::cout << "Getting devices CPU returned: " << platform.getDevices(CL_DEVICE_TYPE_CPU, &devs) << std::endl;
	else if (!strcmp(args_info.hardware_arg, "acc"))
		std::cout << "Getting devices ACCELERATOR returned: " << platform.getDevices(CL_DEVICE_TYPE_ACCELERATOR, &devs) << std::endl;
	else
		std::cout << "Getting devices DEFAULT returned: " << platform.getDevices(CL_DEVICE_TYPE_DEFAULT, &devs) << std::endl;

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

	// Read and process patient Dicom CT data. ITK has origin at upper left
	std::string dicom_path = args_info.dir_arg; // to help debugging

	std::cout << dicom_path << std::endl;
	// Initialize source protons with arrays of energy (T), position (pos), direction (dir) and weight (weight) of each proton.
	// Position and direction should be defined in Dicom CT coordinate.
	
	size_t N_dicom = 0;
	if (args_info.plan_given){
		size_t size = 1;
		if (args_info.Nplans_given)
			size = args_info.Nplans_arg;

		for (size_t i = 0; i < size; i++){
			std::cout << args_info.plan_arg[i] << "\nReading dicom RT plan..." << std::endl;
			N_dicom += GetNFromDicom(args_info.plan_arg[i]);

			if (N_dicom == 69) return -1; // yes, 69 is a joke, but 69 is a highly improbable and specific number of spots.

			if (args_info.verbose_flag)
				std::cout << " Dicom plan readable! " << N_dicom << " spots will be simulated." << std::endl;
		}
	}

	if (args_info.verbose_flag){
		std::cout << "Context created! Initialising physics... " << "\n";

		std::cout << "Using LookUpTable path: " << (!args_info.lut_arg ? "../lut" : std::string(args_info.lut_arg)) << std::endl;
	}
	// Read and process physics data.
	mcEngine.initializePhysics(!args_info.lut_arg ? "../lut" : std::string(args_info.lut_arg)); // Look Up Tables path, relative path accepted

	if (args_info.verbose_flag)
		std::cout << "Physics initialised! Now reading dicom... " << std::endl;

	try{
		mcEngine.initializePhantom(dicom_path); // "090737"); // "directoryToDicomData");
	}
	catch (const std::exception& e) {
		std::cout << "Error reading dicom: " << e.what() << "\n";
		std::cout << "Ususally means the implemented writer didn't give a compatible dicom image!" << std::endl;
		return -1;
	}
	if (args_info.verbose_flag){
		std::cout << "Dicom images read!";
	}

	std::tuple<std::vector<cl_float>, std::vector<cl_float>> dose_tuple = simulate(mcEngine, args_info, N_dicom);

	std::vector<cl_float> doseMean = std::get<0>(dose_tuple);
	std::vector<cl_float> doseStd = std::get<1>(dose_tuple);

	if (args_info.verbose_flag){
		// Calculate mean of mean
		std::cout << "doseMean size: " << doseMean.size();
		print_sum_and_mean(doseMean);

		// Calculate mean of SD
		std::cout << "doseStd size: " << doseStd.size();
		print_sum_and_mean(doseStd);
	}

	if (args_info.batch_given && platform.getInfo<CL_PLATFORM_NAME>().compare("NVIDIA") <= 0){
		for (size_t i = 0; i < args_info.batch_arg; i++){
			std::tuple<std::vector<cl_float>, std::vector<cl_float>> dose_tuple_tmp = simulate(mcEngine, args_info, N_dicom);
			std::vector<cl_float> doseMean_tmp = std::get<0>(dose_tuple_tmp);
			std::vector<cl_float> doseStd_tmp = std::get<1>(dose_tuple_tmp);

#pragma omp parallel for
			for (int j = 0; j < doseMean.size(); j++){
				doseMean[j] += doseMean_tmp[j];
				// doseStd[j] += laplace(doseMean_tmp,j)^2 * doseStd_tmp[j]^2; (laplace() should be the derivative of the function (3d-matrix) at j)
			}
			// sqrt(doseStd)
		}
	}

	// Do something with doseMean and doseStd //
	if (args_info.verbose_flag){
		// Calculate mean of mean
		std::cout << "doseMean size: " << doseMean.size();
		print_sum_and_mean(doseMean);

		// Calculate mean of SD
		std::cout << "doseStd size: " << doseStd.size();
		print_sum_and_mean(doseStd);
	}
	// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% //
	write_dose_to_mha(doseMean, args_info); // IO function

	// Clear the scoring counters in previous simulation runs.
	mcEngine.clearCounter();

	return 0;
}

void print_sum_and_mean(std::vector<cl_float> dose){
	double sum = 0;

#pragma omp parallel for reduction(+:sum)
	for (int i = 0; i < dose.size(); i++)
		sum += dose[i];

	printf(" sum: %.10e", sum);
	printf(" mean: %.10e\n", sum / dose.size());
}



double inverse_cdf(double random){
	for (size_t i = 0; i < (sizeof(cdf_xs) / sizeof(double)); i++)
		if (random <= cdf_xs[i]) // cdf_xs goes from 0 to 1
			return cdf_val[i];
	return cdf_val[sizeof(cdf_val) / sizeof(double)]; // if we have not returned yet, return the last element (impossible in theory, but compiler doesn't know that)
}

size_t initFromSpotScanning(gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq,
	cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight, size_t& total_spots){

	typedef itk::Euler3DTransform< double > TransformType;
	TransformType::ParametersType fixedParam(3); //rotation center
	fixedParam.put(0, 0);
	fixedParam.put(1, 0);
	fixedParam.put(2, 0);

	const double halfC = M_PI / 180.0;

	for (unsigned int i = 0; i < beam_seq->GetNumberOfItems(); ++i) { // LOOPING BEAMS, assumes there are at least one beam
		const gdcm::Item &it_beams = beam_seq->GetItem(i + 1);
		gdcm::Attribute<0x300a, 0xC6> at_rt_type;
		at_rt_type.SetFromDataElement(it_beams.GetDataElement(at_rt_type.GetTag()));
		const std::string rt_type(at_rt_type.GetValue());

		float one_over_charge = 1.0; // Energy is float
		if (rt_type.find("ION") != std::string::npos)
		{
			gdcm::Attribute<0x300a, 0x306> at_charge;
			at_charge.SetFromDataElement(it_beams.GetDataElement(at_charge.GetTag()));
			one_over_charge = 1.0 / at_charge.GetValue();
			std::cout << "Warning: not proton plan! -> dividing energy by charge for approximation: " << at_charge.GetValue() << std::endl;
		}
		const gdcm::DataElement &cp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a8));
		gdcm::SmartPointer<gdcm::SequenceOfItems> control_point_seq = cp_seq_tag.GetValueAsSQ();

		const gdcm::Item &first_control_points = control_point_seq->GetItem(1);
		// While we're on the first control point we get the data that may only be defined here
		gdcm::Attribute<0x300a, 0x11e> at_gantry;
		at_gantry.SetFromDataElement(first_control_points.GetDataElement(at_gantry.GetTag()));
		const double gantry = at_gantry.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x122> at_couch;
		at_couch.SetFromDataElement(first_control_points.GetDataElement(at_couch.GetTag()));
		const double couch = at_couch.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x12c> at_isocenter;
		at_isocenter.SetFromDataElement(first_control_points.GetDataElement(at_isocenter.GetTag()));
		const double* isocenter = at_isocenter.GetValues();

		const gdcm::DataElement &range_shifter_seq_tag = first_control_points.GetDataElement(gdcm::Tag(0x300a, 0x360));
		gdcm::SmartPointer<gdcm::SequenceOfItems> range_shifter_seq = range_shifter_seq_tag.GetValueAsSQ();
		const gdcm::Item &it_rng_shifter = range_shifter_seq->GetItem(1);
		gdcm::Attribute<0x300a, 0x364> at_sid; // Isocenter to Rangeshifter Distance
		at_sid.SetFromDataElement(it_rng_shifter.GetDataElement(at_sid.GetTag()));
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
		//transform->SetRotation(0, couch, gantry);
		//      ASSUMING:         3., 2., 1.
		transform->SetRotation(couch,  0, gantry);
		transform->SetFixedParameters(fixedParam); //Center of the Transform
		const itk::Matrix<double, 3U, 3U> A = transform->GetMatrix();

		const vnl_vector_fixed<double, 3> beam_offset(
			-isocenter[0] - direction[0] * sid,
			-isocenter[1] - direction[1] * sid,
			isocenter[2] - direction[2] * sid
			);

		for (unsigned int j = 0; j < control_point_seq->GetNumberOfItems(); ++j) { // LOOPING-CONTROL POINTS, assumes there are at least one control point per beam.
			const gdcm::Item &it_control_points = control_point_seq->GetItem(j + 1);
			gdcm::Attribute<0x300a, 0x114> at_nom_beam_energy;
			at_nom_beam_energy.SetFromDataElement(it_control_points.GetDataElement(at_nom_beam_energy.GetTag()));
			const float nom_beam_energy = at_nom_beam_energy.GetValue() * one_over_charge; // gPMC uses float and it's not used for anything else.

			gdcm::Attribute< 0x300a, 0x392> at_n_spots;
			at_n_spots.SetFromDataElement(it_control_points.GetDataElement(at_n_spots.GetTag()));
			const unsigned int n_spots = at_n_spots.GetValue(); // value can't be signed, ignore warning

			gdcm::Attribute<0x300a, 0x394> at_spot_pos_map; // x, y - positions of spots
			at_spot_pos_map.SetFromDataElement(it_control_points.GetDataElement(at_spot_pos_map.GetTag()));
			const float* p_spot_pos_map = at_spot_pos_map.GetValues(); // gPMC uses float and it's not used for anything else.
			const std::vector<float> spot_pos_map{ p_spot_pos_map, p_spot_pos_map + at_spot_pos_map.GetNumberOfValues() };

			gdcm::Attribute<0x300a, 0x396> at_spot_meterset_w;
			at_spot_meterset_w.SetFromDataElement(it_control_points.GetDataElement(at_spot_meterset_w.GetTag()));
			const float* p_spot_w = at_spot_meterset_w.GetValues(); // gPMC uses float and it's not used for anything else.
			const std::vector<float> spot_meterset_w{ p_spot_w, p_spot_w + at_spot_meterset_w.GetNumberOfValues() };

			gdcm::Attribute<0x300a, 0x398> at_spot_size;
			at_spot_size.SetFromDataElement(it_control_points.GetDataElement(at_spot_size.GetTag()));
			const float* spot_size = at_spot_size.GetValues();

			unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
			std::minstd_rand0 g1(seed1);  // minstd_rand0 is a standard linear_congruential_engine

#pragma omp parallel for
			for (int i_spot_xy = 0; i_spot_xy < (2 * n_spots); i_spot_xy += 2)
			{
				for (size_t k = 0; k < 100; k++)
				{
					const unsigned int i_local = total_spots + (i_spot_xy / 2 * N_per_spot) + k;
					// Dir(ection) is the opposite (unit)vector of that pointing from iso->source
					// Theta is the angle from z->phi and phi is angle from x->y
					// theta is gantry if z is towards anterior and phi is couch if y is cranial and x is lateral right
					dir[i_local].s[0] = direction[0]; // x-lateral   = -sin(theta)*cos(phi)
					dir[i_local].s[1] = direction[1]; // y-AP        = -cos(theta)
					dir[i_local].s[2] = direction[2]; // z-CC        = -sin(theta)*sin(phi)

					const vnl_vector_fixed<double, 3> point( // Lefthanded x-y-z y pointing caudal for gantry=0 || couch=0, z towards isocenter, but z is 0, so we can treat it as right handed during rotation IF we invert x-axis dduring rotation
						spot_pos_map[i_spot_xy] + 0.1 * spot_size[0] * inverse_cdf(float(g1()) / g1.max()),
						0.0,
						-spot_pos_map[i_spot_xy + 1] + 0.1 * spot_size[1] * inverse_cdf(float(g1()) / g1.max())
						);

					const vnl_vector_fixed<double, 3> point_after_transform = (A * point) + beam_offset;
					pos[i_local].s[0] = point_after_transform[0]; // x-lateral
					pos[i_local].s[1] = point_after_transform[1]; // y-AP
					pos[i_local].s[2] = point_after_transform[2]; // z-CC

					T[i_local] = nom_beam_energy;
					weight[i_local] = spot_meterset_w[i_spot_xy / 2];
				}
			} // end for and omp parallel for
			total_spots += n_spots * N_per_spot; // for enabeling full parallelism

		}
	}

	return total_spots;
}


size_t initFromPassiveScatter(gdcm::SmartPointer<gdcm::SequenceOfItems> beam_seq,
	cl_float * T, cl_float3 * pos, cl_float3 * dir, cl_float * weight, size_t& n_total_spots){

	typedef itk::Euler3DTransform< double > TransformType;
	TransformType::ParametersType fixedParam(3); //rotation center
	fixedParam.put(0, 0);
	fixedParam.put(1, 0);
	fixedParam.put(2, 0);

	const double halfC = M_PI / 180.0;
	for (unsigned int i = 0; i < beam_seq->GetNumberOfItems(); ++i) { // LOOPING BEAMS, assumes there are at least one beam
		const gdcm::Item &it_beams = beam_seq->GetItem(i + 1);
		gdcm::Attribute<0x300a, 0xC6> at_rt_type;
		at_rt_type.SetFromDataElement(it_beams.GetDataElement(at_rt_type.GetTag()));
		const std::string rt_type(at_rt_type.GetValue());

		float one_over_charge = 1.0; // Energy is float
		if (rt_type.find("ION") != std::string::npos)
		{
			gdcm::Attribute<0x300a, 0x306> at_charge;
			at_charge.SetFromDataElement(it_beams.GetDataElement(at_charge.GetTag()));
			one_over_charge = 1.0 / at_charge.GetValue();
			std::cout << "Warning: not proton plan! -> dividing energy by charge for approximation: " << at_charge.GetValue() << std::endl;
		}
		// Range compensators:
		const gdcm::DataElement &rng_comp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x2ea));
		gdcm::SmartPointer<gdcm::SequenceOfItems> rng_comp_seq = rng_comp_seq_tag.GetValueAsSQ();

		const gdcm::Item &it_rng_comp = rng_comp_seq->GetItem(1);

		//while (it_rng_comp < rng_comp_seq->End()) // can there be more than one?
		//{

		gdcm::Attribute< 0x300a, 0x2e4> at_cid; // compensator isocenter distance in mm
		at_cid.SetFromDataElement(it_rng_comp.GetDataElement(at_cid.GetTag()));
		const float sid = at_cid.GetValue(); // mm
		// for source geometry

		gdcm::Attribute< 0x300a, 0x2e0> at_div; // compensator divergence (ABSENT or PRESENT)
		at_div.SetFromDataElement(it_rng_comp.GetDataElement(at_div.GetTag()));
		const char* div = at_div.GetValue();
		if (!strcmp(div, "PRESENT"))
			std::cout << "WARNING: COMPENSATOR DIVERGENCE PRESENT FOR BUT NOT ACCOUNTED FOR!" << std::endl;

		// gdcm::Attribute< 0x300a, 0x2e1> at_mnt; // compensator mounting position (PATIENT_ or SOURCE_SIDE, or DOUBLE_SIDED)
		// at_mnt.SetFromDataElement(it_rng_comp->GetDataElement(at_mnt.GetTag()));
		// const char* mnt = at_mnt.GetValue();
		// should be taken into account?

		gdcm::Attribute< 0x300a, 0xe7> at_n_rows; // number of rows in rng comp. -> x-direction.
		at_n_rows.SetFromDataElement(it_rng_comp.GetDataElement(at_n_rows.GetTag()));
		const unsigned int n_rows = at_n_rows.GetValue(); // value can't be signed, ignore warning
		gdcm::Attribute< 0x300a, 0xe8> at_n_cols; // number of columns in rng comp. -> y-direction.
		at_n_cols.SetFromDataElement(it_rng_comp.GetDataElement(at_n_cols.GetTag()));
		const unsigned int n_cols = at_n_cols.GetValue(); // value can't be signed, ignore warning

		size_t n_spots_max = n_rows*n_cols;
		cl_float2* spot_map = new cl_float2[n_spots_max];

		gdcm::Attribute< 0x300a, 0xe9> at_spacing; // spacing of rows and columns in rng comp.
		at_spacing.SetFromDataElement(it_rng_comp.GetDataElement(at_spacing.GetTag()));
		const double* spacing = at_spacing.GetValues(); // mm
		gdcm::Attribute< 0x300a, 0xea> at_offset; // position of 0,0 in rng comp. -> upper lefthand.
		at_offset.SetFromDataElement(it_rng_comp.GetDataElement(at_offset.GetTag()));
		const double* offset = at_offset.GetValues(); // mm


		// Block:
		const gdcm::DataElement &block_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a6));
		gdcm::SmartPointer<gdcm::SequenceOfItems> block_seq = block_seq_tag.GetValueAsSQ();
		const gdcm::Item &it_block = block_seq->GetItem(1);

		gdcm::Attribute< 0x300a, 0xf8> at_block_type;
		at_block_type.SetFromDataElement(it_block.GetDataElement(at_block_type.GetTag()));
		const char* block_type = at_block_type.GetValue(); // APERTURE or SHIELDING
		if (!strcmp(block_type, "SHIELDING"))
			std::cout << "WARNING: BLOCK TAGGED " << block_type << ", ONLY APERTURE MODE IMPLEMENTED!" << std::endl;

		gdcm::Attribute< 0x300a, 0xe1> at_block_mat;
		at_block_mat.SetFromDataElement(it_block.GetDataElement(at_block_mat.GetTag()));
		const char* block_mat = at_block_mat.GetValue(); // Brass probably
		gdcm::Attribute< 0x300a, 0x100> at_block_thicc;
		at_block_thicc.SetFromDataElement(it_block.GetDataElement(at_block_thicc.GetTag()));
		const double block_thicc = at_block_thicc.GetValue(); // thickness in mm
		std::cout << "Block material: " << block_mat << " of " << block_thicc << " mm thickness";
		std::cout << " assumed to stop beam completely." << std::endl;

		gdcm::Attribute< 0x300a, 0x104> at_block_n_points;
		at_block_n_points.SetFromDataElement(it_block.GetDataElement(at_block_n_points.GetTag()));
		const unsigned int block_n_points = at_block_n_points.GetValue(); // number of x-y pairs defining block

		gdcm::Attribute< 0x300a, 0x106> at_block_points;
		at_block_points.SetFromDataElement(it_block.GetDataElement(at_block_points.GetTag()));
		const double* block_points = at_block_points.GetValues(); // x-y pairs defining block
		float* block_points_x = new float[block_n_points + 1];
		float* block_points_y = new float[block_n_points + 1];

#pragma omp parallel for
		for (int j = 0; j < block_n_points; j++){
			block_points_x[j] = block_points[j * 2] + (n_rows / 2 * spacing[0] - offset[0]);
			block_points_y[j] = block_points[j * 2 + 1] + (n_cols / 2 * spacing[1] - offset[1]);
		}
		block_points_x[block_n_points] = block_points[0] + (n_rows / 2 * spacing[0] - offset[0]);; //close the loop
		block_points_y[block_n_points] = block_points[0 + 1] + (n_cols / 2 * spacing[1] - offset[1]);;


		size_t n_spots_in_aperture = 0;
		for (int j = 0; j < n_cols; j++){
			const cl_float cur_col_pos = j * spacing[1] - offset[1];
			for (int k = 0; k < n_rows; k++){ // k + j*n_rows
				spot_map[n_spots_in_aperture].x = k * spacing[0] - offset[0]; //row pos // unsure if + or - offset
				spot_map[n_spots_in_aperture].y = cur_col_pos; //col pos
				if (1 == point_in_aperture(block_n_points, block_points_x, block_points_y,
					spot_map[n_spots_in_aperture].x, spot_map[n_spots_in_aperture].y))
					n_spots_in_aperture++;
			}
		}

		delete[] block_points_x;
		delete[] block_points_y;

		gdcm::Attribute< 0x300a, 0x2e5> at_hex_offset; // column offset of in rng comp. -> only applicaple for hexagonal compensators.
		at_hex_offset.SetFromDataElement(it_rng_comp.GetDataElement(at_hex_offset.GetTag()));
		const double hex_offset = at_hex_offset.GetValue();
		if (hex_offset != 0)
			std::cout << "WARNING: HEXAGONAL COMPENSATOR DETECTED FOR BUT NOT IMPLEMENTED PROPERLY! expect errors." << std::endl;

		gdcm::Attribute< 0x300a, 0xec> at_thicc; // thickness of rng comp. in mm
		at_thicc.SetFromDataElement(it_rng_comp.GetDataElement(at_thicc.GetTag()));
		const double* thicc = at_thicc.GetValues();
		// DOUBLE SIDED triggers "isocenter to compensator distances" 300a,0x2e6


		gdcm::Attribute< 0x300a, 0x2e7> at_dedx; // Compensator Linear Stopping Power Ratio, ...
		// relative to water, at the beam energy specified by the Nominal Beam Energy (300A,0114) ...
		// of the first Control Point of the Ion Control Point Sequence (300A,03A8).
		at_dedx.SetFromDataElement(it_rng_comp.GetDataElement(at_dedx.GetTag()));
		const double dedx_rel = at_dedx.GetValue(); //relative to water


		gdcm::Attribute< 0x300a, 0xec> at_milling_diameter; // Compensator Milling Tool Diameter
		at_milling_diameter.SetFromDataElement(it_rng_comp.GetDataElement(at_milling_diameter.GetTag()));
		const double mill = at_milling_diameter.GetValue();

		// Range modulator:
		const gdcm::DataElement &rng_mod_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x342));
		gdcm::SmartPointer<gdcm::SequenceOfItems> rng_mod_seq = rng_mod_seq_tag.GetValueAsSQ();

		const gdcm::Item &it_rng_mod = rng_mod_seq->GetItem(1);

		gdcm::Attribute< 0x300a, 0x346> at_mod_ID; // range modulator type
		at_mod_ID.SetFromDataElement(it_rng_mod.GetDataElement(at_mod_ID.GetTag()));
		const char* mod_ID = at_mod_ID.GetValue();
		if (!isdigit(mod_ID[3])){
			std::cout << "\a" << "Non-standard format of range modulator id: " << mod_ID << " source code edit may be necessary." << std::endl;
			return 0;
		}

		const size_t ID_num = mod_ID[3] - '0'; // apparently the best way to convert from char to int

		//         range modulator angles    range modulator heights   material     lead or alu angles        lead or alu heights
		std::tuple<std::vector<const double>, std::vector<const double>, const char*, std::vector<const double>, std::vector<const double>>
			range_modulator = get_range_modulator(ID_num);

		gdcm::Attribute< 0x300a, 0x348> at_mod_type; // range modulator type
		at_mod_type.SetFromDataElement(it_rng_mod.GetDataElement(at_mod_type.GetTag()));
		const char* mod_type = at_mod_type.GetValue();
		if (strcmp(mod_type, "WHL_MODWEIGHTS"))
			std::cout << "\a" << "WARNING: WHL_MODWEIGHTS NOT DEFINED EXPECT FATAL ERRORS!" << std::endl;

		//control points:
		const gdcm::DataElement &cp_seq_tag = it_beams.GetDataElement(gdcm::Tag(0x300a, 0x3a8));
		gdcm::SmartPointer<gdcm::SequenceOfItems> control_point_seq = cp_seq_tag.GetValueAsSQ();
		const gdcm::Item &it_control_points = control_point_seq->GetItem(1);
		// While we're on the first control point we get the data that may only be defined here

		// Range modulator settings:
		const gdcm::DataElement &rng_mod_set_seq_tag = it_control_points.GetDataElement(gdcm::Tag(0x300a, 0x380));
		gdcm::SmartPointer<gdcm::SequenceOfItems> rng_mod_set_seq = rng_mod_set_seq_tag.GetValueAsSQ();

		const gdcm::Item &it_rng_mod_set = rng_mod_set_seq->GetItem(1);

		gdcm::Attribute< 0x300a, 0x382> at_mod_start; // range modulator start index
		at_mod_start.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_start.GetTag()));
		const unsigned int mod_start = at_mod_start.GetValue();

		gdcm::Attribute< 0x300a, 0x384> at_mod_end; // range modulator end index
		at_mod_end.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_end.GetTag()));
		const unsigned int mod_end = at_mod_end.GetValue();

		/* unused
		gdcm::Attribute< 0x300a, 0x386> at_mod_start_thicc; // range modulator thickness at first angle
		at_mod_start_thicc.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_start_thicc.GetTag()));
		const double mod_start_thicc = at_mod_start_thicc.GetValue();

		gdcm::Attribute< 0x300a, 0x388> at_mod_end_thicc; // range modulator thickness at end angle
		at_mod_end_thicc.SetFromDataElement(it_rng_mod_set.GetDataElement(at_mod_end_thicc.GetTag()));
		const double mod_end_thicc = at_mod_end_thicc.GetValue();
		*/

		// True method:
		size_t n_lead_alu_only = 0;
		for (size_t j = 0; j < std::get<3>(range_modulator).size(); j++)
			if (std::get<3>(range_modulator)[j] < std::get<0>(range_modulator)[0] && std::get<3>(range_modulator)[j] > mod_start)
				n_lead_alu_only++;

		const size_t n_rng_mod = std::get<0>(range_modulator).size() + n_lead_alu_only;
		double* mod_energy_diff = new double[n_rng_mod];
		double* mod_energy_weigth = new double[n_rng_mod];
		// const double total_angle = mod_end - mod_start;

		gdcm::Attribute<0x300a, 0x11e> at_gantry;
		at_gantry.SetFromDataElement(it_control_points.GetDataElement(at_gantry.GetTag()));
		const double gantry = at_gantry.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x122> at_couch;
		at_couch.SetFromDataElement(it_control_points.GetDataElement(at_couch.GetTag()));
		const double couch = at_couch.GetValue() * halfC;

		gdcm::Attribute<0x300a, 0x12c> at_isocenter;
		at_isocenter.SetFromDataElement(it_control_points.GetDataElement(at_isocenter.GetTag()));
		const double* isocenter = at_isocenter.GetValues();

		const std::vector<double> direction = {
			std::sin(gantry)*std::cos(couch),
			std::cos(gantry),
			std::sin(gantry)*std::sin(couch)
		};

		//            ( x )                              ( x' )       ( x )             [a b c]
		// IF point = ( y ) THEN point after rotation is ( y' ) = A * ( y ) , WHERE A = [d e f]
		//            ( z )                              ( z' )       ( z )             [g h i]
		TransformType::Pointer transform = TransformType::New();
		//transform->SetRotation(0, couch, gantry);
		//      ASSUMING:         3., 2., 1.
		transform->SetRotation(couch, 0, gantry);
		transform->SetFixedParameters(fixedParam); //Center of the Transform
		const itk::Matrix<double, 3U, 3U> A = transform->GetMatrix();

		const vnl_vector_fixed<double, 3> beam_offset(
			-isocenter[0] - direction[0] * sid,
			-isocenter[1] - direction[1] * sid,
			isocenter[2] - direction[2] * sid
			);

		gdcm::Attribute<0x300a, 0x114> at_nom_beam_energy;
		at_nom_beam_energy.SetFromDataElement(it_control_points.GetDataElement(at_nom_beam_energy.GetTag()));
		const float nom_beam_energy = at_nom_beam_energy.GetValue() * one_over_charge; // gPMC uses float and it's not used for anything else.

		size_t j_W = 0;
		// Range modulator energy difference
		for (size_t j = 0; j < n_lead_alu_only; j++){
			if (std::get<3>(range_modulator)[j] > mod_start){
				mod_energy_diff[j_W] = ((ID_num != 9) ? bethe_lead(nom_beam_energy) : bethe_aluminium(nom_beam_energy))
					* std::get<4>(range_modulator)[j];
			}
			if (j < (n_lead_alu_only - 1)){
				if (mod_start <= std::get<3>(range_modulator)[j]){
					if (std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j + 1])
						mod_energy_weigth[j_W] = (std::get<3>(range_modulator)[j + 1] - std::get<3>(range_modulator)[j]);
					else if (std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j])
						mod_energy_weigth[j_W] = (std::get<0>(range_modulator)[0] - std::get<3>(range_modulator)[j]);
					else
						mod_energy_weigth[j_W] = 0.0;
				}
				else if ((std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j + 1]) && (mod_start < std::get<3>(range_modulator)[j + 1]))
					mod_energy_weigth[j_W] = (std::get<3>(range_modulator)[j + 1] - mod_start);
				else
					mod_energy_weigth[j_W] = 0.0;
			}
			else if ((mod_start <= std::get<3>(range_modulator)[j]) && (std::get<0>(range_modulator)[0] >= std::get<3>(range_modulator)[j]))
				mod_energy_weigth[j_W] = (std::get<0>(range_modulator)[0] - std::get<3>(range_modulator)[j]);
			else
				mod_energy_weigth[j_W] = 0.0;
			if (mod_energy_weigth[j_W] != 0.0)
				j_W++;
		}

		for (int j = n_lead_alu_only; j < n_rng_mod; j++){
			size_t j_RM = j - n_lead_alu_only;
			if (std::get<1>(range_modulator)[j_RM] == NULL)
				mod_energy_diff[j_W] = -nom_beam_energy * 0.99; // kill this point with 1% uncertainty
			else if (!strcmp(std::get<2>(range_modulator), "carbon"))
				mod_energy_diff[j_W] = -bethe_carbon(nom_beam_energy) * std::get<1>(range_modulator)[j_RM];
			else if (!strcmp(std::get<2>(range_modulator), "carbon2"))
				mod_energy_diff[j_W] = -bethe_carbon2(nom_beam_energy) * std::get<1>(range_modulator)[j_RM];
			else if (!strcmp(std::get<2>(range_modulator), "lexan"))
				mod_energy_diff[j_W] = -bethe_lexan(nom_beam_energy) * std::get<1>(range_modulator)[j_RM];

			for (size_t k = 0; k < std::get<3>(range_modulator).size(); k++){ // angle lead or alu foil
				if (std::get<3>(range_modulator)[k] == std::get<0>(range_modulator)[j_RM]){
					mod_energy_diff[j_W] -= ((ID_num != 9) ? bethe_lead(nom_beam_energy) : bethe_aluminium(nom_beam_energy))
						* std::get<4>(range_modulator)[k];
				}	
			}
			if (j < (n_rng_mod - 1)){
				if (mod_start <= std::get<0>(range_modulator)[j_RM]){
					if (mod_end >= std::get<0>(range_modulator)[j_RM + 1])
						mod_energy_weigth[j_W] = (std::get<0>(range_modulator)[j_RM + 1] - std::get<0>(range_modulator)[j_RM]);
					else if (mod_end >= std::get<0>(range_modulator)[j_RM])
						mod_energy_weigth[j_W] = (mod_end - std::get<0>(range_modulator)[j_RM]);
					else
						mod_energy_weigth[j_W] = 0.0;
				}
				else if ((mod_end >= std::get<0>(range_modulator)[j_RM + 1]) && (mod_start < std::get<0>(range_modulator)[j_RM + 1]))
					mod_energy_weigth[j_W] = (std::get<0>(range_modulator)[j_RM + 1] - mod_start);
				else
					mod_energy_weigth[j_W] = 0.0;
			}
			else
			{
				if (mod_end < std::get<0>(range_modulator)[j_RM])
					mod_energy_weigth[j_W] = 0.0;
				else if (mod_start <= std::get<0>(range_modulator)[j_RM])
					mod_energy_weigth[j_W] = (mod_end - std::get<0>(range_modulator)[j_RM]);
			}
			if (mod_energy_weigth[j_W] != 0.0)
				j_W++;
		}


		// increment to second control point to get weight of beam:
		const gdcm::Item &second_control_point = control_point_seq->GetItem(2);

		gdcm::Attribute<0x300a, 0x134> at_weight;
		at_weight.SetFromDataElement(second_control_point.GetDataElement(at_weight.GetTag()));
		const double meterset_weight_per_spot = at_weight.GetValue() / N_per_spot;

#pragma omp parallel for
		for (int i_spot_xy = 0; i_spot_xy < n_spots_in_aperture; i_spot_xy++)
		{

			const cl_float cur_thicc = thicc[i_spot_xy]; //nom_spot_energy = nom_beam_energy - de[i_spot_xy];
			const unsigned int i_group = n_total_spots + (i_spot_xy * j_W * N_per_spot);

			unsigned seed1 = std::chrono::system_clock::now().time_since_epoch().count();
			std::minstd_rand0 g1(seed1);  // minstd_rand0 is a standard linear_congruential_engine

			for (size_t j_energy_step = 0; j_energy_step < j_W; j_energy_step++){
				for (size_t k = 0; k < N_per_spot; k++){
					const unsigned int i_local = i_group + (j_energy_step * N_per_spot) + k;
					// Dir(ection) is the opposite (unit)vector of that pointing from iso->source
					// Theta is the angle from z->phi and phi is angle from x->y
					// theta is gantry if z is towards anterior and phi is couch if y is cranial and x is lateral right
					dir[i_local].s[0] = direction[0]; // x-lateral   = -sin(theta)*cos(phi)
					dir[i_local].s[1] = direction[1]; // y-AP        = -cos(theta)
					dir[i_local].s[2] = direction[2]; // z-CC        = -sin(theta)*sin(phi)

					const vnl_vector_fixed<double, 3> point(
						spot_map[i_spot_xy].x + 0.1*(mill * float(g1()) / g1.max()) - 0.05*mill,
						0.0,
						-spot_map[i_spot_xy].y + 0.1*(mill * float(g1()) / g1.max()) - 0.05*mill
						);
					const vnl_vector_fixed<double, 3> point_after_transform = (A * point) + beam_offset;
					pos[i_local].s[0] = point_after_transform[0]; // x-lateral
					pos[i_local].s[1] = point_after_transform[1]; // y-AP
					pos[i_local].s[2] = point_after_transform[2]; // z-CC

					const cl_float de_comp = -bethe(nom_beam_energy + mod_energy_diff[j_energy_step]) * dedx_rel * cur_thicc;

					T[i_local] = nom_beam_energy + mod_energy_diff[j_energy_step] + de_comp; // mod_energy_diff already negated
					weight[i_local] = meterset_weight_per_spot * mod_energy_weigth[j_energy_step];
					// printf("E: %.3e, W: %.3e\n", T[i_local], weight[i_local]);
				}
			}
		} // end for and omp parallel for

		delete[] mod_energy_diff;
		delete[] mod_energy_weigth;
		delete[] spot_map;

		n_total_spots += n_spots_in_aperture * j_W * N_per_spot; // for enabeling full parallelism
	}

	return n_total_spots;
}


/** \brief Create 3D image from gengetopt specifications. COPIED FROM RTK
*
* This function sets a ConstantImageSource object from command line options stored in ggo struct.
*  The image is not buffered to allow streaming. The image is filled with zeros.
*  The required options in the ggo struct are:
*     - dimension: image size in pixels
*     - spacing: image spacing in coordinate units
*     - origin: image origin in coordinate units
*
* \author Simon Rit
*
* \ingroup Functions
*/

void
SetConstantImageSourceFromGgo(itk::Image<float, 3>::Pointer source, const args_info_gPMC &args_info)
{
	typedef itk::Image<float, 3> ImageType;
	const unsigned int Dimension = ImageType::GetImageDimension();

	ImageType::SizeType imageDimension;
	imageDimension.Fill(args_info.dimension_arg[0]);
	if (args_info.dimension_given){
		imageDimension[0] = args_info.dimension_arg[0] / 2;
		imageDimension[1] = args_info.dimension_arg[1] / 2;
		imageDimension[2] = args_info.dimension_arg[2];
	}

	ImageType::SpacingType imageSpacing;
	imageSpacing.Fill(1.0);
	if (args_info.spacing_given){
		imageSpacing[0] = args_info.spacing_arg[0] * 2.0;
		imageSpacing[1] = args_info.spacing_arg[1] * 2.0;
		imageSpacing[2] = args_info.spacing_arg[2];
	}

	ImageType::PointType imageOrigin;
	for (unsigned int i = 0; i<Dimension; i++)
		imageOrigin[i] = imageSpacing[i] * (imageDimension[i] - 1) * -0.5;
	for (unsigned int i = 0; i<vnl_math_min(args_info.origin_given, Dimension); i++)
		imageOrigin[i] = args_info.origin_arg[i];

	ImageType::DirectionType imageDirection;
	if (args_info.direction_given)
		for (unsigned int i = 0; i<Dimension; i++)
			for (unsigned int j = 0; j<Dimension; j++)
				imageDirection[i][j] = args_info.direction_arg[i*Dimension + j];
	else
		imageDirection.SetIdentity();

	source->SetOrigin(imageOrigin);
	source->SetSpacing(imageSpacing);
	source->SetDirection(imageDirection);
	source->SetRegions(imageDimension);
	// source->SetConstant(0.);

	source->UpdateOutputInformation();
	source->Allocate();
}



void write_dose_to_mha(std::vector<cl_float> dose, args_info_gPMC &args_info){
	typedef itk::Image<float, 3> OutputImageType;
	OutputImageType::Pointer doseImage = OutputImageType::New();

	SetConstantImageSourceFromGgo(doseImage, args_info);

	//doseImage->Allocate();

	unsigned int i = 0;
	itk::ImageRegionIterator<OutputImageType> imIter(doseImage, doseImage->GetLargestPossibleRegion());
	while (!imIter.IsAtEnd())
	{
		imIter.Set(dose[i]);
		++imIter;
		++i;
	}

	if (args_info.verbose_flag)
		std::cout << "Writing output... " << std::endl;

	typedef  itk::ImageFileWriter<OutputImageType> WriterType;
	WriterType::Pointer outputWriter = WriterType::New();
	outputWriter->SetFileName(args_info.output_arg);
	outputWriter->SetInput(doseImage);
	outputWriter->Update();
}
