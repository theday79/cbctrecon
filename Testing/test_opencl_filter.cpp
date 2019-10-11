// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include "OpenCL/ImageFilters.hpp"

#include <chrono>
#include <cmath>
#include <iostream>

#include <QDir>

#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMedianImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "StructureSet.h"
#include "cbctrecon_io.h"

#include "OpenCL/err_code.hpp"

#define REQ_OPENCL_VER static_cast<cl_uint>(CBCTRECON_OPENCL_VERSION)

// From line integral to raw intensity
class LineInt2Intensity {
public:
  LineInt2Intensity() = default;
  ~LineInt2Intensity() = default;
  float operator()(const float val) const {
    float intensityVal = std::exp(-val) /* I_0=1 */;
    return intensityVal;
  }
};
// From raw intensity to line integral
class Intensity2LineInt {
public:
  Intensity2LineInt() = default;
  ~Intensity2LineInt() = default;
  float operator()(const float val) const {
    // mu = ln(I_0/I) OR mu = ln(I/I0)
    float lineintVal = std::numeric_limits<float>::max();
    if (val > 0) {
      lineintVal = /* log(I_0=1) = 0 */ -std::log(val);
    }
    return lineintVal;
  }
};

// log(ushort_max / val)
class LogInvFunctor {
public:
  LogInvFunctor() = default;
  ~LogInvFunctor() = default;
  float operator()(const unsigned short val) const {
    // mu = ln(I_0/I) OR mu = ln(I/I0)
    const auto log_ushrt_max =
        std::log(std::numeric_limits<unsigned short>::max());
    return log_ushrt_max - std::log(val);
  }
};

template <typename T, size_t DIM> auto GenerateImage(const T init_val = 1) {
  using ImageType = itk::Image<T, DIM>;
  auto image = ImageType::New();
  typename ImageType::IndexType origin;
  origin[0] = 0.0;
  origin[1] = 0.0;
  if (DIM == 3) {
    origin[2] = 0.0;
  }
  typename ImageType::SizeType size;
  size[0] = 1024;
  size[1] = 512 + 256;
  if (DIM == 3) {
    size[2] = 600;
  }
  typename ImageType::RegionType region;
  region.SetIndex(origin);
  region.SetSize(size);

  image->SetRegions(region);
  image->Allocate();

  itk::ImageRegionIterator<ImageType> ImageIter(image, region);
  while (!ImageIter.IsAtEnd()) {
    ImageIter.Set(init_val);
    ++ImageIter;
  }

  return image;
}

template <typename T, size_t DIM>
auto GenerateRandImage(const int seed = 69, const T min_val = 0,
                       const T max_val = std::numeric_limits<T>::max()) {
  using ImageType = itk::Image<T, DIM>;
  auto image = ImageType::New();
  typename ImageType::IndexType origin;
  origin[0] = 0.0;
  origin[1] = 0.0;
  if (DIM == 3) {
    origin[2] = 0.0;
  }
  typename ImageType::SizeType size;
  size[0] = 1024;
  size[1] = 512 + 256;
  if (DIM == 3) {
    size[2] = 600;
  }
  typename ImageType::RegionType region;
  region.SetIndex(origin);
  region.SetSize(size);

  image->SetRegions(region);
  image->Allocate();

  std::srand(seed); // RNG seed
  const auto norm_factor = (max_val - min_val) / static_cast<double>(RAND_MAX);

  itk::ImageRegionIterator<ImageType> ImageIter(image, region);
  while (!ImageIter.IsAtEnd()) {
    // rand is between 0 and RAND_MAX, which is implementation dependend
    ImageIter.Set(static_cast<T>(std::rand() * norm_factor + min_val));
    ++ImageIter;
  }

  return image;
}

UShortImageType::Pointer read_dicom_image(const QString &dcm_dir) {

  auto dir = QDir(dcm_dir);
  const auto filenamelist = get_dcm_image_files(dir);

  ShortImageType::Pointer spShortImg;

  if (!filenamelist.empty()) {
    using dcm_reader_type = itk::ImageSeriesReader<ShortImageType>;
    auto dcm_reader = dcm_reader_type::New();
    const auto dicom_io = itk::GDCMImageIO::New();
    dcm_reader->SetImageIO(dicom_io);
    dcm_reader->SetFileNames(filenamelist);
    dcm_reader->Update();
    spShortImg = dcm_reader->GetOutput();
  }

  auto imageCalculatorFilter =
      itk::StatisticsImageFilter<ShortImageType>::New();

  // Thresholding
  auto thresholdFilter = itk::ThresholdImageFilter<ShortImageType>::New();
  thresholdFilter->SetInput(spShortImg);
  thresholdFilter->ThresholdOutside(-1024, 3072); //--> 0 ~ 4095
  thresholdFilter->SetOutsideValue(-1024);
  thresholdFilter->Update();

  imageCalculatorFilter->SetInput(thresholdFilter->GetOutput());
  imageCalculatorFilter->Update();

  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  const auto outputMinVal = static_cast<unsigned short>(minVal + 1024);
  const auto outputMaxVal = static_cast<unsigned short>(maxVal + 1024);

  auto spRescaleFilter =
      itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();

  return spRescaleFilter->GetOutput();
}

template <typename ImageType>
auto CheckImage(typename ImageType::Pointer image_test,
                typename ImageType::Pointer image_ref) {
  itk::ImageRegionIterator<ImageType> ImageIter_test(
      image_test, image_test->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> ImageIter_ref(
      image_ref, image_ref->GetLargestPossibleRegion());

  auto n_err = 0;
  while (!ImageIter_ref.IsAtEnd()) {
    if (fabs(ImageIter_ref.Get() - ImageIter_test.Get()) > 0.00001f) {
      ++n_err;
    }
    ++ImageIter_ref;
    ++ImageIter_test;
  }
  if (n_err != 0) {
    std::cerr << "Error larger than 0.00001 " << n_err
              << " times, between OpenCL and ITK implementation!\n";
    return -n_err;
  }
  return 0;
}

// Extracts version number so 1.2 -> 120 etc.
cl_uint getOpenCLVersion(const std::string &versionInfo) {
  auto highVersion = 0;
  auto lowVersion = 0;
  auto index = 7;
  while (versionInfo[index] != '.') {
    highVersion *= 10;
    highVersion += versionInfo[index] - '0';
    ++index;
  }
  ++index;
  while (versionInfo[index] != ' ' && versionInfo[index] != '\0') {
    lowVersion *= 10;
    lowVersion += versionInfo[index] - '0';
    ++index;
  }
  // Add a zero after lowVersion if necessary
  if (lowVersion < 10) {
    lowVersion *= 10;
  }
  // Assume lowVersion only contains two digits
  return highVersion * 100 + lowVersion;
}

cl_uint getPlatformOpenCLVersion(cl::Platform &platform) {
  const auto versionInfo = platform.getInfo<CL_PLATFORM_VERSION>();

  return getOpenCLVersion(versionInfo);
}

int main(const int argc, char **argv) {

  if (argc < 2) {
    std::cerr << "Usage:\n" << argv[0] << " opencl_filter\n";
    return -1;
  }

  const auto filter_str = std::string(argv[1]);
  std::cerr << filter_str << "\n";

  // Cuts off about 150 ms per filter:
  auto defines = std::string("");
  OpenCL_initialize(0, defines);

  if (filter_str == "add_const_filter") {

    auto image_ocl = GenerateImage<float, 3U>();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    OpenCL_AddConst_InPlace(image_ocl->GetBufferPointer(),
                            image_ocl->GetLargestPossibleRegion().GetSize(),
                            17.0f);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    auto image_cpu = GenerateImage<float, 3U>();
    const auto start_itk_time = std::chrono::steady_clock::now();
    using ImageType = itk::Image<float, 3U>;
    auto add_filter =
        itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(image_cpu);
    add_filter->SetConstant2(17.0f);
    add_filter->Update();
    image_cpu = add_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_itk_time - start_itk_time)
                     .count()
              << " ms\n";

    const auto result = CheckImage<itk::Image<float, 3U>>(image_ocl, image_cpu);
    if (result != 0) {
      return -2;
    }

  } else if (filter_str == "add_const_2d_filter") {

    auto image_ocl = GenerateImage<float, 2U>();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    OpenCL_AddConst_InPlace_2D(image_ocl->GetBufferPointer(),
                               image_ocl->GetLargestPossibleRegion().GetSize(),
                               17.0f);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    auto image_cpu = GenerateImage<float, 2U>();
    const auto start_itk_time = std::chrono::steady_clock::now();
    using ImageType = itk::Image<float, 2U>;
    auto add_filter =
        itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(image_cpu);
    add_filter->SetConstant2(17.0f);
    add_filter->Update();
    image_cpu = add_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_itk_time - start_itk_time)
                     .count()
              << " ms\n";

    const auto result = CheckImage<itk::Image<float, 2U>>(image_ocl, image_cpu);
    if (result != 0) {
      return -2;
    }
  } else if (filter_str == "add_mul_const_filter") {

    auto image_ocl = GenerateImage<float, 3U>();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    OpenCL_AddConst_MulConst_InPlace(
        image_ocl->GetBufferPointer(),
        image_ocl->GetLargestPossibleRegion().GetSize(), 17.0f, 2.0f);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    auto image_cpu = GenerateImage<float, 3U>();
    const auto start_itk_time = std::chrono::steady_clock::now();
    using ImageType = itk::Image<float, 3U>;
    auto add_filter =
        itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(image_cpu);
    add_filter->SetConstant2(17.0f);

    auto mul_filter =
        itk::MultiplyImageFilter<ImageType, ImageType, ImageType>::New();
    mul_filter->SetInput(add_filter->GetOutput());
    mul_filter->SetConstant(2.0f);
    mul_filter->Update();
    image_cpu = mul_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_itk_time - start_itk_time)
                     .count()
              << " ms\n";

    const auto result = CheckImage<itk::Image<float, 3U>>(image_ocl, image_cpu);
    if (result != 0) {
      return -2;
    }
  } else if (filter_str == "min_max_filter") {

    auto image_in = GenerateRandImage<float, 3U>();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto minmax_ocl =
        OpenCL_min_max_3D(image_in->GetBufferPointer(),
                          image_in->GetLargestPossibleRegion().GetSize());
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    const auto start_itk_time = std::chrono::steady_clock::now();
    using ImageType = itk::Image<float, 3U>;
    auto stat_filter = itk::StatisticsImageFilter<ImageType>::New();
    stat_filter->SetInput(image_in);
    stat_filter->Update();
    const auto max_itk = stat_filter->GetMaximum();
    const auto min_itk = stat_filter->GetMinimum();
    const auto end_itk_time = std::chrono::steady_clock::now();
    std::cerr << "ITKCPU:   "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_itk_time - start_itk_time)
                     .count()
              << " ms\n";

    if (fabs(min_itk - minmax_ocl.x) > 0.01f) {
      std::cerr << "Mininums were different: itk: " << min_itk
                << " ocl: " << minmax_ocl.x << "\n";
      return -2;
    }
    if (fabs((max_itk / minmax_ocl.y) - 1.f) > 0.01f) {
      std::cerr << "Maximums were different: itk: " << max_itk
                << " ocl: " << minmax_ocl.y << "\n";
      return -3;
    }

  } else if (filter_str == "min_max_2d_filter") {

    auto image_in = GenerateRandImage<float, 2U>();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto minmax_ocl =
        OpenCL_min_max_2D(image_in->GetBufferPointer(),
                          image_in->GetLargestPossibleRegion().GetSize());
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    const auto start_itk_time = std::chrono::steady_clock::now();
    using ImageType = itk::Image<float, 2U>;
    auto stat_filter = itk::StatisticsImageFilter<ImageType>::New();
    stat_filter->SetInput(image_in);
    stat_filter->Update();
    const auto max_itk = stat_filter->GetMaximum();
    const auto min_itk = stat_filter->GetMinimum();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_itk_time - start_itk_time)
                     .count()
              << " ms\n";

    if (fabs(min_itk - minmax_ocl.x) > 0.01f) {
      std::cerr << "Mininums were different: itk: " << min_itk
                << " ocl: " << minmax_ocl.x << "\n";
      return -2;
    }
    if (fabs((max_itk / minmax_ocl.y) - 1.f) > 0.01f) {
      std::cerr << "Maximums were different: itk: " << max_itk
                << " ocl: " << minmax_ocl.y << "\n";
      return -3;
    }

  } else if (filter_str == "divide_3Dby3D_loginv_filter") {

    using ImageType = itk::Image<unsigned short, 3U>;
    const auto image_in1 = GenerateRandImage<unsigned short, 3U>(69, 0, 4096);
    const auto image_in2 = GenerateRandImage<unsigned short, 3U>(69, 0, 1024);

    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto image_ocl =
        OpenCL_divide3Dby3D_loginv_OutOfPlace(image_in1, image_in2);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    const auto start_itk_time = std::chrono::steady_clock::now();
    using loginv_filter_type =
        itk::UnaryFunctorImageFilter<ImageType, FloatImageType, LogInvFunctor>;
    auto loginv_filter = loginv_filter_type::New();
    loginv_filter->SetInput(image_in1);
    auto loginv_filter_2 = loginv_filter_type::New();
    loginv_filter_2->SetInput(image_in2);

    auto div_filter = itk::DivideImageFilter<FloatImageType, FloatImageType,
                                             FloatImageType>::New();
    div_filter->SetInput1(loginv_filter->GetOutput());
    div_filter->SetInput2(loginv_filter_2->GetOutput());
    div_filter->Update();
    const auto image_itk = div_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_itk_time - start_itk_time)
                     .count()
              << " ms\n";

    const auto result = CheckImage<itk::Image<float, 3U>>(image_ocl, image_itk);
    if (result != 0) {
      return -2;
    }
  } else if (filter_str == "padding_filter") {
  } else if (filter_str == "subtract_2Dfrom3D_filter") {
  } else if (filter_str == "crop_by_struct_filter") {
    if (argc < 3) {
      std::cerr << "Crop filter requires a dicom directory as the 3rd input\n";
      return -1;
    }
    const auto dcmdir_str =
        QString(argv[2]).split(".", QString::SkipEmptyParts).at(0);
    const auto structures = load_rtstruct(
        dcmdir_str +
        "/RS.1.2.246.352.71.4.453824782.282736.20120706180259.dcm");

    const auto body_struct = structures->get_roi_ref_by_name("BODY");
    auto image = read_dicom_image(dcmdir_str);

    const auto start_ocl_time = std::chrono::steady_clock::now();
    OpenCL_crop_by_struct_InPlace(image, body_struct);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    UShortImageType::IndexType index;
    index.SetElement(0, 237);
    index.SetElement(1, 438);
    index.SetElement(2, 99);
    if (image->GetPixel(index) != 0) {
      std::cerr << "Test pixel was not 0, it was: " << image->GetPixel(index)
                << "\n";
      return -2;
    }

    index.SetElement(0, 189);
    index.SetElement(1, 366);
    index.SetElement(2, 99);
    if (image->GetPixel(index) != 1059) {
      std::cerr << "Test pixel was not 1059, it was: " << image->GetPixel(index)
                << "\n";
      return -3;
    }

  } else if (filter_str == "ItoLogI_subtract_median_filter") {
    using ImageType = itk::Image<float, 2U>;
    const auto median_radius = 3U;

    const auto proj_raw = GenerateRandImage<float, 2U>(69, 1.0f, 5.0f);

    const auto proj_scatter = GenerateRandImage<float, 2U>(6969, 5.0f, 5.3f);

    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto proj_corr = OpenCL_LogItoI_subtract_median_ItoLogI(
        proj_raw, proj_scatter, median_radius);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    const auto start_itk_time = std::chrono::steady_clock::now();
    auto convert_filter =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     LineInt2Intensity>::New();
    convert_filter->SetInput(proj_raw);
    auto convert_filter_2 =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     LineInt2Intensity>::New();
    convert_filter_2->SetInput(proj_scatter);

    auto subtract_filter =
        itk::SubtractImageFilter<ImageType, ImageType>::New();
    subtract_filter->SetInput1(convert_filter->GetOutput());
    subtract_filter->SetInput2(convert_filter_2->GetOutput());

    auto median_filter = itk::MedianImageFilter<ImageType, ImageType>::New();
    median_filter->SetInput(subtract_filter->GetOutput());
    median_filter->SetRadius(median_radius);

    auto convert_back_filter =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     Intensity2LineInt>::New();
    convert_back_filter->SetInput(median_filter->GetOutput());
    convert_back_filter->Update();
    const auto itk_proj_corr = convert_back_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_itk_time - start_itk_time)
                     .count()
              << " ms\n";

    const auto result = CheckImage<ImageType>(proj_corr, itk_proj_corr);
    auto border_size = 0;
    for (int i = 0; i < median_radius; ++i) {
      border_size += 2 * (1024 - 2 * i + 768 - 2 * i - 2);
    }
    if (result > border_size) {
      return -2;
    } else {
      std::cerr << "ITK handles the border differently."
                   "So we only fail if more than the border is wrong\n";
    }
  } else if (filter_str == "ItoLogI_subtract_median_gaussian_filter") {
    using ImageType = itk::Image<float, 2U>;

    const auto proj_prim = GenerateRandImage<float, 2U>(69, 1.5f, 5.0f);
    const auto proj_rand = GenerateRandImage<float, 2U>(6969, 0.0f, -0.5f);

    auto add_filter =
        itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(proj_prim);
    add_filter->SetInput2(proj_rand);
    add_filter->Update();
    const auto proj_raw = add_filter->GetOutput();

    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto proj_scatter = OpenCL_LogItoI_subtract_median_gaussian_ItoLogI(
        proj_raw, proj_prim, 3, 1.5);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    const auto start_itk_time = std::chrono::steady_clock::now();
    auto convert_filter =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     LineInt2Intensity>::New();
    convert_filter->SetInput(proj_raw);

    auto convert_filter_2 =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     LineInt2Intensity>::New();
    convert_filter_2->SetInput(proj_prim);

    auto subtract_filter =
        itk::SubtractImageFilter<ImageType, ImageType>::New();
    subtract_filter->SetInput1(convert_filter->GetOutput());
    subtract_filter->SetInput2(convert_filter_2->GetOutput());

    auto median_filter = itk::MedianImageFilter<ImageType, ImageType>::New();
    median_filter->SetInput(subtract_filter->GetOutput());
    median_filter->SetRadius(3);

    auto gaussian_filter =
        itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>::New();
    gaussian_filter->SetInput(median_filter->GetOutput());

    itk::SmoothingRecursiveGaussianImageFilter<
        ImageType, ImageType>::SigmaArrayType gauss_sigma;
    gauss_sigma[0] = 1.5;
    gauss_sigma[1] = gauss_sigma[0] * 0.75;
    gaussian_filter->SetSigmaArray(gauss_sigma);

    auto convert_back_filter =
        itk::UnaryFunctorImageFilter<ImageType, ImageType,
                                     Intensity2LineInt>::New();
    convert_back_filter->SetInput(gaussian_filter->GetOutput());
    convert_back_filter->Update();
    auto itk_proj_scatter = convert_back_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_itk_time - start_itk_time)
                     .count()
              << " ms\n";

    saveImageAsMHA<ImageType>(itk_proj_scatter, "itk_submedgauss.mha");
    saveImageAsMHA<ImageType>(proj_scatter, "ocl_submedgauss.mha");
    const auto result = CheckImage<ImageType>(proj_scatter, itk_proj_scatter);
    if (result != 0) {
      return -2;
    }

  } else {
    std::cerr << "This filter does not exists: " << filter_str << "\n";
    return -1;
  }

  return 0;
}
