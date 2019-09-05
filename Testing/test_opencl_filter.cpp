// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include "OpenCL/ImageFilters.h"

#include <chrono>
#include <iostream>

#include <QDir>

#include "itkGDCMImageIO.h"
#include "itkImageSeriesReader.h"
#include "itkRescaleIntensityImageFilter.h"

#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkStatisticsImageFilter.h"

#include "StructureSet.h"
#include "cbctrecon_io.h"

#include "OpenCL/err_code.h"

#define REQ_OPENCL_VER static_cast<cl_uint>(CBCTRECON_OPENCL_VERSION)

template <typename T, size_t DIM> auto GenerateImage() {
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
    ImageIter.Set(static_cast<T>(1));
    ++ImageIter;
  }

  return image;
}

template <typename T, size_t DIM> auto GenerateRandImage() {
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

  std::srand(69); // RNG seed

  itk::ImageRegionIterator<ImageType> ImageIter(image, region);
  while (!ImageIter.IsAtEnd()) {
    // rand is between 0 and RAND_MAX, which is implementation dependend
    ImageIter.Set(static_cast<T>(std::rand()));
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
    if (fabs(ImageIter_ref.Get() - ImageIter_test.Get()) > 0.01f) {
      ++n_err;
    }
    ++ImageIter_ref;
    ++ImageIter_test;
  }
  if (n_err != 0) {
    std::cerr << "Error larger than 0.01 " << n_err
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
      return result;
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
      return result;
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
      return result;
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

  } else if (filter_str == "divide_3Dby3D_filter") {

    using ImageType = itk::Image<unsigned short, 3U>;
    auto image_in1 = GenerateRandImage<unsigned short, 3U>();
    const auto image_in2 = GenerateRandImage<unsigned short, 3U>();
    auto add_filter =
        itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(image_in1);
    add_filter->SetConstant2(17);
    add_filter->Update();
    image_in1 = add_filter->GetOutput();

    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto image_ocl = OpenCL_divide3Dby3D_OutOfPlace(image_in1, image_in2);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(
                     end_ocl_time - start_ocl_time)
                     .count()
              << " ms\n";

    const auto start_itk_time = std::chrono::steady_clock::now();
    auto div_filter = itk::DivideImageFilter<ImageType, ImageType,
                                             itk::Image<float, 3U>>::New();
    div_filter->SetInput1(image_in1);
    div_filter->SetInput2(image_in2);
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
      return result;
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

  } else {
    std::cerr << "This filter does not exists: " << filter_str << "\n";
    return -1;
  }

  return 0;
}
