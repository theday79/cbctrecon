// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include "OpenCL/ImageFilters.h"

#include <chrono>
#include <iostream>

#include "itkAddImageFilter.h"
#include "itkMultiplyImageFilter.h"
#include "itkStatisticsImageFilter.h"
#include "itkDivideImageFilter.h"

template<typename T, size_t DIM>
auto GenerateImage(){
  using ImageType = itk::Image<T, DIM>;
  auto image = ImageType::New();
  typename ImageType::IndexType origin;
  origin[0] = -256.0;
  origin[1] = -256.0;
  if (DIM == 3){
    origin[2] = -100.0;
  }
  typename ImageType::SizeType size;
  size[0] = 512;
  size[1] = 512;
  if (DIM == 3){
    size[2] = 200;
  }
  typename ImageType::RegionType region;
  region.SetIndex(origin);
  region.SetSize(size);

  image->SetRegions(region);
  image->Allocate();

  itk::ImageRegionIterator<ImageType> ImageIter(image, region);
  while(!ImageIter.IsAtEnd()){
    ImageIter.Set(static_cast<T>(1));
    ++ImageIter;
  }

  return image;
}

template<typename T, size_t DIM>
auto GenerateRandImage(){
  using ImageType = itk::Image<T, DIM>;
  auto image = ImageType::New();
  typename ImageType::IndexType origin;
  origin[0] = -256.0;
  origin[1] = -256.0;
  if (DIM == 3){
    origin[2] = -100.0;
  }
  typename ImageType::SizeType size;
  size[0] = 512;
  size[1] = 512;
  if (DIM == 3){
    size[2] = 200;
  }
  typename ImageType::RegionType region;
  region.SetIndex(origin);
  region.SetSize(size);

  image->SetRegions(region);
  image->Allocate();

  std::srand(69); // RNG seed

  itk::ImageRegionIterator<ImageType> ImageIter(image, region);
  while(!ImageIter.IsAtEnd()){
    // rand is between 0 and RAND_MAX, which is implementation dependend
    ImageIter.Set(static_cast<T>(std::rand()));
    ++ImageIter;
  }

  return image;
}

template<typename ImageType>
auto CheckImage(typename ImageType::Pointer image_test, typename ImageType::Pointer image_ref){
  itk::ImageRegionIterator<ImageType> ImageIter_test(image_test, image_test->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> ImageIter_ref(image_ref, image_ref->GetLargestPossibleRegion());

  auto n_err = 0;
  while (!ImageIter_ref.IsAtEnd()){
    if (fabs(ImageIter_ref.Get() - ImageIter_test.Get()) > 0.01f){
      ++n_err;
    }
    ++ImageIter_ref;
    ++ImageIter_test;
  }
  if (n_err != 0){
    std::cerr << "Error larger than 0.01 " << n_err << " times, between OpenCL and ITK implementation!\n";
    return -n_err;
  }
  return 0;
}


int main(int argc, char** argv){

  if (argc < 2) {
    std::cerr << "Usage:\n"
              << argv[0] << " opencl_filter\n";
    return -1;
  }

  const auto filter_str = std::string(argv[1]);
  std::cerr << filter_str << "\n";

  // Cuts off about 150 ms per filter:
  OpenCL_initialize(512*512*200*sizeof(cl_float));

  if (filter_str == "add_const_filter"){

    auto image_ocl = GenerateImage<float, 3U>();
    auto *image_buffer = image_ocl->GetBufferPointer();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    OpenCL_AddConst_InPlace(image_ocl->GetBufferPointer(),
            image_ocl->GetLargestPossibleRegion().GetSize(), 17.0f);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_ocl_time - start_ocl_time).count()
        << " ms\n";

    auto image_cpu = GenerateImage<float, 3U>();
    const auto start_itk_time = std::chrono::steady_clock::now();
    using ImageType = itk::Image<float, 3U>;
    auto add_filter = itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(image_cpu);
    add_filter->SetConstant2(17.0f);
    add_filter->Update();
    image_cpu = add_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_itk_time - start_itk_time).count()
        << " ms\n";

    return CheckImage<itk::Image<float, 3U>>(
            image_ocl, image_cpu);

  }
  else if (filter_str == "add_const_2d_filter"){

    auto image_ocl = GenerateImage<float, 2U>();
    auto *image_buffer = image_ocl->GetBufferPointer();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    OpenCL_AddConst_InPlace_2D(image_ocl->GetBufferPointer(),
            image_ocl->GetLargestPossibleRegion().GetSize(), 17.0f);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_ocl_time - start_ocl_time).count()
        << " ms\n";

    auto image_cpu = GenerateImage<float, 2U>();
    const auto start_itk_time = std::chrono::steady_clock::now();
    using ImageType = itk::Image<float, 2U>;
    auto add_filter = itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(image_cpu);
    add_filter->SetConstant2(17.0f);
    add_filter->Update();
    image_cpu = add_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_itk_time - start_itk_time).count()
        << " ms\n";

    return CheckImage<itk::Image<float, 2U>>(
            image_ocl, image_cpu);
  }
  else if (filter_str == "add_mul_const_filter"){

    auto image_ocl = GenerateImage<float, 3U>();
    auto *image_buffer = image_ocl->GetBufferPointer();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    OpenCL_AddConst_MulConst_InPlace(image_ocl->GetBufferPointer(),
            image_ocl->GetLargestPossibleRegion().GetSize(), 17.0f, 2.0f);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_ocl_time - start_ocl_time).count()
        << " ms\n";

    auto image_cpu = GenerateImage<float, 3U>();
    const auto start_itk_time = std::chrono::steady_clock::now();
    using ImageType = itk::Image<float, 3U>;
    auto add_filter = itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(image_cpu);
    add_filter->SetConstant2(17.0f);

    auto mul_filter = itk::MultiplyImageFilter<ImageType, ImageType, ImageType>::New();
    mul_filter->SetInput(add_filter->GetOutput());
    mul_filter->SetConstant(2.0f);
    mul_filter->Update();
    image_cpu = mul_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_itk_time - start_itk_time).count()
        << " ms\n";

    return CheckImage<itk::Image<float, 3U>>(
            image_ocl, image_cpu);
  }
  else if (filter_str == "min_max_filter"){

    auto image_in = GenerateRandImage<float, 3U>();
    auto *image_buffer = image_in->GetBufferPointer();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto minmax_ocl = OpenCL_min_max_3D(image_in->GetBufferPointer(),
            image_in->GetLargestPossibleRegion().GetSize());
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_ocl_time - start_ocl_time).count()
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
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_itk_time - start_itk_time).count()
        << " ms\n";

    const auto start_itk_iter_time = std::chrono::steady_clock::now();
    auto iter = itk::ImageRegionIterator<ImageType>(image_in, image_in->GetLargestPossibleRegion());
    auto it_min = std::numeric_limits<float>::max();
    auto it_max = std::numeric_limits<float>::min();
    while (!iter.IsAtEnd()){
      if (iter.Get() < it_min){
        it_min = iter.Get();
      }
      if (iter.Get() > it_max){
        it_max = iter.Get();
      }
      ++iter;
    }
    const auto end_itk_iter_time = std::chrono::steady_clock::now();

    std::cerr << "IterCPU:  "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_itk_iter_time - start_itk_iter_time).count()
        << " ms\n";

    if (fabs(min_itk - it_min) > 0.01f){
        std::cerr << "Mininums were different: itk: " << min_itk << " iter: " << it_min << "\n";
        return -2;
    }
    if (fabs((max_itk / it_max) - 1.f) > 0.01f){
        std::cerr << "Maximums were different: itk: " << max_itk << " iter: " << it_max << "\n";
        return -3;
    }

    if (fabs(min_itk - minmax_ocl.x) > 0.01f){
        std::cerr << "Mininums were different: itk: " << min_itk << " ocl: " << minmax_ocl.x << "\n";
        return -2;
    }
    if (fabs((max_itk / minmax_ocl.y) - 1.f) > 0.01f){
        std::cerr << "Maximums were different: itk: " << max_itk << " ocl: " << minmax_ocl.y << "\n";
        return -3;
    }

    return 0;
  }
  else if (filter_str == "min_max_2d_filter"){

    auto image_in = GenerateRandImage<float, 2U>();
    auto *image_buffer = image_in->GetBufferPointer();
    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto minmax_ocl = OpenCL_min_max_2D(image_in->GetBufferPointer(),
            image_in->GetLargestPossibleRegion().GetSize());
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_ocl_time - start_ocl_time).count()
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
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_itk_time - start_itk_time).count()
        << " ms\n";

    if (fabs(min_itk - minmax_ocl.x) > 0.01f){
        std::cerr << "Mininums were different: itk: " << min_itk << " ocl: " << minmax_ocl.x << "\n";
        return -2;
    }
    if (fabs((max_itk / minmax_ocl.y) - 1.f) > 0.01f){
        std::cerr << "Maximums were different: itk: " << max_itk << " ocl: " << minmax_ocl.y << "\n";
        return -3;
    }

    return 0;

  }
  else if (filter_str == "divide_3Dby3D_filter"){

    using ImageType = itk::Image<unsigned short, 3U>;
    auto image_in1 = GenerateRandImage<unsigned short, 3U>();
    auto image_in2 = GenerateRandImage<unsigned short, 3U>();
    auto add_filter = itk::AddImageFilter<ImageType, ImageType, ImageType>::New();
    add_filter->SetInput1(image_in1);
    add_filter->SetConstant2(17);
    add_filter->Update();
    image_in1 = add_filter->GetOutput();

    const auto start_ocl_time = std::chrono::steady_clock::now();
    const auto image_ocl = OpenCL_divide3Dby3D_OutOfPlace(image_in1, image_in2);
    const auto end_ocl_time = std::chrono::steady_clock::now();

    std::cerr << "OpenCL: "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_ocl_time - start_ocl_time).count()
        << " ms\n";

    const auto start_itk_time = std::chrono::steady_clock::now();
    auto div_filter = itk::DivideImageFilter<ImageType, ImageType, itk::Image<float, 3U>>::New();
    div_filter->SetInput1(image_in1);
    div_filter->SetInput2(image_in2);
    div_filter->Update();
    const auto image_itk = div_filter->GetOutput();
    const auto end_itk_time = std::chrono::steady_clock::now();

    std::cerr << "ITKCPU:   "
        << std::chrono::duration_cast<std::chrono::milliseconds>(end_itk_time - start_itk_time).count()
        << " ms\n";

    return CheckImage<itk::Image<float, 3U>>(image_ocl, image_itk);
  }
  else if (filter_str == "padding_filter"){}
  else if (filter_str == "subtract_2Dfrom3D_filter"){}
  else { std::cerr << "This filter does not exists: " << filter_str << "\n"; }


  return 0;
}
