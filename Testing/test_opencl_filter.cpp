// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include "OpenCL/ImageFilters.h"

#include <chrono>
#include <iostream>

#include "itkAddImageFilter.h"

template<typename T, size_t DIM>
auto GenerateImage(){
  using ImageType = itk::Image<T, DIM>;
  auto image = ImageType::New();
  typename ImageType::IndexType origin;
  origin[0] = -256.0;
  origin[1] = -256.0;
  origin[2] = -100.0;
  typename ImageType::SizeType size;
  size[0] = 512;
  size[1] = 512;
  size[2] = 200;
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

template<typename ImageType>
auto CheckImage(typename ImageType::Pointer image_test, typename ImageType::Pointer image_ref){
  itk::ImageRegionIterator<ImageType> ImageIter_test(image_test, image_test->GetLargestPossibleRegion());
  itk::ImageRegionIterator<ImageType> ImageIter_ref(image_ref, image_ref->GetLargestPossibleRegion());

  auto n_err = 0;
  while (!ImageIter_ref.IsAtEnd()){
    if (ImageIter_ref.Get() - ImageIter_test.Get() > 0.01f){
      ++n_err;
      std::cerr << "Ref: " << ImageIter_ref.Get() << " Test: " << ImageIter_test.Get() << "\n";
      break;
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
  if (filter_str == "add_const_filter"){
    std::cerr << filter_str << "\n";

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


  return 0;
}
