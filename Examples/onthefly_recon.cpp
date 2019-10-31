/*=========================================================================
 *
 *  Copyright RTK Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include <chrono>
#include <mutex>
#include <thread>

#include <rtkConfiguration.h>
#include <rtkGgoFunctions.h>

#include "onthefly_recon_ggo.h"

#include <rtkDisplacedDetectorImageFilter.h>
#include <rtkFDKConeBeamReconstructionFilter.h>
#include <rtkParkerShortScanImageFilter.h>
#include <rtkProjectionsReader.h>
#include <rtkThreeDCircularProjectionGeometryXMLFile.h>
#ifdef USE_CUDA
#include <rtkCudaFDKConeBeamReconstructionFilter.h>
#endif
#ifdef RTK_USE_OPENCL
#include "rtkOpenCLFDKConeBeamReconstructionFilter.h"
#endif

#include <itkImageFileWriter.h>

// Pass projection name, projection parameters, last
struct ThreadInfoStruct {
  std::mutex mutex;
  args_info_onthefly_recon *args_info;
  bool stop;
  unsigned int nproj;
  double radius;
  double sid;
  double sdd;
  double gantryAngle;
  double projOffsetX;
  double projOffsetY;
  double outOfPlaneAngle;
  double inPlaneAngle;
  double sourceOffsetX;
  double sourceOffsetY;
  double collimationUInf;
  double collimationUSup;
  double collimationVInf;
  double collimationVSup;
  double minimumOffsetX; // Used for Wang weighting
  double maximumOffsetX;
  std::string fileName;
};

void computeOffsetsFromGeometry(
    rtk::ThreeDCircularProjectionGeometry::Pointer geometry, double *minOffset,
    double *maxOffset) {
  double min = std::numeric_limits<double>::max();
  double max = std::numeric_limits<double>::min();

  for (unsigned int i = 0; i < geometry->GetProjectionOffsetsX().size(); i++) {
    min = std::min(min, geometry->GetProjectionOffsetsX()[i]);
    max = std::max(max, geometry->GetProjectionOffsetsX()[i]);
  }
  *minOffset = min;
  *maxOffset = max;
}

// This thread reads in a geometry file and a sequence of projection file names
// and communicates them one by one to the other thread via a ThreadinfoStruct.
void AcquisitionCallback(ThreadInfoStruct &threadInfo) {
  double minOffset, maxOffset;

  threadInfo.mutex.lock();

  // Get file names
  std::vector<std::string> names =
      rtk::GetProjectionsFileNamesFromGgo(*(threadInfo.args_info));

  // Geometry
  if (threadInfo.args_info->verbose_flag)
    std::cout << "Reading geometry information from "
              << threadInfo.args_info->geometry_arg << "..." << std::endl;
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(threadInfo.args_info->geometry_arg);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation());

  // Computes the minimum and maximum offsets from Geometry
  computeOffsetsFromGeometry(geometryReader->GetOutputObject(), &minOffset,
                             &maxOffset);
  std::cout << " main :" << minOffset << " " << maxOffset << std::endl;

  // Set the radius of the cylindrical detector (0 means flat)
  threadInfo.radius =
      geometryReader->GetOutputObject()->GetRadiusCylindricalDetector();

  threadInfo.mutex.unlock();

  // Mock an inline acquisition
  unsigned int nproj = geometryReader->GetOutputObject()->GetMatrices().size();
  rtk::ThreeDCircularProjectionGeometry *geometry =
      geometryReader->GetOutputObject();
  for (unsigned int i = 0; i < nproj; i++) {
    threadInfo.mutex.lock();
    threadInfo.sdd = geometry->GetSourceToDetectorDistances()[i];
    threadInfo.sid = geometry->GetSourceToIsocenterDistances()[i];
    threadInfo.gantryAngle = geometry->GetGantryAngles()[i];
    threadInfo.sourceOffsetX = geometry->GetSourceOffsetsX()[i];
    threadInfo.sourceOffsetY = geometry->GetSourceOffsetsY()[i];
    threadInfo.inPlaneAngle = geometry->GetInPlaneAngles()[i];
    threadInfo.outOfPlaneAngle = geometry->GetOutOfPlaneAngles()[i];
    threadInfo.collimationUInf = geometry->GetCollimationUInf()[i];
    threadInfo.collimationUSup = geometry->GetCollimationUSup()[i];
    threadInfo.collimationVInf = geometry->GetCollimationVInf()[i];
    threadInfo.collimationVSup = geometry->GetCollimationVSup()[i];
    threadInfo.minimumOffsetX = minOffset;
    threadInfo.maximumOffsetX = maxOffset;
    threadInfo.fileName = names[std::min(i, (unsigned int)names.size() - 1)];
    threadInfo.nproj = i + 1;
    threadInfo.stop = (i == nproj - 1);
    if (threadInfo.args_info->verbose_flag)
      std::cout
          << std::endl
          << "AcquisitionCallback has simulated the acquisition of projection #"
          << i << std::endl;
    threadInfo.mutex.unlock();

    std::this_thread::sleep_for(std::chrono::milliseconds(200));
  }

  return;
}

enum enDeviceType { CUDA_DEVT, CPU_DEVT, OPENCL_DEVT };

template <enDeviceType Tdev, typename pixel_type, unsigned int dim>
struct fdk_type {
  using ImageType = itk::Image<pixel_type, dim>;
  using type = rtk::FDKConeBeamReconstructionFilter<ImageType>;
};

#ifdef USE_CUDA
template <typename pixel_type, unsigned int dim>
struct fdk_type<CUDA_DEVT, pixel_type, dim> {
  using ImageType = itk::CudaImage<pixel_type, dim>;
  using type = rtk::FDKConeBeamReconstructionFilter<ImageType>;
};
#endif

#ifdef RTK_USE_OPENCL
template <> struct fdk_type<OPENCL_DEVT, float, 3U> {
  using ImageType = itk::Image<float, 3U>;
  using type = rtk::OpenCLFDKConeBeamReconstructionFilter;
};
#endif

// This thread receives information of each projection (one-by-one) and process
// directly the projections for which it has enough information. This thread
// currently assumes that the projections are sequentially sent with
// increasing gantry angles. Specific management with a queue must be
// implemented if the projections are not exactly sequential. Short
// scans has not been implemented yet because this filter currently
// require the full geometry of the acquisition. Management with a mock
// geometry file would be possible but it is still to be implemented.
template <enDeviceType device>
void InlineThreadCallback(ThreadInfoStruct &threadInfo) {

  threadInfo.mutex.lock();
  typedef float OutputPixelType;
  const unsigned int Dimension = 3;
  using CPUOutputImageType =
      fdk_type<CPU_DEVT, OutputPixelType, Dimension>::ImageType;

  using OutputImageType =
      typename fdk_type<device, OutputPixelType, Dimension>::ImageType;

  rtk::ThreeDCircularProjectionGeometry::Pointer geometry =
      rtk::ThreeDCircularProjectionGeometry::New();
  std::vector<std::string> fileNames;

  // Projections reader
  auto reader = rtk::ProjectionsReader<OutputImageType>::New();

  // Create reconstructed image
  auto constantImageSource = rtk::ConstantImageSource<OutputImageType>::New();
  rtk::SetConstantImageSourceFromGgo<rtk::ConstantImageSource<OutputImageType>,
                                     args_info_onthefly_recon>(
      constantImageSource, *(threadInfo.args_info));

  // Extract filter to process one projection at a time
  using ExtractFilterType =
      itk::ExtractImageFilter<OutputImageType, OutputImageType>;
  auto extract = ExtractFilterType::New();
  extract->SetInput(reader->GetOutput());
  typename ExtractFilterType::InputImageRegionType subsetRegion;

  // Displaced detector weighting
  auto ddf = rtk::DisplacedDetectorImageFilter<OutputImageType>::New();
  ddf->SetInput(extract->GetOutput());
  ddf->SetGeometry(geometry);
  ddf->SetDisable(threadInfo.args_info->nodisplaced_flag);

  // Short scan image filter
  //  typedef rtk::ParkerShortScanImageFilter< OutputImageType > PSSFType;
  //  PSSFType::Pointer pssf = PSSFType::New();
  //  pssf->SetInput( ddf->GetOutput() );
  //  pssf->SetGeometry( geometryReader->GetOutputObject() );
  //  pssf->InPlaceOff();

  // FDK reconstruction filtering
  using FDKType = typename fdk_type<device, OutputPixelType, Dimension>::type;
  typename FDKType::Pointer fdk;
  fdk->SetInput(0, constantImageSource->GetOutput());
  fdk->SetInput(1, ddf->GetOutput());
  fdk->SetGeometry(geometry);
  fdk->GetRampFilter()->SetTruncationCorrection(threadInfo.args_info->pad_arg);
  fdk->GetRampFilter()->SetHannCutFrequency(threadInfo.args_info->hann_arg);

  // Writer
  typedef itk::ImageFileWriter<CPUOutputImageType> WriterType;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(threadInfo.args_info->output_arg);

  // Set the cylindrical detector's radius
  geometry->SetRadiusCylindricalDetector(threadInfo.radius);

  threadInfo.mutex.unlock();

  // Inline loop
  std::cout << "Reconstruction thread has entered in the processing loop"
            << std::endl;
  for (;;) {
    threadInfo.mutex.lock();

    if (geometry->GetMatrices().size() < threadInfo.nproj) {
      if (threadInfo.args_info->verbose_flag)
        std::cerr << "InlineThreadCallback has received projection #"
                  << threadInfo.nproj - 1 << std::endl;

      if (threadInfo.fileName != "" &&
          (fileNames.size() == 0 || fileNames.back() != threadInfo.fileName))
        fileNames.push_back(threadInfo.fileName);

      geometry->AddProjectionInRadians(
          threadInfo.sid, threadInfo.sdd, threadInfo.gantryAngle,
          threadInfo.projOffsetX, threadInfo.projOffsetY,
          threadInfo.outOfPlaneAngle, threadInfo.inPlaneAngle,
          threadInfo.sourceOffsetX, threadInfo.sourceOffsetY);
      geometry->SetCollimationOfLastProjection(
          threadInfo.collimationUInf, threadInfo.collimationUSup,
          threadInfo.collimationVInf, threadInfo.collimationVSup);

      std::cout << "Geometry size : " << geometry->GetMatrices().size()
                << std::endl;

      if (geometry->GetMatrices().size() != threadInfo.nproj) {
        std::cerr << "Missed one projection in InlineThreadCallback"
                  << std::endl;
        exit(EXIT_FAILURE);
      }
      if (geometry->GetMatrices().size() < 3) {
        threadInfo.mutex.unlock();
        continue;
      }

      typename OutputImageType::RegionType region =
          reader->GetOutput()->GetLargestPossibleRegion();
      std::cout << "Reader size : " << region.GetSize()[0] << " "
                << region.GetSize()[1] << " " << region.GetSize()[2]
                << std::endl;
      std::cout << "Reader index : " << region.GetIndex()[0] << " "
                << region.GetIndex()[1] << " " << region.GetIndex()[2]
                << std::endl;

      reader->SetFileNames(fileNames);
      TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->UpdateOutputInformation())
      subsetRegion = reader->GetOutput()->GetLargestPossibleRegion();
      subsetRegion.SetIndex(Dimension - 1, geometry->GetMatrices().size() - 2);
      subsetRegion.SetSize(Dimension - 1, 1);
      extract->SetExtractionRegion(subsetRegion);

      std::cout << "Region size : " << subsetRegion.GetSize()[0] << " "
                << subsetRegion.GetSize()[1] << " " << subsetRegion.GetSize()[2]
                << std::endl;
      std::cout << "Region index : " << subsetRegion.GetIndex()[0] << " "
                << subsetRegion.GetIndex()[1] << " "
                << subsetRegion.GetIndex()[2] << std::endl;

      typename ExtractFilterType::InputImageRegionType extractRegion =
          extract->GetOutput()->GetLargestPossibleRegion();

      std::cout << "Extract region size : " << extractRegion.GetSize()[0] << " "
                << extractRegion.GetSize()[1] << " "
                << extractRegion.GetSize()[2] << std::endl;
      std::cout << "Extract region index : " << extractRegion.GetIndex()[0]
                << " " << extractRegion.GetIndex()[1] << " "
                << extractRegion.GetIndex()[2] << std::endl;

      ddf->SetOffsets(threadInfo.minimumOffsetX, threadInfo.maximumOffsetX);

      TRY_AND_EXIT_ON_ITK_EXCEPTION(fdk->Update())
      typename OutputImageType::Pointer pimg = fdk->GetOutput();
      pimg->DisconnectPipeline();
      fdk->SetInput(pimg);
      TRY_AND_EXIT_ON_ITK_EXCEPTION(fdk->GetOutput()->UpdateOutputInformation())
      TRY_AND_EXIT_ON_ITK_EXCEPTION(
          fdk->GetOutput()->PropagateRequestedRegion());

      if (threadInfo.args_info->verbose_flag)
        std::cout << "Projection #" << subsetRegion.GetIndex(Dimension - 1)
                  << " has been processed in reconstruction." << std::endl;

      if (threadInfo.stop) {
        // Process first projection
        subsetRegion.SetIndex(Dimension - 1, 0);
        extract->SetExtractionRegion(subsetRegion);
        TRY_AND_EXIT_ON_ITK_EXCEPTION(fdk->Update())
        typename OutputImageType::Pointer pimg = fdk->GetOutput();
        pimg->DisconnectPipeline();
        fdk->SetInput(pimg);
        TRY_AND_EXIT_ON_ITK_EXCEPTION(
            fdk->GetOutput()->UpdateOutputInformation())
        TRY_AND_EXIT_ON_ITK_EXCEPTION(
            fdk->GetOutput()->PropagateRequestedRegion());
        if (threadInfo.args_info->verbose_flag)
          std::cout << "Projection #" << subsetRegion.GetIndex(Dimension - 1)
                    << " has been processed in reconstruction." << std::endl;

        // Process last projection
        subsetRegion.SetIndex(Dimension - 1,
                              geometry->GetMatrices().size() - 1);
        extract->SetExtractionRegion(subsetRegion);
        TRY_AND_EXIT_ON_ITK_EXCEPTION(fdk->Update())
        writer->SetInput(fdk->GetOutput());

        if (threadInfo.args_info->verbose_flag)
          std::cout << "Projection #" << subsetRegion.GetIndex(Dimension - 1)
                    << " has been processed in reconstruction." << std::endl;

        // Write to disk and exit
        TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update())
        exit(EXIT_SUCCESS);
      }
    }

    threadInfo.mutex.unlock();
  }

  return;
}

int main(int argc, char *argv[]) {
  GGO(onthefly_recon, args_info);

  // Launch threads, one for acquisition, one for reconstruction with inline
  // processing
  ThreadInfoStruct threadInfo;
  threadInfo.args_info = &args_info;
  threadInfo.nproj = 0;
  threadInfo.minimumOffsetX = 0.0;
  threadInfo.maximumOffsetX = 0.0;

  std::thread thread_acq(AcquisitionCallback, std::ref(threadInfo));

  std::thread thread_fdk;

  if (!strcmp(args_info.hardware_arg, "cuda")) {
#ifndef USE_CUDA
    std::cerr << "Program was not compiled with cuda option\n";
    return EXIT_FAILURE;
#else
    thread_fdk =
        std::thread(InlineThreadCallback<CUDA_DEVT>, std::ref(threadInfo));
#endif
  } else if (!strcmp(args_info.hardware_arg, "opencl")) {
#ifndef RTK_USE_OPENCL
    std::cerr << "Program was not compiled with OpenCL option\n";
    return EXIT_FAILURE;
#else
    thread_fdk =
        std::thread(InlineThreadCallback<OPENCL_DEVT>, std::ref(threadInfo));
#endif
  } else {
    thread_fdk =
        std::thread(InlineThreadCallback<CPU_DEVT>, std::ref(threadInfo));
  }

  thread_acq.join();
  thread_fdk.join();
  return EXIT_SUCCESS;
}

