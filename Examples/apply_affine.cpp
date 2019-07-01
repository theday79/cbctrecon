// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

// For testing CbctRecon

#include "cbctrecon.h"

#include <iostream>

#include <QDir>

#include "itkEuler3DTransform.h"


constexpr auto deg2rad(const double deg) { return deg / 180.0 * itk::Math::pi; }

int main(const int argc, char *argv[]) {

  if (argc < 6) {
    std::cerr << "Usage:\n"
              << argv[0]
              << " input_file output_filename translation rotation imglike\n"
              << "Where translation[mm] and rotation[deg] is comma-seperated "
                 "triples in the eclipse coordinate system,\n"
                 "imglike is an itk image of the wished output geometry\n";
    return -1;
  }

  std::cerr << "Running BiGART script!\n"
            << " " << argv[0] << " " << argv[1] << " " << argv[2] << " "
            << argv[3] << " " << argv[4] << " " << argv[5] << "\n";

  using ImageType = itk::Image<uint32_t, 3>;

  /* Read dicom and structures */
  const auto input_file = QString(argv[1]);
  auto reader = itk::ImageFileReader<ImageType>::New();
  reader->SetFileName(input_file.toStdString());
  reader->Update();

  const auto imglike_file = QString(argv[5]);
  auto reader_like = itk::ImageFileReader<ImageType>::New();
  reader_like->SetFileName(imglike_file.toStdString());
  reader_like->UpdateOutputInformation();

  auto ct_img = reader->GetOutput();

  if (!ct_img) {
    std::cerr << "Manual Rigid CT was NULL -> Dicom dir was not read!\n";
  }

  auto resampler =
      itk::ResampleImageFilter<ImageType, ImageType>::New();

  std::string::size_type sz;
  auto translation_str = QString(argv[3]);
  const auto translation_vec = DoubleVector{
      std::stod(translation_str.split(",").at(0).toStdString(), &sz),
      std::stod(translation_str.split(",").at(1).toStdString(), &sz),
      std::stod(translation_str.split(",").at(2).toStdString(), &sz)};
  auto rotation_str = QString(argv[4]);
  const auto rotation_vec =
      DoubleVector{std::stod(rotation_str.split(",").at(0).toStdString(), &sz),
                   std::stod(rotation_str.split(",").at(1).toStdString(), &sz),
                   std::stod(rotation_str.split(",").at(2).toStdString(), &sz)};

  auto transform = itk::Euler3DTransform<double>::New();
  transform->SetRotation(deg2rad(-rotation_vec.x), deg2rad(-rotation_vec.y),
                         deg2rad(-rotation_vec.z));

  auto translation = itk::Euler3DTransform<double>::InputVectorType();
  translation.SetElement(
      0, -translation_vec.x); // 5.5185);  // -X in eclipse -> X itk
  translation.SetElement(
      1, -translation_vec.y); // -2.6872); // -Y in eclipse -> Y itk
  translation.SetElement(
      2, -translation_vec.z); // -6.4281); // -Z in eclipse -> Z itk
  transform->SetTranslation(translation);

  resampler->SetInput(ct_img);
  auto interpolator =
      itk::LinearInterpolateImageFunction<ImageType, double>::New();
  resampler->SetInterpolator(interpolator);
  resampler->SetSize(reader_like->GetOutput()->GetLargestPossibleRegion().GetSize());
  resampler->SetOutputSpacing(reader_like->GetOutput()->GetSpacing());
  resampler->SetOutputOrigin(reader_like->GetOutput()->GetOrigin());
  resampler->SetOutputDirection(reader_like->GetOutput()->GetDirection());
  resampler->SetTransform(transform);
  resampler->Update();

  const auto output_file = QString(argv[2]);
  auto writer = itk::ImageFileWriter<ImageType>::New();
  writer->SetInput(resampler->GetOutput());
  writer->SetFileName(output_file.toStdString());
  writer->Update();

  return 0;
}

