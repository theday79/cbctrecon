// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

// For testing CbctRecon

#include "cbctrecon.h"

#include <iostream>
#include <memory>
#include <numeric>
#include <string>

#include <QDir>

#include "itkEuler3DTransform.h"
#include "vnl_vector_fixed.h"

#include "cbctrecon_io.h"
#include "cbctrecon_test.hpp"

constexpr auto deg2rad(double deg) { return deg / 180.0 * itk::Math::pi; }

UShortImageType::Pointer get_image_from_dicom(
    const QString &dir, const UShortImageType::SpacingType &new_spacing,
    const UShortImageType::SizeType &new_size,
    const UShortImageType::PointType &new_origin,
    const UShortImageType::DirectionType &new_direction,
    const DoubleVector &translation_vec, const DoubleVector &rotation_vec) {
  auto cbctrecon_test = std::make_unique<CbctReconTest>();

  auto dcm_dir = QDir(dir);
  auto dcm_path = dcm_dir.absolutePath();
  if (!dcm_dir.exists()) {
    std::cerr << "Directory didn't exist: " << dcm_path.toStdString() << "\n";
  }
  if (dcm_dir.isEmpty(QDir::AllEntries | QDir::NoDotAndDotDot)) {
    std::cerr << "Directory was empty: " << dcm_path.toStdString() << "\n";
  }

  cbctrecon_test->m_cbctrecon->m_strPathDirDefault = dcm_path;
  cbctrecon_test->test_LoadDICOMdir();
  auto &ct_img = cbctrecon_test->m_cbctrecon->m_spManualRigidCT;

  if (ct_img.IsNull()) {
    std::cerr << "Manual Rigid CT was NULL -> Dicom dir was not read!\n";
  }
  /* Some debug info: */
  UShortImageType::PointType first_point;
  UShortImageType::IndexType index{0, 0, 0};
  ct_img->TransformIndexToPhysicalPoint(index, first_point);
  std::cerr << "First pixel point: " << first_point[0] << ", " << first_point[1]
            << ", " << first_point[2] << "\n";

  auto size = ct_img->GetLargestPossibleRegion().GetSize();
  index.SetElement(0, size[0]);
  index.SetElement(1, size[1]);
  index.SetElement(2, size[2]);
  UShortImageType::PointType last_point;
  ct_img->TransformIndexToPhysicalPoint(index, last_point);
  std::cerr << "Last pixel point: " << last_point[0] << ", " << last_point[1]
            << ", " << last_point[2] << "\n";

  auto resampler =
      itk::ResampleImageFilter<UShortImageType, UShortImageType>::New();

  auto transform = itk::Euler3DTransform<double>::New();
  transform->SetRotation(deg2rad(-rotation_vec.x), deg2rad(-rotation_vec.y),
                         deg2rad(-rotation_vec.z));
  // transform->SetRotation(deg2rad(-0.8843), deg2rad(1.6274),
  // deg2rad(-0.3450));

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
      itk::LinearInterpolateImageFunction<UShortImageType, double>::New();
  resampler->SetInterpolator(interpolator);
  resampler->SetSize(new_size);
  resampler->SetOutputSpacing(new_spacing);
  resampler->SetOutputOrigin(new_origin);
  resampler->SetOutputDirection(new_direction);
  resampler->SetTransform(transform);
  resampler->Update();

  return resampler->GetOutput();
}

int main(const int argc, char *argv[]) {

  if (argc < 6) {
    std::cerr << "Usage:\n"
              << argv[0]
              << " /plan/dicom/dir /recalc/dicom/dir voi_name gantry couch "
                 "translation rotation\n"
              << "Where translation[mm] and rotation[deg] is comma-seperated "
                 "triples in the eclipse coordinate system\n";
    return -1;
  }

  std::cerr << "Running BiGART script!\n"
            << " " << argv[0] << " " << argv[1] << " " << argv[2] << " "
            << argv[3] << " " << argv[4] << " " << argv[5] << " " << argv[6]
            << " " << argv[7] << "\n";

  /* Read dicom and structures */
  auto cbctrecon_test = std::make_unique<CbctReconTest>();
  auto dcm_dir_str = QString(argv[1]);

  auto dcm_dir = QDir(dcm_dir_str);
  auto dcm_path = dcm_dir.absolutePath();
  if (!dcm_dir.exists()) {
    std::cerr << "Directory didn't exist: " << dcm_path.toStdString() << "\n";
    return -2;
  }
  if (dcm_dir.isEmpty(QDir::AllEntries | QDir::NoDotAndDotDot)) {
    std::cerr << "Directory was empty: " << dcm_path.toStdString() << "\n";
    return -3;
  }

  cbctrecon_test->m_cbctrecon->m_strPathDirDefault = dcm_path;
  cbctrecon_test->test_LoadDICOMdir();

  auto &ct_img = cbctrecon_test->m_cbctrecon->m_spManualRigidCT;
  if (ct_img.IsNull()) {
    std::cerr << "Manual Rigid CT was NULL -> Dicom dir was not read!\n";
    return -4;
  }
  std::string::size_type sz;
  auto translation_str = QString(argv[6]);
  auto translation = DoubleVector{
      std::stod(translation_str.split(",").at(0).toStdString(), &sz),
      std::stod(translation_str.split(",").at(1).toStdString(), &sz),
      std::stod(translation_str.split(",").at(2).toStdString(), &sz)};
  auto rotation_str = QString(argv[7]);
  auto rotation =
      DoubleVector{std::stod(rotation_str.split(",").at(0).toStdString(), &sz),
                   std::stod(rotation_str.split(",").at(1).toStdString(), &sz),
                   std::stod(rotation_str.split(",").at(2).toStdString(), &sz)};

  auto recalc_dcm_dir = QString(argv[2]);

  auto recalc_img = get_image_from_dicom(
      recalc_dcm_dir, ct_img->GetSpacing(),
      ct_img->GetLargestPossibleRegion().GetSize(), ct_img->GetOrigin(),
      ct_img->GetDirection(), translation, rotation);

  saveImageAsMHA<UShortImageType>(recalc_img,
                                  recalc_dcm_dir.toStdString() + "/recalc.mha");
  saveImageAsMHA<UShortImageType>(ct_img, dcm_path.toStdString() + "/orig.mha");

  /* Structure test: */
  auto ss = cbctrecon_test->m_cbctrecon->m_structures->get_ss(PLAN_CT);

  auto voi = std::string();

  for (auto &structure : ss->slist) {
    if (structure.name.compare(argv[3]) == 0) { // "GTV 74/37 LB") == 0) {
      voi = structure.name;
    }
  }
  if (voi.empty()) {
    std::cerr << "Couldn't find: \"" << argv[3] << "\"\n";
    return -5;
  }
  std::cerr << voi << "\n";

  /* calculate WEPL coordinates */
  UShortImageType::PointType first_point;
  UShortImageType::IndexType index{0, 0, 0};
  ct_img->TransformIndexToPhysicalPoint(index, first_point);
  std::cerr << "First pixel point: " << first_point[0] << ", " << first_point[1]
            << ", " << first_point[2] << "\n";

  auto size = ct_img->GetLargestPossibleRegion().GetSize();
  index.SetElement(0, size[0]);
  index.SetElement(1, size[1]);
  index.SetElement(2, size[2]);
  UShortImageType::PointType last_point;
  ct_img->TransformIndexToPhysicalPoint(index, last_point);
  std::cerr << "Last pixel point: " << last_point[0] << ", " << last_point[1]
            << ", " << last_point[2] << "\n";

  auto &orig_voi = ss->get_roi_ref_by_name(voi);
  const auto gantry_angle = static_cast<int>(round(std::stod(argv[4], &sz)));
  const auto couch_angle = static_cast<int>(round(std::stod(argv[5], &sz)));
  cbctrecon_test->m_cbctregistration->CalculateWEPLtoVOI(
      voi, gantry_angle, couch_angle, ct_img, recalc_img);

  /* Generate a vector of vectors with distances */
  const auto &wepl_voi = cbctrecon_test->m_cbctregistration->WEPL_voi;
  auto output = std::vector<std::vector<float>>(orig_voi.pslist.size());

  std::transform(
      std::begin(wepl_voi->pslist), std::end(wepl_voi->pslist),
      std::begin(orig_voi.pslist), std::begin(output),
      [](const Rtss_contour_modern &wepl_contour,
         const Rtss_contour_modern &orig_contour) {
        auto out_contour = std::vector<float>(orig_contour.coordinates.size());
        std::transform(
            std::begin(wepl_contour.coordinates),
            std::end(wepl_contour.coordinates),
            std::begin(orig_contour.coordinates), std::begin(out_contour),
            [](const FloatVector &wepl_coord, const FloatVector &orig_coord) {
              return sqrt(pow(wepl_coord.x - orig_coord.x, 2) +
                          pow(wepl_coord.y - orig_coord.y, 2) +
                          pow(wepl_coord.z - orig_coord.z, 2));
            });
        return out_contour;
      });

  /* Write distances to file */
  auto better_name = orig_voi.name;
  std::replace(std::begin(better_name), std::end(better_name), ' ', '_');
  std::replace(std::begin(better_name), std::end(better_name), '/', '-');
  auto output_filename =
      "Dist_" + better_name + "_CT" +
      dcm_path.at(dcm_path.size() - 1).toLatin1() + "_to_CT" +
      recalc_dcm_dir.at(recalc_dcm_dir.size() - 1).toLatin1() + "_at_G" +
      std::to_string(gantry_angle) + "C" + std::to_string(couch_angle) + ".txt";

  std::ofstream f_stream;
  f_stream.open(output_filename);
  if (f_stream.is_open()) {
    f_stream << std::fixed << std::setprecision(5);

    std::for_each(std::begin(output), std::end(output),
                  [&f_stream](const std::vector<float> out_vec) {
                    std::for_each(std::begin(out_vec), std::end(out_vec),
                                  [&f_stream](const float val) {
                                    f_stream << val << "\n";
                                  });
                  });

    f_stream.close();
  } else {
    std::cerr << "Could not open file: " << output_filename << "\n";
  }

  return 0;
}
