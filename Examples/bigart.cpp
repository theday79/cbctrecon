// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

// For testing CbctRecon

#include "cbctrecon.h"

#include <execution>
#include <algorithm>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "itkEuler3DTransform.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkRescaleIntensityImageFilter.h"

#include "WEPL.h"
#include "cbctrecon_io.h"
#include "cbctrecon_test.hpp"

namespace fs = std::filesystem;
using namespace std::literals;

constexpr float ce_distance(const FloatVector& point_a, const FloatVector& point_b) {
  return crl::ce_sqrt<float>(crl::ce_pow(point_a.x - point_b.x, 2) +
                      crl::ce_pow(point_a.y - point_b.y, 2) +
                      crl::ce_pow(point_a.z - point_b.z, 2));
}

template<typename T>
struct hausdorff_result {
  T h_min = std::numeric_limits<T>::max();
  T h_max = std::numeric_limits<T>::max();
  T h_percent = std::numeric_limits<T>::max();
};

template<typename T>
T min_distance(const FloatVector from_point, const Rtss_roi_modern &to_roi) {
  auto to_roi_distances = std::vector<T>();
  for (const auto &to_contour : to_roi.pslist) {
    auto &to_coords = to_contour.coordinates;
    auto to_contour_distances = std::vector<T>(to_coords.size());
    std::transform(std::execution::par_unseq, to_coords.begin(),
                   to_coords.end(), to_contour_distances.begin(),
                   [from_point](const auto to_point) {
                     return ce_distance(from_point, to_point);
                   });
    to_roi_distances.push_back(*std::min_element(to_contour_distances.begin(),
                                                to_contour_distances.end()));
  }
  return *std::min_element(to_roi_distances.begin(), to_roi_distances.end());
}

template<typename T, unsigned char percent>
hausdorff_result<T> calculate_hausdorff(const Rtss_roi_modern &from_roi,
                                        const Rtss_roi_modern &to_roi) {
  static_assert(
      percent <= 100,
      "Percent should be less than 100 as a higher value doesn't make sense.");

  const auto n_total_points =
      std::transform_reduce(std::execution::par_unseq, from_roi.pslist.begin(), from_roi.pslist.end(),
                  static_cast<size_t>(0), std::plus<size_t>(), [](const Rtss_contour_modern &contour) {
                    return contour.coordinates.size();
                  });

  auto distances = std::vector<T>(n_total_points);
  auto dist_iterator = distances.begin();
  for (const auto &from_contour : from_roi.pslist) {
    auto &from_coords = from_contour.coordinates;
    dist_iterator = std::transform(
        std::execution::par_unseq, from_coords.begin(), from_coords.end(),
        dist_iterator, [&to_roi](const auto from_point) {
          return min_distance<T>(from_point, to_roi);
        });
  }

  std::sort(std::execution::par_unseq, distances.begin(), distances.end());

  hausdorff_result<T> hausdorff_output;
  hausdorff_output.h_min = distances.front();
  hausdorff_output.h_max = distances.back();

  const auto percent_index = static_cast<size_t>(
      crl::ce_round((percent / 100.0) * static_cast<double>(distances.size())));

  hausdorff_output.h_percent = distances.at(percent_index);

  return hausdorff_output;
}


constexpr auto deg2rad(const double deg) { return deg / 180.0 * itk::Math::pi; }

UShortImageType::Pointer get_image_from_dicom(
    const fs::path &dcm_dir, const UShortImageType::SpacingType &new_spacing,
    const UShortImageType::SizeType &new_size,
    const UShortImageType::PointType &new_origin,
    const UShortImageType::DirectionType &new_direction,
    const DoubleVector &translation_vec, const DoubleVector &rotation_vec) {
  // auto cbctrecon_test = std::make_unique<CbctReconTest>();

  auto dcm_path = fs::absolute(dcm_dir);
  if (!fs::exists(dcm_path)) {
    std::cerr << "Directory didn't exist: " << dcm_path << "\n";
  }
  if (fs::is_empty(dcm_path)) {
    std::cerr << "Directory was empty: " << dcm_path << "\n";
  }

  // cbctrecon_test->m_cbctrecon->m_strPathDirDefault = dcm_path;
  // cbctrecon_test->test_LoadDICOMdir();
  // auto ct_img = cbctrecon_test->m_cbctrecon->m_spManualRigidCT;
  auto mha_reader = itk::ImageFileReader<ShortImageType>::New();
  mha_reader->SetFileName((dcm_path / "recalc_in.mha").string());
  mha_reader->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  auto imageCalculatorFilter = ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(mha_reader->GetOutput());
  imageCalculatorFilter->Compute();

  const auto minVal0 = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal0 = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Original Min and Max Values are	" << minVal0 << "	"
            << maxVal0 << std::endl;

  auto bNKI = false;
  if (minVal0 > -600) // impossible for normal Short image. IN NKI, always -512.
                      // don't know why
  {
    bNKI = true;
  }

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  auto thresholdFilter = ThresholdImageFilterType::New();

  if (!bNKI) {
    thresholdFilter->SetInput(mha_reader->GetOutput());
    thresholdFilter->ThresholdOutside(-1024, 3072); //--> 0 ~ 4095
    thresholdFilter->SetOutsideValue(-1024);
    thresholdFilter->Update();
  } else {
    thresholdFilter->SetInput(mha_reader->GetOutput());
    thresholdFilter->ThresholdOutside(0, 4095); //--> 0 ~ 4095
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();
  }

  imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
  imageCalculatorFilter->Compute();

  const auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  const auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Current Min and Max Values are	" << minVal << "	"
            << maxVal << std::endl;

  USHORT_PixelType outputMinVal, outputMaxVal;
  if (!bNKI) {
    outputMinVal = static_cast<USHORT_PixelType>(minVal) + 1024U;
    outputMaxVal = static_cast<USHORT_PixelType>(maxVal) + 1024U;
  } else {
    outputMinVal = static_cast<USHORT_PixelType>(minVal);
    outputMaxVal = static_cast<USHORT_PixelType>(maxVal);
  }

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>;
  auto spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();
  const auto ct_img = spRescaleFilter->GetOutput();

  if (!ct_img) {
    std::cerr << "Manual Rigid CT was NULL -> Dicom dir was not read!\n";
    return nullptr;
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

  auto translation = itk::Euler3DTransform<double>::InputVectorType();
  translation.SetElement(
      0, -translation_vec.x); // 5.5185);  // -X in eclipse -> X itk
  translation.SetElement(
      1, -translation_vec.y); // -2.6872); // -Y in eclipse -> Y itk
  translation.SetElement(
      2, -translation_vec.z); // -6.4281); // -Z in eclipse -> Z itk
  transform->SetTranslation(translation);

  resampler->SetInput(ct_img);
  const auto interpolator =
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

// Currently just writes the points to a vector of vectors of QStrings
auto get_signed_difference(const Rtss_roi_modern *wepl_voi,
                           const Rtss_roi_modern *orig_voi/*,
                           const std::array<double, 3> basis*/) {

  auto output = std::vector<std::vector<QString>>(orig_voi->pslist.size());

  std::transform(std::begin(wepl_voi->pslist), std::end(wepl_voi->pslist),
                 std::begin(orig_voi->pslist), std::begin(output),
                 [/*&basis*/](const Rtss_contour_modern &wepl_contour,
                              const Rtss_contour_modern &orig_contour) {
                   auto out_contour =
                       std::vector<QString>(orig_contour.coordinates.size());
                   std::transform(std::begin(wepl_contour.coordinates),
                                  std::end(wepl_contour.coordinates),
                                  /*std::begin(orig_contour.coordinates),*/
                                  std::begin(out_contour),
                                  [/*&basis*/](const FloatVector &wepl_coord /*,
                                               const FloatVector &orig_coord*/) {
                                    return QString("%1,%2,%3\n")
                                        .arg(wepl_coord.x)
                                        .arg(wepl_coord.y)
                                        .arg(wepl_coord.z);
                                    /*return basis.at(0) * (orig_coord.x -
                                       wepl_coord.x) + basis.at(1) *
                                       (orig_coord.y - wepl_coord.y) +
                                           basis.at(2) * (orig_coord.z -
                                       wepl_coord.z);*/
                                  });
                   return out_contour;
                 });

  return output;
}

auto roi_to_distal_only_roi(const Rtss_roi_modern *roi, const double gantry,
                                const double couch) {
  const auto direction = crl::wepl::get_basis_from_angles(gantry, couch);
  auto out_roi = Rtss_roi_modern(*roi);

  for (auto &contour : out_roi.pslist) {
    contour.coordinates = crl::wepl::distal_points_only(contour, direction);
  }
  
  return out_roi;
}

bool save_orig_with_only_distal(const std::string &voi_name,
                                const Rtss_modern *ss,
                                const std::string &suffix,
                                const fs::path &orig_file, const double gantry,
                                const double couch) {
  auto ss_distal = std::make_unique<Rtss_modern>(*ss);

  const auto direction = crl::wepl::get_basis_from_angles(gantry, couch);
  for (auto &voi_distal_it : ss_distal->slist) {
    if (voi_distal_it.name.find(voi_name) != std::string::npos) {
      for (auto &contour : voi_distal_it.pslist) {
        contour.coordinates = crl::wepl::distal_points_only(contour, direction);
      }
      voi_distal_it.name = "WEPLorig" + voi_distal_it.name;
    }
  }

  const auto out_distal_orig_dcm = "RS.orig_distal_structure_"s;
  const fs::path out_distal_orig_dcm_file(out_distal_orig_dcm + voi_name +
                                          "_CT"s + suffix);
  std::cerr << "Writing distal struct to dicom...\n";
  return crl::AlterData_RTStructureSetStorage(orig_file, ss_distal.get(),
                                              out_distal_orig_dcm_file);
}

auto get_ss_from_dir(const fs::path &dir) {
  const auto filenamelist = crl::get_dcm_image_files(dir);

  for (auto &&filename : fs::directory_iterator(dir)) {
    if (!(filename.is_regular_file() ||
          filename.is_symlink())) { // I guess symlinks should be allowed?
      continue;
    }
    if (filename.path().string().find("-hash-stamp") != std::string::npos) {
      continue; // Just so the test data is less annoying.
    }
    const auto fullfilename = filename.path();
    const auto modality = crl::get_dcm_modality(fullfilename);
    switch (modality) {
    case DCM_MODALITY::RTSTRUCT:
      return crl::load_rtstruct(fullfilename);
      break;
    default:
      break;
    }
  }
  auto out = std::make_unique<Rtss_modern>();
  out = nullptr;
  return out;
}

template<typename T>
auto stringify_hausdorff(const hausdorff_result<T> &hd) {
  return crl::stringify(hd.h_max) + ";" + crl::stringify(hd.h_min) + ";" +
         crl::stringify(hd.h_percent);
}

int main(const int argc, char **argv) {

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
  const auto dcm_dir_str = std::string(argv[1]);

  auto dcm_dir = fs::path(dcm_dir_str);
  const auto dcm_path = fs::absolute(dcm_dir);
  if (!fs::exists(dcm_dir)) {
    std::cerr << "Directory didn't exist: " << dcm_path << "\n";
    return -2;
  }
  if (fs::is_empty(dcm_dir)) {
    std::cerr << "Directory was empty: " << dcm_path << "\n";
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
  auto translation_str_list = QString(argv[6]).split(",");
  const auto translation = DoubleVector{
      std::stod(translation_str_list.at(0).toStdString(), &sz),
      std::stod(translation_str_list.at(1).toStdString(), &sz),
      std::stod(translation_str_list.at(2).toStdString(), &sz)};
  auto rotation_str_list = QString(argv[7]).split(",");
  const auto rotation =
      DoubleVector{std::stod(rotation_str_list.at(0).toStdString(), &sz),
                   std::stod(rotation_str_list.at(1).toStdString(), &sz),
                   std::stod(rotation_str_list.at(2).toStdString(), &sz)};

  auto recalc_dcm_dir = fs::path(argv[2]);
  const auto recalc_dcm_path = fs::absolute(recalc_dcm_dir);

  auto recalc_img = get_image_from_dicom(
      recalc_dcm_path, ct_img->GetSpacing(),
      ct_img->GetLargestPossibleRegion().GetSize(), ct_img->GetOrigin(),
      ct_img->GetDirection(), translation, rotation);

  if (!recalc_img) {
    return -5;
  }

  cbctrecon_test->m_cbctrecon->m_spRawReconImg = recalc_img;
  cbctrecon_test->m_dlgRegistration->UpdateListOfComboBox(0);
  cbctrecon_test->m_dlgRegistration->UpdateListOfComboBox(1);
  QString raw_str("RAW_CBCT");
  QString man_str("MANUAL_RIGID_CT");
  cbctrecon_test->m_dlgRegistration->LoadImgFromComboBox(0, raw_str);
  cbctrecon_test->m_dlgRegistration->LoadImgFromComboBox(1, man_str);

  cbctrecon_test->m_dlgRegistration->SLT_MovingImageSelected("MANUAL_RIGID_CT");
  const auto body_index =
      cbctrecon_test->m_dlgRegistration->ui_comboBox_VOItoCropBy->findText(
          "BODY", Qt::MatchStartsWith | Qt::MatchCaseSensitive);

  crl::saveImageAsMHA<UShortImageType>(
      cbctrecon_test->m_cbctrecon->m_spRawReconImg,
      (recalc_dcm_path / "recalc.mha").string());
  crl::saveImageAsMHA<UShortImageType>(
      cbctrecon_test->m_cbctrecon->m_spManualRigidCT,
      (dcm_path / "orig.mha").string());
  ct_img = cbctrecon_test->m_cbctrecon->m_spManualRigidCT;
  recalc_img = cbctrecon_test->m_cbctrecon->m_spRawReconImg;
  std::cout << "Images saved" << std::endl;

  /* Structure test: */
  auto ss = cbctrecon_test->m_cbctrecon->m_structures->get_ss(ctType::RIGID_CT);

  auto voi = std::string();

  for (auto &structure : ss->slist) {
    if (structure.name == argv[3]) {
      voi = structure.name;
    }
  }
  if (voi.empty()) {
    for (auto &structure : ss->slist) {
      if (structure.name.find(argv[3]) != std::string::npos) {
        voi = structure.name;
      }
    }
  }
  if (voi.empty()) {
    std::cerr << "Couldn't find: \"" << argv[3] << "\"\n";
    return -5;
  }
  std::cerr << voi << "\n";

  const auto gantry_angle = std::stod(argv[4], &sz);
  const auto couch_angle = std::stod(argv[5], &sz);
  const auto descriptive_suffix = "_G"s + argv[4] + "_C" + argv[5] + ".dcm";

  if (!save_orig_with_only_distal(argv[3], ss, descriptive_suffix,
                                  cbctrecon_test->m_cbctrecon->m_strPathRS,
                                  gantry_angle, couch_angle)) {
    return -6;
  }

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

  const auto orig_voi = ss->get_roi_by_name(voi);
  // const auto basis = get_basis_from_angles(gantry_angle, couch_angle);

  constexpr auto distal_only = true;
  const auto wepl_voi = crl::wepl::CalculateWEPLtoVOI<distal_only>(
      orig_voi.get(), gantry_angle, couch_angle, ct_img, recalc_img);

  /* Generate a vector of vectors with distances */
  //const auto output =
  //    get_signed_difference(wepl_voi, orig_voi.get()); //, basis);

  for (auto &cur_voi : ss->slist) {
    if (cur_voi.name.find(argv[3]) != std::string::npos) {
      cur_voi = *wepl_voi;
    }
  }
  const auto out_dcm = "RS.wepl_structure_"s;
  // http://www.plastimatch.org/plastimatch.html#plastimatch-dice
  const fs::path out_dcm_file(out_dcm + argv[3] + "_"s +
                              recalc_dcm_dir.filename().string() +
                              descriptive_suffix);
  std::cerr << "Writing WEPL struct to dicom...\n";
  if (!crl::AlterData_RTStructureSetStorage(
          cbctrecon_test->m_cbctrecon->m_strPathRS, ss, out_dcm_file)) {
    std::cerr << "\a"
              << "Could not write dcm\n";
    return -7;
  }

  /* Write distances to file */
  auto better_name = orig_voi->name;
  std::replace(std::begin(better_name), std::end(better_name), ' ', '_');
  std::replace(std::begin(better_name), std::end(better_name), '/', '-');
  const auto output_filename = "Hausdorff_" + better_name + "_" +
                               dcm_dir.filename().string() + "_to_" +
                               recalc_dcm_dir.filename().string() + ".txt";

  const auto rct_ss = get_ss_from_dir(recalc_dcm_dir);
  const auto rct_voi = rct_ss->get_roi_by_name(voi);
  const auto distal_rct_voi = roi_to_distal_only_roi( rct_voi.get(), gantry_angle, couch_angle);
  const auto distal_orig_voi =
      roi_to_distal_only_roi(orig_voi.get(), gantry_angle, couch_angle);
  const auto &cref_wepl_voi = *wepl_voi;

  constexpr auto percent = 95;
  const auto hd_orig_to_wepl =
      calculate_hausdorff<float, percent>(distal_orig_voi, cref_wepl_voi);
  const auto hd_wepl_to_orig =
      calculate_hausdorff<float, percent>(cref_wepl_voi, distal_orig_voi);

  const auto hd_orig_to_rct =
      calculate_hausdorff<float, percent>(distal_orig_voi, distal_rct_voi);
  const auto hd_rct_to_orig =
      calculate_hausdorff<float, percent>(distal_rct_voi, distal_orig_voi);

  const auto hd_rct_to_wepl =
      calculate_hausdorff<float, percent>(distal_rct_voi, cref_wepl_voi);
  const auto hd_wepl_to_rct =
      calculate_hausdorff<float, percent>(cref_wepl_voi, distal_rct_voi);

  const auto str_percent = crl::stringify(percent);

  std::cerr << "Writing hausdorff to txt...\n";
  std::ofstream f_stream;
  f_stream.open(output_filename);
  if (f_stream.is_open()) {
    // f_stream << std::fixed << std::setprecision(5);
    f_stream << "from;to;HD;min;max;percent\n";

    f_stream << "orig;wepl;" + stringify_hausdorff(hd_orig_to_wepl) + "\n";
    f_stream << "orig;wepl;" + stringify_hausdorff(hd_wepl_to_orig) + "\n";

    f_stream << "orig;rct;" + stringify_hausdorff(hd_orig_to_rct) + "\n";
    f_stream << "rct;orig;" + stringify_hausdorff(hd_rct_to_orig) + "\n";

    f_stream << "rct;wepl;" + stringify_hausdorff(hd_rct_to_wepl) + "\n";
    f_stream << "wepl;rct;" + stringify_hausdorff(hd_wepl_to_rct) + "\n";

    f_stream.close();
  } else {
    std::cerr << "Could not open file: " << output_filename << "\n";
    return -8;
  }

  std::cerr << "Done!\n";

  return 0;
}
