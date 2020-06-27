// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

// For testing CbctRecon

#include "cbctrecon.h"

#include <algorithm>
#include <execution>
#include <filesystem>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include <QString>

#include "itkMinimumMaximumImageCalculator.h"
#include "itkQuaternionRigidTransform.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkTransformFileWriter.h"

#include "WEPL.h"
#include "cbctrecon_io.h"
#include "cbctrecon_test.hpp"

namespace fs = std::filesystem;
using namespace std::literals;

constexpr auto deg2rad(const double deg) { return deg / 180.0 * itk::Math::pi; }

auto get_transform(const DoubleVector &rot, const DoubleVector &translation) {
  using TransformType = itk::QuaternionRigidTransform<double>;
  auto transform = TransformType::New();

  auto quart = TransformType::VnlQuaternionType(deg2rad(rot.x), deg2rad(rot.y),
                                                deg2rad(rot.z));

  transform->SetRotation(quart);

  auto translation_vec = TransformType::InputVectorType();
  translation_vec.SetElement(0, translation.x);
  translation_vec.SetElement(1, translation.y);
  translation_vec.SetElement(2, translation.z);
  transform->SetTranslation(translation_vec);

  std::cerr << "Transform Matrix: " << transform->GetMatrix() << "\n";

  auto out_transform = TransformType::New();
  if (!transform->GetInverse(out_transform)) {
    std::cerr << "Could not invert transformation, cringe!?\n";
  }

  return out_transform;
}

std::pair<UShortImageType::Pointer, std::unique_ptr<Rtss_modern>>
get_image_from_dicom(const fs::path &dcm_dir,
                     const UShortImageType::SpacingType &new_spacing,
                     const UShortImageType::SizeType &new_size,
                     const UShortImageType::PointType &new_origin,
                     const UShortImageType::DirectionType &new_direction,
                     const DoubleVector &translation_vec,
                     const DoubleVector &rotation_vec) {

  auto dcm_path = fs::absolute(dcm_dir);
  if (!fs::exists(dcm_path)) {
    std::cerr << "Directory didn't exist: " << dcm_path << "\n";
  }
  if (fs::is_empty(dcm_path)) {
    std::cerr << "Directory was empty: " << dcm_path << "\n";
  }

  auto cbctrecon_test = std::make_unique<CbctReconTest>();
  cbctrecon_test->m_cbctrecon->m_strPathDirDefault = dcm_path;
  cbctrecon_test->test_LoadDICOMdir();
  auto ct_img = cbctrecon_test->m_cbctrecon->m_spManualRigidCT;
  // auto mha_reader = itk::ImageFileReader<ShortImageType>::New();
  // mha_reader->SetFileName((dcm_path / "recalc_in.mha").string());
  // mha_reader->Update();

  if (!ct_img) {
    std::cerr << "Manual Rigid CT was NULL -> Dicom dir was not read!\n";
    return std::make_pair(nullptr, nullptr);
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

  const auto transform = get_transform(rotation_vec, translation_vec);

  // Write transform to file so it can be handled by plastimatch and translated
  // into an appropriate lambda
  auto transform_writer = itk::TransformFileWriterTemplate<double>::New();
  transform_writer->SetInput(transform);
  const auto transform_file = dcm_dir / "QuarternionRigid3D_transform.txt"s;
  transform_writer->SetFileName(transform_file.string());
  transform_writer->Update();

  cbctrecon_test->m_cbctrecon->m_structures->ApplyTransformTo<ctType::PLAN_CT>(
      transform_file);

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

  return std::make_pair(resampler->GetOutput(),
                        std::move(cbctrecon_test->m_cbctrecon->m_structures
                                      ->get_ss<ctType::RIGID_CT>()));
}

// Currently just writes the points to a vector of vectors of QStrings
auto get_signed_difference(const Rtss_roi_modern &wepl_voi,
                           const Rtss_roi_modern &orig_voi/*,
                           const std::array<double, 3> basis*/) {

  auto output = std::vector<std::vector<std::string>>(orig_voi.pslist.size());

  std::transform(std::begin(wepl_voi.pslist), std::end(wepl_voi.pslist),
                 std::begin(orig_voi.pslist), std::begin(output),
                 [/*&basis*/](const Rtss_contour_modern &wepl_contour,
                              const Rtss_contour_modern &orig_contour) {
                   auto out_contour = std::vector<std::string>(
                       orig_contour.coordinates.size());
                   std::transform(std::begin(wepl_contour.coordinates),
                                  std::end(wepl_contour.coordinates),
                                  /*std::begin(orig_contour.coordinates),*/
                                  std::begin(out_contour),
                                  [/*&basis*/](const FloatVector &wepl_coord /*,
                                               const FloatVector &orig_coord*/) {
                                    return crl::make_sep_str<','>(wepl_coord.x, wepl_coord.y, wepl_coord.z) + "\n";
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

bool save_orig_with_only_distal(const std::string_view voi_name,
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
  const fs::path out_distal_orig_dcm_file(
      out_distal_orig_dcm + std::string(voi_name) + "_CT"s + suffix);
  std::cerr << "Writing distal struct to dicom...\n";
  return crl::AlterData_RTStructureSetStorage(orig_file, ss_distal.get(),
                                              out_distal_orig_dcm_file);
}

auto get_rct_voi(Rtss_modern *rct_ss, const std::string &rct_name,
                 const std::string &orig_voi_name) {

  const auto remove_space = orig_voi_name.back() == ' ' ? 1 : 0;
  // length of "pCT" is 3, so:
  const auto voi_postfix = std::string_view(
      orig_voi_name.data() + 3, orig_voi_name.size() - (3 + remove_space));
  const auto rct_voi_name = rct_name + std::string(voi_postfix);
  std::cerr << "\"" << rct_voi_name << "\"\n";

  return rct_ss->get_roi_ref_by_name(rct_voi_name);
}

// max(%);max(min);max(max)
template <typename T, class BinOp>
auto stringify_hausdorff(const crl::hausdorff_result<T> &hd_lhs,
                         const crl::hausdorff_result<T> &hd_rhs, BinOp op) {
  return crl::stringify(op(hd_lhs.h_percent, hd_rhs.h_percent)) + ";" +
         crl::stringify(op(hd_lhs.h_min, hd_rhs.h_min)) + ";" +
         crl::stringify(op(hd_lhs.h_max, hd_rhs.h_max));
}

template <typename T, class BinOp>
bool write_hd_to_file(const fs::path &out_fn, const int dose_level,
                      const double gantry_angle, const std::string &rct,
                      const crl::hausdorff_result<T> &hd_lhs,
                      const crl::hausdorff_result<T> &hd_rhs, BinOp op) {
  std::cerr << "Writing hausdorff to txt...\n";
  std::ofstream f_stream;
  f_stream.open(out_fn, std::ios::out | std::ios::app);
  if (f_stream.is_open()) {
    if (f_stream.tellp() == 0) {
      f_stream << "Dose;gantry;rct;hausdorff_percent;min;max\n";
    }
    f_stream << dose_level << ";" << gantry_angle << ";" << rct << ";"
             << stringify_hausdorff(hd_lhs, hd_rhs, op) << "\n";
    f_stream.close();
  } else {
    std::cerr << "Could not open file: " << out_fn << "\n";
    return false;
  }
  return true;
}

auto sv_to_dvec(std::string_view sv) {
  const auto sv_v = crl::split_string(sv, ",");
  const auto vec = crl::from_sv_v<double>(sv_v);
  return DoubleVector{vec.at(0), vec.at(1), vec.at(2)};
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

  const auto dcm_dir = fs::path(argv[1]);
  const auto recalc_dcm_dir = fs::path(argv[2]);
  const auto voi_str = std::string_view(argv[3]);
  const auto gantry_str = std::string_view(argv[4]);
  const auto couch_str = std::string_view(argv[5]);
  const auto translation = sv_to_dvec(argv[6]);
  const auto rotation = sv_to_dvec(argv[7]);

  /* Read dicom and structures */
  auto cbctrecon_test = std::make_unique<CbctReconTest>();

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

  const auto recalc_dcm_path = fs::absolute(recalc_dcm_dir);

  auto [recalc_img, rct_ss] = get_image_from_dicom(
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
    if (structure.name == voi_str) {
      voi = structure.name;
    }
  }
  if (voi.empty()) {
    for (auto &structure : ss->slist) {
      if (structure.name.find(voi_str) != std::string::npos) {
        voi = structure.name;
      }
    }
  }
  if (voi.empty()) {
    std::cerr << "Couldn't find: \"" << voi_str << "\"\n";
    return -5;
  }
  std::cerr << voi << "\n";

  auto gantry_angle = 0;
  crl::from_sv(gantry_str, gantry_angle);
  auto couch_angle = 0;
  crl::from_sv(couch_str, couch_angle);

  const auto descriptive_suffix = "_G"s + crl::stringify(gantry_angle) + "_C" +
                                  crl::stringify(couch_angle) + ".dcm";

  if (!save_orig_with_only_distal(voi_str, ss, descriptive_suffix,
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

  const auto orig_voi = ss->get_roi_ref_by_name(voi);
  // const auto basis = get_basis_from_angles(gantry_angle, couch_angle);

  constexpr auto distal_only = true;
  const auto wepl_voi = crl::wepl::CalculateWEPLtoVOI<distal_only>(
      &orig_voi, gantry_angle, couch_angle, ct_img, recalc_img);

  /* Generate a vector of vectors with distances */
  // const auto output =
  //    get_signed_difference(wepl_voi, orig_voi.get()); //, basis);

  for (auto &cur_voi : ss->slist) {
    if (cur_voi.name.find(voi_str) != std::string::npos) {
      cur_voi = *wepl_voi;
    }
  }
  const auto out_dcm = "RS.wepl_structure_"s;
  // http://www.plastimatch.org/plastimatch.html#plastimatch-dice
  const fs::path out_dcm_file(out_dcm + std::string(voi_str) + "_"s +
                              recalc_dcm_dir.filename().string() +
                              descriptive_suffix);
  std::cerr << "Writing WEPL struct to dicom...\n";
  if (!crl::AlterData_RTStructureSetStorage(
          cbctrecon_test->m_cbctrecon->m_strPathRS, ss, out_dcm_file)) {
    std::cerr << "\a"
              << "Could not write dcm\n";
    return -7;
  }

  const auto distal_orig_voi =
      crl::roi_to_distal_only_roi(orig_voi, gantry_angle, couch_angle);

  const auto &cref_wepl_voi = *wepl_voi;

  constexpr auto percent = 95;
  const auto hd_orig_to_wepl =
      crl::calculate_hausdorff<float, percent>(distal_orig_voi, cref_wepl_voi);
  const auto hd_wepl_to_orig =
      crl::calculate_hausdorff<float, percent>(cref_wepl_voi, distal_orig_voi);

  // pCT vs rCT
  const auto rct_str = recalc_dcm_dir.filename().string();
  const auto rct_voi = get_rct_voi(rct_ss.get(), rct_str, orig_voi.name);

  const auto distal_rct_voi =
      crl::roi_to_distal_only_roi(rct_voi, gantry_angle, couch_angle);

  const auto hd_orig_to_rct =
      crl::calculate_hausdorff<float, percent>(distal_orig_voi, distal_rct_voi);
  const auto hd_rct_to_orig =
      crl::calculate_hausdorff<float, percent>(distal_rct_voi, distal_orig_voi);

  // rCT vs WEPL
  const auto hd_rct_to_wepl =
      crl::calculate_hausdorff<float, percent>(distal_rct_voi, cref_wepl_voi);
  const auto hd_wepl_to_rct =
      crl::calculate_hausdorff<float, percent>(cref_wepl_voi, distal_rct_voi);

  const auto str_percent = crl::stringify(percent);

  /* Write distances to file */
  const auto pos_dose = orig_voi.name.find_first_of('D');
  const auto pos_pct = orig_voi.name.find_first_of('%');
  const auto sv_dose_level = std::string_view(&orig_voi.name.at(pos_dose + 1),
                                              (pos_pct - 1) - (pos_dose + 1));
  auto dose_level = 0;
  if (!crl::from_sv(sv_dose_level, dose_level)) {
    std::cerr << "Could not determine dose level from: " << sv_dose_level
              << "\n";
    return -8;
  }

  const auto out_filebase =
      "Hausdorff_" + recalc_dcm_path.parent_path().filename().string();

  auto write_hd_max = [dose_level, gantry_angle, rct_str, out_filebase](
                          const std::string &name,
                          const crl::hausdorff_result<float> &hd_lhs,
                          const crl::hausdorff_result<float> &hd_rhs) {
    fs::path out_file = out_filebase + "_" + name + "_max.txt";
    return write_hd_to_file(out_file, dose_level, gantry_angle, rct_str, hd_lhs,
                            hd_rhs,
                            [](auto x, auto y) { return std::max(x, y); });
  };
  auto write_hd_min = [dose_level, gantry_angle, rct_str, out_filebase](
                          const std::string &name,
                          const crl::hausdorff_result<float> &hd_lhs,
                          const crl::hausdorff_result<float> &hd_rhs) {
    fs::path out_file = out_filebase + "_" + name + "_min.txt";
    return write_hd_to_file(out_file, dose_level, gantry_angle, rct_str, hd_lhs,
                            hd_rhs,
                            [](auto x, auto y) { return std::min(x, y); });
  };
  auto write_hd_directed = [dose_level, gantry_angle, rct_str, out_filebase](
                               const std::string &name,
                               const crl::hausdorff_result<float> &hd_lhs) {
    fs::path out_file = out_filebase + "_" + name + ".txt";
    return write_hd_to_file(out_file, dose_level, gantry_angle, rct_str, hd_lhs,
                            hd_lhs, [](auto x, auto y) { return x; });
  };

  auto success = write_hd_max("pct_vs_wepl", hd_orig_to_wepl, hd_wepl_to_orig);
  success =
      success && write_hd_min("pct_vs_wepl", hd_orig_to_wepl, hd_wepl_to_orig);
  success = success && write_hd_directed("pct_vs_wepl", hd_orig_to_wepl);
  success = success && write_hd_directed("wepl_vs_pct", hd_wepl_to_orig);

  success =
      success && write_hd_max("rct_vs_wepl", hd_rct_to_wepl, hd_wepl_to_rct);
  success =
      success && write_hd_min("rct_vs_wepl", hd_rct_to_wepl, hd_wepl_to_rct);
  success = success && write_hd_directed("rct_vs_wepl", hd_rct_to_wepl);
  success = success && write_hd_directed("wepl_vs_rct", hd_wepl_to_rct);

  success =
      success && write_hd_max("pct_vs_rct", hd_orig_to_rct, hd_rct_to_orig);
  success =
      success && write_hd_min("pct_vs_rct", hd_orig_to_rct, hd_rct_to_orig);
  success = success && write_hd_directed("pct_vs_rct", hd_orig_to_rct);
  success = success && write_hd_directed("rct_vs_pct", hd_rct_to_orig);

  if (!success) {
    return -9;
  }

  std::cerr << "Done!\n";

  return 0;
}
