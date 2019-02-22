// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

// For testing CbctRecon

#include "cbctrecon.h"

#include <iostream>
#include <memory>
#include <numeric>

#include <QDir>

#include "cbctrecon_test.hpp"

UShortImageType::Pointer get_image_from_dicom(QString &dir) {
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
  if (cbctrecon_test->m_cbctrecon->m_spManualRigidCT.IsNull()) {
    std::cerr << "Manual Rigid CT was NULL -> Dicom dir was not read!\n";
  }

  return cbctrecon_test->m_cbctrecon->m_spManualRigidCT;
}

int main(const int argc, char *argv[]) {

  if (argc < 3) {
    std::cerr << "Usage:\n"
              << argv[0] << " /plan/dicom/dir /recalc/dicom/dir\n";
    return -1;
  }

  std::cerr << "Running BiGART script!\n";

  auto recalc_dcm_dir = QString(argv[2]);
  auto recalc_img = get_image_from_dicom(recalc_dcm_dir);

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
  if (cbctrecon_test->m_cbctrecon->m_spManualRigidCT.IsNull()) {
    std::cerr << "Manual Rigid CT was NULL -> Dicom dir was not read!\n";
    return -4;
  }

  /* Structure test: */
  auto ss = cbctrecon_test->m_cbctrecon->m_structures->get_ss(PLAN_CT);

  auto voi = std::string();

  for (auto &structure : ss->slist) {
    if (structure.name.compare("GTV 74/37 LB") == 0) {
      voi = structure.name;
    }
  }
  std::cerr << voi << "\n";

  /* calculate WEPL coordinates */
  UShortImageType::PointType point;
  UShortImageType::IndexType index{0, 0, 0};
  cbctrecon_test->m_cbctrecon->m_spManualRigidCT->TransformIndexToPhysicalPoint(
      index, point);
  std::cerr << "First pixel point: " << point[0] << ", " << point[1] << ", "
            << point[2];

  auto size =
      cbctrecon_test->m_cbctrecon->m_spManualRigidCT->GetLargestPossibleRegion()
          .GetSize();
  index.SetElement(0, size[0]);
  index.SetElement(1, size[1]);
  index.SetElement(2, size[2]);
  cbctrecon_test->m_cbctrecon->m_spManualRigidCT->TransformIndexToPhysicalPoint(
      index, point);
  std::cerr << "Last pixel point: " << point[0] << ", " << point[1] << ", "
            << point[2];

  const auto &offset =
      cbctrecon_test->m_cbctrecon->m_spManualRigidCT->GetOrigin();

  auto &orig_voi = ss->get_roi_ref_by_name(voi);
  std::for_each(std::begin(orig_voi.pslist), std::end(orig_voi.pslist),
                [&offset](Rtss_contour_modern &contour) {
                  std::for_each(std::begin(contour.coordinates),
                                std::end(contour.coordinates),
                                [&offset](FloatVector &coord) {
                                  coord.x += offset[0];
                                  coord.y += offset[1];
                                  coord.z += offset[2];
                                });
                });

  cbctrecon_test->m_cbctregistration->CalculateWEPLtoVOI(
      voi, 90, 0, cbctrecon_test->m_cbctrecon->m_spManualRigidCT, recalc_img);

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
  auto pFile = fopen("out.txt", "w");

  std::for_each(std::begin(output), std::end(output),
                [&pFile](const std::vector<float> out_vec) {
                  std::for_each(std::begin(out_vec), std::end(out_vec),
                                [&pFile](const float val) {
                                  fprintf(pFile, "%.3f\n", val);
                                });
                });

  fclose(pFile);

  return 0;
}
