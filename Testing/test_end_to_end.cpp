// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

// For testing CbctRecon

#ifdef USE_TINYREFL
#include <tinyrefl/api.hpp> // MUST BE INCLUDED FIRST, FFS!

#include "cbctrecon.h"
#include "cbctrecon.h.tinyrefl"
#else
#include "cbctrecon.h"
#endif

#include <chrono>
#include <iostream>
#include <memory>

#include <QDir>

#include "cbctrecon_test.hpp"

int main(const int argc, char *argv[]) {

  if (argc < 3) {
    std::cerr << "Usage:\n"
              << argv[0]
              << " ./dicom/directory.tar.gz ./CB_proj/directory.tar.gz\n";
    return -1;
  }

  std::cerr << "Running cbctrecon_test!\n";
  auto cbctrecon_test = std::make_unique<CbctReconTest>();

  const auto dcm_dir_str =
      QString(argv[1]).split(".", QString::SkipEmptyParts).at(0);

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

  /* Load projections (Needs to be uploaded to girder first) */
  const auto cbct_dir_str =
      QString(argv[2]).split(".", QString::SkipEmptyParts).at(0);
  auto cbct_dir = QDir(cbct_dir_str);
  auto cbct_path = cbct_dir.absolutePath();
  if (!cbct_dir.exists()) {
    std::cerr << "Directory didn't exist: " << cbct_path.toStdString() << "\n";
    return -2;
  }
  if (cbct_dir.isEmpty(QDir::AllEntries | QDir::NoDotAndDotDot)) {
    std::cerr << "Directory was empty: " << cbct_path.toStdString() << "\n";
    return -3;
  }

  // Set and guess some member variables from projection directory:
  auto proj_dir = QDir(cbct_path + "/Acquisitions/746879825/");
  if (!proj_dir.exists()) {
    std::cerr << "Projection directory: "
              << proj_dir.absolutePath().toStdString() << " doesn't exists!\n";
    return -2;
  }
  auto proj_path = proj_dir.absolutePath();
  cbctrecon_test->test_SetHisDir(proj_path);

  auto start_time = std::chrono::steady_clock::now();
  if (!cbctrecon_test->test_LoadSelectedProjFiles(proj_path, true)) {
    std::cerr << "Could not load or reconstruct CB projections!"
              << "\n";
    return -4;
  }
  if (cbctrecon_test->m_cbctrecon->m_spRawReconImg == nullptr) {
    std::cerr << "Raw reconstruction was Null!\n";
    return -5;
  }
  auto end_time = std::chrono::steady_clock::now();
  std::cerr << "Proj. was loaded and reconstructed in: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                     start_time)
                   .count()
            << " ms"
            << "\n";

  /* Scatter correction algorithm "Batch" style */
  // Move MovingImage by -7.0, -85.0, -237.0: By emulating
  // DlgRegistration::ImageManualMove

  cbctrecon_test->m_dlgRegistration->SLT_KeyMoving(true);
  cbctrecon_test->m_dlgRegistration->ImageManualMoveOneShot(-7.0f, -85.0f,
                                                            -237.0f);

  cbctrecon_test->m_dlgRegistration->SLT_ConfirmManualRegistration();
  cbctrecon_test->m_dlgRegistration->SLT_DoRegistrationRigid();
  cbctrecon_test->m_dlgRegistration->SLT_DoRegistrationDeform();

  /* WEPL structure test: */
  auto ss = cbctrecon_test->m_cbctrecon->m_structures->get_ss(DEFORM_CT);
  for (auto &structure : ss->slist) {
    std::cerr << structure.name << "\n";
  }

  const auto voi = std::string("CTV1");
  start_time = std::chrono::steady_clock::now();
  cbctrecon_test->m_cbctregistration->CalculateWEPLtoVOI(
      voi, 45, 45, cbctrecon_test->m_cbctrecon->m_spDeformedCT_Final,
      cbctrecon_test->m_cbctrecon->m_spRawReconImg);
  end_time = std::chrono::steady_clock::now();
  std::cerr << "WEPL was calculated in: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                     start_time)
                   .count()
            << " ms"
            << "\n";

  /* Some verification of the WEPL results should go here */

  return 0;
}
