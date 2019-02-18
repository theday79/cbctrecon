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

#include <iostream>
#include <memory>

#include <QDir>

#include "cbctrecon_test.hpp"

const static auto true_voi_names = std::array<std::string, 13>{
    {"BODY", "CTV1", "PTV1", "Airways_Lungs ", "Brain ", "Spinal cord ",
     "Brain Stem", "Mandible", "Parotid L ", "Parotid R ", "Thyroid L ",
     "Thyroid R ", "Hyoid "}};

int main(const int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage:\n" << argv[0] << " /dicom/directory.tar.gz\n";
    return -1;
  }

  std::cerr << "Running test_dcm_reader!\n";
  auto cbctrecon_test = std::make_unique<CbctReconTest>();

  auto dcm_dir_str = QString(argv[1]).split(".", QString::SkipEmptyParts).at(0);

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
  auto true_vois = std::begin(true_voi_names);
  for (auto &structure : ss->slist) {
    if (true_vois->compare(structure.name) == 0) {
      ++true_vois;
    } else {
      std::cerr << "­\"" << structure.name << "\" is not \"" << *true_vois
                << "\"\n";
      return -5;
    }
  }

  return 0;
}
