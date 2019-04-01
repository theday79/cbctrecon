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

#include "itkImageFileReader.h"

#include "cbctrecon_io.h"

#include "cbctrecon_test.hpp"

int main(const int argc, char *argv[]) {

  if (argc < 3) {
    std::cerr << "Usage:\n"
              << argv[0] << " /CB_proj/directory.tar.gz /baseline/CBproj.mha\n";
    return -1;
  }

  std::cerr << "Running test_cb_reader!\n";
  auto cbctrecon_test = std::make_unique<CbctReconTest>();

  /* Load projections */
  const auto cbct_dir_str =
      QString(argv[1]).split(".", QString::SkipEmptyParts).at(0);
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

  const auto start_time = std::chrono::steady_clock::now();
  if (!cbctrecon_test->test_LoadSelectedProjFiles(proj_path, false)) {
    std::cerr << "Could not load or reconstruct CB projections!"
              << "\n";
    return -4;
  }
  if (cbctrecon_test->m_cbctrecon->m_spProjImg3DFloat == nullptr) {
    std::cerr << "Projections was Null!\n";
    return -5;
  }
  const auto end_time = std::chrono::steady_clock::now();
  std::cerr << "Proj. was loaded in: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                     start_time)
                   .count()
            << " ms"
            << "\n";

  auto reader = itk::ImageFileReader<UShortImageType>::New();
  reader->SetFileName(argv[2]);
  reader->Update();
  // Andreas from the past says: You should write m_spProjImgRaw3D to a file and
  // inspect it by eye vs the test data!
  CheckImageQuality<UShortImageType>(
      cbctrecon_test->m_cbctrecon->m_spProjImgRaw3D, reader->GetOutput(), 1e-8,
      100, 2.0);

  return 0;
}
