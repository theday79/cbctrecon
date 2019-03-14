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
#include <numeric>

#include <QDir>

#include "cbctrecon_test.hpp"

static const auto true_voi_names = std::array<std::string, 13>{
    {"BODY", "CTV1", "PTV1", "Airways_Lungs ", "Brain ", "Spinal cord ",
     "Brain Stem", "Mandible", "Parotid L ", "Parotid R ", "Thyroid L ",
     "Thyroid R ", "Hyoid "}};

// Oops, this may depend on the binary reprentation of floats on the given
// system!
static const auto true_super_sums = std::array<std::string, 13>{
    {"81b8a34b", "aa254449", "12079d49", "9a7de24a", "7718a24a", "49957c49",
     "9a027048", "6069a649", "1054e548", "3063a348", "a8fd6b48", "4d333c48",
     "9e1f5a48"}};

template <typename T> std::string float_to_hex(T x) {
  const auto pf = reinterpret_cast<const unsigned char *>(&x);
  std::stringstream ss;
  ss << std::hex << std::setfill('0');
  for (auto i = 0U; i < sizeof(float); ++i) {
    ss << std::setw(2) << static_cast<unsigned>(pf[i]);
  }
  return ss.str();
}

int main(const int argc, char *argv[]) {

  if (argc < 2) {
    std::cerr << "Usage:\n" << argv[0] << " /dicom/directory.tar.gz\n";
    return -1;
  }

  std::cerr << "Running test_dcm_reader!\n";
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

  /* Structure test: */
  auto ss = cbctrecon_test->m_cbctrecon->m_structures->get_ss(PLAN_CT);
  auto true_vois = std::begin(true_voi_names);
  auto true_super_sum = std::begin(true_super_sums);
  for (auto &structure : ss->slist) {
    if (*true_vois == structure.name) {
      ++true_vois;
    } else {
      std::cerr << "­\"" << structure.name << "\" is not \"" << *true_vois
                << "\"\n";
      return -5;
    }
    auto super_sum = 0.0f;
    for (auto &contour : structure.pslist) {
      const auto start = FloatVector{.0f, .0f, .0f};
      const auto sum = std::accumulate(
          std::begin(contour.coordinates), std::end(contour.coordinates), start,
          [](auto acc_val, auto val) {
            return FloatVector{val.x + acc_val.x, val.y + acc_val.y,
                               val.z + acc_val.z};
          });
      super_sum += sum.x + sum.y + sum.z;
    }

    if (*true_super_sum == float_to_hex(super_sum)) {
      ++true_super_sum;
    } else {
      std::cerr << "\"" << float_to_hex(super_sum) << "\" != \""
                << *true_super_sum << "\"\n";
    }
  }

  return 0;
}
