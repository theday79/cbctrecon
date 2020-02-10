#ifndef FREE_FUNCTIONS_H
#define FREE_FUNCTIONS_H

#include <algorithm>
#include <charconv>
#include <execution>
#include <filesystem>
#include <limits>
#include <optional>
#include <string_view>
#include <vector>

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

class YK16GrayImage;
class CbctRecon;

// CbctReconLib
namespace crl {

CBCTRECON_API
std::string HexStr2IntStr(std::string_view str_hex);

CBCTRECON_API
std::unique_ptr<YK16GrayImage>
ApplyCalibrationMaps(CbctRecon *p_cr, YK16GrayImage *const &rawImg,
                     const bool DarkCorr, const bool GainCorr,
                     const bool DefectCorr);

CBCTRECON_API
std::filesystem::path CorrectSingleFile(CbctRecon *p_cr,
                                        const std::filesystem::path &filePath,
                                        const bool DarkCorr,
                                        const bool GainCorr,
                                        const bool DefectCorr);

CBCTRECON_API
void CorrectSingleFile(CbctRecon *p_cr, YK16GrayImage *pYKRawImg,
                       const bool DarkCorr, const bool GainCorr,
                       const bool DefectCorr);

CBCTRECON_API
void LoadBadPixelMap(std::vector<BADPIXELMAP> &vPixelReplMap,
                     const std::filesystem::path &filePath);

CBCTRECON_API
std::unique_ptr<YK16GrayImage>
BadPixReplacement(std::vector<BADPIXELMAP> &vPixelReplMap,
                  std::unique_ptr<YK16GrayImage> targetImg);

CBCTRECON_API
std::vector<std::string> GetProjFileNames(
    CbctRecon *p_cr,
    std::filesystem::path &dirPath); // main loading fuction for projections

CBCTRECON_API
std::vector<size_t> GetExcludeProjFiles(
    const rtk::ThreeDCircularProjectionGeometry::Pointer spFullGeometry,
    const bool bScanDirectionCW, const bool bManAngleGap,
    const double gantryAngleInterval);

CBCTRECON_API
void GetSelectedIndices(const std::vector<double> &vFullAngles,
                        std::vector<double> &vNormAngles,
                        std::vector<size_t> &vTargetIdx, const bool bCW,
                        std::vector<size_t> &vExcludingIdx);

CBCTRECON_API
bool IsFileNameOrderCorrect(std::vector<std::string> &vFileNames);

CBCTRECON_API
void CopyDictionary(itk::MetaDataDictionary &fromDict,
                    itk::MetaDataDictionary &toDict); // NOT COMPLETED YET!!
                                                      // Export DICOM without
                                                      // Source DICOM is not
                                                      // possible

/// Template Meta Functions:

template <typename T> bool from_sv(std::string_view sv, T &val) {
  if (auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), val);
      ec == std::errc()) {
    return true;
  }
  return false;
}

// code inspired by
// https://marcoarena.wordpress.com/2017/01/03/string_view-odi-et-amo and B.
// Filipek's book C++17 in Detail
std::vector<std::string_view> split_string(std::string_view strv,
                                           std::string_view delims = " ") {
  std::vector<std::string_view> output;
  auto first = strv.begin();

  while (first != strv.end()) {
    const auto second = std::find_first_of(
        first, std::cend(strv), std::cbegin(delims), std::cend(delims));

    if (first != second) {
      output.emplace_back(strv.substr(std::distance(strv.begin(), first),
                                      std::distance(first, second)));
    }

    if (second == strv.end())
      break;

    first = std::next(second);
  }

  return output;
}

template <bool spaces_only = false>
constexpr std::string_view trim_string(const std::string_view string_to_trim) {
  const auto whitespaces = spaces_only ? " " : "\t\n\v\f\r ";
  const auto new_start = string_to_trim.find_first_not_of(whitespaces);
  const auto new_end = string_to_trim.find_last_not_of(whitespaces);
  return string_to_trim.substr(new_start, new_end - new_start + 1);
}

template <typename T>
constexpr std::optional<T> from_string(const std::string &number) {
  T out_var;
  if (auto [p, ec] = std::from_chars(number.c_str(),
                                     number.c_str() + number.length(), out_var);
      ec != std::errc()) {
    return out_var;
  }
  return std::nullopt;
}

template <typename T> constexpr T ce_pow(T val, size_t exponent) {
  T out_val = 1;
  for (auto i_exp = 0; i_exp < exponent; ++i_exp) {
    out_val *= val;
  }
  return out_val;
}

template <typename T>
constexpr T BeamHardModel(const T val, const T a, const T b, const T c,
                          const T d) {
  return a * ce_pow(val, 3) + b * ce_pow(val, 2) + c * val + d;
}

template <enProjFormat PF, typename T> constexpr T BeamHardModel(const T val);

template <typename T>
constexpr T BeamHardModel<enProjFormat::HND_FORMAT, T>(const T val) {
  // a * x^3 + b * x^2 + c * x + d
  return 6.0e-08 * ce_pow(val, 3) - 1.0e-08 * ce_pow(val, 2) - 5.0e-07 * val +
         8.0e-01;
}

template <typename T>
constexpr T BeamHardModel<enProjFormat::HIS_FORMAT, T>(const T val) {
  // a * x^3 + b * x^2 + c * x + d
  return 9.321e-05 * ce_pow(val, 3) - 2.609e-03 * ce_pow(val, 2) +
         3.374e-02 * val + 9.691e-01;
}

template <typename T>
constexpr T BeamHardModel<enProjFormat::XIM_FORMAT, T>(
    const T val) { // a * x^3 + b * x^2 + c * x + d
  return 6.0e-08 * ce_pow(val, 3) + 9.0e-5 * ce_pow(val, 2) + 1.0e-2 * val +
         0.8;
}

template <enProjFormat PF, typename T>
void BeamHardening(T *pBuffer, const int nPix) {
  std::transform(std::execution::par_unseq, &pBuffer[0], &pBuffer[0] + nPix,
                 &pBuffer[0], [](const T val) {
                   return (val < 1.189 ? val
                                       : val * BeamHardModel<PF, T>(val)) +
                          (PF == enProjFormat::HIS_FORMAT ? 0 : 1.47);
                 });
}

template <typename T>
constexpr int ce_sgn(const T val) { return (T(0) < val) - (val < T(0)); }

template <typename T>
constexpr T ce_heaviside(const T x) { return static_cast<T>(.5 * ce_sgn(x) + 0.5); }

template <typename T> constexpr T ce_abs(const T x) {
  return x < T(0) ? -x : x;
}
 
template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0 >
T constexpr sqrtNewtonRaphson(T x, T curr, T prev) {
  return curr == prev ? curr
                      : sqrtNewtonRaphson(x, 0.5 * (curr + x / curr), curr);
}

/*
 * Constexpr version of the square root
 * Return value:
 *	- For a finite and non-negative value of "x", returns an approximation
 *for the square root of "x"
 *   - Otherwise, returns NaN
 */
template <typename T, typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
T constexpr ce_sqrt(T x) {
  return x >= 0 && x < std::numeric_limits<T>::infinity()
             ? sqrtNewtonRaphson<T>(x, x, 0)
             : std::numeric_limits<T>::quiet_NaN();
}

template<typename T, typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
constexpr
T fullFan_subFunction(const T a, const T b,
                                  const T c, const T d,
                                  const T x) {
  return c - ce_sqrt(ce_abs(ce_pow(a, 2) - ce_pow(x * d - b, 2))) *
                 ce_heaviside(x * d - b + a) * ce_heaviside(-(x * d - b - a));
}

template<typename T, typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
inline T fullFan_Function(const T a, const T b, const T c,
                               const T d, const T e, const T x) {
  auto es = std::array<int, 7>();
  std::iota(es.begin(), es.end(), -3); //iota not constexpr until C++20
  return std::accumulate(es.begin(), es.end(), T(0), [=](auto i_e) {
    return fullFan_subFunction(a, b, c, d, x - i_e * e);
  }) / es.size();
}

/// Functors:

// From line integral to raw intensity
class LineInt2Intensity {
public:
  LineInt2Intensity() = default;
  ~LineInt2Intensity() = default;
  float operator()(const float val) const {
    float intensityVal = std::exp(-val) /* I_0=1 */;
    return intensityVal;
  }
};
// From raw intensity to line integral
class Intensity2LineInt {
public:
  Intensity2LineInt() = default;
  ~Intensity2LineInt() = default;
  float operator()(const float val) const {
    // mu = ln(I_0/I) OR mu = ln(I/I0)
    constexpr auto I_0_div_I_air = 1.0 / (0.1541 * 1.225e-3 * 10.0);
    float lineintVal = std::log(I_0_div_I_air); // 10 cm of air

    if (val > 0) {
      lineintVal = /* log(I_0=1) = 0 */ -std::log(val);
    }
    return lineintVal;
  }
};

} // namespace crl

#endif // FREE_FUNCTIONS_H