#ifndef FREE_FUNCTIONS_H
#define FREE_FUNCTIONS_H

#include <algorithm>
#include <charconv>
#include <execution>
#include <filesystem>
#include <limits>
#include <optional>
#include <string_view>
#include <type_traits>
#include <variant>
#include <vector>

#if defined(__GNUC__)
#include "absl/strings/charconv.h"
#define float_from_chars absl::from_chars
#else
#define float_from_chars std::from_chars
#endif

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

class YK16GrayImage;
class CbctRecon;

// CbctReconLib
namespace crl {

CBCTRECON_API
std::vector<std::string_view> split_string(std::string_view strv,
                                           std::string_view delims = " ");

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

template <typename T, std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
bool is_scan_direction_CW(const std::vector<T> &angles) {

  std::vector<T> vTempConvAngles;

  const auto itBegin = angles.begin();
  const auto itEnd = angles.end();

  for (auto it = itBegin; it != itEnd; ++it) {
    auto tmpAngle = *it;

    if (tmpAngle > T(180.0)) {
      tmpAngle = tmpAngle - T(360.0);
    }

    vTempConvAngles.push_back(tmpAngle);
  }
  const auto geoDataSize = angles.size();
  // compare 2 points in the middle of the angle list
  const auto iLowerIdx = static_cast<size_t>(geoDataSize * T(1.0) / T(3.0));
  const auto iUpperIdx = static_cast<size_t>(geoDataSize * T(2.0) / T(3.0));
  auto bScanDirectionCW = false;
  if (vTempConvAngles.at(iLowerIdx) <
      vTempConvAngles.at(iUpperIdx)) // ascending
  {
    bScanDirectionCW = true;
    std::cout << "The scan direction is CW" << std::endl;
  } else {
    std::cout << "The scan direction is CCW" << std::endl;
  }
  return bScanDirectionCW;
}

/// Template Meta Functions:

template <typename T> constexpr bool from_sv(std::string_view sv, T &val) {
  bool success = false;
  if constexpr (std::is_floating_point_v<T>) {
    if (auto [ptr, ec] =
            float_from_chars(sv.data(), sv.data() + sv.size(), val);
        ec == std::errc()) {
      success = true;
    }
  } else {
    if (auto [ptr, ec] = std::from_chars(sv.data(), sv.data() + sv.size(), val);
        ec == std::errc()) {
      success = true;
    }
  }
  return success;
}

template <>
constexpr bool
from_sv<std::variant<int, double>>(std::string_view sv,
                                   std::variant<int, double> &val) {
  bool success = false;
  if (sv.find_first_of('.') != std::string::npos) {
    double tmp_val = 0.0;
    if (from_sv(sv, tmp_val)) {
      val = tmp_val;
      success = true;
    }
  } else {
    int tmp_val = 0;
    if (from_sv(sv, tmp_val)) {
      val = tmp_val;
      success = true;
    }
  }
  return success;
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
  if (from_sv(number, out_var)) {
    return out_var;
  }
  return std::nullopt;
}

template <typename T> std::string stringify(T arg) {
  if constexpr (std::is_floating_point_v<T> || std::is_integral_v<T>) {
    return std::to_string(arg);
  } else if constexpr (std::is_convertible_v<T, std::string>) {
    return arg;
  } else if constexpr (std::is_constructible_v<std::string, T>) {
    return std::string(arg);
  } else if constexpr (std::is_same_v<T, std::filesystem::path>) {
    return arg.string();
  }
  return "UNHANDLED_ERROR";
}

template <char SEP, typename... Args> std::string make_sep_str(Args &&... arg) {
  auto tmp_str = ((stringify(std::forward<Args>(arg)) + SEP) + ...);
  return tmp_str.substr(0, tmp_str.size() - 1);
}

template <typename T>
auto from_sv_v(const std::vector<std::string_view> &sv_v) {
  auto v_val = std::vector<T>();
  std::transform(sv_v.begin(), sv_v.end(), std::back_inserter(v_val),
                 [](std::string_view sv) {
                   T val;
                   crl::from_sv(sv, val);
                   return val;
                 });
  return v_val;
}

template <typename T> constexpr T ce_pow(T val, size_t exponent) {
  T out_val = 1;
  for (size_t i_exp = 0; i_exp < exponent; ++i_exp) {
    out_val *= val;
  }
  return out_val;
}

template <typename T>
constexpr T BeamHardModelCustom(const T val, const T a, const T b, const T c,
                                const T d) {
  return a * ce_pow(val, 3) + b * ce_pow(val, 2) + c * val + d;
}

template <enProjFormat PF, typename T> constexpr T BeamHardModel(const T val) {

  static_assert(PF != enProjFormat::HND_FORMAT ||
                    PF != enProjFormat::HIS_FORMAT ||
                    PF != enProjFormat::XIM_FORMAT,
                "Projection format not recognised!");

  if constexpr (PF == enProjFormat::HND_FORMAT) {
    // a * x^3 + b * x^2 + c * x + d
    return 6.0e-08 * ce_pow(val, 3) - 1.0e-08 * ce_pow(val, 2) - 5.0e-07 * val +
           8.0e-01;
  } else if constexpr (PF == enProjFormat::HIS_FORMAT) {
    // a * x^3 + b * x^2 + c * x + d
    return 9.321e-05 * ce_pow(val, 3) - 2.609e-03 * ce_pow(val, 2) +
           3.374e-02 * val + 9.691e-01;
  } else if constexpr (PF == enProjFormat::XIM_FORMAT) {
    // a * x^3 + b * x^2 + c * x + d
    return 6.0e-08 * ce_pow(val, 3) + 9.0e-5 * ce_pow(val, 2) + 1.0e-2 * val +
           0.8;
  }
}

template <enProjFormat PF, typename T> T BeamHardening(T val) {
  return (val < 1.189 ? val : val * BeamHardModel<PF>(val)) +
         (PF == enProjFormat::HIS_FORMAT ? 0 : 1.47);
}

template <enProjFormat PF, typename T>
void BeamHardening(T *pBuffer, const size_t nPix) {
  std::transform(std::execution::par_unseq, &pBuffer[0], &pBuffer[0] + nPix,
                 &pBuffer[0], [](auto val) { return BeamHardening<PF>(val); });
}

template <typename T> constexpr int ce_round(const T val) {
  static_assert(std::is_floating_point_v<T>,
                "Rounding not necessary for non-float values");
  return (val >= 0.0) ? int(val + 0.5) : int(val - 0.5);
}

template <typename T> constexpr int ce_sgn(const T val) {
  return (T(0) < val) - (val < T(0));
}

template <typename T> constexpr T ce_heaviside(const T x) {
  return static_cast<T>(.5 * ce_sgn(x) + 0.5);
}

template <typename T> constexpr T ce_abs(const T x) {
  return x < T(0) ? -x : x;
}

template <typename T,
          typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
T constexpr sqrtNewtonRaphson(T x, T curr, T prev) {
  return curr == prev ? curr
                      : sqrtNewtonRaphson<T>(x, 0.5 * (curr + x / curr), curr);
}

/*
 * Constexpr version of the square root
 * Return value:
 *	- For a finite and non-negative value of "x", returns an approximation
 *for the square root of "x"
 *   - Otherwise, returns NaN
 */
template <typename T,
          typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
T constexpr ce_sqrt(T x) {
  return x >= 0 && x < std::numeric_limits<T>::infinity()
             ? sqrtNewtonRaphson<T>(x, x, 0)
             : std::numeric_limits<T>::quiet_NaN();
}

template <typename T,
          typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
constexpr T fullFan_subFunction(const T a, const T b, const T c, const T d,
                                const T x) {
  return c - ce_sqrt(ce_abs(ce_pow(a, 2) - ce_pow(x * d - b, 2))) *
                 ce_heaviside(x * d - b + a) * ce_heaviside(-(x * d - b - a));
}

template <typename T,
          typename std::enable_if_t<std::is_floating_point_v<T>, int> = 0>
inline T fullFan_Function(const T a, const T b, const T c, const T d, const T e,
                          const T x) {
  auto es = std::array<int, 7>();
  std::iota(es.begin(), es.end(), -3); // iota not constexpr until C++20
  return std::accumulate(es.begin(), es.end(), T(0),
                         [=](T sum, auto i_e) {
                           return sum +
                                  fullFan_subFunction(a, b, c, d, x - i_e * e);
                         }) /
         es.size();
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

class LineInt2Intensity_ushort {
public:
  LineInt2Intensity_ushort() = default;
  ~LineInt2Intensity_ushort() = default;
  float operator()(const float val) const {
    constexpr auto max_ushort = std::numeric_limits<unsigned short>::max();
    float intensityVal = std::exp(static_cast<double>(val) * -1.0) *
                         static_cast<double>(max_ushort);

    if (intensityVal <= 1.0) {
      intensityVal = 1.0;
    }
    if (intensityVal >= (max_ushort - 1)) {
      intensityVal = static_cast<double>(max_ushort - 1);
    }

    return static_cast<unsigned short>(intensityVal);
  }
};

template <typename Tinput> class Intensity2LineInt_ushort {
public:
  Intensity2LineInt_ushort() = default;
  ~Intensity2LineInt_ushort() = default;
  float operator()(const Tinput val) const {
    constexpr auto max_ushort = std::numeric_limits<unsigned short>::max();
    // mu = ln(I_0/I) OR mu = ln(I/I0)
    const float mu_t_val =
        log(static_cast<double>(max_ushort) / static_cast<double>(val));

    return mu_t_val;
  }
};

} // namespace crl

#endif // FREE_FUNCTIONS_H
