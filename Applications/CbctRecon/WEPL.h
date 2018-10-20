#ifndef WEPL_H
#define WEPL_H
#include <array>
#include <vector>

#include "PlmWrapper.h" // for FloatVector, Rtss_contour_modern
#include "cbctrecon.h"  // for FloatImageType, UShortImageType

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

struct WEPLVector {
  double WEPL;
  FloatVector point;
};

struct DoubleVector {
  double x;
  double y;
  double z;
};

struct IntVector {
  int x;
  int y;
  int z;
};

double WEPL_from_point(std::array<size_t, 3> cur_point_id,
                       std::array<double, 3> vec_basis,
                       std::array<double, 3> vec_cubesize,
                       std::array<size_t, 3> cubedim,
                       const FloatImageType::Pointer &wepl_cube);

std::array<double, 3> get_basis_from_angles(double gantry, double couch);

std::vector<double> WEPL_trace_from_point(
    std::array<size_t, 3> cur_point_id, std::array<double, 3> vec_basis,
    std::array<double, 3> vec_cubesize, std::array<size_t, 3> cubedim,
    const FloatImageType::Pointer &wepl_cube);

std::vector<WEPLVector>
WEPLContourFromRtssContour(const Rtss_contour_modern &rt_contour,
                           std::array<double, 3> vec_basis,
                           const FloatImageType::Pointer &wepl_cube);

DoubleVector point_from_WEPL(DoubleVector start_point, const double fWEPL,
                             std::array<double, 3> vec_basis,
                             const FloatImageType::Pointer &wepl_cube);

FloatVector NewPoint_from_WEPLVector(WEPLVector vwepl,
                                     std::array<double, 3> vec_basis,
                                     const FloatImageType::Pointer &wepl_cube);

FloatImageType::Pointer
ConvertUshort2WeplFloat(UShortImageType::Pointer &spImgUshort);

#endif
