#ifndef WEPL_H
#define WEPL_H
#include <array>
#include <vector>

#include "cbctrecon_types.h" // for FloatImageType, UShortImageType

class Rtss_contour_modern;

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }

double WEPL_from_point(const std::array<size_t, 3> &cur_point_id,
                       const std::array<double, 3> &vec_basis,
                       const std::array<double, 3> &vec_cubesize,
                       const std::array<size_t, 3> &cubedim,
                       const FloatImageType::Pointer &wepl_cube);

std::array<double, 3> get_basis_from_angles(double gantry, double couch);

std::vector<double> WEPL_trace_from_point(
    const std::array<size_t, 3> &cur_point_id, const std::array<double, 3> &vec_basis,
    const std::array<double, 3> &vec_cubesize, const std::array<size_t, 3> &cubedim,
    const FloatImageType::Pointer &wepl_cube);

std::vector<WEPLVector>
WEPLContourFromRtssContour(const Rtss_contour_modern &rt_contour,
                           const std::array<double, 3> &vec_basis,
                           const FloatImageType::Pointer &wepl_cube);

DoubleVector point_from_WEPL(DoubleVector &start_point, double fWEPL,
                             const std::array<double, 3> &vec_basis,
                             const FloatImageType::Pointer &wepl_cube);

FloatVector NewPoint_from_WEPLVector(const WEPLVector &vwepl,
                                     const std::array<double, 3> &vec_basis,
                                     const FloatImageType::Pointer &wepl_cube);

FloatImageType::Pointer
ConvertUshort2WeplFloat(UShortImageType::Pointer &spImgUshort);

#endif
