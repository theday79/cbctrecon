#ifndef WEPL_H
#define WEPL_H
#include <array>
#include <vector>

#include "cbctrecon_config.h"
#include "cbctrecon_types.h" // for FloatImageType, UShortImageType

#include <itkPoint.h>

struct Rtss_roi_modern;
struct Rtss_contour_modern;

namespace crl {

namespace wepl {

double WEPL_from_point(const std::array<size_t, 3> &cur_point_id,
                       const std::array<double, 3> &vec_basis,
                       const std::array<double, 3> &vec_voxelsize,
                       const FloatImageType::Pointer &wepl_cube);

CBCTRECON_API std::array<double, 3> get_basis_from_angles(double gantry,
                                                          double couch);

std::vector<double>
WEPL_trace_from_point(const std::array<size_t, 3> &cur_point_id,
                      const std::array<double, 3> &vec_basis,
                      const std::array<double, 3> &vec_cubesize,
                      const FloatImageType::Pointer &wepl_cube);

CBCTRECON_API
std::vector<FloatVector>
distal_points_only(const Rtss_contour_modern &points,
                   const std::array<double, 3> &direction);

std::vector<WEPLVector>
WEPLContourFromRtssContour(const Rtss_contour_modern &rt_contour,
                           const std::array<double, 3> &vec_basis,
                           const FloatImageType::Pointer &wepl_cube);

std::vector<WEPLVector>
DistalWEPLContourFromRtssContour(const Rtss_contour_modern &rt_contour,
                           const std::array<double, 3> &vec_basis,
                           const FloatImageType::Pointer &wepl_cube);

FloatImageType::PointType
point_from_WEPL(const vnl_vector_fixed<double, 3> &start_point, double fWEPL,
                const vnl_vector_fixed<double, 3> &vec_basis,
                const FloatImageType::Pointer &wepl_cube);

FloatVector NewPoint_from_WEPLVector(const WEPLVector &vwepl,
                                     const std::array<double, 3> &arr_basis,
                                     const FloatImageType::Pointer &wepl_cube);

FloatImageType::Pointer
ConvertUshort2WeplFloat(const UShortImageType::Pointer &spImgUshort);


template<bool DISTAL_ONLY=false>
Rtss_roi_modern *CalculateWEPLtoVOI(const Rtss_roi_modern *voi,
                                    const double gantry_angle,
                                    const double couch_angle,
                                    const UShortImageType::Pointer &spMoving,
                                    const UShortImageType::Pointer &spFixed) {

  // Get basis from angles
  const auto vec_basis = get_basis_from_angles(gantry_angle, couch_angle);

  // Get Fixed and Moving
  // Tranlate fixed and moving to dEdx
  const auto wepl_cube = ConvertUshort2WeplFloat(spMoving);
  const auto wepl_cube_fixed = ConvertUshort2WeplFloat(spFixed);

  // Initialize WEPL contour
  auto WEPL_voi = std::make_unique<Rtss_roi_modern>();
  WEPL_voi->name = "WEPL" + voi->name;
  WEPL_voi->color = "255 0 0";
  WEPL_voi->id = voi->id;   /* Used for import/export (must be >= 1) */
  WEPL_voi->bit = voi->bit; /* Used for ss-img (-1 for no bit) */
  WEPL_voi->pslist.resize(voi->pslist.size());

  auto i = 0U;
  // Calculate WEPL
  for (const auto &contour : voi->pslist) {
    auto WEPL_contour = std::make_unique<Rtss_contour_modern>(contour);
    WEPL_contour->ct_slice_uid = contour.ct_slice_uid;
    WEPL_contour->slice_no = contour.slice_no;
    WEPL_contour->coordinates.clear();

    // Actually calculate WEPL on spMoving
    auto WEPL_points =
        DISTAL_ONLY
            ? DistalWEPLContourFromRtssContour(contour, vec_basis, wepl_cube)
            : WEPLContourFromRtssContour(contour, vec_basis, wepl_cube);

    // Inversely calc WEPL on spFixed
    // And put WEPL point in contour
    std::transform(std::begin(WEPL_points), std::end(WEPL_points),
                   std::back_inserter(WEPL_contour->coordinates),
                   [&vec_basis, &wepl_cube_fixed](const WEPLVector &val) {
                     return NewPoint_from_WEPLVector(val, vec_basis,
                                                     wepl_cube_fixed);
                   });
    WEPL_voi->pslist.at(i++) = *WEPL_contour.release();
  }
  return WEPL_voi.release();
}

} // namespace wepl

} // namespace crl

#endif
