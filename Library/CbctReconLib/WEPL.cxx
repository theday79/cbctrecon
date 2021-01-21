// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include <algorithm> // for find_if, min_element, transform
#include <array>
#include <chrono>   // for chrono
#include <cmath>    // for sqrt, pow, sin, cos
#include <iostream> // for operator<<, endl, basic_ostream, cerr, ostream
#include <numeric>  // for accumulate, adjacent_difference
#include <utility>  // for pair
#include <valarray> // for valarray

#include <itkFixedArray.h> // for FixedArray
#include <itkImage.h> // for Image<>::Pointer, Image, Image<>::IndexType, Image<>::PointType
#include <itkImageFileWriter.h> // for ImageFileWriter
#include <itkLinearInterpolateImageFunction.h>
#include <itkMacro.h> // for CastImageFilter::New, ImageFileWriter::New
#include <itkUnaryFunctorImageFilter.h>
#include <vnl/vnl_vector_fixed.h> // for vnl_vector_fixed

#include "StructureSet.h"
#include "WEPL.h"
#include "cbctrecon_io.h"
#include "cbctrecon_types.h"
#include "free_functions.h"

namespace crl {
namespace wepl {

double lin_interpolate_old(const std::array<int, 3> &point_id,
                           const std::array<double, 3> &point_id_pos,
                           const int idx_2, const int idy_2, const int idz_2,
                           const FloatImageType::Pointer &wepl_cube) {

  std::array<double, 8> weights{};
  // x                    xyz
  weights.at(0) = std::sqrt( // 000 =
      std::pow(point_id.at(0) - point_id_pos.at(0), 2) +
      std::pow(point_id.at(1) - point_id_pos.at(1), 2) +
      std::pow(point_id.at(2) - point_id_pos.at(2), 2));
  weights.at(1) = std::sqrt( // 100 =
      std::pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
      std::pow(point_id.at(1) - point_id_pos.at(1), 2) +
      std::pow(point_id.at(2) - point_id_pos.at(2), 2));
  // y
  weights.at(2) = std::sqrt( // 010 =
      std::pow(point_id.at(0) - point_id_pos.at(0), 2) +
      std::pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
      std::pow(point_id.at(2) - point_id_pos.at(2), 2));
  weights.at(3) = std::sqrt( // 110 =
      std::pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
      std::pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
      std::pow(point_id.at(2) - point_id_pos.at(2), 2));
  // z
  weights.at(4) = std::sqrt( // 001 =
      std::pow(point_id.at(0) - point_id_pos.at(0), 2) +
      std::pow(point_id.at(1) - point_id_pos.at(1), 2) +
      std::pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
  weights.at(5) = std::sqrt( // 101 =
      std::pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
      std::pow(point_id.at(1) - point_id_pos.at(1), 2) +
      std::pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
  weights.at(6) = std::sqrt( // 011 =
      std::pow(point_id.at(0) - point_id_pos.at(0), 2) +
      std::pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
      std::pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));
  weights.at(7) = std::sqrt( // 111 =
      std::pow(point_id.at(0) + idx_2 - point_id_pos.at(0), 2) +
      std::pow(point_id.at(1) + idy_2 - point_id_pos.at(1), 2) +
      std::pow(point_id.at(2) + idz_2 - point_id_pos.at(2), 2));

  auto sum_weights = 1.0 / std::accumulate(weights.begin(), weights.end(), 0.0);

  const auto divide_by_weight = [&sum_weights](double &val) -> void {
    val = 1 - val * sum_weights;
  };

  std::for_each(std::begin(weights), std::end(weights), divide_by_weight);

  auto out_val = 0.0;
  for (auto i = 0; i < 8; i++) {
    // convert point_id to cube_ids
    FloatImageType::IndexType cube_id{};
    cube_id[0] = point_id.at(0) + idx_2 * (i % 2);     // x= 0,1,0,1,0,1,0,1
    cube_id[1] = point_id.at(1) + idy_2 * (i / 2 % 2); // y= 0,0,1,1,0,0,1,1
    cube_id[2] = point_id.at(2) + idz_2 * (i / 4);     // z= 0,0,0,0,1,1,1,1

    const auto val = wepl_cube->GetPixel(cube_id);
    out_val += val * weights.at(i);
  }
  return out_val;
}

double WEPL_from_point(const std::array<size_t, 3> &cur_point_id,
                       const std::array<double, 3> &vec_basis,
                       const std::array<double, 3> &vec_voxelsize,
                       const FloatImageType::Pointer &wepl_cube) {
  auto interpolator =
      itk::LinearInterpolateImageFunction<FloatImageType, double>::New();
  interpolator->SetInputImage(wepl_cube);
  const auto step_length = 0.1;
  const auto step =
      DoubleVector{vec_basis.at(0) * step_length, vec_basis.at(1) * step_length,
                   vec_basis.at(2) * step_length};

  const auto inv_voxelsize =
      DoubleVector{1.0 / vec_voxelsize.at(0), 1.0 / vec_voxelsize.at(1),
                   1.0 / vec_voxelsize.at(2)};

  auto point = DoubleVector{
      static_cast<double>(cur_point_id.at(0)) * vec_voxelsize.at(0),
      static_cast<double>(cur_point_id.at(1)) * vec_voxelsize.at(1),
      static_cast<double>(cur_point_id.at(2)) * vec_voxelsize.at(2)};

  auto out = 0.0;

  while (true) {
    auto index_pos = itk::ContinuousIndex<double, 3>();
    index_pos[0] = point.x * inv_voxelsize.x;
    index_pos[1] = point.y * inv_voxelsize.y;
    index_pos[2] = point.z * inv_voxelsize.z;

    if (!wepl_cube->GetBufferedRegion().IsInside(index_pos)) {
      break;
    }
    const auto val = interpolator->EvaluateAtContinuousIndex(index_pos);
    if (!std::isnan(val)) { // Hopefully this is only happening at the edge!?
      out += val;
    }

    // point = point - step
    point.x -= step.x;
    point.y -= step.y;
    point.z -= step.z;
  }

  return out * step_length;
}

std::array<double, 3> get_basis_from_angles(double gantry, double couch) {
  gantry += 180.0;
  gantry *= itk::Math::pi / 180.0;
  couch *= itk::Math::pi / 180.0;

  const std::array<double, 3> basis = {
      {sin(gantry) * cos(couch), -cos(gantry), sin(couch) * sin(gantry)}};
  return basis;
}

std::vector<double>
WEPL_trace_from_point(const std::array<size_t, 3> &cur_point_id,
                      const std::array<double, 3> &vec_basis,
                      const std::array<double, 3> &vec_cubesize,
                      const FloatImageType::Pointer &wepl_cube) {
  auto interpolator =
      itk::LinearInterpolateImageFunction<FloatImageType, double>::New();
  interpolator->SetInputImage(wepl_cube);

  std::vector<double> cumWEPL; // cumulative WEPL

  const auto step_length = 0.1;
  const auto step =
      DoubleVector{vec_basis.at(0) * step_length, vec_basis.at(1) * step_length,
                   vec_basis.at(2) * step_length};

  const auto inv_voxelsize =
      DoubleVector{1.0 / vec_cubesize.at(0), 1.0 / vec_cubesize.at(1),
                   1.0 / vec_cubesize.at(2)};

  auto point = DoubleVector{
      static_cast<double>(cur_point_id.at(0)) * vec_cubesize.at(0),
      static_cast<double>(cur_point_id.at(1)) * vec_cubesize.at(1),
      static_cast<double>(cur_point_id.at(2)) * vec_cubesize.at(2)};

  auto out_wepl = 0.0;
  while (true) {

    auto index_pos = itk::ContinuousIndex<double, 3>();
    index_pos[0] = point.x * inv_voxelsize.x;
    index_pos[1] = point.y * inv_voxelsize.y;
    index_pos[2] = point.z * inv_voxelsize.z;

    if (!wepl_cube->GetBufferedRegion().IsInside(index_pos)) {
      break;
    }
    const auto val = interpolator->EvaluateAtContinuousIndex(index_pos);
    if (!std::isnan(val)) { // Hopefully this is only happening at the edge!?
      out_wepl += val;
    }

    cumWEPL.push_back(out_wepl);

    // point = point - step
    point.x -= step.x;
    point.y -= step.y;
    point.z -= step.z;
  }
  std::vector<double> vdiff(cumWEPL.size());
  std::adjacent_difference(std::begin(cumWEPL), std::end(cumWEPL),
                           std::begin(vdiff));

  std::vector<double> out_vec;
  auto rev_wepl = cumWEPL.at(cumWEPL.size() - 1);

  for (auto val : vdiff) {
    rev_wepl -= val;
    out_vec.push_back(rev_wepl * step_length);
  }

  return out_vec;
}

std::vector<WEPLVector>
WEPLContourFromRtssContour(const Rtss_contour_modern &rt_contour,
                           const std::array<double, 3> &vec_basis,
                           const FloatImageType::Pointer &wepl_cube) {

  const std::array<double, 3> pixel_size = {{wepl_cube->GetSpacing()[0],
                                             wepl_cube->GetSpacing()[1],
                                             wepl_cube->GetSpacing()[2]}};
  const auto img_size = wepl_cube->GetLargestPossibleRegion().GetSize();

  std::vector<WEPLVector> WEPL_contour;

  for (auto &point : rt_contour.coordinates) {
    FloatImageType::PointType p(&point[0]);
    FloatImageType::IndexType cur_idx{};
    if (!wepl_cube->TransformPhysicalPointToIndex(p, cur_idx)) {
      const FloatImageType::IndexType first_idx = {0, 0, 0};
      FloatImageType::IndexType last_idx;
      last_idx.SetElement(0, img_size[0]);
      last_idx.SetElement(1, img_size[1]);
      last_idx.SetElement(2, img_size[2]);

      FloatImageType::PointType first_p;
      FloatImageType::PointType last_p;
      wepl_cube->TransformIndexToPhysicalPoint(first_idx, first_p);
      wepl_cube->TransformIndexToPhysicalPoint(last_idx, last_p);

      std::cerr << "Physical point {" << p.GetElement(0) << ", "
                << p.GetElement(1) << ", " << p.GetElement(2)
                << "} out of bounds: {" << first_p.GetElement(0) << ", "
                << first_p.GetElement(1) << ", " << first_p.GetElement(2)
                << "} -> {" << last_p.GetElement(0) << ", "
                << last_p.GetElement(1) << ", " << last_p.GetElement(2) << "} "
                << "in WEPL calculation"
                << "\n";
      return WEPL_contour;
    }
    const std::array<size_t, 3> point_id = {{static_cast<size_t>(cur_idx[0]),
                                             static_cast<size_t>(cur_idx[1]),
                                             static_cast<size_t>(cur_idx[2])}};
    const auto wepl =
        WEPL_from_point(point_id, vec_basis, pixel_size, wepl_cube);
    WEPL_contour.push_back({wepl, point});
  }

  return WEPL_contour;
}

std::vector<FloatVector>
distal_points_only(const Rtss_contour_modern &points,
                   const std::array<double, 3> &direction) {
  const auto eps_scale = 1.0;
  const auto eps_dir =
      FloatVector{static_cast<float>(eps_scale * std::get<0>(direction)),
                  static_cast<float>(eps_scale * std::get<1>(direction)),
                  static_cast<float>(eps_scale * std::get<2>(direction))};

  auto out = std::vector<FloatVector>();
  std::copy_if(points.coordinates.begin(), points.coordinates.end(),
               std::back_inserter(out),
               [eps_dir, &points](const FloatVector point) {
                 return points.is_distal(point, points.get_centre(), eps_dir);
                 // return points.is_inside(point - eps_dir);
               });

  return out;
}

std::vector<WEPLVector>
DistalWEPLContourFromRtssContour(const Rtss_contour_modern &rt_contour,
                                 const std::array<double, 3> &vec_basis,
                                 const FloatImageType::Pointer &wepl_cube) {
  auto tmp_rt_contour = rt_contour;
  tmp_rt_contour.coordinates = distal_points_only(rt_contour, vec_basis);
  return WEPLContourFromRtssContour(tmp_rt_contour, vec_basis, wepl_cube);
}

FloatImageType::PointType
point_from_WEPL(const FloatVector &start_point /* Physical point */,
                const double fWEPL, const FloatVector &vec_basis,
                const FloatImageType::Pointer &wepl_cube) {

  using precision = float;
  auto interpolator =
      itk::LinearInterpolateImageFunction<FloatImageType, precision>::New();
  interpolator->SetInputImage(wepl_cube);

  const auto start_point_phys = FloatImageType::PointType(&start_point[0]);
  auto start_idx = itk::ContinuousIndex<double, 3>();
  if (!wepl_cube->TransformPhysicalPointToContinuousIndex(start_point_phys,
                                                          start_idx)) {
    std::cerr << "Start point {" << start_point[0] << ", " << start_point[1]
              << ", " << start_point[2]
              << "} for reverse WEPL calc was not in image!\n";
    return {};
  }

  const auto step_length = 0.1f;
  const auto step = vec_basis * step_length;

  const auto pixel_size = static_vec_cast<precision, double, 3>(
      wepl_cube->GetSpacing().GetVnlVector());
  const auto ones = FloatVector(1.0);
  const auto inv_pixel_size = element_quotient(ones, pixel_size);

  auto point = element_product(
      static_vec_cast<precision, double, 3>(start_idx.GetVnlVector()),
      pixel_size);

  // Acumulate WEPL until fWEPL is reached
  auto accumWEPL = 0.0;
  while (accumWEPL < fWEPL) {
    const auto index_pos_vec = element_quotient(point, pixel_size);
    const auto index_pos =
        itk::ContinuousIndex<precision, 3>(&index_pos_vec[0]);

    if (!wepl_cube->GetBufferedRegion().IsInside(index_pos)) {
      const auto size = wepl_cube->GetBufferedRegion().GetSize();
      std::cerr << "Hit {" << index_pos.GetElement(0) << ", "
                << index_pos.GetElement(1) << ", " << index_pos.GetElement(2)
                << "} outside image boundary: {" << size[0] << ", " << size[1]
                << ", " << size[2] << "} on WEPL re-calc\n";
      break;
    }
    const auto val = interpolator->EvaluateAtContinuousIndex(index_pos);
    accumWEPL += val * step_length;

    // point = point + step (Reverse of WEPL_from_point)
    point += step;
  }

  const auto point_idx_vec = element_product(point, inv_pixel_size);
  const auto point_idx = itk::ContinuousIndex<precision, 3>(&point_idx_vec[0]);

  auto out_point = FloatImageType::PointType();

  wepl_cube->TransformContinuousIndexToPhysicalPoint(point_idx, out_point);

  return out_point;
}

FloatVector NewPoint_from_WEPLVector(const WEPLVector &vwepl,
                                     const std::array<double, 3> &arr_basis,
                                     const FloatImageType::Pointer &wepl_cube) {
  /* Find point of intersection with edge of wepl cube
   * Then step into wepl cube from point until vwepl.WEPL is reached.
   */

  // VNL should use SSEX.Y when available (?)
  using VectorType = vnl_vector_fixed<float, 3>;
  const auto vec_basis =
      VectorType(arr_basis.at(0), arr_basis.at(1), arr_basis.at(2));
  // Find direction of itk coordinate system vs. index coordinate system:
  auto first_idx = FloatImageType::IndexType();
  first_idx.SetElement(0, 0);
  first_idx.SetElement(1, 0);
  first_idx.SetElement(2, 0);

  auto first_phys_point = FloatImageType::PointType();
  wepl_cube->TransformIndexToPhysicalPoint(first_idx, first_phys_point);
  const auto sign_vec =
      IntVector{crl::ce_sgn(first_phys_point.GetElement(0) - vwepl.point[0]),
                crl::ce_sgn(first_phys_point.GetElement(1) - vwepl.point[1]),
                crl::ce_sgn(first_phys_point.GetElement(2) - vwepl.point[2])};
  // The sign of the first minus any point in cube should yield the sign
  // transformation we want on vec_basis

  const auto img_size = wepl_cube->GetLargestPossibleRegion().GetSize();
  const IntVector cubedim = {static_cast<int>(img_size[0]),
                             static_cast<int>(img_size[1]),
                             static_cast<int>(img_size[2])};

  /* Intersection of line with plane:
   * https://en.wikipedia.org/wiki/Line-plane_intersection
   * d = (p_0 - l_0) dot n / ( l dot n)
   * where p_0 is a point in the plane,
   * l_0 a point in the line, <- vwepl.point
   * n is normal vector to plane
   * l is a vector in direction of line <- vec_basis
   * => p = d * l + l_0
   */
  const VectorType l_0 = vwepl.point;
  const VectorType l(sign_vec.x * vec_basis.get(0),
                     sign_vec.y * vec_basis.get(1),
                     sign_vec.z * vec_basis.get(2));

  /* We will only need to check three planes
   * Depending on the sign of vec_basis:
   *    _________         __________
   *   /   y+   /|       /|        /|
   *  /_______ / |      /_|______ / |
   * |        |x+|     |x-|  z-  |  |
   * |   z+   |  |     |  |______|__|
   * |        | /      | /   y-  | /
   * |________|/       |/________|/
   */
  // Find p_0 and n of the three planes:
  // p_0 can be the corner where the three planes intersect:
  const FloatImageType::IndexType itk_p0_idx = {
      {l.get(0) < 0.0 ? 0 : cubedim.x, l.get(1) < 0.0 ? 0 : cubedim.y,
       l.get(2) < 0.0 ? 0 : cubedim.z}};

  FloatImageType::PointType itk_p0;
  wepl_cube->TransformIndexToPhysicalPoint(itk_p0_idx, itk_p0);
  const VectorType p_0 =
      static_vec_cast<float, double, 3>(itk_p0.GetVnlVector());

  std::array<VectorType, 3U> n = {{
      // Normal vectors to the three planes
      VectorType(1.0, 0.0, 0.0),
      VectorType(0.0, 1.0, 0.0),
      VectorType(0.0, 0.0, 1.0),
  }};

  // d = (p_0 - l_0) dot n / ( l dot n)
  const VectorType d(dot_product(p_0 - l_0, n.at(0)) / dot_product(l, n.at(0)),
                     dot_product(p_0 - l_0, n.at(1)) / dot_product(l, n.at(1)),
                     dot_product(p_0 - l_0, n.at(2)) / dot_product(l, n.at(2)));

  auto p = std::array<VectorType, 3U>{
      {d.get(0) * l + l_0, d.get(1) * l + l_0, d.get(2) * l + l_0}};

  // Now find the plane closest to l_0:
  auto p_dist = std::array<double, 3U>{
      {std::sqrt(dot_product(l_0 - p.at(0), l_0 - p.at(0))),
       std::sqrt(dot_product(l_0 - p.at(1), l_0 - p.at(1))),
       std::sqrt(dot_product(l_0 - p.at(2), l_0 - p.at(2)))}};
  const auto it_min_dist =
      std::min_element(std::begin(p_dist), std::end(p_dist));
  // index of min:
  const auto min_dist_plane = std::distance(std::begin(p_dist), it_min_dist);
  // Point of intersection:
  const auto p_intersect = p.at(min_dist_plane);

  const VectorType pixel_size =
      static_vec_cast<float, double, 3>(wepl_cube->GetSpacing().GetVnlVector());
  const auto start_point = p_intersect + element_product(pixel_size, vec_basis);

  const auto out_point =
      point_from_WEPL(start_point, vwepl.WEPL, vec_basis, wepl_cube);

  return static_vec_cast<float, double, 3>(out_point.GetVnlVector());
}

class hu_to_dedx_functor {
public:
  /*{ // plastimatch data:
          std::pair<float, double>(NLMIN(float), 0),
                  std::pair<float, double>(-1000, 0.00106),
                  std::pair<float, double>(0, 1.0),
                  std::pair<float, double>(41.46, 1.048674),
                  std::pair<float, double>(NLMAX(float), 0.005011) // wtf?
  };*/
  /*// Data from TRiP: 19990218.hlut
           std::pair<short, double>(std::numeric_limits<short>::min(), 0.0),
           std::pair<short, double>(-1000, 0.041),
           std::pair<short, double>(-798, 0.244),
           std::pair<short, double>(-750, 0.297),
           std::pair<short, double>(-108, 0.943),
           std::pair<short, double>(-75, 0.977),
           std::pair<short, double>(0, 1.0),
           std::pair<short, double>(40, 1.042),
           std::pair<short, double>(55, 1.049),
           std::pair<short, double>(67, 1.065),
           std::pair<short, double>(262, 1.095),
           std::pair<short, double>(1033, 1.468),
           std::pair<short, double>(1432, 1.634),
           std::pair<short, double>(1974, 1.778),
           std::pair<short, double>(3000, 2.051),
           std::pair<short, double>(std::numeric_limits<short>::max(),
                                    2.051)}
                                    */

  hu_to_dedx_functor() = default;
  ~hu_to_dedx_functor() = default;
  // Linear interpolator (as first class function)
  float operator()(short val) const {
    /**** Convert CT to dEdx ****/
    constexpr static const auto lookup =
        std::array<std::pair<short, double>, 16>{
            // DCPT Dual energy CT, "TB Large"
            {std::pair<short, double>(std::numeric_limits<short>::min(), 0.0),
             // std::pair<short, double>(-1024, 0.0f),
             std::pair<short, double>(-1000, 0.0),
             std::pair<short, double>(-850, 0.108),
             std::pair<short, double>(-266, 0.710),
             std::pair<short, double>(-166, 0.879),
             std::pair<short, double>(147, 1.117),
             std::pair<short, double>(172, 1.117),
             std::pair<short, double>(2500, 2.473),
             std::pair<short, double>(2950, 3.165),
             // std::pair<short, double>(3300, 3.165f),
             std::pair<short, double>(std::numeric_limits<short>::max(),
                                      3.165)}};

    const auto lookup_cond = [&val](const std::pair<short, double> cur_pair) {
      return val < cur_pair.first;
    };
    // Find first index in lookup that satisfies "val < lookup[i].first" :
    const auto lookup_upper_ptr =
        std::find_if(std::begin(lookup), std::end(lookup), lookup_cond);

    const auto lookup_upper = *lookup_upper_ptr;

    // Get the previous index:
    const auto lookup_lower = lookup.at(static_cast<size_t>(
        std::max(0, static_cast<int>(lookup_upper_ptr - lookup.begin() - 1))));

    // Do linear interpolation between upper and lower data point:
    const auto a = (lookup_upper.second - lookup_lower.second) /
                   static_cast<double>(lookup_upper.first - lookup_lower.first);
    const auto b =
        lookup_upper.second - a * static_cast<double>(lookup_upper.first);

    return static_cast<float>(a * static_cast<double>(val) + b);
  }
};

FloatImageType::Pointer
ConvertUshort2WeplFloat(const UShortImageType::Pointer &spImgUshort) {
  ShortImageType::Pointer hu_image_tmp;
  ConvertUshort2Short(spImgUshort, hu_image_tmp);

  auto hu_to_dedx_filter =
      itk::UnaryFunctorImageFilter<ShortImageType, FloatImageType,
                                   hu_to_dedx_functor>::New();
  hu_to_dedx_filter->SetInput(hu_image_tmp);
  hu_to_dedx_filter->Update();
  const auto wepl_image = hu_to_dedx_filter->GetOutput();

  return wepl_image;
}

} // namespace wepl

} // namespace crl
