#include <cmath>                          // for sqrt, pow, sin, cos
#include <algorithm>                      // for find_if, min_element, transform
#include <iostream>                       // for operator<<, endl, basic_ostream, cerr, ostream
#include <numeric>                        // for accumulate, adjacent_difference
#include <utility>                        // for pair
#include <valarray>                       // for valarray

#include <vnl_vector_fixed.h>             // for vnl_vector_fixed
#include <itkCastImageFilter.h>           // for CastImageFilter
#include <itkFixedArray.h>                // for FixedArray
#include <itkImage.h>                     // for Image<>::Pointer, Image, Image<>::IndexType, Image<>::PointType
#include <itkImageFileWriter.h>           // for ImageFileWriter
#include <itkMacro.h>                     // for CastImageFilter::New, ImageFileWriter::New

#include "WEPL.h"
#include "cbctrecon.h"                    // for CbctRecon
#include "plm_math.h"                     // for M_PI, NLMAX, NLMIN

double lin_interpolate(std::array<int, 3> point_id,
                       std::array<double, 3> point_id_pos, const int idx_2,
                       const int idy_2, const int idz_2,
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

  const double sum_weights =
      std::accumulate(weights.begin(), weights.end(), 0.0);

  std::transform(
      weights.begin(), weights.end(), weights.begin(),
      [&sum_weights](double val) -> double { return val / sum_weights; });

  double out_val = 0.0;
  for (int i = 0; i < 8; i++) {
    // convert point_id to cube_ids
    FloatImageType::IndexType cube_id{};
    cube_id[0] = point_id.at(0) + (idx_2 * (i % 2));       // x= 0,1,0,1,0,1,0,1
    cube_id[1] = point_id.at(1) + (idy_2 * ((i / 2) % 2)); // y= 0,0,1,1,0,0,1,1
    cube_id[2] = point_id.at(2) + (idz_2 * (i / 4));       // z= 0,0,0,0,1,1,1,1

    out_val += wepl_cube->GetPixel(cube_id) * weights.at(i);
  }
  return out_val;
}

double WEPL_from_point(const std::array<size_t, 3> cur_point_id,
                       const std::array<double, 3> vec_basis,
                       const std::array<double, 3> vec_cubesize,
                       const std::array<size_t, 3> cubedim,
                       const FloatImageType::Pointer &wepl_cube) {
  const double step_length = 0.1;
  const std::array<double, 3> step = {{vec_basis.at(0) * step_length,
                                       vec_basis.at(1) * step_length,
                                       vec_basis.at(2) * step_length}};

  const std::array<double, 3> inv_cubesize = {{1.0 / vec_cubesize.at(0),
                                               1.0 / vec_cubesize.at(1),
                                               1.0 / vec_cubesize.at(2)}};

  std::array<double, 3> point = {
      {static_cast<double>(cur_point_id.at(0)) * vec_cubesize.at(0),
       static_cast<double>(cur_point_id.at(1)) * vec_cubesize.at(1),
       static_cast<double>(cur_point_id.at(2)) * vec_cubesize.at(2)}};

  double out = 0.0;

  while (true) {
    // point_id = point / cube_size
    const std::array<int, 3> point_id = {
        {static_cast<int>(round(point.at(0) * inv_cubesize.at(0))),
         static_cast<int>(round(point.at(1) * inv_cubesize.at(1))),
         static_cast<int>(round(point.at(2) * inv_cubesize.at(2)))}};

    if (point_id.at(0) < 0.0 ||
        point_id.at(0) >= static_cast<int>(cubedim.at(0)) ||
        point_id.at(1) < 0.0 ||
        point_id.at(1) >= static_cast<int>(cubedim.at(1)) ||
        point_id.at(2) < 0.0 ||
        point_id.at(2) >= static_cast<int>(cubedim.at(2))) {
      break;
    }

    // get nearest neighbors:
    const std::array<double, 3> point_id_pos = {
        {point.at(0) * inv_cubesize.at(0), point.at(1) * inv_cubesize.at(1),
         point.at(2) * inv_cubesize.at(2)}};

    int idx_2 = -1;
    if (point_id.at(0) < (point_id_pos.at(0))) {
      idx_2 = 1;
    }

    int idy_2 = -1;
    if (point_id.at(1) < (point_id_pos.at(1))) {
      idx_2 = 1;
    }

    int idz_2 = -1;
    if (point_id.at(2) < (point_id_pos.at(2))) {
      idz_2 = 1;
    }

    if ((point_id.at(0) + idx_2) < 0.0 ||
        (point_id.at(0) + idx_2) >= static_cast<int>(cubedim.at(0)) ||
        (point_id.at(1) + idy_2) < 0.0 ||
        (point_id.at(1) + idy_2) >= static_cast<int>(cubedim.at(1)) ||
        (point_id.at(2) + idz_2) < 0.0 ||
        (point_id.at(2) + idz_2) >= static_cast<int>(cubedim.at(2))) {
      break;
    }

    out +=
        lin_interpolate(point_id, point_id_pos, idx_2, idy_2, idz_2, wepl_cube);

    // point = point - step
    point.at(0) -= step.at(0);
    point.at(1) -= step.at(1);
    point.at(2) -= step.at(2);
  }

  return out * step_length;
}

std::array<double, 3> get_basis_from_angles(double gantry, double couch) {
  gantry += 180.0;
  gantry *= M_PI / 180.0;
  couch *= M_PI / 180.0;

  std::array<double, 3> basis = {
      {sin(gantry) * cos(couch), -cos(gantry), sin(couch) * sin(gantry)}};
  return basis;
}

std::vector<double>
WEPL_trace_from_point(const std::array<size_t, 3> cur_point_id,
                      const std::array<double, 3> vec_basis,
                      const std::array<double, 3> vec_cubesize,
                      const std::array<size_t, 3> cubedim,
                      const FloatImageType::Pointer &wepl_cube) {

  std::vector<double> cumWEPL; // cumulative WEPL

  const double step_length = 0.1;
  const std::array<double, 3> step = {{vec_basis.at(0) * step_length,
                                       vec_basis.at(1) * step_length,
                                       vec_basis.at(2) * step_length}};

  const std::array<double, 3> inv_cubesize = {{1.0 / vec_cubesize.at(0),
                                               1.0 / vec_cubesize.at(1),
                                               1.0 / vec_cubesize.at(2)}};

  std::array<double, 3> point = {
      {static_cast<double>(cur_point_id.at(0)) * vec_cubesize.at(0),
       static_cast<double>(cur_point_id.at(1)) * vec_cubesize.at(1),
       static_cast<double>(cur_point_id.at(2)) * vec_cubesize.at(2)}};

  double out_wepl = 0.0;
  while (true) {
    // point_id = point / cube_size
    const std::array<int, 3> point_id = {
        {static_cast<int>(round(point.at(0) * inv_cubesize.at(0))),
         static_cast<int>(round(point.at(1) * inv_cubesize.at(1))),
         static_cast<int>(round(point.at(2) * inv_cubesize.at(2)))}};

    if (point_id.at(0) < 0.0 ||
        point_id.at(0) >= static_cast<int>(cubedim.at(0)) ||
        point_id.at(1) < 0.0 ||
        point_id.at(1) >= static_cast<int>(cubedim.at(1)) ||
        point_id.at(2) < 0.0 ||
        point_id.at(2) >= static_cast<int>(cubedim.at(2))) {
      break;
    }

    // get nearest neighbors:
    const std::array<double, 3> point_id_pos = {
        {point.at(0) * inv_cubesize.at(0), point.at(1) * inv_cubesize.at(1),
         point.at(2) * inv_cubesize.at(2)}};

    int idx_2 = -1;
    if (point_id.at(0) < (point_id_pos.at(0))) {
      idx_2 = 1;
    }

    int idy_2 = -1;
    if (point_id.at(1) < (point_id_pos.at(1))) {
      idx_2 = 1;
    }

    int idz_2 = -1;
    if (point_id.at(2) < (point_id_pos.at(2))) {
      idz_2 = 1;
    }

    if ((point_id.at(0) + idx_2) < 0.0 ||
        (point_id.at(0) + idx_2) >= static_cast<int>(cubedim.at(0)) ||
        (point_id.at(1) + idy_2) < 0.0 ||
        (point_id.at(1) + idy_2) >= static_cast<int>(cubedim.at(1)) ||
        (point_id.at(2) + idz_2) < 0.0 ||
        (point_id.at(2) + idz_2) >= static_cast<int>(cubedim.at(2))) {
      break;
    }

    out_wepl +=
        lin_interpolate(point_id, point_id_pos, idx_2, idy_2, idz_2, wepl_cube);

    cumWEPL.push_back(out_wepl);

    // point = point - step
    point.at(0) -= step.at(0);
    point.at(1) -= step.at(1);
    point.at(2) -= step.at(2);
  }
  std::valarray<double> vdiff(cumWEPL.size());
  std::adjacent_difference(cumWEPL.begin(), cumWEPL.end(), &vdiff[0]);

  std::vector<double> out_vec;
  double revWEPL = cumWEPL.at(cumWEPL.size() - 1);

  for (auto val : vdiff) {
    revWEPL -= val;
    out_vec.push_back(revWEPL * step_length);
  }

  return out_vec;
}

std::vector<WEPLVector>
WEPLContourFromRtssContour(Rtss_contour_modern rt_contour,
                           const std::array<double, 3> vec_basis,
                           const FloatImageType::Pointer &wepl_cube) {

  const std::array<double, 3> pixel_size = {{wepl_cube->GetSpacing()[0],
                                             wepl_cube->GetSpacing()[1],
                                             wepl_cube->GetSpacing()[2]}};

  const std::array<size_t, 3> cubedim = {
      {wepl_cube->GetLargestPossibleRegion().GetSize()[0],
       wepl_cube->GetLargestPossibleRegion().GetSize()[1],
       wepl_cube->GetLargestPossibleRegion().GetSize()[2]}};

  std::vector<WEPLVector> WEPL_contour;

  for (auto point : rt_contour.coordinates) {
    FloatImageType::PointType p;
    p.SetElement(0, point.x);
    p.SetElement(1, point.y);
    p.SetElement(2, point.z);
    FloatImageType::IndexType cur_idx{};
    wepl_cube->TransformPhysicalPointToIndex(p, cur_idx);
    const std::array<size_t, 3> point_id = {{static_cast<size_t>(cur_idx[0]),
                                             static_cast<size_t>(cur_idx[1]),
                                             static_cast<size_t>(cur_idx[2])}};
    double wepl =
        WEPL_from_point(point_id, vec_basis, pixel_size, cubedim, wepl_cube);
    WEPL_contour.push_back({wepl, point});
  }

  return WEPL_contour;
}

DoubleVector point_from_WEPL(const DoubleVector start_point, const double fWEPL,
                             const std::array<double, 3> vec_basis,
                             const FloatImageType::Pointer &wepl_cube) {

  const double step_length = 0.1;
  const DoubleVector step = {vec_basis.at(0) * step_length,
                             vec_basis.at(1) * step_length,
                             vec_basis.at(2) * step_length};
  const IntVector cubedim = {
      static_cast<int>(wepl_cube->GetLargestPossibleRegion().GetSize()[0]),
      static_cast<int>(wepl_cube->GetLargestPossibleRegion().GetSize()[1]),
      static_cast<int>(wepl_cube->GetLargestPossibleRegion().GetSize()[2])};

  const DoubleVector pixel_size = {wepl_cube->GetSpacing()[0],
                                   wepl_cube->GetSpacing()[1],
                                   wepl_cube->GetSpacing()[2]};

  const DoubleVector inv_pixel_size = {1.0 / pixel_size.x, 1.0 / pixel_size.y,
                                       1.0 / pixel_size.z};

  DoubleVector point = {start_point.x * pixel_size.x,
                        start_point.y * pixel_size.y,
                        start_point.z * pixel_size.z};

  // Acumulate WEPL until fWEPL is reached
  auto accumWEPL = 0.0;
  while (accumWEPL < fWEPL) {
    // point_id = point / cube_size
    const std::array<int, 3> point_id = {
        {static_cast<int>(round(point.x * inv_pixel_size.x)),
         static_cast<int>(round(point.y * inv_pixel_size.y)),
         static_cast<int>(round(point.z * inv_pixel_size.z))}};

    // get nearest neighbors:
    const std::array<double, 3> point_id_pos = {{point.x * inv_pixel_size.x,
                                                 point.y * inv_pixel_size.y,
                                                 point.z * inv_pixel_size.z}};
    int idx_2 = -1;
    if (point_id.at(0) < (point_id_pos.at(0))) {
      idx_2 = 1;
    }

    int idy_2 = -1;
    if (point_id.at(1) < (point_id_pos.at(1))) {
      idx_2 = 1;
    }

    int idz_2 = -1;
    if (point_id.at(2) < (point_id_pos.at(2))) {
      idz_2 = 1;
    }

    // Check we are still in cube:
    if ((point_id.at(0) + idx_2) < 0.0 ||
        (point_id.at(0) + idx_2) >= cubedim.x ||
        (point_id.at(1) + idy_2) < 0.0 ||
        (point_id.at(1) + idy_2) >= cubedim.y ||
        (point_id.at(2) + idz_2) < 0.0 ||
        (point_id.at(2) + idz_2) >= cubedim.z) {
      std::cerr << "Image boundary was reached on WEPL calc!" << std::endl;
      break;
    }

    accumWEPL +=
        lin_interpolate(point_id, point_id_pos, idx_2, idy_2, idz_2, wepl_cube);
    // point = point + step (Reverse of WEPL_from_point)
    point.x += step.x;
    point.y += step.y;
    point.z += step.z;
  }

  return point;
}

FloatVector NewPoint_from_WEPLVector(const WEPLVector vwepl,
                                     const std::array<double, 3> vec_basis,
                                     const FloatImageType::Pointer &wepl_cube) {
  /* Find point of intersection with edge of wepl cube
   * Then step into wepl cube from point until vwepl.WEPL is reached.
   */

  // VNL should use SSEX.Y when available (?)
  using VectorType = vnl_vector_fixed<double, 3>;

  const IntVector cubedim = {
      static_cast<int>(wepl_cube->GetLargestPossibleRegion().GetSize()[0]),
      static_cast<int>(wepl_cube->GetLargestPossibleRegion().GetSize()[1]),
      static_cast<int>(wepl_cube->GetLargestPossibleRegion().GetSize()[2])};

  /* Intersection of line with plane:
   * https://en.wikipedia.org/wiki/Line-plane_intersection
   * d = (p_0 - l_0) dot n / ( l dot n)
   * where p_0 is a point in the plane,
   * l_0 a point in the line, <- vwepl.point
   * n is normal vector to plane
   * l is a vector in direction of line <- vec_basis
   * => p = d * l + l_0
   */
  VectorType l_0(vwepl.point.x, vwepl.point.y, vwepl.point.z);
  VectorType l(vec_basis.at(0), vec_basis.at(1), vec_basis.at(2));

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
  FloatImageType::IndexType itk_p0_idx = {{0 ? l.get(0) < 0.0 : cubedim.x,
                                           0 ? l.get(1) < 0.0 : cubedim.y,
                                           0 ? l.get(2) < 0.0 : cubedim.z}};

  FloatImageType::PointType itk_p0;
  wepl_cube->TransformIndexToPhysicalPoint(itk_p0_idx, itk_p0);
  VectorType p_0 = itk_p0.Get_vnl_vector();

  std::array<VectorType, 3U> n = {{
      VectorType(sgn(vec_basis.at(0)), 0.0, 0.0),
      VectorType(0.0, sgn(vec_basis.at(1)), 0.0),
      VectorType(0.0, 0.0, sgn(vec_basis.at(3))),
  }};

  // d = (p_0 - l_0) dot n / ( l dot n)
  VectorType d(dot_product(p_0 - l, n.at(0)) / dot_product(l, n.at(0)),
               dot_product(p_0 - l, n.at(1)) / dot_product(l, n.at(1)),
               dot_product(p_0 - l, n.at(2)) / dot_product(l, n.at(2)));

  std::array<VectorType, 3U> p = {
      {d.get(0) * l + l_0, d.get(1) * l + l_0, d.get(2) * l + l_0}};

  // Now find the plane closest to l_0:
  std::array<double, 3U> p_dist = {
      {std::sqrt(dot_product(l_0 - p.at(0), l_0 - p.at(0))),
       std::sqrt(dot_product(l_0 - p.at(1), l_0 - p.at(1))),
       std::sqrt(dot_product(l_0 - p.at(2), l_0 - p.at(2)))}};
  auto it_min_dist = std::min_element(p_dist.begin(), p_dist.end());
  // index of min:
  auto min_dist_plane = std::distance(p_dist.begin(), it_min_dist);
  // Point of intersection:
  DoubleVector intersect = {p.at(min_dist_plane).get(0),
                            p.at(min_dist_plane).get(1),
                            p.at(min_dist_plane).get(2)};

  auto out_point = point_from_WEPL(intersect, vwepl.WEPL, vec_basis, wepl_cube);

  return FloatVector{static_cast<float>(out_point.x),
                     static_cast<float>(out_point.y),
                     static_cast<float>(out_point.z)};
}

/**** Convert CT to dEdx ****/

using double_pair_list = std::array<std::pair<float, double>, 16>;
const double_pair_list lookup = {
    {// Data from TRiP: 19990218.hlut
     std::pair<float, double>(NLMIN(float), 0),
     std::pair<float, double>(-1000.0f, 0.041),
     std::pair<float, double>(-798.0f, 0.244),
     std::pair<float, double>(-750.0f, 0.297),
     std::pair<float, double>(-108.0f, 0.943),
     std::pair<float, double>(-75.0f, 0.977),
     std::pair<float, double>(0.0f, 1.0),
     std::pair<float, double>(40.0f, 1.042),
     std::pair<float, double>(55.0f, 1.049),
     std::pair<float, double>(67.0f, 1.065),
     std::pair<float, double>(262.0f, 1.095),
     std::pair<float, double>(1033.0f, 1.468),
     std::pair<float, double>(1432.0f, 1.634),
     std::pair<float, double>(1974.0f, 1.778),
     std::pair<float, double>(3000.0f, 2.051),
     std::pair<float, double>(NLMAX(float), 2.051)}};
/*{ // plastimatch data:
        std::pair<float, double>(NLMIN(float), 0),
                std::pair<float, double>(-1000, 0.00106),
                std::pair<float, double>(0, 1.0),
                std::pair<float, double>(41.46, 1.048674),
                std::pair<float, double>(NLMAX(float), 0.005011) // wtf?
};*/

// Linear interpolator (as first class function)
float hu_to_dEdx(float val) {
  // Find first index in lookup that satisfies "val < lookup[i].first" :
  auto lookup_upper_ptr = std::find_if(
      lookup.begin(), lookup.end(), [val](std::pair<float, double> cur_pair) {
        return val < cur_pair.first;
      });
  const auto lookup_upper = *lookup_upper_ptr;

  // Get the previous index:
  const auto lookup_lower = lookup.at((lookup_upper_ptr - lookup.begin()) - 1);

  // Do linear interpolation between upper and lower data point:
  const double a = (lookup_upper.second - lookup_lower.second) /
                   static_cast<double>(lookup_upper.first - lookup_lower.first);
  const double b =
      lookup_upper.second - a * static_cast<double>(lookup_upper.first);

  return a * val + b;
};

FloatImageType::Pointer
ConvertUshort2WeplFloat(UShortImageType::Pointer &spImgUshort) {
  ShortImageType::Pointer hu_image_tmp;
  CbctRecon::ConvertUshort2Short(spImgUshort, hu_image_tmp);

  using CastFilterType = itk::CastImageFilter<ShortImageType, FloatImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(hu_image_tmp);
  castFilter->Update();

  FloatImageType::Pointer wepl_image = castFilter->GetOutput();

  itk::ImageRegionIterator<FloatImageType> it(
      wepl_image, wepl_image->GetLargestPossibleRegion());

  // std::transform(it.Begin(), it.End(), hu_to_dEdx); //ITK iterators doesn't
  // support <algorithm> (yet?)
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    it.Set(hu_to_dEdx(it.Get()));
  }

  using WriterType = itk::ImageFileWriter<FloatImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(wepl_image);
  writer->SetFileName("wepl_image.mha");
  writer->Update();

  return wepl_image;
}
