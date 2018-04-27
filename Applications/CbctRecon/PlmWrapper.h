#ifndef PLMWRAPPER_H
#define PLMWRAPPER_H

#include <array>
#include <functional>
#include <memory>
#include <vector>

#include <direction_cosines.h>
#include <rtss.h>
#include <rtss_contour.h>
#include <rtss_roi.h>

using ItkVectorType = itk::Vector<float, 3U>;
using VectorFieldType = itk::Image<ItkVectorType, 3U>;
using PointType = itk::Point<double, 3U>;
// Sorry, I can't control myself, I just love std::function
using TransformType = std::function<PointType(PointType)>;

struct FloatVector {
  float x;
  float y;
  float z;
};

class Plm_image_friend : public Plm_image {
public:
  VectorFieldType::Pointer friend_convert_to_itk(Volume *vol);
};

class Rtss_contour_modern : public Rtss_contour {
public:
  Rtss_contour_modern() = default;
  Rtss_contour_modern(const Rtss_contour *old_contour);
  Rtss_contour_modern(const Rtss_contour_modern *old_contour);
  std::vector<FloatVector> coordinates;
};

class Rtss_roi_modern : public Rtss_roi {
public:
  Rtss_roi_modern() = default;
  Rtss_roi_modern(const Rtss_roi *old_roi);
  Rtss_roi_modern(const Rtss_roi_modern *old_roi);

  std::vector<Rtss_contour_modern> pslist;
};

class Rtss_modern : public Rtss {
public:
  Rtss_modern() = default;
  // Unique pointer, to make sure it's killed by its destructor after copy
  Rtss_modern(std::unique_ptr<Rtss> old_rtss);
  Rtss_modern(const Rtss *old_rtss);
  Rtss_modern(const Rtss_modern *old_rtss);

  /* Output geometry */
  std::array<plm_long, 3> m_dim{};
  std::array<float, 3> m_spacing{};
  std::array<float, 3> m_offset{};
  /* Rasterization geometry */
  std::array<plm_long, 3> rast_dim{};
  std::array<float, 3> rast_spacing{};
  std::array<float, 3> rast_offset{};
  std::unique_ptr<Direction_cosines> rast_dc;
  /* Structures */
  std::vector<Rtss_roi_modern> slist;
};

#endif
