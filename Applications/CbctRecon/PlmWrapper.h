#ifndef PLMWRAPPER_H
#define PLMWRAPPER_H

#include <array>
#include <memory>

#include "plm_image.h"
#include "rtss.h"
#include "rtss_contour.h"
#include "rtss_roi.h"

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

class Plm_image_friend : public Plm_image {
public:
  VectorFieldType::Pointer friend_convert_to_itk(Volume *vol);
};

class CBCTRECON_API Rtss_contour_modern : public Rtss_contour {
public:
  Rtss_contour_modern() = default;
  ~Rtss_contour_modern() = default;
  Rtss_contour_modern(const Rtss_contour *old_contour);
  Rtss_contour_modern(const Rtss_contour_modern &old_contour);
  std::vector<FloatVector> coordinates;
};

class CBCTRECON_API Rtss_roi_modern : public Rtss_roi {
public:
  Rtss_roi_modern() = default;
  ~Rtss_roi_modern() = default;
  Rtss_roi_modern(const Rtss_roi *old_roi);
  Rtss_roi_modern(const Rtss_roi_modern &old_roi);

  std::vector<Rtss_contour_modern> pslist;
};

class CBCTRECON_API Rtss_modern : public Rtss {
public:
  Rtss_modern() = default;
  ~Rtss_modern() = default;
  // Unique pointer, to make sure it's killed by its destructor after copy
  Rtss_modern(std::unique_ptr<Rtss> old_rtss);
  Rtss_modern(const Rtss *old_rtss);
  Rtss_modern(const Rtss_modern &old_rtss);

  std::unique_ptr<Rtss_roi_modern> get_roi_by_name(std::string &name);
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
