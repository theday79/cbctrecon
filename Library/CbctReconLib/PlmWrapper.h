#ifndef PLMWRAPPER_H
#define PLMWRAPPER_H

#include <array>
#include <memory>

#undef TIMEOUT
#undef CUDA_FOUND
#include "plm_image.h"

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

class Plm_image_friend : public Plm_image {
public:
  VectorFieldType::Pointer friend_convert_to_itk(Volume *vol);
};

struct CBCTRECON_API Rtss_contour_modern { // : public Rtss_contour {
  std::vector<FloatVector> coordinates;
  /* Plastimatch specific */
  int slice_no = -1;
  std::string ct_slice_uid = "";
  size_t num_vertices = 0U;
};

struct CBCTRECON_API Rtss_roi_modern { // : public Rtss_roi {
  std::vector<Rtss_contour_modern> pslist;
  std::string name = "";
  std::string color = "255 0 0";
  /* Plastimatch specific */
  size_t id = 1; /* Used for import/export (must be >= 1) */
  int bit = -1;  /* Used for ss-img (-1 for no bit) */
  size_t num_contours = 0;
};

struct CBCTRECON_API Rtss_modern { // : public Rtss {
  Rtss_modern() = default;
  ~Rtss_modern() = default;
  Rtss_modern(const Rtss_modern &old);
  std::unique_ptr<Rtss_roi_modern> get_roi_by_name(const std::string &name);
  Rtss_roi_modern &get_roi_ref_by_name(const std::string &name);
  bool wait();
  /* Output geometry */
  std::array<plm_long, 3> m_dim{};
  std::array<float, 3> m_spacing{};
  std::array<float, 3> m_offset{};
  /* Rasterization geometry */
  std::array<plm_long, 3> rast_dim{};
  std::array<float, 3> rast_spacing{};
  std::array<float, 3> rast_offset{};
  Direction_cosines *rast_dc = nullptr;
  /* Structures */
  std::vector<Rtss_roi_modern> slist;
  /* Plastimatch specific */
  bool have_geometry = false;
  size_t num_structures = 0;

  bool ready = false;
  std::thread thread_obj;
};

#endif
