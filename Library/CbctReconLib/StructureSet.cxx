// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include <memory>

#undef TIMEOUT
#include "plm_image_header.h"
#include "volume_header.h"
#include "xform.h"

#include "StructureSet.h"

#define USE_THREADING true

class Volume_header_private {
public:
  plm_long m_dim[3]{};
  float m_origin[3]{};
  float m_spacing[3]{};
  Direction_cosines m_direction_cosines;

  Volume_header_private() { m_direction_cosines.set_identity(); }
};

StructureSet::StructureSet() = default;

StructureSet::~StructureSet() = default;

void StructureSet::set_planCT_ss(std::unique_ptr<Rtss_modern> &&struct_set) {
  m_plan_ss = std::move(struct_set);
}
void StructureSet::set_rigidCT_ss(std::unique_ptr<Rtss_modern> &&struct_set) {
  m_rigid_ss = std::move(struct_set);
}
void StructureSet::set_deformCT_ss(std::unique_ptr<Rtss_modern> &&struct_set) {
  m_deform_ss = std::move(struct_set);
}

Rtss_modern *StructureSet::get_ss(const ctType struct_set) const {
  switch (struct_set) {
  case PLAN_CT:
    m_plan_ss->wait();
    return m_plan_ss.get();
  case RIGID_CT:
    m_rigid_ss->wait();
    return m_rigid_ss.get();
  case DEFORM_CT:
    m_deform_ss->wait();
    return m_deform_ss.get();
  }
  std::cerr << "Invalid CT type" << std::endl;
  return nullptr;
}

void StructureSet::transform_by_vector(
    const ctType struct_set, const FloatVector &vec,
    std::unique_ptr<Rtss_modern> &out_ss) const {
  const auto ss = get_ss(struct_set);
  out_ss = std::make_unique<Rtss_modern>(*ss);
  out_ss->ready = false;

  const auto trn_by_vec = [&out_ss](const FloatVector &vec) {
    for (auto &roi : out_ss->slist) {
      for (auto &contour : roi.pslist) {
        for (auto &coord : contour.coordinates) { // should SIMD
          coord.x += vec.x;
          coord.y += vec.y;
          coord.z += vec.z;
        }
      }
    }
  };

#if USE_THREADING
  out_ss->thread_obj = std::thread(trn_by_vec, vec);
#else
  trn_by_vec(vec);
  out_ss->ready = true;
#endif
}

void StructureSet::transform_by_vectorField(
    const ctType ct_type, const VectorFieldType::Pointer &vf,
    std::unique_ptr<Rtss_modern> &out_ss) const {

  const auto ss = get_ss(ct_type);
  out_ss = std::make_unique<Rtss_modern>(*ss);
  out_ss->ready = false;

  const auto trn_by_vf = [&out_ss](const VectorFieldType::Pointer &vf) {
    for (auto &roi : out_ss->slist) {
      for (auto &contour : roi.pslist) {
        for (auto &coord : contour.coordinates) {
          VectorFieldType::PointType physIndex;
          physIndex[0] = coord.x;
          physIndex[1] = coord.y;
          physIndex[2] = coord.z;

          VectorFieldType::IndexType index{};
          if (!vf->TransformPhysicalPointToIndex(physIndex, index)) {
            std::cerr << "Index: " << index << " out of bounds: " << physIndex
                      << std::endl;
            continue;
          }

          coord.x += vf->GetPixel(index)[0];
          coord.y += vf->GetPixel(index)[1];
          coord.z += vf->GetPixel(index)[2];
        }
      }
    }
  };
#if USE_THREADING
  out_ss->thread_obj = std::thread(trn_by_vf, vf);
#else
  trn_by_vf(vf);
  out_ss->ready = true;
#endif
}

void StructureSet::transform_by_Lambda(
    const ctType ct_type, const TransformType &transform_function,
    std::unique_ptr<Rtss_modern> &out_ss) const {

  const auto ss = get_ss(ct_type);
  out_ss = std::make_unique<Rtss_modern>(*ss);
  out_ss->ready = false;

  const auto trn_by_fun = [&out_ss, &transform_function]() {
    for (auto &roi : out_ss->slist) {
      for (auto &contour : roi.pslist) {
        for (auto &coord : contour.coordinates) {
          VectorFieldType::PointType physIndex;
          physIndex[0] = coord.x;
          physIndex[1] = coord.y;
          physIndex[2] = coord.z;

          auto new_point = transform_function(physIndex);

          coord.x = new_point[0];
          coord.y = new_point[1];
          coord.z = new_point[2];
        }
      }
    }
  };
#if USE_THREADING
  out_ss->thread_obj = std::thread(trn_by_fun);
#else
  trn_by_fun();
  out_ss->ready = true;
#endif
}

std::pair<TransformType, VectorFieldType::Pointer>
StructureSet::get_transform_function(const Xform::Pointer &xform) const {
  const auto xform_type = xform->get_type();

  TransformType transform;
  VectorFieldType::Pointer vf;

  switch (xform_type) {
  case XFORM_ITK_TRANSLATION:
    transform = [&xform](const itk::Point<double, 3U> point) {
      return xform->get_trn()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_VERSOR:
    transform = [&xform](const itk::Point<double, 3U> point) {
      return xform->get_vrs()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_QUATERNION:
    transform = [&xform](const itk::Point<double, 3U> point) {
      return xform->get_quat()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_AFFINE:
    transform = [&xform](const itk::Point<double, 3U> point) {
      return xform->get_aff()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_BSPLINE:
    transform = [&xform](const itk::Point<double, 3U> point) {
      return xform->get_itk_bsp()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_TPS:
    transform = [&xform](const itk::Point<double, 3U> point) {
      return xform->get_itk_tps()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_SIMILARITY:
    transform = [&xform](const itk::Point<double, 3U> point) {
      return xform->get_similarity()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_VECTOR_FIELD: {
    vf = xform->get_itk_vf();
    break;
  }
  case XFORM_GPUIT_BSPLINE: {
    auto plm_header = xform->get_plm_image_header();
    plm_header.print();
    const auto xform_itk = xform_to_itk_bsp(xform, &plm_header, nullptr);
    const auto bsp_transform = xform_itk->get_itk_bsp();
    bsp_transform->GetValidRegion().Print(std::cerr);
    transform = [bsp_transform](const itk::Point<double, 3U> point) {
      return bsp_transform->TransformPoint(point);
    };
    break;
  }
  case XFORM_GPUIT_VECTOR_FIELD: {
    auto plm_img = std::make_unique<Plm_image_friend>();
    const auto vol = xform->get_gpuit_vf().get();
    // The below function is expensive compared to the other options
    // consider methods to avoid using this format in DlgRegistration.
    vf = plm_img->friend_convert_to_itk(vol);
    break;
  }
  case XFORM_NONE:
    break;
  }

  return std::make_pair(transform, vf);
}
