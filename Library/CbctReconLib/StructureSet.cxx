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
  if (m_plan_ss == nullptr) {
    return nullptr;
  }
  switch (struct_set) {
  case ctType::PLAN_CT:
    m_plan_ss->wait();
    return m_plan_ss.get();
  case ctType::RIGID_CT:
    if (m_rigid_ss == nullptr) {
      std::cerr << "Rigid reg. structs not ready, falling back to plan CT!\n";
      m_plan_ss->wait();
      return m_plan_ss.get();
    }
    m_rigid_ss->wait();
    return m_rigid_ss.get();
  case ctType::DEFORM_CT:
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
  if (!ss) {
    out_ss = nullptr;
    return;
  }
  out_ss = std::make_unique<Rtss_modern>(*ss);
  out_ss->ready = false;

  const auto trn_by_vec = [&out_ss](const FloatVector &vec) {
    for (auto &roi : out_ss->slist) {
      for (auto &contour : roi.pslist) {
        for (auto &coord : contour.coordinates) { // should SIMD
          coord += vec;
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
  if (!ss) {
    out_ss = nullptr;
    return;
  }
  out_ss = std::make_unique<Rtss_modern>(*ss);
  out_ss->ready = false;

  const auto trn_by_vf = [&out_ss](const VectorFieldType::Pointer &vf) {
    for (auto &roi : out_ss->slist) {
      for (auto &contour : roi.pslist) {
        for (auto &coord : contour.coordinates) {
          VectorFieldType::PointType physIndex(&coord[0]);

          VectorFieldType::IndexType index{};
          if (!vf->TransformPhysicalPointToIndex(physIndex, index)) {
            std::cerr << "Index: " << index << " out of bounds: " << physIndex
                      << std::endl;
            continue;
          }

          coord += vf->GetPixel(index).GetVnlVector();
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
  if (!ss) {
    out_ss = nullptr;
    return;
  }
  out_ss = std::make_unique<Rtss_modern>(*ss);
  out_ss->ready = false;

  const auto trn_by_fun = [&out_ss, transform_function]() {
    for (auto &roi : out_ss->slist) {
      for (auto &contour : roi.pslist) {
        for (auto &coord : contour.coordinates) {
          const auto tmp_point = coord;

          const auto new_point = transform_function(tmp_point);

          coord = new_point;
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
  case XFORM_ITK_TRANSLATION: {
    const auto &trnsl = xform->get_trn();
    const auto offset = -(trnsl->GetOffset().GetVnlVector());
    transform = [offset](const FloatVector &point) {
      return point + static_vec_cast<float, double, 3>(offset);
    };
    break;
  }
  case XFORM_ITK_VERSOR: {
    const auto &versor = xform->get_vrs();
    const auto invtransform =
        itk::MatrixOffsetTransformBase<double, 3, 3>::New();
    if (!versor->GetInverse(invtransform)) {
      std::cerr << "Could not get inverse versor transform structures won't be "
                   "registered\n";
    }
    const auto invmatrix = invtransform->GetMatrix().GetVnlMatrix();
    const auto invoffset = invtransform->GetOffset().GetVnlVector();
    transform = [invmatrix, invoffset](const FloatVector &point) {
      return static_vec_cast<float, double, 3>(
          invmatrix * static_vec_cast<double, float, 3>(point) + invoffset);
    };
    break;
  }
  case XFORM_ITK_QUATERNION: {
    const auto &quarternion = xform->get_quat();
    const auto invtransform =
        itk::MatrixOffsetTransformBase<double, 3, 3>::New();
    if (!quarternion->GetInverse(invtransform)) {
      std::cerr << "Could not get inverse versor transform structures won't be "
                   "registered\n";
    }
    const auto invmatrix = invtransform->GetMatrix().GetVnlMatrix();
    const auto invoffset = invtransform->GetOffset().GetVnlVector();
    transform = [invmatrix, invoffset](const FloatVector &point) {
      return static_vec_cast<float, double, 3>(
          invmatrix * static_vec_cast<double, float, 3>(point) + invoffset);
    };
    break;
  }
  case XFORM_ITK_AFFINE: {
    const auto &affine = xform->get_aff();
    const auto invtransform = itk::AffineTransform<double, 3>::New();
    if (!affine->GetInverse(invtransform)) {
      std::cerr << "Could not get inverse affine transform structures won't be "
                   "registered\n";
    }
    const auto invmatrix = invtransform->GetMatrix().GetVnlMatrix();
    const auto invoffset = invtransform->GetOffset().GetVnlVector();
    transform = [invmatrix, invoffset](const FloatVector &point) {
      return static_vec_cast<float, double, 3>(
          invmatrix * static_vec_cast<double, float, 3>(point) + invoffset);
    };
    break;
  }
  case XFORM_ITK_BSPLINE: {
    const auto bspline = xform->get_itk_bsp();
    transform = [bspline](const FloatVector &point) {
      return static_vec_cast<float, double, 3>(
          bspline->TransformPoint(&point[0]).GetVnlVector());
    };
    break;
  }
  case XFORM_ITK_TPS: {
    const auto tps = xform->get_itk_tps();
    transform = [tps](const FloatVector &point) {
      return static_vec_cast<float, double, 3>(
          tps->TransformPoint(&point[0]).GetVnlVector());
    };
    break;
  }
  case XFORM_ITK_SIMILARITY: {
    const auto similarity = xform->get_similarity();
    transform = [similarity](const FloatVector &point) {
      return static_vec_cast<float, double, 3>(
          similarity->TransformPoint(&point[0]).GetVnlVector());
    };
    break;
  }
  case XFORM_ITK_VECTOR_FIELD: {
    vf = xform->get_itk_vf();
    break;
  }
  case XFORM_GPUIT_BSPLINE: {
    auto plm_header = xform->get_plm_image_header();
    const auto xform_itk = xform_to_itk_bsp(xform, &plm_header, nullptr);
    const auto bsp_transform = xform_itk->get_itk_bsp();
    bsp_transform->GetValidRegion().Print(std::cerr);
    transform = [bsp_transform](const FloatVector &point) {
      return static_vec_cast<float, double, 3>(
          bsp_transform->TransformPoint(&point[0]).GetVnlVector());
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
