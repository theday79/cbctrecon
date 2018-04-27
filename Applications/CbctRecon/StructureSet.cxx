#include <memory>

#include "StructureSet.h"

StructureSet::StructureSet() = default;

StructureSet::~StructureSet() = default;

void StructureSet::set_planCT_ss(std::unique_ptr<Rtss> struct_set) {
  m_plan_ss = std::make_unique<Rtss_modern>(std::move(struct_set));
}
void StructureSet::set_rigidCT_ss(std::unique_ptr<Rtss> struct_set) {
  m_rigid_ss = std::make_unique<Rtss_modern>(std::move(struct_set));
}
void StructureSet::set_deformCT_ss(std::unique_ptr<Rtss> struct_set) {
  m_deform_ss = std::make_unique<Rtss_modern>(std::move(struct_set));
}

void StructureSet::set_planCT_ss(const Rtss *struct_set) {
  m_plan_ss = std::make_unique<Rtss_modern>(struct_set);
}
void StructureSet::set_rigidCT_ss(const Rtss *struct_set) {
  m_rigid_ss = std::make_unique<Rtss_modern>(struct_set);
}
void StructureSet::set_deformCT_ss(const Rtss *struct_set) {
  m_deform_ss = std::make_unique<Rtss_modern>(struct_set);
}

Rtss_modern *StructureSet::get_ss(ctType struct_set) {
  switch (struct_set) {
  case PLAN_CT:
    return m_plan_ss.get();
  case RIGID_CT:
    return m_rigid_ss.get();
  case DEFORM_CT:
    return m_deform_ss.get();
  }
  
  return nullptr;
}

std::unique_ptr<Rtss_modern>
StructureSet::transform_by_vector(const ctType struct_set,
                                  const FloatVector vec) {
  auto out_ss = std::make_unique<Rtss_modern>(get_ss(struct_set));

  for (auto roi : out_ss->slist) {
    for (auto contour : roi.pslist) {
      for (auto coord : contour.coordinates) { // should SIMD
        coord.x += vec.x;
        coord.y += vec.y;
        coord.z += vec.z;
      }
    }
  }
  return out_ss;
}

std::unique_ptr<Rtss_modern>
StructureSet::transform_by_vectorField(const ctType struct_set,
                                       const VectorFieldType::Pointer &vf) {

  auto out_ss = std::make_unique<Rtss_modern>(get_ss(struct_set));

  for (auto roi : out_ss->slist) {
    for (auto contour : roi.pslist) {
      for (auto coord : contour.coordinates) {
        VectorFieldType::PointType physIndex;
        physIndex[0] = coord.x;
        physIndex[1] = coord.y;
        physIndex[2] = coord.z;

        VectorFieldType::IndexType index{};
        vf->TransformPhysicalPointToIndex(physIndex, index);

        coord.x += vf->GetPixel(index)[0];
        coord.y += vf->GetPixel(index)[1];
        coord.z += vf->GetPixel(index)[2];
      }
    }
  }
  return out_ss;
}

std::unique_ptr<Rtss_modern>
StructureSet::transform_by_Lambda(const ctType struct_set,
                                  const TransformType &transform_function) {

  auto out_ss = std::make_unique<Rtss_modern>(get_ss(struct_set));

  for (auto roi : out_ss->slist) {
    for (auto contour : roi.pslist) {
      for (auto coord : contour.coordinates) {
        VectorFieldType::PointType physIndex;
        physIndex[0] = coord.x;
        physIndex[1] = coord.y;
        physIndex[2] = coord.z;

        PointType new_point = transform_function(physIndex);

        coord.x = new_point[0];
        coord.y = new_point[1];
        coord.z = new_point[2];
      }
    }
  }
  return out_ss;
}

bool StructureSet::ApplyRigidTransformToPlan(QFile rigid_transform_file) {
  auto xform = std::make_unique<Xform>();
  xform->load(rigid_transform_file.fileName().toStdString());

  auto xform_type = xform->get_type();
  // First we make sure, it was a rigid-transform.
  switch (xform_type) {
  case XFORM_ITK_TRANSLATION:
  case XFORM_ITK_VERSOR:
  case XFORM_ITK_QUATERNION:
  case XFORM_ITK_AFFINE:
  case XFORM_ITK_SIMILARITY:
    std::cout << "Warning, rotation (and scaling) not yet implemented"
              << std::endl;
    break;
  default:
    std::cerr << "\a"
              << "You are passing the wrong type of transform file to the "
                 "rigid transform"
              << std::endl;
    return false;
  }

  // Now we get the translation of the transform
  // params is equivalent to itk::Array<T>
  itk::TransformBase::ParametersType params;

  switch (xform_type) {
  case XFORM_ITK_TRANSLATION: {
    auto trn = xform->get_trn();
    params = trn->GetParameters();
    break;
  }
  case XFORM_ITK_VERSOR: {
    auto trn = xform->get_vrs();
    params = trn->GetParameters();
    break;
  }
  case XFORM_ITK_QUATERNION: {
    auto trn = xform->get_quat();
    params = trn->GetParameters();
    break;
  }
  case XFORM_ITK_AFFINE: {
    auto trn = xform->get_aff();
    params = trn->GetParameters();
    break;
  }
  case XFORM_ITK_SIMILARITY: {
    auto trn = xform->get_similarity();
    params = trn->GetParameters();
    break;
  }
  default:
    std::cerr << "we should never get here" << std::endl;
    return false;
  }

  FloatVector trn_vec{static_cast<float>(params.x()),
                      static_cast<float>(params.y()),
                      static_cast<float>(params.z())};

  m_rigid_ss = transform_by_vector(PLAN_CT, trn_vec);

  return true;
}

bool StructureSet::ApplyDeformTransformToRigid(QFile deform_transform_file) {
  auto xform = std::make_unique<Xform>();
  xform->load(deform_transform_file.fileName().toStdString());

  auto xform_type = xform->get_type();

  TransformType transform;
  VectorFieldType::Pointer vf;

  switch (xform_type) {
  case XFORM_ITK_TRANSLATION:
    transform = [&xform](itk::Point<double, 3U> point) {
      return xform->get_trn()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_VERSOR:
    transform = [&xform](itk::Point<double, 3U> point) {
      return xform->get_vrs()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_QUATERNION:
    transform = [&xform](itk::Point<double, 3U> point) {
      return xform->get_quat()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_AFFINE:
    transform = [&xform](itk::Point<double, 3U> point) {
      return xform->get_aff()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_BSPLINE:
    transform = [&xform](itk::Point<double, 3U> point) {
      return xform->get_itk_bsp()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_TPS:
    transform = [&xform](itk::Point<double, 3U> point) {
      return xform->get_itk_tps()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_SIMILARITY:
    transform = [&xform](itk::Point<double, 3U> point) {
      return xform->get_similarity()->TransformPoint(point);
    };
    break;
  case XFORM_ITK_VECTOR_FIELD: {
    vf = xform->get_itk_vf();
    break;
  }
  case XFORM_GPUIT_BSPLINE: {
    auto xform_converter = std::make_unique<Xform_convert>();
    xform_converter->set_input_xform(
        std::move(xform)); // we dont need xform anymore after this anyway
    xform_converter->m_xf_out_type = XFORM_ITK_VECTOR_FIELD;
    xform_converter->run();
    auto out_xform =
        std::unique_ptr<Xform>(xform_converter->get_output_xform().get());
    vf = out_xform->get_itk_vf();
    break;
  }
  case XFORM_GPUIT_VECTOR_FIELD: {
    auto plm_img = std::make_unique<Plm_image_friend>();
    Volume *vol = xform->get_gpuit_vf().get();
    // The below function is expensive compared to the other options
    // consider methods to avoid using this format in DlgRegistration.
    vf = plm_img->friend_convert_to_itk(vol);
    break;
  }
  default:
    std::cerr << "\a"
              << "You are passing the wrong type of transform file to the "
                 "deform transform"
              << std::endl;
    return false;
  }

  if (transform != nullptr) {
    m_deform_ss = transform_by_Lambda(RIGID_CT, transform);
    return true;
  }
  if (vf.IsNotNull()) {
    m_deform_ss = transform_by_vectorField(RIGID_CT, vf);
    return true;
  }

  std::cerr
      << "\a"
      << "Transform function were not created and no vector field were applied"
      << std::endl;
  return false;
}
