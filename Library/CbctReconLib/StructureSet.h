#ifndef STRUCTURESET_H
#define STRUCTURESET_H

#include <filesystem>

#undef TIMEOUT
#undef CUDA_FOUND
#include "xform.h"

#include "PlmWrapper.h"

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

namespace fs = std::filesystem;

class CBCTRECON_API StructureSet {
public:
  StructureSet();
  ~StructureSet();
  StructureSet(const StructureSet &) =
      delete; // Maybe these should be implemented?
  void operator=(const StructureSet &) = delete;
  StructureSet(StructureSet &&) = delete;
  void operator=(StructureSet &&) = delete;

  void set_planCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);
  void set_rigidCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);
  void set_deformCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);

  Rtss_modern *get_ss(ctType struct_set) const;

  template <ctType CT_TYPE> bool is_ss_null() const {
    switch (CT_TYPE) {
    case ctType::PLAN_CT:
      return m_plan_ss == nullptr;
    case ctType::RIGID_CT:
      return m_rigid_ss == nullptr;
    case ctType::DEFORM_CT:
      return m_deform_ss == nullptr;
    }
    return true;
  }

  template <ctType CT_TYPE>
  void ApplyVectorTransformTo(const FloatVector &vec) {
    transform_by_vector(CT_TYPE, vec, m_rigid_ss);
  }

  template <ctType CT_TYPE>
  void ApplyVectorTransform_InPlace(const FloatVector &vec) {
    switch (CT_TYPE) {
    case ctType::PLAN_CT:
      transform_by_vector(CT_TYPE, vec, m_plan_ss);
      break;
    case ctType::RIGID_CT:
      transform_by_vector(CT_TYPE, vec, m_rigid_ss);
      break;
    case ctType::DEFORM_CT:
      transform_by_vector(CT_TYPE, vec, m_deform_ss);
      break;
    }
  }

  template <ctType CT_TYPE> bool ApplyTransformTo(const fs::path &transform_file) {
    switch (CT_TYPE) {
    case ctType::PLAN_CT:
      if (m_plan_ss == nullptr) {
        return false;
      }
      break;
    case ctType::RIGID_CT:
      if (m_rigid_ss == nullptr) {
        return false;
      }
      break;
    case ctType::DEFORM_CT:
      if (m_deform_ss == nullptr) {
        return false;
      }
      break;
    }

    auto xform = Xform::New();
    xform->load(fs::absolute(transform_file).string());

    auto transform_pair = this->get_transform_function(xform);
    auto &transform = transform_pair.first;
    auto &vf = transform_pair.second;

    if (transform != nullptr) {
      switch (CT_TYPE) {
      case ctType::PLAN_CT:
        transform_by_Lambda(CT_TYPE, transform, m_rigid_ss);
        break;
      case ctType::RIGID_CT:
        [[fallthrough]];
      case ctType::DEFORM_CT:
        transform_by_Lambda(CT_TYPE, transform, m_deform_ss);
        break;
      }
      return true;
    }
    if (vf.IsNotNull()) {
      transform_by_vectorField(CT_TYPE, vf, m_deform_ss);
      return true;
    }

    std::cerr << "\a"
              << "Transform function were not created and no vector field were "
                 "applied\n";
    return false;
  }

private:
  std::pair<TransformType, VectorFieldType::Pointer>
  get_transform_function(const Xform::Pointer &xform) const;

  void transform_by_vector(ctType struct_set, const FloatVector &vec,
                           std::unique_ptr<Rtss_modern> &out_ss) const;

  void transform_by_Lambda(ctType ct_type,
                           const TransformType &transform_function,
                           std::unique_ptr<Rtss_modern> &out_ss) const;

  void transform_by_vectorField(ctType ct_type,
                                const VectorFieldType::Pointer &vf,
                                std::unique_ptr<Rtss_modern> &out_ss) const;

  std::unique_ptr<Rtss_modern> m_plan_ss;
  std::unique_ptr<Rtss_modern> m_rigid_ss;
  std::unique_ptr<Rtss_modern> m_deform_ss;
};

#endif
