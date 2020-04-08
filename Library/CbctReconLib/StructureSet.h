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

  void set_planCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);
  void set_rigidCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);
  void set_deformCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);

  Rtss_modern *get_ss(ctType struct_set) const;

  template <ctType CT_TYPE> constexpr Rtss_modern *get_ss() {
    if (m_plan_ss == nullptr) {
      return nullptr;
    }
    if constexpr (CT_TYPE == ctType::PLAN_CT) {
      m_plan_ss->wait();
      return m_plan_ss.get();
    }
    if constexpr (CT_TYPE == ctType::RIGID_CT) {
      if (m_rigid_ss == nullptr) {
        std::cerr << "Rigid reg. structs not ready, falling back to plan CT!\n";
        m_plan_ss->wait();
        return m_plan_ss.get();
      }
      m_rigid_ss->wait();
      return m_rigid_ss.get();
    }
    if constexpr (CT_TYPE == ctType::DEFORM_CT) {
      m_deform_ss->wait();
      return m_deform_ss.get();
    }
    std::cerr << "Invalid CT type" << std::endl;
    return nullptr;
  }

  template <ctType CT_TYPE> bool is_ss_null() const {
    bool is_null = true;

    if constexpr (CT_TYPE == ctType::PLAN_CT) {
      is_null = this->m_plan_ss == nullptr;
    }
    else if constexpr (CT_TYPE == ctType::RIGID_CT) {
      is_null = this->m_rigid_ss == nullptr;
    }
    else if constexpr (CT_TYPE == ctType::DEFORM_CT) {
      is_null = this->m_deform_ss == nullptr;
    }

    return is_null;
  }

  template <ctType CT_TYPE>
  void ApplyVectorTransformTo(const FloatVector &vec) {
    transform_by_vector(CT_TYPE, vec, m_rigid_ss);
  }

  template <ctType CT_TYPE>
  void ApplyVectorTransform_InPlace(const FloatVector &vec) {
    assert(!is_ss_null<CT_TYPE>());
    transform_by_vector(CT_TYPE, vec, get_ss<CT_TYPE>());
  }

  template <ctType CT_TYPE> bool ApplyTransformTo(const fs::path &transform_file) {

    if (is_ss_null<CT_TYPE>()) {
      return false;
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
      default:
      //case ctType::RIGID_CT:
      //case ctType::DEFORM_CT:
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
