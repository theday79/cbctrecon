#ifndef STRUCTURESET_H
#define STRUCTURESET_H

#if __has_include(<oneapi/dpl/execution>)
#include <oneapi/dpl/algorithm>
#include <oneapi/dpl/execution>
namespace execution = oneapi::dpl::execution;
#else
#include <algorithm>
#include <execution>
namespace execution = std::execution;
#endif

#include <filesystem>

#undef TIMEOUT
#undef CUDA_FOUND
#include "xform.h"

#include "PlmWrapper.h"

#include "WEPL.h"
#include "cbctrecon_config.h"
#include "cbctrecon_types.h"
#include "free_functions.h"

namespace fs = std::filesystem;

class CBCTRECON_API StructureSet {
public:
  void set_planCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);
  void set_rigidCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);
  void set_deformCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);

  Rtss_modern *get_ss(ctType struct_set) const;

  template <ctType CT_TYPE> constexpr auto &get_ss() {
    if (m_plan_ss == nullptr) {
      // return & unique nullptr:
      return m_plan_ss;
    }
    if constexpr (CT_TYPE == ctType::PLAN_CT) {
      m_plan_ss->wait();
      return m_plan_ss;
    }
    if constexpr (CT_TYPE == ctType::RIGID_CT) {
      if (m_rigid_ss == nullptr) {
        std::cerr << "Rigid reg. structs not ready, falling back to plan CT!\n";
        m_plan_ss->wait();
        return m_plan_ss;
      }
      m_rigid_ss->wait();
      return m_rigid_ss;
    }
    if constexpr (CT_TYPE == ctType::DEFORM_CT) {
      m_deform_ss->wait();
      return m_deform_ss;
    }
    static_assert(CT_TYPE == ctType::PLAN_CT || CT_TYPE == ctType::RIGID_CT ||
                      CT_TYPE == ctType::DEFORM_CT,
                  "Invalid CT type");
    return m_plan_ss;
  }

  template <ctType CT_TYPE> bool is_ss_null() const {
    bool is_null = true;

    if constexpr (CT_TYPE == ctType::PLAN_CT) {
      is_null = this->m_plan_ss == nullptr;
    } else if constexpr (CT_TYPE == ctType::RIGID_CT) {
      is_null = this->m_rigid_ss == nullptr;
    } else if constexpr (CT_TYPE == ctType::DEFORM_CT) {
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

  template <ctType CT_TYPE>
  bool ApplyTransformTo(const fs::path &transform_file) {

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
        // case ctType::RIGID_CT:
        // case ctType::DEFORM_CT:
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

namespace crl {

template <typename T>
auto roi_to_distal_only_roi(const Rtss_roi_modern &roi, const T gantry,
                            const T couch) {
  const auto direction = crl::wepl::get_basis_from_angles(gantry, couch);
  auto out_roi = Rtss_roi_modern(roi);

  for (auto &contour : out_roi.pslist) {
    contour.coordinates = crl::wepl::distal_points_only(contour, direction);
  }

  return out_roi;
}

template <typename T> struct hausdorff_result {
  T h_min = std::numeric_limits<T>::max();
  T h_max = std::numeric_limits<T>::max();
  T h_percent = std::numeric_limits<T>::max();
  std::vector<T> h_all;
};

template <typename T>
T min_distance(const FloatVector from_point, const Rtss_roi_modern &to_roi) {
  auto to_roi_distances = std::vector<T>();
  for (const auto &to_contour : to_roi.pslist) {
    auto &to_coords = to_contour.coordinates;
    if (to_coords.size() == 0) {
      continue;
    }
    auto to_contour_distances = std::vector<T>(to_coords.size());
    std::transform(execution::par_unseq, to_coords.begin(), to_coords.end(),
                   to_contour_distances.begin(),
                   [from_point](const auto to_point) {
                     return crl::ce_distance<T>(from_point, to_point);
                   });
    to_roi_distances.push_back(*std::min_element(to_contour_distances.begin(),
                                                 to_contour_distances.end()));
  }
  return *std::min_element(to_roi_distances.begin(), to_roi_distances.end());
}

template <typename T, unsigned char percent>
auto calculate_hausdorff(const Rtss_roi_modern &from_roi,
                         const Rtss_roi_modern &to_roi) {
  static_assert(
      percent <= 100,
      "Percent should be less than 100 as a higher value doesn't make sense.");

  const auto n_total_points = std::transform_reduce(
      execution::par_unseq, from_roi.pslist.begin(), from_roi.pslist.end(),
      static_cast<size_t>(0), std::plus<size_t>(),
      [](const Rtss_contour_modern &contour) {
        return contour.coordinates.size();
      });

  auto distances = std::vector<T>(n_total_points);
  auto dist_iterator = distances.begin();
  for (const auto &from_contour : from_roi.pslist) {
    auto &from_coords = from_contour.coordinates;
    dist_iterator = std::transform(execution::par_unseq, from_coords.begin(),
                                   from_coords.end(), dist_iterator,
                                   [&to_roi](const auto from_point) {
                                     return min_distance<T>(from_point, to_roi);
                                   });
  }

  std::sort(execution::par_unseq, distances.begin(), distances.end());

  hausdorff_result<T> hausdorff_output;
  hausdorff_output.h_min = distances.front();
  hausdorff_output.h_max = distances.back();

  const auto percent_index = static_cast<size_t>(
      crl::ce_round((percent / 100.0) * static_cast<double>(distances.size())));

  hausdorff_output.h_percent = distances.at(percent_index);
  hausdorff_output.h_all = distances;
  return hausdorff_output;
}

template <typename T, unsigned char percent>
auto calculate_hausdorff_and_top5(const Rtss_roi_modern &from_roi,
                                  const Rtss_roi_modern &to_roi) {
  const auto hausdorff_output =
      calculate_hausdorff<T, percent>(from_roi, to_roi);

  // Return also a roi for displaying the top 5% points
  auto top5 = std::make_unique<Rtss_roi_modern>();
  top5->color = "0/127/127";
  top5->name = "Top5%";
  for (const auto &from_contour : from_roi.pslist) {
    const auto &from_coords = from_contour.coordinates;
    auto top5_contour = Rtss_contour_modern();
    std::copy_if(from_coords.begin(), from_coords.end(),
                 std::back_inserter(top5_contour.coordinates),
                 [&to_roi, hausdorff_output](const auto &from_point) {
                   return min_distance<T>(from_point, to_roi) >
                          hausdorff_output.h_percent;
                 });
    top5->pslist.push_back(top5_contour);
  }

  return std::make_tuple(hausdorff_output, std::move(top5));
}
} // namespace crl

#endif
