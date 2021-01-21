#ifndef ITK_SYCL_WRAPPER_HPP
#define ITK_SYCL_WRAPPER_HPP

#include "itkImage.h"

namespace crl {
namespace sycl {
namespace util {

template <typename T, size_t D>
size_t get_total_size_of_itk_image(typename itk::Image<T, D>::Pointer &img) {
  auto img_size = img->GetLargestPossibleRegion().size();
  size_t total_size = img_size[0];
  for (auto i = 1; i < D; ++i) {
    total_size *= img_size[i];
  }
  return total_size;
}

template <typename T, size_t D> class itk_image_buffer {
  using ItkImageType = itk::Image<T, D>;

private:
  T *m_data = nullptr;
  size_t m_size = 0;

public:
  explicit itk_image_buffer(typename ItkImageType::Pointer &img)
      : m_data{img->GetBufferPointer()}, m_size{get_total_size_of_itk_image(
                                             img)} {}
};

} // namespace util
} // namespace sycl
} // namespace crl

#endif // ITK_SYCL_WRAPPER_HPP
