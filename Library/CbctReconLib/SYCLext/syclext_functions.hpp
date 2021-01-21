#ifndef SYCLEXT_FUNCTIONS_HPP
#define SYCLEXT_FUNCTIONS_HPP

// #include "itk_sycl_wrapper.hpp"

#include <CL/sycl.hpp>

#include "syclext_utils.hpp"

namespace crl {
namespace sycl {
namespace util {
template <typename T> class image_buffer {

private:
  T *m_data = nullptr;
  size_t m_size = 0;
};
} // namespace util

template <typename T> class add;

template <typename T, size_t D>
void add_const_inplace(util::image_buffer<T> &image, const T value) {
  using namespace cl::sycl;

  auto size = image.size();

  queue defaultQueue;

  buffer<T, 1> imageBuf(image.data(), range<1>(size));

  defaultQueue.submit([&](handler &cgh) {
    auto imageAcc = imageBuf.template get_access<access::mode::read_write>(cgh);

    cgh.parallel_for<add<T>>(range<1>(size),
                             [=](id<1> i) { imageAcc[i] += value; });
  });
}

} // namespace sycl
} // namespace crl

#endif // SYCLEXT_FUNCTIONS_HPP
