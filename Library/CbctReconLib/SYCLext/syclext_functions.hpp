#ifndef SYCLEXT_FUNCTIONS_HPP
#define SYCLEXT_FUNCTIONS_HPP

#include "itk_sycl_wrapper.hpp"

#include "syclext_utils.hpp"

namespace crl {
namespace sycl {

template <typename T> class add;

template <typename T, size_t D>
void add_const_inplace(util::itk_image_buffer<T, D> &image, const T value) {
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
