// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include <memory>

#include "PlmWrapper.h"
#include <itkImageRegionIterator.h>

VectorFieldType::Pointer Plm_image_friend::friend_convert_to_itk(Volume *vol) {
  /*return this->convert_gpuit_to_itk<VectorFieldType::Pointer, float(*)[3]>(
    vol);*/

  using ImageType = VectorFieldType;
  const auto img = reinterpret_cast<FloatVector *>(vol->img);
  ImageType::SizeType sz;
  ImageType::IndexType st;
  ImageType::RegionType rg;
  ImageType::PointType og;
  ImageType::SpacingType sp;
  ImageType::DirectionType dc;

  /* Copy header & allocate data for itk */
  for (auto d1 = 0U; d1 < 3U; d1++) {
    st[d1] = 0U;
    sz[d1] = static_cast<ImageType::SizeValueType>(vol->dim[d1]);
    sp[d1] = static_cast<ImageType::SpacingValueType>(vol->spacing[d1]);
    og[d1] = static_cast<ImageType::PointValueType>(vol->origin[d1]);
    for (auto d2 = 0U; d2 < 3U; d2++) {
      dc[d1][d2] = static_cast<ImageType::DirectionType::ValueType>(
          vol->direction_cosines[d1 * 3 + d2]);
    }
  }
  rg.SetSize(sz);
  rg.SetIndex(st);

  auto itk_img = ImageType::New();
  itk_img->SetRegions(rg);
  itk_img->SetOrigin(og);
  itk_img->SetSpacing(sp);
  itk_img->SetDirection(dc);

  itk_img->Allocate();

  /* Copy data into itk */
  using IteratorType = itk::ImageRegionIterator<ImageType>;
  IteratorType it(itk_img, rg);
  size_t i;
  for (it.GoToBegin(), i = 0; !it.IsAtEnd(); ++it, ++i) {
    auto vec = itk::Vector<float, 3>();
    vec.SetElement(0, img[i].x);
    vec.SetElement(1, img[i].y);
    vec.SetElement(2, img[i].z);
    it.Set(vec);
  }

  /* Free input volume */
  this->free_volume();

  /* Return the new image; caller will assign to correct member
  and set type */
  return itk_img;
}

Rtss_modern::Rtss_modern(const Rtss_modern &old)
    : slist(old.slist), have_geometry(old.have_geometry),
      ready(old.ready) { /* thread_obj shouldn't need initialization */
}

std::unique_ptr<Rtss_roi_modern>
Rtss_modern::get_roi_by_name(const std::string &name) {
  if (!ready) {
    thread_obj.join();
    ready = true;
  }
  for (auto &roi : slist) {
    if (roi.name.find(name) != std::string::npos) {
      return std::make_unique<Rtss_roi_modern>(roi);
    }
  }
  std::cerr << "VOI name was not in structure set !?" << std::endl;
  return nullptr;
}

Rtss_roi_modern &Rtss_modern::get_roi_ref_by_name(const std::string &name) {
  if (!ready) {
    thread_obj.join();
    ready = true;
  }
  for (auto &roi : slist) {
    if (roi.name.find(name) != std::string::npos) {
      return roi;
    }
  }
  std::cerr << "VOI name was not in structure set !?\n";
  std::cerr << "Warning: Returning " << slist.at(0).name << " instead of "
            << name << "\n";
  return slist.at(0);
}

bool Rtss_modern::wait() {
  if (!this->ready || this->thread_obj.joinable()) {
    this->thread_obj.join();
    this->ready = true;
  }
  return this->ready;
}

// Copyright 2000 softSurfer, 2012 Dan Sunday
// This code may be freely used and modified for any purpose
// providing that this copyright notice is included with it.
// SoftSurfer makes no warranty for this code, and cannot be held
// liable for any real or imagined damage resulting from its use.
// Users of this code must verify correctness for their application.

// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
inline int cn_PnPoly(const FloatVector &P, const std::vector<FloatVector> &V) {
  const auto n = V.size();
  int cn = 0; // the  crossing number counter

  // loop through all edges of the polygon
  for (size_t i = 0; i < n; i++) {              // edge from V[i]  to V[i+1]
    if (((V[i].y <= P.y) && (V[i + 1].y > P.y)) // an upward crossing
        || ((V[i].y > P.y) && (V[i + 1].y <= P.y))) { // a downward crossing
      // compute  the actual edge-ray intersect x-coordinate
      float vt = (float)(P.y - V[i].y) / (V[i + 1].y - V[i].y);
      if (P.x < V[i].x + vt * (V[i + 1].x - V[i].x)) // P.x < intersect
        ++cn; // a valid crossing of y=P.y right of P.x
    }
  }
  return (cn & 1); // 0 if even (out), and 1 if  odd (in)
}

bool Rtss_contour_modern::is_inside(const FloatVector point) const {
  return cn_PnPoly(point, this->coordinates) == 1;
}
