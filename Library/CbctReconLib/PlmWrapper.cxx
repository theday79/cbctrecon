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
  int i;
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

std::unique_ptr<Rtss_roi_modern>
Rtss_modern::get_roi_by_name(const std::string &name) {
  for (auto &roi : slist) {
    if (roi.name == name) {
      return std::make_unique<Rtss_roi_modern>(roi);
    }
  }
  std::cerr << "VOI name was not in structure set !?" << std::endl;
  return nullptr;
}

Rtss_roi_modern &Rtss_modern::get_roi_ref_by_name(const std::string &name) {
  for (auto &roi : slist) {
    if (roi.name == name) {
      return roi;
    }
  }
  std::cerr << "VOI name was not in structure set !?\n";
  std::cerr << "Warning: Returning " << slist.at(0).name << " instead of "
            << name << "\n";
  return slist.at(0);
}
