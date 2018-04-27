#include <memory>

#include "PlmWrapper.h"
#include <itkImageRegionIterator.h>

VectorFieldType::Pointer Plm_image_friend::friend_convert_to_itk(Volume *vol) {
  /*return this->convert_gpuit_to_itk<VectorFieldType::Pointer, float(*)[3]>(
    vol);*/

  using ImageType = VectorFieldType;
  int i, d1, d2;
  auto img = reinterpret_cast<FloatVector *>(vol->img);
  ImageType::SizeType sz;
  ImageType::IndexType st;
  ImageType::RegionType rg;
  ImageType::PointType og;
  ImageType::SpacingType sp;
  ImageType::DirectionType dc;

  /* Copy header & allocate data for itk */
  for (d1 = 0; d1 < 3; d1++) {
    st[d1] = 0;
    sz[d1] = vol->dim[d1];
    sp[d1] = vol->spacing[d1];
    og[d1] = vol->origin[d1];
    for (d2 = 0; d2 < 3; d2++) {
      dc[d1][d2] = vol->direction_cosines[d1 * 3 + d2];
    }
  }
  rg.SetSize(sz);
  rg.SetIndex(st);

  VectorFieldType::Pointer itk_img = ImageType::New();
  itk_img->SetRegions(rg);
  itk_img->SetOrigin(og);
  itk_img->SetSpacing(sp);
  itk_img->SetDirection(dc);

  itk_img->Allocate();

  /* Copy data into itk */
  using IteratorType = itk::ImageRegionIterator<ImageType>;
  IteratorType it(itk_img, rg);
  for (it.GoToBegin(), i = 0; !it.IsAtEnd(); ++it, ++i) {
    float vec[] = {img[i].x, img[i].y, img[i].z};
    it.Set(vec);
  }

  /* Free input volume */
  this->free_volume();

  /* Return the new image; caller will assign to correct member
  and set type */
  return itk_img;
}

Rtss_contour_modern::Rtss_contour_modern(const Rtss_contour *old_contour) {
  ct_slice_uid = old_contour->ct_slice_uid;
  slice_no = old_contour->slice_no;
  num_vertices = old_contour->num_vertices;

  int i = -1;
  std::generate_n(coordinates.begin(), num_vertices, [&i, &old_contour]() {
    i++;
    return FloatVector{old_contour->x[i], old_contour->y[i], old_contour->z[i]};
  });
};

Rtss_contour_modern::Rtss_contour_modern(
    const Rtss_contour_modern *old_contour) {
  ct_slice_uid = old_contour->ct_slice_uid;
  slice_no = old_contour->slice_no;
  num_vertices = old_contour->num_vertices;
  coordinates.resize(num_vertices);
  std::copy(old_contour->coordinates.begin(), old_contour->coordinates.end(),
            coordinates.begin());
};

Rtss_roi_modern::Rtss_roi_modern(const Rtss_roi *old_roi) {
  name = old_roi->name;
  color = old_roi->color;
  id = old_roi->id;   /* Used for import/export (must be >= 1) */
  bit = old_roi->bit; /* Used for ss-img (-1 for no bit) */
  num_contours = old_roi->num_contours;
  std::copy_n(&old_roi->pslist[0], num_contours, pslist.begin());
};

Rtss_roi_modern::Rtss_roi_modern(const Rtss_roi_modern *old_roi) {
  name = old_roi->name;
  color = old_roi->color;
  id = old_roi->id;   /* Used for import/export (must be >= 1) */
  bit = old_roi->bit; /* Used for ss-img (-1 for no bit) */
  num_contours = old_roi->num_contours;
  pslist.resize(num_contours);
  std::copy(old_roi->pslist.begin(), old_roi->pslist.end(), pslist.begin());
};

Rtss_modern::Rtss_modern(std::unique_ptr<Rtss> old_rtss) {
  m_dim = {{old_rtss->m_dim[0], old_rtss->m_dim[1], old_rtss->m_dim[2]}};
  m_spacing = {
      {old_rtss->m_spacing[0], old_rtss->m_spacing[1], old_rtss->m_spacing[2]}};
  m_offset = {
      {old_rtss->m_offset[0], old_rtss->m_offset[1], old_rtss->m_offset[2]}};
  rast_dim = {
      {old_rtss->rast_dim[0], old_rtss->rast_dim[1], old_rtss->rast_dim[2]}};
  rast_spacing = {{old_rtss->rast_spacing[0], old_rtss->rast_spacing[1],
                   old_rtss->rast_spacing[2]}};
  rast_offset = {{old_rtss->rast_offset[0], old_rtss->rast_offset[1],
                  old_rtss->rast_offset[2]}};

  have_geometry = old_rtss->have_geometry;

  auto unique_rast_dc =
      std::make_unique<Direction_cosines>(old_rtss->rast_dc.get_matrix());
  rast_dc = std::move(unique_rast_dc);

  num_structures = old_rtss->num_structures;

  std::copy_n(&old_rtss->slist[0], num_structures, slist.begin());
};

Rtss_modern::Rtss_modern(const Rtss *old_rtss) {
  m_dim = {{old_rtss->m_dim[0], old_rtss->m_dim[1], old_rtss->m_dim[2]}};
  m_spacing = {
      {old_rtss->m_spacing[0], old_rtss->m_spacing[1], old_rtss->m_spacing[2]}};
  m_offset = {
      {old_rtss->m_offset[0], old_rtss->m_offset[1], old_rtss->m_offset[2]}};
  rast_dim = {
      {old_rtss->rast_dim[0], old_rtss->rast_dim[1], old_rtss->rast_dim[2]}};
  rast_spacing = {{old_rtss->rast_spacing[0], old_rtss->rast_spacing[1],
                   old_rtss->rast_spacing[2]}};
  rast_offset = {{old_rtss->rast_offset[0], old_rtss->rast_offset[1],
                  old_rtss->rast_offset[2]}};

  have_geometry = old_rtss->have_geometry;

  auto unique_rast_dc =
      std::make_unique<Direction_cosines>(old_rtss->rast_dc.get_matrix());
  rast_dc = std::move(unique_rast_dc);

  num_structures = old_rtss->num_structures;

  std::copy_n(&old_rtss->slist[0], num_structures, slist.begin());
};

Rtss_modern::Rtss_modern(const Rtss_modern *old_rtss) {
  m_dim = old_rtss->m_dim;
  m_spacing = old_rtss->m_spacing;
  m_offset = old_rtss->m_offset;
  rast_dim = old_rtss->rast_dim;
  rast_spacing = old_rtss->rast_spacing;
  rast_offset = old_rtss->rast_offset;
  have_geometry = old_rtss->have_geometry;

  auto unique_rast_dc =
      std::make_unique<Direction_cosines>(old_rtss->rast_dc->get_matrix());
  rast_dc = std::move(unique_rast_dc);

  num_structures = old_rtss->num_structures;
  slist.resize(num_structures);
  std::copy(old_rtss->slist.begin(), old_rtss->slist.end(), slist.begin());
};
