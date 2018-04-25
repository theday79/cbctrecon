#include "StructureSet.h"
#include <itkImage.h>
#include <itkVector.h>

StructureSet::StructureSet() {
  m_plan_ss = new Rtss();
  m_rigid_ss = new Rtss();
  m_deform_ss = new Rtss();
}

StructureSet::~StructureSet() {
  delete m_plan_ss;
  delete m_rigid_ss;
  delete m_deform_ss;
}

void StructureSet::set_planCT_ss(Rtss::Pointer &struct_set) {
  m_plan_ss = struct_set.get();
}
void StructureSet::set_rigidCT_ss(Rtss::Pointer &struct_set) {
  m_rigid_ss = struct_set.get();
}
void StructureSet::set_deformCT_ss(Rtss::Pointer &struct_set) {
  m_deform_ss = struct_set.get();
}

// This function is very C, because Rtss is used.
void clone_rtss(const Rtss *orig_ss, Rtss *new_ss) {
  auto n_structs = orig_ss->num_structures;
  for (auto i = 0; i < n_structs; i++) {
    new_ss->add_structure( // reallocs slist
        orig_ss->slist[i]->name, orig_ss->slist[i]->color,
        orig_ss->slist[i]->id, orig_ss->slist[i]->bit);
    auto n_contours = orig_ss->slist[i]->num_contours;

    for (auto j = 0; j < n_contours; j++) {
      auto n_verts = orig_ss->slist[i]->pslist[j]->num_vertices;
      new_ss->slist[i]->add_polyline(n_verts); // mallocs x, y, z

      new_ss->slist[i]->pslist[j]->slice_no =
          orig_ss->slist[i]->pslist[j]->slice_no;
      new_ss->slist[i]->pslist[j]->ct_slice_uid =
          orig_ss->slist[i]->pslist[j]->ct_slice_uid;

      for (auto k = 0; k < n_verts; k++) {
        new_ss->slist[i]->pslist[j]->x[k] = orig_ss->slist[i]->pslist[j]->x[k];
        new_ss->slist[i]->pslist[j]->y[k] = orig_ss->slist[i]->pslist[j]->y[k];
        new_ss->slist[i]->pslist[j]->z[k] = orig_ss->slist[i]->pslist[j]->z[k];
      }
    }
  }
}

void StructureSet::transform_by_vector(const ctType struct_set,
                                       const FloatVector vec, Rtss *out_ss) {
  if (out_ss->num_structures != 0) {
    out_ss->clear();
  }

  switch (struct_set) {
  case PLAN_CT:
    clone_rtss(m_plan_ss, out_ss);
    break;
  case RIGID_CT:
    clone_rtss(m_rigid_ss, out_ss);
    break;
  case DEFORM_CT:
    clone_rtss(m_deform_ss, out_ss);
    break;
  }

  auto n_structs = m_plan_ss->num_structures;
#pragma omp parallel for
  for (auto i = 0; i < n_structs; i++) {
    Rtss_roi *cur_s = m_plan_ss->slist[i];
    auto n_contours = cur_s->num_contours;
    for (auto j = 0; j < n_contours; j++) {
      Rtss_contour *cur_c = cur_s->pslist[j];
      auto n_points = cur_c->num_vertices;
      for (auto k = 0; k < n_points; k++) {
        cur_c->x[k] += vec.x;
        cur_c->y[k] += vec.y;
        cur_c->z[k] += vec.z;
      }
    }
  }
}
