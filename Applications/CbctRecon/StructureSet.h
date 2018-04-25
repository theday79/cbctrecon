#ifndef STRUCTURESET_H
#define STRUCTURESET_H
#include "cbctrecon.h"

#include <rtss.h>
#include <rtss_contour.h>
#include <rtss_roi.h>

struct FloatVector {
  float x;
  float y;
  float z;
};

enum ctType {
  PLAN_CT = 0,
  RIGID_CT = 1,
  DEFORM_CT = 2,
};

class StructureSet {
public:
  StructureSet();
  ~StructureSet();

  void set_planCT_ss(Rtss::Pointer &struct_set);
  void set_rigidCT_ss(Rtss::Pointer &struct_set);
  void set_deformCT_ss(Rtss::Pointer &struct_set);
  void transform_by_vector(ctType struct_set, FloatVector vec, Rtss *out_ss);

private:
  Rtss *m_plan_ss;
  Rtss *m_rigid_ss;
  Rtss *m_deform_ss;
};

#endif
