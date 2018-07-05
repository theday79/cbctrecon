#ifndef STRUCTURESET_H
#define STRUCTURESET_H
#include "PlmWrapper.h"
#include "cbctrecon.h"

#include <bspline_xform.h>
#include <plm_image.h>
#include <xform.h>
#include <xform_convert.h>

enum ctType {
  PLAN_CT = 0,
  RIGID_CT = 1,
  DEFORM_CT = 2,
};

class StructureSet {
public:
  StructureSet();
  ~StructureSet();

  void set_planCT_ss(std::unique_ptr<Rtss> struct_set);
  void set_rigidCT_ss(std::unique_ptr<Rtss> struct_set);
  void set_deformCT_ss(std::unique_ptr<Rtss> struct_set);
  void set_planCT_ss(Rtss *struct_set);
  void set_rigidCT_ss(Rtss *struct_set);
  void set_deformCT_ss(Rtss *struct_set);

  std::unique_ptr<Rtss_modern> transform_by_vector(ctType struct_set,
                                                   FloatVector vec);

  std::unique_ptr<Rtss_modern>
  transform_by_vectorField(ctType struct_set,
                           const VectorFieldType::Pointer &vf);

  std::unique_ptr<Rtss_modern>
  transform_by_Lambda(ctType struct_set,
                      const TransformType &transform_function);

  Rtss_modern *get_ss(ctType struct_set);

  bool ApplyRigidTransformToPlan(QFile rigid_transform_file);
  bool ApplyDeformTransformToRigid(QFile deform_transform_file);

private:
  std::unique_ptr<Rtss_modern> m_plan_ss;
  std::unique_ptr<Rtss_modern> m_rigid_ss;
  std::unique_ptr<Rtss_modern> m_deform_ss;
};

#endif
