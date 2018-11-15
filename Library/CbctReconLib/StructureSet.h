#ifndef STRUCTURESET_H
#define STRUCTURESET_H

#include "PlmWrapper.h"

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

class QFile;

class CBCTRECON_API StructureSet {
public:
  StructureSet();
  ~StructureSet();
  StructureSet(const StructureSet &) =
      delete; // Maybe these should be implemented?
  void operator=(const StructureSet &) = delete;
  StructureSet(StructureSet &&) = delete;
  void operator=(StructureSet &&) = delete;

  void set_planCT_ss(std::unique_ptr<Rtss> struct_set);
  void set_planCT_ss(std::unique_ptr<Rtss_modern> &&struct_set);
  void set_rigidCT_ss(std::unique_ptr<Rtss> struct_set);
  void set_deformCT_ss(std::unique_ptr<Rtss> struct_set);
  void set_planCT_ss(Rtss *struct_set);
  void set_rigidCT_ss(Rtss *struct_set);
  void set_deformCT_ss(Rtss *struct_set);

  std::unique_ptr<Rtss_modern> transform_by_vector(ctType struct_set,
                                                   FloatVector vec) const;

  std::unique_ptr<Rtss_modern>
  transform_by_vectorField(ctType struct_set,
                           const VectorFieldType::Pointer &vf) const;

  std::unique_ptr<Rtss_modern>
  transform_by_Lambda(ctType struct_set,
                      const TransformType &transform_function) const;

  Rtss_modern *get_ss(ctType struct_set) const;

  bool ApplyRigidTransformToPlan(QFile rigid_transform_file);
  bool ApplyDeformTransformToRigid(QFile deform_transform_file);

private:
  std::unique_ptr<Rtss_modern> m_plan_ss;
  std::unique_ptr<Rtss_modern> m_rigid_ss;
  std::unique_ptr<Rtss_modern> m_deform_ss;
};

#endif
