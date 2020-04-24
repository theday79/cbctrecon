#ifndef CBCTREGISTRATION_H
#define CBCTREGISTRATION_H

#include <filesystem>
#include <string_view>

#include <itkImage.h>

#undef TIMEOUT
#undef CUDA_FOUND
#include "itk_mask.h"

#include "cbctrecon_config.h"

#include "AG17RGBAImage.h"
#include "StructureSet.h"
#include "YK16GrayImage.h"
#include "cbctrecon_types.h"

class CbctRecon;
class Dmap_parms;
class Pcmd_threshold;
class Mask_parms;
class Plm_image_header;
class Dcmtk_rt_study;

enum class enDevice {
  CPU_DEV,
  GPU_DEV,
};

enum class enViewArrange {
  AXIAL_FRONTAL_SAGITTAL = 0,
  FRONTAL_SAGITTAL_AXIAL = 1,
  SAGITTAL_AXIAL_FRONTAL = 2,
};

enum class enRegisterOption {
  PLAST_RIGID = 0,
  PLAST_GRADIENT,
  PLAST_AFFINE,
  PLAST_BSPLINE,
};

namespace fs = std::filesystem;

class CBCTRECON_API CbctRegistration {

public:
  explicit CbctRegistration(CbctRecon *parent);
  ~CbctRegistration();
  CbctRegistration(const CbctRegistration &) = delete;
  void operator=(const CbctRegistration &) = delete;
  CbctRegistration(CbctRegistration &&) = delete;
  void operator=(CbctRegistration &&) = delete;

  UShortImageType::Pointer &
  get_image_from_combotext(const std::string_view ct_type) const;
  void GenPlastiRegisterCommandFile(
      const fs::path &strPathCommandFile, const fs::path &strPathFixedImg,
      const fs::path &strPathMovingImg, const fs::path &strPathOutImg,
      const fs::path &strPathXformOut, enRegisterOption regiOption,
      const std::string &strStageOption1, const std::string &strStageOption2,
      const std::string &strStageOption3, const fs::path &strPathFixedMask,
      bool optim_mse, bool use_cuda, std::string &GradOptionStr) const;

  // get val mm
  static VEC3D GetShiftValueFromGradientXForm(const fs::path &file_path,
                                              bool b_inverse);

  bool PreprocessCT(UShortImageType::Pointer &ct_img, int iAirThresholdShort,
                    const Rtss_modern *rt_structs, const std::string &strRSName,
                    bool fill_bubble, int iBubbleFillingVal,
                    int iAirFillValShort);
  static void autoPreprocessCT(int iAirThresholdShort,
                               UShortImageType::Pointer &spFixed,
                               UShortImageType::Pointer &spMoving);
  float *ManualMoveByDCM() const;
  static UShortImageType::Pointer
  MoveByEclRegistration(const DoubleVector &translation_vec,
                        const DoubleVector &rotation_vec,
                        const UShortImageType::Pointer &ct_img);
  void LoadRTPlan(const fs::path &strDCMPath);

  static void CallingPLMCommand(const fs::path &command_filepath);
  static DoubleVector3DType
  CallingPLMCommandXForm(const fs::path &command_filepath);
  bool CallingGPMCcommand(enDevice device, int n_sims, int n_plans,
                          const fs::path &comma_sep_planfilepath,
                          UShortImageType::Pointer &spFixed,
                          UShortImageType::Pointer &spMoving,
                          UShortImageType::Pointer &spFixedDose,
                          UShortImageType::Pointer &spMovingDose);

  void plm_dmap_main(const fs::path &img_in_fn, const fs::path &img_out_fn) const;
  void plm_threshold_main(std::string &strRange, const fs::path &img_in_fn,
                          const fs::path &img_out_fn) const;
  void plm_mask_main(Mask_operation mask_option, const fs::path &input_fn,
                     const fs::path &mask_fn, const fs::path &output_fn,
                     float mask_value) const;

  void plm_mask_img(Mask_operation mask_option, const fs::path &mask_fn,
                    float mask_value, Plm_image::Pointer &img) const;

  void plm_expansion_contract_msk(const fs::path &strPath_msk,
                                  const fs::path &strPath_msk_exp_cont,
                                  double fExpVal) const;

  void plm_synth_trans_xf(const fs::path &strPath_fixed, const fs::path &strPath_out_xf,
                          double transX, double transY, double transZ) const;

  fs::path gen_and_expand_skinmask_plm(const fs::path &mskSkinCT_manRegi_path,
                                      const fs::path &XFAutoRigid_path,
                                      const fs::path &rawCBCT_path,
                                      double skinExp);
  fs::path gen_bubble_mask_plm(float bubble_thresh, float bubble_fill,
                              const fs::path &strPathOutputCBCT);

  void ProcessCBCT_beforeDeformRegi(
      const fs::path &strPathRawCBCT, const fs::path &strPath_mskSkinCT_manRegi,
      fs::path &strPathOutputCBCT, const fs::path &strPathXFAutoRigid,
      bool bBubbleFilling, bool bPrepareMaskOnly, double skinExp,
      int bubbleThresh,
      int bubbleFill); // 8 mm skin cut + fill air bubbles inside CBCT
  void SetPlmOutputDir(std::string &endFix);

  void
  PostSkinRemovingCBCT(UShortImageType::Pointer &spCBCT,
                       const std::string &voi_name) const; // this function
                                                           // will be called
                                                           // from main Dlg.

  void ThermoMaskRemovingCBCT(UShortImageType::Pointer &spCBCTraw,
                              UShortImageType::Pointer &spCBCTcor,
                              int diffThreshold,
                              int noTouchThreshold /*= 1100*/,
                              double innerMargin, double outerMargin) const;

  void GenShellMask(const fs::path &strPathInputMask, const fs::path &strPathOutputMask,
                    double fInnerMargin, double fOuterMargin) const;

  VEC3D GetIsocenterDCM_FromRTPlan(const fs::path &strFilePath) const;

  // still public:
  CbctRecon *m_pParent{};            // to pull 3D images
  YK16GrayImage m_YKImgFixed[3];     // CBCT in this study
  YK16GrayImage m_YKImgMoving[3];    // CBCT in this study
  YK16GrayImage m_YKDisp[3];         // CBCT in this study
  AG17RGBAImage m_DoseImgFixed[3];   // CBCT in this study
  AG17RGBAImage m_DoseImgMoving[3];  // CBCT in this study
  AG17RGBAImage m_AGDisp_Overlay[3]; // CBCT in this study
  std::unique_ptr<Rtss_roi_modern> WEPL_voi;
  std::unique_ptr<Rtss_roi_modern> extra_voi;
  std::unique_ptr<Rtss_roi_modern> cur_voi;
  bool dose_loaded = false;

  fs::path m_strPathPlastimatch;     // full path
  fs::path m_strPathCTSkin;          // shared data among functions
  fs::path m_strPathCTSkin_autoRegi; // shared data among functions w/ margin (10
                                    // mm), registered to CBCT
  fs::path m_strPathCTSkin_deformRegi; // this is for WEPL calculation
  fs::path m_strPathXFAutoRigid;       // Updated during SLT_DoRegistrationRigid
  fs::path m_strPathMskCBCTBubble;     // Updated during SLT_DoRegistrationRigid

  std::unique_ptr<Dcmtk_rt_study> m_pDcmStudyPlan;
};

#endif // CBCTREGISTRATION_H
