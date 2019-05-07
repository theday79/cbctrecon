#ifndef CbctRegistrationTest_H
#define CbctRegistrationTest_H

#include <itkImage.h>

#include "MyQComboBox.hxx"
#include "StructureSet.h"
#include "cbctrecon.h"
#include "cbctregistration.h"

class qyklabel;
class QDialog;
class QString;
class CbctReconTest;

class CbctRegistrationTest {

public:
  CbctReconTest *m_pParent{}; // to pull 3D images
  std::unique_ptr<CbctRegistration> m_cbctregistration;
  UShortImageType::Pointer m_spFixed;  // pointer only, for display
  UShortImageType::Pointer m_spMoving; // pointer only, for display
  std::unique_ptr<MyQComboBox> ui_comboBox_VOItoCropBy;

private: // Just pointers to m_cbctregistration members, for convienience
  YK16GrayImage *m_YKImgFixed;
  YK16GrayImage *m_YKImgMoving;
  YK16GrayImage *m_YKDisp;
  AG17RGBAImage *m_DoseImgFixed;
  AG17RGBAImage *m_DoseImgMoving;
  AG17RGBAImage *m_AGDisp_Overlay;
  UShortImageType::Pointer m_spFixedDose;  // pointer only, for display
  UShortImageType::Pointer m_spMovingDose; // pointer only, for display
  // YK16GrayImage m_YKImgMoving[3];    //RefCT

  QPoint m_ptWindowLevelStart; // data point
  QPoint m_ptPanStart;         // data point
  QPoint m_ptTmpOriginalDataOffset;

  /*UI simulation:*/
  bool ui_checkBoxKeyMoving = false;
  bool ui_checkBoxCropBkgroundCT = true;
  bool ui_checkBoxUseROIForRigid = true;
  bool ui_checkBoxUseROIForDIR = true;
  bool ui_checkBoxFillBubbleCT = false;
  bool ui_checkBoxFillBubbleCBCT = false;
  QString ui_lineEditFOVPos = "0.0,0.0,190.0";
  QString ui_lineEditGradOption = "0.7,0.7,0.7";
  QString ui_lineEditCBCTSkinCropBfRegid = "10.0";
  QString ui_lineEditCBCTSkinCropBfDIR = "8.0";
  int ui_spinBoxBkFillCT = 0;
  int ui_spinBoxBubFillCT = 1000;
  int ui_spinBoxBkDetectCT = 500;
  QString ui_comboBoxImToCropFill = "Moving";
  QString ui_lineEditArgument1 = "2,2,1,30,0.00001,0.005,5";
  QString ui_lineEditArgument2 = "";
  QString ui_lineEditArgument3 = "";
  std::unique_ptr<MyQComboBox> ui_comboBoxImgMoving;
  std::unique_ptr<MyQComboBox> ui_comboBoxImgFixed;
  std::unique_ptr<MyQComboBox> ui_comboBox_VOI;
  bool ui_radioButton_mse = true;
  bool ui_radioButton_UseCUDA = false;
  bool ui_radioButton_UseCPU = true;

  // experimental:
  bool ui_checkBoxRemoveMaskAfterCor = false;
  QString ui_lineEditRawCorThre = "100";
  QString ui_lineEditiNoTouchThreshold = "1100";
  QString ui_lineEditThermoInner = "5.0";
  QString ui_lineEditThermoOuter = "10.0";

public:
  CbctRegistrationTest();
  explicit CbctRegistrationTest(CbctReconTest *parent);
  // ~CbctRegistrationTest() = default;
  void UpdateVOICombobox(ctType ct_type) const;
  void UpdateListOfComboBox(int idx) const;
  void SelectComboExternal(int idx, enREGI_IMAGES iImage);
  void LoadImgFromComboBox(int idx, QString &strSelectedComboTxt);
  void initCbctRegistrationTest(QString &strDCMUID);

  void ImageManualMove(int direction, double resol) const;
  void ImageManualMoveOneShot(float shiftX, float shiftY,
                              float shiftZ) const; // DICOM coordinate
  void
  AddImageToCombo(int comboIdx,
                  enREGI_IMAGES option) const; // comboIdx 0: fixed, 1: moving
  void LoadVOIFromComboBox(int idx, QString &strSelectedComboTxt) const;

  void SLT_RestoreImageSingle() const;
  void SLT_RestoreImageAll() const;

  void SLT_DoRegistrationRigid();
  void SLT_DoRegistrationDeform();
  void SLT_DoRegistrationGradient();
  void SLT_ManualMoveByDCMPlan();
  void SLT_ManualMoveByDCMPlanOpen(QString &filePath);
  void SLT_gPMCrecalc(std::vector<QString> &dcm_plans, size_t n_sims);
  void SLT_WEPLcalc(int gantry_angle, int couch_angle) const;
  void SLT_DoEclRegistration(const DoubleVector &translation_vec,
                             const DoubleVector &rotation_vec);
  void SLT_ResetEclRegistration();

  void SLT_KeyMoving(bool bChecked);

  // Image selection event from Combobox
  void SLT_FixedImageSelected(
      QString selText); // here, when fixed_image_loaded function will be called
  void SLT_MovingImageSelected(QString selText); // here, when
                                                 // mvoing_image_loaded function
                                                 // will be called

  void SLT_RestoreMovingImg(); // copy refCT to Manual Moving Image.

  void SLT_PreProcessCT();

  /* After manual move: this will trigger skin cropping and uncheck the key
   * moving*/
  void SLT_ConfirmManualRegistration();
  void SLT_IntensityNormCBCT(float fROI_Radius);
  void SLT_DoLowerMaskIntensity(); // button
  void SLT_Override(int sliderPosIdxX, int sliderPosIdxY, int sliderPosIdxZ,
                    size_t radius, int new_value,
                    const QString &img_to_override) const;
};

#endif
