#ifndef DLGREGISTRATION_H
#define DLGREGISTRATION_H

// #include <QDialog>
// #include <QString>

#include <itkImage.h>

#include "StructureSet.h"
#include "cbctrecon.h"
#include "cbctregistration.h"
#include "qyklabel.h"
#include "ui_DlgRegistration.h"

class QDialog;
class QString;
class CbctReconWidget;

class DlgRegistration : public QDialog, public Ui::DlgRegistrationClass {
  Q_OBJECT

public:
  CbctReconWidget *m_pParent{}; // to pull 3D images
  std::unique_ptr<CbctRegistration> m_cbctregistration;
  UShortImageType::Pointer m_spFixed;  // pointer only, for display
  UShortImageType::Pointer m_spMoving; // pointer only, for display

private: // Just pointers to m_cbctregistration members, for convienience
  YK16GrayImage *m_YKImgFixed;
  YK16GrayImage *m_YKImgMoving;
  YK16GrayImage *m_YKDisp;
  AG17RGBAImage *m_DoseImgFixed;
  AG17RGBAImage *m_DoseImgMoving;
  AG17RGBAImage *m_AGDisp_Overlay;
  UShortImageType::Pointer m_spFixedDose;  // pointer only, for display
  UShortImageType::Pointer m_spMovingDose; // pointer only, for display
  int m_enViewArrange{};
  // YK16GrayImage m_YKImgMoving[3];    //RefCT

  bool m_bPressedLeft[3]{}; // Left Mouse Pressed but not released
  bool m_bPressedRight[3]{};
  QPoint m_ptWindowLevelStart; // data point
  QPoint m_ptPanStart;         // data point
  QPoint m_ptTmpOriginalDataOffset;
  int m_iTmpOriginalW{};
  int m_iTmpOriginalL{};

public:
  DlgRegistration();
  DlgRegistration(CbctReconWidget *parent);
  // ~DlgRegistration() = default;
  void UpdateVOICombobox(ctType ct_type);
  void UpdateListOfComboBox(int idx);
  void SelectComboExternal(int idx, enREGI_IMAGES iImage);
  void LoadImgFromComboBox(int idx, QString &strSelectedComboTxt);
  void initDlgRegistration(QString &strDCMUID);

private:
  void whenFixedImgLoaded(); // should be called by comboBox
  void whenMovingImgLoaded();
  void initOverlapWndSize();
  void shiftSliceSlider();
  void updateSliceLabel();
  void UpdateSplit(int viewIdx, qyklabel *pOverlapWnd);
  void MousePressedRight(int wndIdx, qyklabel *pWnd);

  void childEvent(QChildEvent *event) override;
  bool eventFilter(QObject *target, QEvent *event) override;

  void ImageManualMove(int direction, double resol);
  void ImageManualMoveOneShot(float shiftX, float shiftY,
                              float shiftZ); // DICOM coordinate
  void AddImageToCombo(int comboIdx,
                       enREGI_IMAGES option); // comboIdx 0: fixed, 1: moving
  void LoadVOIFromComboBox(int idx, QString &strSelectedComboTxt);

public slots:
  void SLT_CrntPosGo();
  void SLT_DrawImageWhenSliceChange(); // upper level drawing: big calculation
  void SLT_DrawImageInFixedSlice();    // lower level Drawing func.

  void SLT_UpdateSplit1(); // lower level Drawing func. //Mouse Move even
  void SLT_UpdateSplit2(); // lower level Drawing func.//Mouse Move even
  void SLT_UpdateSplit3(); // lower level Drawing func.//Mouse Move even

  void SLT_CancelMouseAction();

  void SLT_MouseWheelUpdate1();
  void SLT_MouseWheelUpdate2();
  void SLT_MouseWheelUpdate3();

  void SLT_MousePressedLeft1();
  void SLT_MousePressedLeft2();
  void SLT_MousePressedLeft3();

  void SLT_MousePressedRight1();
  void SLT_MousePressedRight2();
  void SLT_MousePressedRight3();

  void SLT_MouseReleasedLeft1();
  void SLT_MouseReleasedLeft2();
  void SLT_MouseReleasedLeft3();

  void SLT_MouseReleasedRight1();
  void SLT_MouseReleasedRight2();
  void SLT_MouseReleasedRight3();
  void SLT_ChangeView(); // 3 toggle button

  void SLT_RestoreImageSingle();
  void SLT_RestoreImageAll();

  void SLT_DoRegistrationRigid();
  void SLT_DoRegistrationDeform();
  void SLT_DoRegistrationGradient();
  void SLT_ManualMoveByDCMPlan();
  void SLT_ManualMoveByDCMPlanOpen();
  void SLT_gPMCrecalc();
  void SLT_WEPLcalc();
  void SLT_BringFocusToEnableArrow(bool bChecked);

  void SLT_KeyMoving(bool bChecked);

  // Image selection event from Combobox
  void SLT_FixedImageSelected(
      QString selText); // here, when fixed_image_loaded function will be called
  void SLT_MovingImageSelected(QString selText); // here, when
                                                 // mvoing_image_loaded function
                                                 // will be called

  void SLT_RestoreMovingImg(); // copy refCT to Manual Moving Image.

  void SLT_PreProcessCT();

  void SLT_PassFixedImgForAnalysis();
  void SLT_PassMovingImgForAnalysis();

  void SLT_ExchangeRawRef();

  /* After manual move: this will trigger skin cropping and uncheck the key
   * moving*/
  void SLT_ConfirmManualRegistration();
  void SLT_IntensityNormCBCT();
  void SLT_DoLowerMaskIntensity(); // button
  void SLT_Override();

public:
  Ui::DlgRegistrationClass ui{};
};

#endif
