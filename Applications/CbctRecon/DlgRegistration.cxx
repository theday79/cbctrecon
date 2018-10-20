#include "DlgRegistration.h"

#include <QDialog>
#include <QFileDialog>
#include <QMessageBox>
#include <QProcess>
#include <QString>
#include <qcombobox.h>

// configs
// #include <itkConfigure.h>
// #include <plm_config.h>

#include "StructureSet.h"
#include "cbctrecon.h"
#include "cbctrecon_compute.h"
#include "cbctrecon_mainwidget.h"
#include "cbctregistration.h"

#define FIXME_BACKGROUND_MAX (-1200)

// using namespace std;

//#include "itkImageSliceConstIteratorWithIndex.h"
//#include "itkImageSliceIteratorWithIndex.h"

DlgRegistration::DlgRegistration() {
  /* Sets up the GUI */
  ui.setupUi(this);
  m_YKImgFixed = nullptr;
  m_YKImgMoving = nullptr;
  m_YKDisp = nullptr;

  m_DoseImgFixed = nullptr;
  m_DoseImgMoving = nullptr;
  m_AGDisp_Overlay = nullptr;
}

DlgRegistration::DlgRegistration(CbctReconWidget *parent) : QDialog(parent) {
  /* Sets up the GUI */
  ui.setupUi(this);
  m_pParent = parent;
  m_cbctregistration =
      std::make_unique<CbctRegistration>(parent->m_cbctrecon.get());

  m_YKImgFixed = m_cbctregistration->m_YKImgFixed;
  m_YKImgMoving = m_cbctregistration->m_YKImgMoving;
  m_YKDisp = m_cbctregistration->m_YKDisp;

  m_DoseImgFixed = m_cbctregistration->m_DoseImgFixed;
  m_DoseImgMoving = m_cbctregistration->m_DoseImgMoving;
  m_AGDisp_Overlay = m_cbctregistration->m_AGDisp_Overlay;

  // m_spFixed = m_cbctregistration->m_spFixed;
  // m_spMoving = m_cbctregistration->m_spMoving;
  // m_spFixedDose = m_cbctregistration->m_spFixedDose;
  // m_spMovingDose = m_cbctregistration->m_spMovingDose;

  connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Move()), this,
          SLOT(SLT_UpdateSplit1())); // added
  connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Move()), this,
          SLOT(SLT_UpdateSplit2())); // added
  connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Move()), this,
          SLOT(SLT_UpdateSplit3())); // added

  connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Pressed_Left()), this,
          SLOT(SLT_MousePressedLeft1())); // added
  connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Pressed_Left()), this,
          SLOT(SLT_MousePressedLeft2())); // added
  connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Pressed_Left()), this,
          SLOT(SLT_MousePressedLeft3())); // added
  connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Pressed_Right()), this,
          SLOT(SLT_MousePressedRight1())); // added
  connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Pressed_Right()), this,
          SLOT(SLT_MousePressedRight2())); // added
  connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Pressed_Right()), this,
          SLOT(SLT_MousePressedRight3())); // added

  connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Released_Left()), this,
          SLOT(SLT_MouseReleasedLeft1())); // added
  connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Released_Left()), this,
          SLOT(SLT_MouseReleasedLeft2())); // added
  connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Released_Left()), this,
          SLOT(SLT_MouseReleasedLeft3())); // added
  connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Released_Right()), this,
          SLOT(SLT_MouseReleasedRight1())); // added
  connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Released_Right()), this,
          SLOT(SLT_MouseReleasedRight2())); // added
  connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Released_Right()), this,
          SLOT(SLT_MouseReleasedRight3())); // added

  connect(ui.labelOverlapWnd1, SIGNAL(FocusOut()), this,
          SLOT(SLT_CancelMouseAction())); // added
  connect(ui.labelOverlapWnd2, SIGNAL(FocusOut()), this,
          SLOT(SLT_CancelMouseAction())); // added
  connect(ui.labelOverlapWnd3, SIGNAL(FocusOut()), this,
          SLOT(SLT_CancelMouseAction())); // added

  connect(ui.labelOverlapWnd1, SIGNAL(Mouse_Wheel()), this,
          SLOT(SLT_MouseWheelUpdate1())); // added
  connect(ui.labelOverlapWnd2, SIGNAL(Mouse_Wheel()), this,
          SLOT(SLT_MouseWheelUpdate2())); // added
  connect(ui.labelOverlapWnd3, SIGNAL(Mouse_Wheel()), this,
          SLOT(SLT_MouseWheelUpdate3())); // added

  SLT_CancelMouseAction();
}

void DlgRegistration::initDlgRegistration(QString &strDCMUID) {
  m_cbctregistration->SetPlmOutputDir(strDCMUID);

  UShortImageType::Pointer spNull;
  // unlink all of the pointers
  // m_pParent->m_spReconImg->Delete(); //fixed image // ID: RawCBCT
  m_cbctregistration->m_pParent->m_spRefCTImg = spNull;
  m_cbctregistration->m_pParent->m_spManualRigidCT =
      spNull; // copied from RefCTImg; ID: RefCT --> Moving Img, cloned
  m_cbctregistration->m_pParent->m_spAutoRigidCT = spNull; // ID: AutoRigidCT
  m_cbctregistration->m_pParent->m_spDeformedCT1 =
      spNull; // Deformmation will be carried out based
              // on Moving IMage of GUI //AutoDeformCT1
  m_cbctregistration->m_pParent->m_spDeformedCT2 = spNull;      // AutoDeformCT2
  m_cbctregistration->m_pParent->m_spDeformedCT3 = spNull;      // AutoDeformCT3
  m_cbctregistration->m_pParent->m_spDeformedCT_Final = spNull; // AutoDeformCT3

  ui.checkBoxKeyMoving->setChecked(false);
  ui.lineEditOriginChanged->setText("");

  // show();

  UpdateListOfComboBox(0);
  UpdateListOfComboBox(1);
  // if not found, just skip
  SelectComboExternal(0, REGISTER_RAW_CBCT);     // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID); // WILL BE IGNORED
}

void DlgRegistration::SLT_CrntPosGo() {
  // m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_AXIAL, 0.0,
  // m_dspDlgRegi1);  m_pParent->GetValueFrom3DImageUshort(5, 5, 5,
  // m0_pParent->m_spReconImg);  CbctRecon* pParent=
  // (CbctRecon*)(this->parent());
  /*whenFixedImgLoaded();
  SLT_DrawImageWhenSliceChange();    	    */
  if (m_spFixed == nullptr) {
    return;
  }

  // DICOMN position, mm
  double curDCMPosX = ui.lineEditCurPosX->text().toDouble();
  double curDCMPosY = ui.lineEditCurPosY->text().toDouble();
  double curDCMPosZ = ui.lineEditCurPosZ->text().toDouble();

  // UShortImageType::SizeType imgSize =
  // m_spFixed->GetRequestedRegion().GetSize(); //1016x1016 x z
  UShortImageType::PointType imgOrigin = m_spFixed->GetOrigin();
  UShortImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();

  // double curPhysPos[3];
  // curPhysPos[0] = imgOrigin[2] + sliderPosIdxZ*imgSpacing[2] ; //Z in default
  // setting  curPhysPos[1] = imgOrigin[1] + sliderPosIdxY*imgSpacing[1]; //Y
  // curPhysPos[2] = imgOrigin[0] + sliderPosIdxX*imgSpacing[0]; //Z

  int iSliderPosIdxZ =
      qRound((curDCMPosZ - imgOrigin[2]) / static_cast<double>(imgSpacing[2]));
  int iSliderPosIdxY =
      qRound((curDCMPosY - imgOrigin[1]) / static_cast<double>(imgSpacing[1]));
  int iSliderPosIdxX =
      qRound((curDCMPosX - imgOrigin[0]) / static_cast<double>(imgSpacing[0]));

  ui.sliderPosDisp1->setValue(iSliderPosIdxZ);
  ui.sliderPosDisp2->setValue(iSliderPosIdxY);
  ui.sliderPosDisp3->setValue(iSliderPosIdxX);
}

void DlgRegistration::SLT_DrawImageWhenSliceChange() {
  /*if (m_pParent == NULL)
      return;	*/
  if (m_spFixed == nullptr) {
    return;
  }

  int sliderPosIdxZ, sliderPosIdxY, sliderPosIdxX;

  switch (m_enViewArrange) {
  case AXIAL_FRONTAL_SAGITTAL:
    sliderPosIdxZ =
        ui.sliderPosDisp1
            ->value(); // Z corresponds to axial, Y to frontal, X to sagittal
    sliderPosIdxY = ui.sliderPosDisp2->value();
    sliderPosIdxX = ui.sliderPosDisp3->value();
    break;
  case FRONTAL_SAGITTAL_AXIAL:
    sliderPosIdxY = ui.sliderPosDisp1->value();
    sliderPosIdxX = ui.sliderPosDisp2->value();
    sliderPosIdxZ = ui.sliderPosDisp3->value();

    break;
  case SAGITTAL_AXIAL_FRONTAL:
    sliderPosIdxX = ui.sliderPosDisp1->value();
    sliderPosIdxZ = ui.sliderPosDisp2->value();
    sliderPosIdxY = ui.sliderPosDisp3->value();
    break;
  default:
    sliderPosIdxZ =
        ui.sliderPosDisp1
            ->value(); // Z corresponds to axial, Y to frontal, X to sagittal
    sliderPosIdxY = ui.sliderPosDisp2->value();
    sliderPosIdxX = ui.sliderPosDisp3->value();
    break;
  }

  UShortImageType::SizeType imgSize =
      m_spFixed->GetRequestedRegion().GetSize(); // 1016x1016 x z
  UShortImageType::PointType imgOrigin = m_spFixed->GetOrigin();
  UShortImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();

  double curPhysPos[3];
  curPhysPos[0] =
      imgOrigin[2] + sliderPosIdxZ * imgSpacing[2]; // Z in default setting
  curPhysPos[1] = imgOrigin[1] + sliderPosIdxY * imgSpacing[1]; // Y
  curPhysPos[2] = imgOrigin[0] + sliderPosIdxX * imgSpacing[0]; // Z

  // This caused the problem!!!
  /*m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_AXIAL, curPhysPos1,
  m_dspDlgRegi1); m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg,
  PLANE_FRONTAL, curPhysPos2, m_dspDlgRegi2);
  m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_SAGITTAL, curPhysPos3,
  m_dspDlgRegi3);*/

  /*m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_AXIAL, curPhysPos1,
  m_dspDlgRegi1); m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg,
  PLANE_FRONTAL, curPhysPos2, m_dspDlgRegi2);
  m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_SAGITTAL, curPhysPos3,
  m_dspDlgRegi3);        */

  int refIdx = 3 - m_enViewArrange;

  if (ui.checkBoxDrawCrosshair->isChecked()) {
    m_YKDisp[refIdx % 3].m_bDrawCrosshair = true;
    m_YKDisp[(refIdx + 1) % 3].m_bDrawCrosshair = true;
    m_YKDisp[(refIdx + 2) % 3].m_bDrawCrosshair = true;

    // m_YKDisp[0]// Left Top image, the largest

    m_YKDisp[refIdx % 3].m_ptCrosshair.setX(sliderPosIdxX); // axial
    m_YKDisp[refIdx % 3].m_ptCrosshair.setY(sliderPosIdxY);

    m_YKDisp[(refIdx + 1) % 3].m_ptCrosshair.setX(sliderPosIdxX); // Frontal
    m_YKDisp[(refIdx + 1) % 3].m_ptCrosshair.setY(static_cast<int>(imgSize[2]) -
                                                  sliderPosIdxZ - 1);

    m_YKDisp[(refIdx + 2) % 3].m_ptCrosshair.setX(sliderPosIdxY); // Sagittal
    m_YKDisp[(refIdx + 2) % 3].m_ptCrosshair.setY(static_cast<int>(imgSize[2]) -
                                                  sliderPosIdxZ - 1);

    m_YKImgFixed[0].m_bDrawCrosshair = true;
    m_YKImgFixed[1].m_bDrawCrosshair = true;
    m_YKImgFixed[2].m_bDrawCrosshair = true;

    m_YKImgFixed[0].m_ptCrosshair.setX(sliderPosIdxX); // sagittal slider
    m_YKImgFixed[0].m_ptCrosshair.setY(sliderPosIdxY);

    m_YKImgFixed[1].m_ptCrosshair.setX(sliderPosIdxX); // sagittal slider
    m_YKImgFixed[1].m_ptCrosshair.setY(static_cast<int>(imgSize[2]) -
                                       sliderPosIdxZ - 1);

    m_YKImgFixed[2].m_ptCrosshair.setX(sliderPosIdxY); // sagittal slider
    m_YKImgFixed[2].m_ptCrosshair.setY(static_cast<int>(imgSize[2]) -
                                       sliderPosIdxZ - 1);
  } else {
    m_YKDisp[0].m_bDrawCrosshair = false;
    m_YKDisp[1].m_bDrawCrosshair = false;
    m_YKDisp[2].m_bDrawCrosshair = false;

    m_YKImgFixed[0].m_bDrawCrosshair = false;
    m_YKImgFixed[1].m_bDrawCrosshair = false;
    m_YKImgFixed[2].m_bDrawCrosshair = false;
  }

  // center of the split value is passed by m_YKImgFixed;
  // m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_AXIAL,
  // curPhysPos[0], m_YKImgFixed[0]);
  // m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_FRONTAL,
  // curPhysPos[1], m_YKImgFixed[1]);
  // m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_SAGITTAL,
  // curPhysPos[2], m_YKImgFixed[2]);

  if (m_spMoving != nullptr) {
    m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
        m_spFixed, m_spMoving, PLANE_AXIAL, curPhysPos[0], m_YKImgFixed[0],
        m_YKImgMoving[0]);
    m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
        m_spFixed, m_spMoving, PLANE_FRONTAL, curPhysPos[1], m_YKImgFixed[1],
        m_YKImgMoving[1]);
    m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
        m_spFixed, m_spMoving, PLANE_SAGITTAL, curPhysPos[2], m_YKImgFixed[2],
        m_YKImgMoving[2]);
    if (m_cbctregistration->dose_loaded) {
      m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
          m_spFixedDose, m_spMovingDose, PLANE_AXIAL, curPhysPos[0],
          m_DoseImgFixed[0], m_DoseImgMoving[0]);
      m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
          m_spFixedDose, m_spMovingDose, PLANE_FRONTAL, curPhysPos[1],
          m_DoseImgFixed[1], m_DoseImgMoving[1]);
      m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
          m_spFixedDose, m_spMovingDose, PLANE_SAGITTAL, curPhysPos[2],
          m_DoseImgFixed[2], m_DoseImgMoving[2]);
    }
  } else {
    m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
        m_spFixed, m_spFixed, PLANE_AXIAL, curPhysPos[0], m_YKImgFixed[0],
        m_YKImgMoving[0]);
    m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
        m_spFixed, m_spFixed, PLANE_FRONTAL, curPhysPos[1], m_YKImgFixed[1],
        m_YKImgMoving[1]);
    m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
        m_spFixed, m_spFixed, PLANE_SAGITTAL, curPhysPos[2], m_YKImgFixed[2],
        m_YKImgMoving[2]);
    if (m_cbctregistration->dose_loaded) {
      m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
          m_spFixedDose, m_spFixedDose, PLANE_AXIAL, curPhysPos[0],
          m_DoseImgFixed[0], m_DoseImgMoving[0]);
      m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
          m_spFixedDose, m_spFixedDose, PLANE_FRONTAL, curPhysPos[1],
          m_DoseImgFixed[1], m_DoseImgMoving[1]);
      m_pParent->m_cbctrecon->Draw2DFrom3DDouble(
          m_spFixedDose, m_spFixedDose, PLANE_SAGITTAL, curPhysPos[2],
          m_DoseImgFixed[2], m_DoseImgMoving[2]);
    }
  }

  // Update position lineEdit
  QString strPos1, strPos2, strPos3;
  strPos1.sprintf("%3.1f", curPhysPos[0]);
  strPos2.sprintf("%3.1f", curPhysPos[1]);
  strPos3.sprintf("%3.1f", curPhysPos[2]);

  ui.lineEditCurPosX->setText(strPos3);
  ui.lineEditCurPosY->setText(strPos2);
  ui.lineEditCurPosZ->setText(strPos1);

  ////Update Origin text box
  UShortImageType::PointType imgOriginFixed = m_spFixed->GetOrigin();
  QString strOriFixed;
  strOriFixed.sprintf("%3.4f, %3.4f, %3.4f", imgOriginFixed[0],
                      imgOriginFixed[1], imgOriginFixed[2]);
  ui.lineEditOriginFixed->setText(strOriFixed);

  if (m_spMoving != nullptr) {
    UShortImageType::PointType imgOriginMoving = m_spMoving->GetOrigin();
    QString strOriMoving;
    strOriMoving.sprintf("%3.4f, %3.4f, %3.4f", imgOriginMoving[0],
                         imgOriginMoving[1], imgOriginMoving[2]);
    ui.lineEditOriginMoving->setText(strOriMoving);
  }

  if (m_cbctregistration->cur_voi != nullptr) {
    auto *Wnd1_contour = &ui.labelOverlapWnd1->m_vPt;
    auto *Wnd2_contour = &ui.labelOverlapWnd2->m_vPt;
    auto *Wnd3_contour = &ui.labelOverlapWnd3->m_vPt;
    Wnd1_contour->clear();
    Wnd2_contour->clear();
    Wnd3_contour->clear();

    for (auto contour : m_cbctregistration->cur_voi->pslist) {
      if (contour.coordinates.size() == 0) {
        continue;
      }
      auto first_point = contour.coordinates.at(0);
      // Axial
      if (first_point.z > curPhysPos[0] - imgSpacing[2] &&
          first_point.z < curPhysPos[0] + imgSpacing[2]) {
        for (auto point : contour.coordinates) {
          Wnd1_contour->push_back(QPoint(point.x, point.y));
        }
      }
      for (auto point : contour.coordinates) {
        // Frontal
        if (point.y > curPhysPos[1] - imgSpacing[1] &&
            point.y < curPhysPos[1] + imgSpacing[1]) {
          Wnd2_contour->push_back(QPoint(point.x, point.z));
        }
        // Sagittal
        if (point.x > curPhysPos[2] - imgSpacing[0] &&
            point.x < curPhysPos[2] + imgSpacing[0]) {
          Wnd3_contour->push_back(QPoint(point.y, point.z));
        }
      }
    }
    // Get contour for axial, sagittal and frontal
    // create plotable Qt objects from the contours
    // plot Qt objects on ui.labelOverlapWnd*
    ui.labelOverlapWnd1->m_bDrawPoints = true;
    ui.labelOverlapWnd2->m_bDrawPoints = true;
    ui.labelOverlapWnd3->m_bDrawPoints = true;
  }
  /*qDebug() << strOriFixed;
  qDebug() << strOriMoving;*/
  //

  SLT_DrawImageInFixedSlice();
}
// Display is not included here
void DlgRegistration::whenFixedImgLoaded() {
  if (m_spFixed == nullptr) {
    return;
  }
  /*if (!m_pParent->m_spReconImg)
      return;

      m_spFixed = m_pParent->m_spReconImg;*/

  UShortImageType::SizeType imgSize =
      m_spFixed->GetRequestedRegion().GetSize(); // 1016x1016 x z
  // UShortImageType::PointType imgOrigin = m_spFixed->GetOrigin();
  // UShortImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();

  // to avoid first unnecessary action.
  disconnect(ui.sliderPosDisp1, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageWhenSliceChange()));
  disconnect(ui.sliderPosDisp2, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageWhenSliceChange()));
  disconnect(ui.sliderPosDisp3, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageWhenSliceChange()));

  disconnect(ui.sliderFixedW, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageInFixedSlice()));
  disconnect(ui.sliderFixedL, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageInFixedSlice()));
  disconnect(ui.sliderMovingW, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageInFixedSlice()));
  disconnect(ui.sliderMovingL, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageInFixedSlice()));

  initOverlapWndSize();

  ui.sliderPosDisp1->setMinimum(0);
  ui.sliderPosDisp1->setMaximum(static_cast<int>(imgSize[2] - 1));
  int curPosZ = static_cast<int>(imgSize[2] / 2);
  ui.sliderPosDisp1->setValue(curPosZ);

  ui.sliderPosDisp2->setMinimum(0);
  ui.sliderPosDisp2->setMaximum(static_cast<int>(imgSize[1] - 1));
  int curPosY = static_cast<int>(imgSize[1] / 2);
  ui.sliderPosDisp2->setValue(curPosY);

  ui.sliderPosDisp3->setMinimum(0);
  ui.sliderPosDisp3->setMaximum(static_cast<int>(imgSize[0] - 1));
  int curPosX = static_cast<int>(imgSize[0] / 2);
  ui.sliderPosDisp3->setValue(curPosX);

  QPoint x_split = QPoint(static_cast<int>(imgSize[0] / 2),
                          static_cast<int>(imgSize[1] / 2));
  QPoint y_split = QPoint(static_cast<int>(imgSize[0] / 2),
                          static_cast<int>(imgSize[2] / 2));
  QPoint z_split = QPoint(static_cast<int>(imgSize[1] / 2),
                          static_cast<int>(imgSize[2] / 2));

  m_YKDisp[0].SetSplitCenter(x_split);
  m_YKDisp[1].SetSplitCenter(y_split);
  m_YKDisp[2].SetSplitCenter(z_split);

  // Temporarily load the Disp Image to calculate the window level
  // m_pParent->Draw2DFrom3D(m_pParent->m_spReconImg, PLANE_AXIAL, 0.0,
  // m_YKDisp[0]);

  ////Window Level based on half square
  // m_YKImgFixed[0].setROI((int)(m_YKImgFixed[0].m_iWidth/4.0),
  //						(int)(m_YKImgFixed[0].m_iHeight/4.0),
  //						(int)(m_YKImgFixed[0].m_iWidth*3.0/4.0),
  //						(int)(m_YKImgFixed[0].m_iHeight*3.0/4.0));
  // m_YKImgFixed[0].CalcImageInfo_ROI();
  //
  // int fixedImgWinMid = (int)m_YKImgFixed[0].m_fPixelMean_ROI;
  // int fixedImgWinWidth = (int)(m_YKImgFixed[0].m_fPixelSD_ROI*2.0);

  // int iSliderMin = fixedImgWinMid - fixedImgWinWidth/2.0;
  // int iSliderMax = fixedImgWinMid + fixedImgWinWidth/2.0;

  // if (iSliderMin <= 0)
  //  iSliderMin = 0;
  // if (iSliderMax >= 65535)
  //  iSliderMax = 10000;

  int iSliderW = 2000;
  int iSliderL = 1024;

  ui.sliderFixedW->setValue(iSliderW);
  ui.sliderMovingW->setValue(iSliderW);

  ui.sliderFixedL->setValue(iSliderL);
  ui.sliderMovingL->setValue(iSliderL);

  connect(ui.sliderPosDisp1, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageWhenSliceChange()));
  connect(ui.sliderPosDisp2, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageWhenSliceChange()));
  connect(ui.sliderPosDisp3, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageWhenSliceChange()));

  connect(ui.sliderFixedW, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageInFixedSlice()));
  connect(ui.sliderFixedL, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageInFixedSlice()));
  connect(ui.sliderMovingW, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageInFixedSlice()));
  connect(ui.sliderMovingL, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageInFixedSlice()));
}

void DlgRegistration::whenMovingImgLoaded() {
  // do nothing so far
  /*if (!m_pParent->m_spRefCTImg)
        return;

  m_spMoving = m_pParent->m_spRefCTImg;*/
}

void DlgRegistration::SLT_DrawImageInFixedSlice() // Display Swap here!
{
  // Constitute m_YKDisp from Fixed and Moving

  if (ui.checkBoxDrawSplit->isChecked()) {
    for (int i = 0; i < 3; i++) {
      int idxAdd = m_enViewArrange; // m_iViewArrange = 0,1,2
      if (idxAdd + i >= 3) {
        idxAdd = idxAdd - 3;
      }

      m_YKDisp[i].SetSpacing(m_YKImgFixed[i + idxAdd].m_fSpacingX,
                             m_YKImgFixed[i + idxAdd].m_fSpacingY);

      m_YKDisp[i].SetSplitOption(PRI_LEFT_TOP);
      // m_YKDisp[i].SetSplitCenter(QPoint dataPt);//From mouse event
      if (!m_YKDisp[i].ConstituteFromTwo(m_YKImgFixed[i + idxAdd],
                                         m_YKImgMoving[i + idxAdd])) {
        std::cout << "Image error " << i + 1 << " th view" << std::endl;
      }
    }
  } else {
    for (int i = 0; i < 3; i++) {
      int addedViewIdx = m_enViewArrange;
      if (i + addedViewIdx >= 3) {
        addedViewIdx = addedViewIdx - 3;
      }

      // m_YKDisp[i].CreateImage(m_YKImgFixed[i].m_iWidth,
      // m_YKImgFixed[i].m_iHeight, 0);
      // m_YKDisp[i].CopyFromBuffer(m_YKImgFixed[i].m_pData,
      // m_YKImgFixed[i].m_iWidth, m_YKImgFixed[i].m_iHeight);  std::cout <<
      // "width" << m_YKImgFixed[i].m_iWidth << std::endl;  std::cout <<
      // "height"
      // << m_YKImgFixed[i].m_iHeight << std::endl;
      m_YKDisp[i].CloneImage(m_YKImgFixed[i + addedViewIdx]);
    }
  }

  // For dose overlay
  if (m_cbctregistration->dose_loaded) {
    if (ui.checkBoxDrawSplit->isChecked()) {
      for (int i = 0; i < 3; i++) {
        int idxAdd = m_enViewArrange; // m_iViewArrange = 0,1,2
        if (idxAdd + i >= 3) {
          idxAdd = idxAdd - 3;
        }

        m_AGDisp_Overlay[i].SetSpacing(m_DoseImgFixed[i + idxAdd].m_fSpacingX,
                                       m_DoseImgFixed[i + idxAdd].m_fSpacingY);

        m_AGDisp_Overlay[i].SetSplitOption(PRI_LEFT_TOP);
        if (!m_AGDisp_Overlay[i].ConstituteFromTwo(
                m_DoseImgFixed[i + idxAdd], m_DoseImgMoving[i + idxAdd])) {
          std::cout << "Dose Image error " << i + 1 << " th view" << std::endl;
        }
      }
    } else {
      for (int i = 0; i < 3; i++) {
        int addedViewIdx = m_enViewArrange;
        if (i + addedViewIdx >= 3) {
          addedViewIdx = addedViewIdx - 3;
        }

        m_AGDisp_Overlay[i].CloneImage(m_DoseImgFixed[i + addedViewIdx]);
      }
    }
  }

  int sliderW1 = ui.sliderFixedW->value();
  int sliderW2 = ui.sliderMovingW->value();

  int sliderL1 = ui.sliderFixedL->value();
  int sliderL2 = ui.sliderMovingL->value();

  m_YKDisp[0].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);
  m_YKDisp[1].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);
  m_YKDisp[2].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);

  ui.labelOverlapWnd1->SetBaseImage(&m_YKDisp[0]);
  ui.labelOverlapWnd2->SetBaseImage(&m_YKDisp[1]);
  ui.labelOverlapWnd3->SetBaseImage(&m_YKDisp[2]);

  // here gPMC results could be checked for and displayed, possibly with
  // modification to the qyklabel class /AGA 02/08/2017
  if (m_cbctregistration->dose_loaded) {
    m_AGDisp_Overlay[0].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);
    m_AGDisp_Overlay[1].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);
    m_AGDisp_Overlay[2].FillPixMapDual(sliderL1, sliderL2, sliderW1, sliderW2);

    ui.labelOverlapWnd1->SetOverlayImage(&m_AGDisp_Overlay[0]);
    ui.labelOverlapWnd2->SetOverlayImage(&m_AGDisp_Overlay[1]);
    ui.labelOverlapWnd3->SetOverlayImage(&m_AGDisp_Overlay[2]);
  }

  ui.labelOverlapWnd1->update();
  ui.labelOverlapWnd2->update();
  ui.labelOverlapWnd3->update();
}

void DlgRegistration::SLT_UpdateSplit1() // Mouse Move event
{
  int idx = 0;
  UpdateSplit(idx, ui.labelOverlapWnd1);
}

void DlgRegistration::SLT_UpdateSplit2() {
  int idx = 1;
  UpdateSplit(idx, ui.labelOverlapWnd2);
}
void DlgRegistration::SLT_UpdateSplit3() {
  int idx = 2;
  UpdateSplit(idx, ui.labelOverlapWnd3);
}

void DlgRegistration::UpdateSplit(int viewIdx, qyklabel *pOverlapWnd) {
  int idx = viewIdx;

  if (pOverlapWnd == nullptr) {
    return;
  }

  if (m_YKDisp[idx].IsEmpty()) {
    return;
  }

  if (!m_bPressedLeft[idx] && !m_bPressedRight[idx]) {
    return;
  }

  double dspWidth = pOverlapWnd->width();
  double dspHeight = pOverlapWnd->height();

  int dataWidth = m_YKDisp[idx].m_iWidth;
  int dataHeight = m_YKDisp[idx].m_iHeight;
  if (dataWidth * dataHeight == 0) {
    return;
  }

  int dataX = pOverlapWnd->GetDataPtFromMousePos().x();
  int dataY = pOverlapWnd->GetDataPtFromMousePos().y();

  // only works when while the left Mouse is being clicked
  if (m_bPressedLeft[idx]) {
    QPoint xy_point = QPoint(dataX, dataY);
    m_YKDisp[idx].SetSplitCenter(xy_point);
    /* if (cond) {
     * calculateAngularWEPL(xy_point);
     * SLT_DrawWEPLContour(); // Radial (Spiderweb) plot
     * }
     */
    SLT_DrawImageInFixedSlice();
  } else if (m_bPressedRight[idx] && ui.checkBoxPan->isChecked()) {
    ////Update offset information of dispImage

    // GetOriginalDataPos (PanStart)
    // offset should be 0.. only relative distance matters. offset is in
    // realtime changing
    QPoint ptDataPanStartRel = pOverlapWnd->View2DataExt(
        m_ptPanStart, static_cast<int>(dspWidth), static_cast<int>(dspHeight),
        dataWidth, dataHeight, QPoint(0, 0), m_YKDisp[idx].m_fZoom);

    QPoint ptDataPanEndRel = pOverlapWnd->View2DataExt(
        QPoint(pOverlapWnd->x, pOverlapWnd->y), static_cast<int>(dspWidth),
        static_cast<int>(dspHeight), dataWidth, dataHeight, QPoint(0, 0),
        m_YKDisp[idx].m_fZoom);

    // int dspOffsetX = pOverlapWnd->x - m_ptPanStart.x();
    // int dspOffsetY = m_ptPanStart.y() - pOverlapWnd->y;

    /*QPoint ptDataStart= pOverlapWnd->GetDataPtFromViewPt(m_ptPanStart.x(),
    m_ptPanStart.y()); QPoint ptDataEnd=
    pOverlapWnd->GetDataPtFromViewPt(pOverlapWnd->x, pOverlapWnd->y);*/

    int curOffsetX = ptDataPanEndRel.x() - ptDataPanStartRel.x();
    int curOffsetY = ptDataPanEndRel.y() - ptDataPanStartRel.y();

    int prevOffsetX = m_ptTmpOriginalDataOffset.x();
    int prevOffsetY = m_ptTmpOriginalDataOffset.y();

    // double fZoom = m_YKDisp[idx].m_fZoom;

    m_YKDisp[idx].SetOffset(prevOffsetX - curOffsetX, prevOffsetY - curOffsetY);

    SLT_DrawImageInFixedSlice();
  } else if (m_bPressedRight[idx] &&
             !ui.checkBoxPan->isChecked()) // Window Level
  {
    double wWidth = 2.0;
    double wLevel = 2.0;

    auto iAddedWidth =
        static_cast<int>((m_ptWindowLevelStart.y() - pOverlapWnd->y) * wWidth);
    auto iAddedLevel =
        static_cast<int>((pOverlapWnd->x - m_ptWindowLevelStart.x()) * wLevel);

    // Which image is clicked first??
    QPoint crntDataPt = pOverlapWnd->GetDataPtFromViewPt(
        m_ptWindowLevelStart.x(), m_ptWindowLevelStart.y());

    if (pOverlapWnd->m_pYK16Image != nullptr) {
      if (m_YKDisp[idx].isPtInFirstImage(crntDataPt.x(), crntDataPt.y())) {
        // ui.sliderFixedW->setValue(ui.sliderFixedW->value() + iAddedWidth);
        // //SLT_DrawImageInFixedSlice will be called
        // ui.sliderFixedL->setValue(ui.sliderFixedL->value() + iAddedLevel);
        ui.sliderFixedW->setValue(
            m_iTmpOriginalW +
            iAddedWidth); // SLT_DrawImageInFixedSlice will be called
        ui.sliderFixedL->setValue(m_iTmpOriginalL + iAddedLevel);
      } else {
        ui.sliderMovingW->setValue(m_iTmpOriginalW + iAddedWidth);
        ui.sliderMovingL->setValue(m_iTmpOriginalL + iAddedLevel);
      }
    }
  }
}

// Slide change by scrolling
void DlgRegistration::SLT_MouseWheelUpdate1() {
  if (ui.checkBoxZoom->isChecked()) {
    double oldZoom = ui.labelOverlapWnd1->m_pYK16Image->m_fZoom;

    double fWeighting = 0.2;

    ui.labelOverlapWnd1->m_pYK16Image->SetZoom(
        oldZoom + ui.labelOverlapWnd1->m_iMouseWheelDelta * fWeighting);
    this->SLT_DrawImageInFixedSlice();
  } else {
    ui.sliderPosDisp1->setValue(ui.sliderPosDisp1->value() +
                                ui.labelOverlapWnd1->m_iMouseWheelDelta);
  }
}

void DlgRegistration::SLT_MouseWheelUpdate2() {
  if (ui.checkBoxZoom->isChecked()) {
    double oldZoom = ui.labelOverlapWnd2->m_pYK16Image->m_fZoom;
    double fWeighting = 0.2;

    ui.labelOverlapWnd2->m_pYK16Image->SetZoom(
        oldZoom + ui.labelOverlapWnd2->m_iMouseWheelDelta * fWeighting);
    this->SLT_DrawImageInFixedSlice();

  } else {
    ui.sliderPosDisp2->setValue(ui.sliderPosDisp2->value() +
                                ui.labelOverlapWnd2->m_iMouseWheelDelta);
  }
}

void DlgRegistration::SLT_MouseWheelUpdate3() {
  if (ui.checkBoxZoom->isChecked()) {
    double oldZoom = ui.labelOverlapWnd3->m_pYK16Image->m_fZoom;
    double fWeighting = 0.2;

    ui.labelOverlapWnd3->m_pYK16Image->SetZoom(
        oldZoom + ui.labelOverlapWnd3->m_iMouseWheelDelta * fWeighting);
    this->SLT_DrawImageInFixedSlice();

  } else {
    ui.sliderPosDisp3->setValue(ui.sliderPosDisp3->value() +
                                ui.labelOverlapWnd3->m_iMouseWheelDelta);
  }
}

// release everything
void DlgRegistration::SLT_CancelMouseAction() {
  for (int i = 0; i < 3; i++) {
    m_bPressedLeft[i] = false; // Left Mouse Pressed but not released
    m_bPressedRight[i] = false;
  }
}

void DlgRegistration::SLT_MousePressedLeft1() { m_bPressedLeft[0] = true; }

void DlgRegistration::SLT_MousePressedLeft2() { m_bPressedLeft[1] = true; }

void DlgRegistration::SLT_MousePressedLeft3() { m_bPressedLeft[2] = true; }

void DlgRegistration::MousePressedRight(int wndIdx, qyklabel *pWnd) {
  m_bPressedRight[wndIdx] = true;

  if (ui.checkBoxPan->isChecked()) {
    m_ptPanStart.setX(pWnd->x);
    m_ptPanStart.setY(pWnd->y);

    m_ptTmpOriginalDataOffset.setX(m_YKDisp[wndIdx].m_iOffsetX);
    m_ptTmpOriginalDataOffset.setY(m_YKDisp[wndIdx].m_iOffsetY);
  } else {
    m_ptWindowLevelStart.setX(pWnd->x);
    m_ptWindowLevelStart.setY(pWnd->y);

    QPoint crntDataPt = pWnd->GetDataPtFromViewPt(m_ptWindowLevelStart.x(),
                                                  m_ptWindowLevelStart.y());

    if (m_YKDisp[wndIdx].isPtInFirstImage(crntDataPt.x(), crntDataPt.y())) {
      m_iTmpOriginalL = ui.sliderFixedL->value();
      m_iTmpOriginalW = ui.sliderFixedW->value();
    } else {
      m_iTmpOriginalL = ui.sliderMovingL->value();
      m_iTmpOriginalW = ui.sliderMovingW->value();
    }
  }
}

void DlgRegistration::SLT_MousePressedRight1() {
  int idx = 0;
  MousePressedRight(idx, ui.labelOverlapWnd1);
}

void DlgRegistration::SLT_MousePressedRight2() {
  int idx = 1;
  MousePressedRight(idx, ui.labelOverlapWnd2);
}

void DlgRegistration::SLT_MousePressedRight3() {
  int idx = 2;
  MousePressedRight(idx, ui.labelOverlapWnd3);
}

void DlgRegistration::SLT_MouseReleasedLeft1() { m_bPressedLeft[0] = false; }

void DlgRegistration::SLT_MouseReleasedLeft2() { m_bPressedLeft[1] = false; }

void DlgRegistration::SLT_MouseReleasedLeft3() { m_bPressedLeft[2] = false; }

void DlgRegistration::SLT_MouseReleasedRight1() { m_bPressedRight[0] = false; }

void DlgRegistration::SLT_MouseReleasedRight2() { m_bPressedRight[1] = false; }

void DlgRegistration::SLT_MouseReleasedRight3() { m_bPressedRight[2] = false; }

void DlgRegistration::SLT_ChangeView() {
  m_enViewArrange = (m_enViewArrange + 1) % 3;

  YK16GrayImage tmpBufYK[3];

  for (int i = 0; i < 3; i++) {
    tmpBufYK[i].CloneImage(m_YKDisp[i]);
  }

  for (int i = 0; i < 3; i++) {
    int nextIdx = (i + 1) % 3;
    m_YKDisp[i].CloneImage(tmpBufYK[nextIdx]);
  }
  // /*if (m_enViewArrange > 2)
  // m_enViewArrange = m_enViewArrange-3;*/

  // //Split center should be transfered as well.
  // for (int i = 0 ; i<3 ; i++)
  // {
  // int nextIdx = (i+1)%3;
  ///*if (nextIdx >= 3)
  //  nextIdx = nextIdx -3;*/
  // m_YKDisp[i].SetSplitCenter(m_YKDisp[i+1].m_ptSplitCenter);
  ////no need of data copy, but other display information (zoom and others
  /// should be copied here.
  // }
  initOverlapWndSize();
  shiftSliceSlider();
  updateSliceLabel();

  SLT_DrawImageWhenSliceChange(); // only image data will be updated. Zoom and
                                  // other things are not.
  // SLT_DrawImageInFixedSlice();
}

void DlgRegistration::initOverlapWndSize() {
  ui.labelOverlapWnd1->setFixedWidth(DEFAULT_LABEL_SIZE1);
  ui.labelOverlapWnd1->setFixedHeight(DEFAULT_LABEL_SIZE1);

  ui.labelOverlapWnd2->setFixedWidth(DEFAULT_LABEL_SIZE2);
  ui.labelOverlapWnd2->setFixedHeight(DEFAULT_LABEL_SIZE2);

  ui.labelOverlapWnd3->setFixedWidth(DEFAULT_LABEL_SIZE2);
  ui.labelOverlapWnd3->setFixedHeight(DEFAULT_LABEL_SIZE2);
}

void DlgRegistration::shiftSliceSlider() // shift one slice slider information
{
  // to avoid first unnecessary action.
  disconnect(ui.sliderPosDisp1, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageWhenSliceChange()));
  disconnect(ui.sliderPosDisp2, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageWhenSliceChange()));
  disconnect(ui.sliderPosDisp3, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawImageWhenSliceChange()));

  int crntMin[3];
  int crntMax[3];
  int crntValue[3];

  int newMin[3];
  int newMax[3];
  int newValue[3];

  crntMin[0] = ui.sliderPosDisp1->minimum();
  crntMin[1] = ui.sliderPosDisp2->minimum();
  crntMin[2] = ui.sliderPosDisp3->minimum();

  crntMax[0] = ui.sliderPosDisp1->maximum();
  crntMax[1] = ui.sliderPosDisp2->maximum();
  crntMax[2] = ui.sliderPosDisp3->maximum();

  crntValue[0] = ui.sliderPosDisp1->value();
  crntValue[1] = ui.sliderPosDisp2->value();
  crntValue[2] = ui.sliderPosDisp3->value();

  for (int i = 0; i < 3; i++) {
    int newIdx = (i + 1) % 3;
    newMin[i] = crntMin[newIdx];
    newMax[i] = crntMax[newIdx];
    newValue[i] = crntValue[newIdx];
  }

  ui.sliderPosDisp1->setMinimum(newMin[0]);
  ui.sliderPosDisp2->setMinimum(newMin[1]);
  ui.sliderPosDisp3->setMinimum(newMin[2]);

  ui.sliderPosDisp1->setMaximum(newMax[0]);
  ui.sliderPosDisp2->setMaximum(newMax[1]);
  ui.sliderPosDisp3->setMaximum(newMax[2]);

  ui.sliderPosDisp1->setValue(newValue[0]);
  ui.sliderPosDisp2->setValue(newValue[1]);
  ui.sliderPosDisp3->setValue(newValue[2]);

  connect(ui.sliderPosDisp1, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageWhenSliceChange()));
  connect(ui.sliderPosDisp2, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageWhenSliceChange()));
  connect(ui.sliderPosDisp3, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawImageWhenSliceChange()));
}

void DlgRegistration::updateSliceLabel() {

  switch (m_enViewArrange) {
  case AXIAL_FRONTAL_SAGITTAL:
    ui.labelDisp1->setText("AXIAL");
    ui.labelDisp2->setText("FRONTAL");
    ui.labelDisp3->setText("SAGITTAL");
    break;

  case FRONTAL_SAGITTAL_AXIAL:
    ui.labelDisp1->setText("FRONTAL");
    ui.labelDisp2->setText("SAGITTAL");
    ui.labelDisp3->setText("AXIAL");
    break;

  case SAGITTAL_AXIAL_FRONTAL:
    ui.labelDisp1->setText("SAGITTAL");
    ui.labelDisp2->setText("AXIAL");
    ui.labelDisp3->setText("FRONTAL");
    break;
  default:
    std::cerr << "WTF!?" << std::endl;
    break;
  }
}

void DlgRegistration::LoadImgFromComboBox(
    int idx,
    QString
        &strSelectedComboTxt) // -->when fixed image loaded will be called here!
{
  // std::cout << "LoadImgFromComboBox " << "index " << idx << "text " <<
  // strSelectedComboTxt.toLocal8Bit().constData() << std::endl;

  UShortImageType::Pointer spTmpImg;
  if (strSelectedComboTxt.compare(QString("RAW_CBCT"), Qt::CaseSensitive) ==
      0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spRawReconImg;
  } else if (strSelectedComboTxt.compare(QString("REF_CT"),
                                         Qt::CaseSensitive) == 0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spRefCTImg;
  } else if (strSelectedComboTxt.compare(QString("MANUAL_RIGID_CT"),
                                         Qt::CaseSensitive) == 0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spManualRigidCT;
  } else if (strSelectedComboTxt.compare(QString("AUTO_RIGID_CT"),
                                         Qt::CaseSensitive) == 0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spAutoRigidCT;
  } else if (strSelectedComboTxt.compare(QString("DEFORMED_CT1"),
                                         Qt::CaseSensitive) == 0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spDeformedCT1;
  } else if (strSelectedComboTxt.compare(QString("DEFORMED_CT2"),
                                         Qt::CaseSensitive) == 0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spDeformedCT2;
  } else if (strSelectedComboTxt.compare(QString("DEFORMED_CT3"),
                                         Qt::CaseSensitive) == 0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spDeformedCT3;
  } else if (strSelectedComboTxt.compare(QString("DEFORMED_CT_FINAL"),
                                         Qt::CaseSensitive) == 0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spDeformedCT_Final;
  } else if (strSelectedComboTxt.compare(QString("COR_CBCT"),
                                         Qt::CaseSensitive) == 0) {
    spTmpImg = m_cbctregistration->m_pParent->m_spScatCorrReconImg;
  }

  if (spTmpImg == nullptr) {
    // std::cout << "Selected image is not ready: " <<
    // strSelectedComboTxt.toLocal8Bit().constData() << std::endl;
    return;
  }

  if (idx == 0) {
    m_spFixed = spTmpImg;

    whenFixedImgLoaded();
  } else if (idx == 1) {
    m_spMoving = spTmpImg;
    // std::cout << "idx: " << idx << "m_spMoving"  << m_spMoving << std::endl;
    whenMovingImgLoaded();
  }

  SLT_DrawImageWhenSliceChange();
}

void DlgRegistration::LoadVOIFromComboBox(int idx,
                                          QString &strSelectedComboTxt) {

  ctType ct_type = PLAN_CT;
  auto ct = ui.comboBoxImgMoving->currentText().toStdString();
  if (ct == std::string("REF_CT")) {
    // ct_type = PLAN_CT;
  } else if (ct == std::string("AUTO_RIGID_CT")) {
    ct_type = RIGID_CT;
  } else if (ct == std::string("DEFORMED_CT_FINAL")) {
    ct_type = DEFORM_CT;
  } else {
    std::cout << "This moving image does not own any VOIs" << std::endl;
    return;
  }

  auto struct_set =
      m_cbctregistration->m_pParent->m_structures->get_ss(ct_type);
  if (struct_set == nullptr) {
    return;
  }
  if (struct_set->slist.empty()) {
    std::cerr << "Structures not initialized yet" << std::endl;
    return;
  }

  for (auto voi : struct_set->slist) {
    if (strSelectedComboTxt.compare(voi.name.c_str()) == 0) {
      m_cbctregistration->cur_voi = std::make_unique<Rtss_roi_modern>(voi);
    }
  }

  SLT_DrawImageWhenSliceChange();
}

// search  for the  main data, if there  is, add  the predefined name to the
// combobox
void DlgRegistration::UpdateListOfComboBox(int idx) {
  QComboBox *crntCombo;

  if (idx == 0) {
    crntCombo = ui.comboBoxImgFixed;
  } else {
    crntCombo = ui.comboBoxImgMoving;
  }

  // remove all the list
  crntCombo->clear();

  if (m_cbctregistration->m_pParent->m_spRawReconImg != nullptr) {
    crntCombo->addItem("RAW_CBCT");
  }

  if (m_cbctregistration->m_pParent->m_spRefCTImg != nullptr) {
    crntCombo->addItem("REF_CT");
  }

  if (m_cbctregistration->m_pParent->m_spManualRigidCT != nullptr) {
    crntCombo->addItem("MANUAL_RIGID_CT");
  }

  if (m_cbctregistration->m_pParent->m_spAutoRigidCT != nullptr) {
    crntCombo->addItem("AUTO_RIGID_CT");
  }

  if (m_cbctregistration->m_pParent->m_spDeformedCT1 != nullptr) {
    crntCombo->addItem("DEFORMED_CT1");
  }

  if (m_cbctregistration->m_pParent->m_spDeformedCT2 != nullptr) {
    crntCombo->addItem("DEFORMED_CT2");
  }

  if (m_cbctregistration->m_pParent->m_spDeformedCT3 != nullptr) {
    crntCombo->addItem("DEFORMED_CT3");
  }

  if (m_cbctregistration->m_pParent->m_spDeformedCT_Final != nullptr) {
    crntCombo->addItem("DEFORMED_CT_FINAL");
  }

  if (m_cbctregistration->m_pParent->m_spScatCorrReconImg != nullptr) {
    crntCombo->addItem("COR_CBCT");
  }
}

void DlgRegistration::SLT_RestoreImageSingle() {
  int mainWndIdx = 0;
  m_YKDisp[mainWndIdx].SetZoom(1.0);
  m_YKDisp[mainWndIdx].SetOffset(0, 0);

  SLT_DrawImageInFixedSlice();
}

void DlgRegistration::SLT_RestoreImageAll() {
  for (auto &i : m_cbctregistration->m_YKDisp) {
    i.SetZoom(1.0);
    i.SetOffset(0, 0);

    SLT_DrawImageInFixedSlice();
  }
}
// Auto registration
void DlgRegistration::SLT_DoRegistrationRigid() // plastimatch auto registration
{
  // 1) Save current image files
  if ((m_spFixed == nullptr) || (m_spMoving == nullptr)) {
    return;
  }

  // Before the registration, remove background and fill the air-bubble of the
  // cbct  1) calculate manual shift value  PreProcess_CBCT(); // remove skin
  // and fill bubble before registration

  if ((m_pParent->m_cbctrecon->m_spRefCTImg == nullptr) ||
      (m_pParent->m_cbctrecon->m_spManualRigidCT == nullptr)) {
    return;
  }

  if (ui.checkBoxKeyMoving->isChecked()) {
    SLT_KeyMoving(false);
  }

  /*- Make a synthetic std::vector field according to the translation
        plastimatch synth-vf --fixed [msk_skin.mha] --output
  [xf_manual_trans.mha] --xf-trans "[origin diff (raw - regi)]" plastimatch
  synth-vf --fixed E:\PlastimatchData\DicomEg\OLD\msk_skin.mha --output
  E:\PlastimatchData\DicomEg\OLD\xf_manual_trans.mha --xf-trans "21 -29 1"

        - Move the skin contour according to the std::vector field
        plastimatch warp --input [msk_skin.mha] --output-img
  [msk_skin_manRegi.mha] --xf [xf_manual_trans.mha] plastimatch warp --input
  E:\PlastimatchData\DicomEg\OLD\msk_skin.mha --output-img
  E:\PlastimatchData\DicomEg\OLD\msk_skin_manRegi.mha --xf
  E:\PlastimatchData\DicomEg\OLD\xf_manual_trans.mha*/

  bool bPrepareMaskOnly; // prepare mask but not apply

  bPrepareMaskOnly = !ui.checkBoxCropBkgroundCBCT->isChecked();

  std::cout << "1: writing temporary files" << std::endl;

  // Both image type: Unsigned Short
  QString filePathFixed =
      m_cbctregistration->m_strPathPlastimatch + "/" + "fixed_rigid.mha";
  QString filePathMoving =
      m_cbctregistration->m_strPathPlastimatch + "/" + "moving_rigid.mha";
  QString filePathOutput =
      m_cbctregistration->m_strPathPlastimatch + "/" + "output_rigid.mha";
  QString filePathXform =
      m_cbctregistration->m_strPathPlastimatch + "/" + "xform_rigid.txt";
  QString filePathROI = m_cbctregistration->m_strPathPlastimatch + "/" +
                        "fixed_roi_rigid.mha"; // optional

  using writerType = itk::ImageFileWriter<UShortImageType>;

  writerType::Pointer writer1 = writerType::New();
  writer1->SetFileName(filePathFixed.toLocal8Bit().constData());
  writer1->SetUseCompression(true);
  writer1->SetInput(m_spFixed);

  writerType::Pointer writer2 = writerType::New();
  writer2->SetFileName(filePathMoving.toLocal8Bit().constData());
  writer2->SetUseCompression(true);
  writer2->SetInput(m_spMoving);

  writer1->Update();
  writer2->Update();

  if (ui.checkBoxUseROIForRigid->isChecked()) {
    std::cout << "Creating a ROI mask for Rigid registration " << std::endl;
    QString strFOVGeom = ui.lineEditFOVPos->text();

    QStringList strListFOV = strFOVGeom.split(",");
    if (strListFOV.count() == 3) {
      const auto FOV_DcmPosX = strListFOV.at(0).toFloat(); // mm
      const auto FOV_DcmPosY = strListFOV.at(1).toFloat();
      const auto FOV_Radius = strListFOV.at(2).toFloat();

      // Create Image using FixedImage sp

      // Image Pointer here
      UShortImageType::Pointer spRoiMask;
      AllocateByRef(m_spFixed, spRoiMask);
      m_pParent->m_cbctrecon->GenerateCylinderMask(spRoiMask, FOV_DcmPosX,
                                                   FOV_DcmPosY, FOV_Radius);

      writerType::Pointer writer3 = writerType::New();
      writer3->SetFileName(filePathROI.toLocal8Bit().constData());
      writer3->SetUseCompression(true);
      writer3->SetInput(spRoiMask);
      writer3->Update();
    }
  }

  // Preprocessing
  // 2) move CT-based skin mask on CBCT based on manual shift
  //  if (m_strPathCTSkin)

  // m_strPathCTSkin is prepared after PreProcessCT()

  // QString filePathFixed_proc = filePathFixed;

  using readerType = itk::ImageFileReader<UShortImageType>;

  QString strPathOriginalCTSkinMask;

  QString filePathFixed_proc =
      m_cbctregistration->m_strPathPlastimatch + "/" +
      "fixed_rigid_proc.mha"; // After autoRigidbody Regi

  QString strPath_mskSkinCT_manRegi_exp =
      m_cbctregistration->m_strPathPlastimatch +
      "/msk_skin_CT_manRegi_exp.mha"; // proof of preproceesing for CBCT

  QFileInfo finfoFixedProc = QFileInfo(filePathFixed_proc);
  QFileInfo finfoManMask = QFileInfo(strPath_mskSkinCT_manRegi_exp);

  // if this file already exists, no need of proprocess for CBCT (skin cropping)
  // is needed.
  // if skin cropping is not done for this case (supposed to be done during
  // confirm manual regi)

  // if (!finfoFixedProc.exists() && !finfoManMask.exists() &&
  // ui.checkBoxCropBkgroundCBCT->isChecked())
  if (!finfoFixedProc.exists() && !finfoManMask.exists()) {
    std::cout << "Preprocessing for CBCT is not done. It is being done here "
                 "before rigid body registration"
              << std::endl;

    UShortImageType::PointType originBefore =
        m_pParent->m_cbctrecon->m_spRefCTImg->GetOrigin();
    UShortImageType::PointType originAfter =
        m_pParent->m_cbctrecon->m_spManualRigidCT->GetOrigin();

    double fShift[3];
    fShift[0] = (originBefore[0] - originAfter[0]); // DICOM
    fShift[1] = (originBefore[1] - originAfter[1]);
    fShift[2] = (originBefore[2] - originAfter[2]);

    QFileInfo finfoSkinFile1 = QFileInfo(m_cbctregistration->m_strPathCTSkin);
    QString strPathAlternateSkin =
        m_cbctregistration->m_strPathPlastimatch + "/" + "msk_skin_CT.mha";
    QFileInfo finfoSkinFile2 = QFileInfo(strPathAlternateSkin);

    //&& ui.checkBoxCropBkgroundCBCT->isChecked()
    double skinExp = ui.lineEditCBCTSkinCropBfRegid->text().toDouble();
    int bkGroundValUshort = ui.lineEditBkFillCBCT->text().toInt(); // 0

    if (finfoSkinFile1.exists()) {
      strPathOriginalCTSkinMask = m_cbctregistration->m_strPathCTSkin;
      // This was OK.
      m_cbctregistration->ProcessCBCT_beforeAutoRigidRegi(
          filePathFixed, strPathOriginalCTSkinMask, filePathFixed_proc, fShift,
          bPrepareMaskOnly, skinExp, bkGroundValUshort);

      if (bPrepareMaskOnly) { // currently, filePathFixed_proc == "";
        filePathFixed_proc = filePathFixed;
      }
    } else if (finfoSkinFile2.exists()) {
      std::cout << "alternative skin file will be used" << std::endl;
      strPathOriginalCTSkinMask = strPathAlternateSkin;
      m_cbctregistration->ProcessCBCT_beforeAutoRigidRegi(
          filePathFixed, strPathOriginalCTSkinMask, filePathFixed_proc, fShift,
          bPrepareMaskOnly, skinExp, bkGroundValUshort);

      if (bPrepareMaskOnly) {
        filePathFixed_proc = filePathFixed;
      }
    } else {
      std::cout << "CT skin file(msk_skin_CT.mha) is not found. preprocessing "
                   "will be skipped"
                << std::endl;
      filePathFixed_proc = filePathFixed;
    }

    QFileInfo fixedProcInfo(filePathFixed_proc);
    if (!fixedProcInfo.exists()) {
      std::cout << "Error! proc file doesn't exist!" << std::endl;
      return;
    }

    readerType::Pointer reader2 = readerType::New();
    reader2->SetFileName(filePathFixed_proc.toLocal8Bit().constData());
    reader2->Update();
    m_pParent->m_cbctrecon->m_spRawReconImg = reader2->GetOutput();

    double tmpSkinMargin = ui.lineEditCBCTSkinCropBfRegid->text().toDouble();
    QString update_message =
        QString("Skin removed CBCT with margin %1 mm").arg(tmpSkinMargin);
    m_pParent->UpdateReconImage(m_pParent->m_cbctrecon->m_spRawReconImg,
                                update_message);
  }

  if (!finfoFixedProc.exists()) {
    filePathFixed_proc = filePathFixed; //"fixed_rigid.mha";
  }

  std::cout << "2: Creating a plastimatch command file" << std::endl;

  QString fnCmdRegisterRigid = "cmd_register_rigid.txt";
  // QString fnCmdRegisterDeform = "cmd_register_deform.txt";
  QString pathCmdRegister =
      m_cbctregistration->m_strPathPlastimatch + "/" + fnCmdRegisterRigid;

  // GenPlastiRegisterCommandFile(pathCmdRegister, filePathFixed_proc,
  // filePathMoving,
  //						filePathOutput, filePathXform,
  // PLAST_RIGID,
  //"","","");

  QString strDummy = "";

  bool mse = ui.radioButton_mse->isChecked();
  bool cuda = m_pParent->ui.radioButton_UseCUDA->isChecked();
  QString GradOptionStr = ui.lineEditGradOption->text();
  // For Cropped patients, FOV mask is applied.
  if (ui.checkBoxUseROIForRigid->isChecked()) {
    std::cout << "2.A:  ROI-based Rigid body registration will be done"
              << std::endl;
    m_cbctregistration->GenPlastiRegisterCommandFile(
        pathCmdRegister, filePathFixed_proc, filePathMoving, filePathOutput,
        filePathXform, PLAST_RIGID, strDummy, strDummy, strDummy, filePathROI,
        mse, cuda, GradOptionStr);
  } else {
    m_cbctregistration->GenPlastiRegisterCommandFile(
        pathCmdRegister, filePathFixed_proc, filePathMoving, filePathOutput,
        filePathXform, PLAST_RIGID, strDummy, strDummy, strDummy, strDummy, mse,
        cuda, GradOptionStr);
  }

  std::string str_command_filepath = pathCmdRegister.toLocal8Bit().constData();

  m_cbctregistration->CallingPLMCommand(str_command_filepath);

  std::cout << "5: Reading output image-CT" << std::endl;

  readerType::Pointer reader = readerType::New();
  reader->SetFileName(filePathOutput.toLocal8Bit().constData());
  reader->Update();
  m_pParent->m_cbctrecon->m_spAutoRigidCT = reader->GetOutput();

  std::cout << "6: Reading is completed" << std::endl;

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_AUTO_RIGID);

  m_cbctregistration->m_strPathXFAutoRigid = filePathXform; // for further use
}

bool DlgRegistration::eventFilter(QObject *target, QEvent *event) {
  if (event->type() == QEvent::Paint) { // eventy->type() = 12 if paint event
    return false;
  }

  // int aa = (int)(event->type());

  // std::cout << target << std::endl;
  // std::cout <<"Event Type " << aa << std::endl;

  if (event->type() == QEvent::ShortcutOverride) // 51
  {
    if (!ui.checkBoxKeyMoving->isChecked()) {
      return false;
    }

    auto *keyEvent = dynamic_cast<QKeyEvent *>(event);

    auto iArrowKey = static_cast<int>(keyEvent->nativeVirtualKey());
    // std::cout << iArrowKey << std::endl;
    double resol = ui.lineEditMovingResol->text().toDouble(); // mm

    switch (iArrowKey) {
      // case Qt::Key_Left: //37
    case 37: // 37
      ImageManualMove(37, resol);
      break;
    case 38:                      // 38
      ImageManualMove(38, resol); // up
      break;
    case 39: // 39
      ImageManualMove(39, resol);
      break;
    case 40: // 40
      ImageManualMove(40, resol);
      break;
    case 33: //
      ImageManualMove(33, resol);
      break;
    case 34: //
      ImageManualMove(34, resol);
      break;
    default:
      break; // do nothing
    }
    // std::cout << keyEvent->->nativeVirtualKey() << std::endl;
    // event->accept();
    return true; // No more event transfering?
  }
  /*else if (event->type() == QEvent::KeyPress)
  {
        std::cout << "key pressed" << std::endl;
  }
  else if (event->type() == QEvent::KeyRelease)
  {
        std::cout << "key Released" << std::endl;
  }*/

  // return true;
  return QDialog::eventFilter(target,
                              event); // This is mandatory to deliver
                                      // the event signal to the children
                                      // components
}

void DlgRegistration::childEvent(QChildEvent *event) {
  // std::cout << "childEvent Called" << std::endl;
  if (event->added()) {
    event->child()->installEventFilter(this);
  }
}

void DlgRegistration::ImageManualMove(int direction, double resol) {
  if (m_spMoving == nullptr) {
    return;
  }

  // USHORT_ImageType::SizeType imgSize =
  // m_spMoving->GetRequestedRegion().GetSize(); //1016x1016 x z
  UShortImageType::PointType imgOrigin = m_spMoving->GetOrigin();
  // USHORT_ImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();

  if (direction == 37) {                 // LEFT
    imgOrigin[0] = imgOrigin[0] - resol; // mm unit
  } else if (direction == 39) {          // RIGHT
    imgOrigin[0] = imgOrigin[0] + resol;
  } else if (direction == 38) { // Up
    imgOrigin[1] = imgOrigin[1] - resol;
  } else if (direction == 40) { // DOWN
    imgOrigin[1] = imgOrigin[1] + resol;
  } else if (direction == 33) { // PageUp
    imgOrigin[2] = imgOrigin[2] + resol;
  } else if (direction == 34) { // PageDown
    imgOrigin[2] = imgOrigin[2] - resol;
  }

  m_spMoving->SetOrigin(imgOrigin);

  SLT_DrawImageWhenSliceChange();

  // Display relative movement
  // Starting point? RefCT image
  // Only Valid when Moving image is the ManualMove
  UShortImageType::PointType imgOriginRef =
      m_pParent->m_cbctrecon->m_spRefCTImg->GetOrigin();

  QString strDelta;
  strDelta.sprintf(
      "delta(mm): %3.1f, %3.1f, %3.1f", (imgOrigin[0] - imgOriginRef[0]),
      (imgOrigin[1] - imgOriginRef[1]), (imgOrigin[2] - imgOriginRef[2]));
  ui.lineEditOriginChanged->setText(strDelta);

  // std::cout << "direction  " << direction << std::endl;
  // std::cout << "resolution " << resol << std::endl;
}

// change origin of moving image by shift value acquired from gradient
// registration
void DlgRegistration::ImageManualMoveOneShot(float shiftX, float shiftY,
                                             float shiftZ) // DICOM coordinate
{
  if (m_spMoving == nullptr) {
    return;
  }

  // USHORT_ImageType::SizeType imgSize =
  // m_spMoving->GetRequestedRegion().GetSize(); //1016x1016 x z
  using USPointType = UShortImageType::PointType;
  USPointType imgOrigin = m_spMoving->GetOrigin();
  // USHORT_ImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();
  imgOrigin[0] = imgOrigin[0] - static_cast<USPointType::ValueType>(shiftX);
  imgOrigin[1] = imgOrigin[1] - static_cast<USPointType::ValueType>(shiftY);
  imgOrigin[2] = imgOrigin[2] - static_cast<USPointType::ValueType>(shiftZ);

  m_spMoving->SetOrigin(imgOrigin);

  SLT_DrawImageWhenSliceChange();

  // Display relative movement
  // Starting point? RefCT image
  // Only Valid when Moving image is the ManualMove
  UShortImageType::PointType imgOriginRef =
      m_pParent->m_cbctrecon->m_spRefCTImg->GetOrigin();

  QString strDelta;
  strDelta.sprintf(
      "delta(mm): %3.1f, %3.1f, %3.1f", (imgOrigin[0] - imgOriginRef[0]),
      (imgOrigin[1] - imgOriginRef[1]), (imgOrigin[2] - imgOriginRef[2]));
  ui.lineEditOriginChanged->setText(strDelta);
}

// Bring Focus
void DlgRegistration::SLT_BringFocusToEnableArrow(bool bChecked) {
  if (bChecked) {
    ui.labelDisp1->setFocus(); // if focus is in label, Key Event will trigger
                               // 51 (override)
  }
}

void DlgRegistration::SLT_KeyMoving(bool bChecked) // Key Moving check box
{
  ui.lineEditMovingResol->setDisabled(bChecked);
  if (bChecked) {
    SelectComboExternal(1,
                        REGISTER_MANUAL_RIGID); // should be
  }
  ui.comboBoxImgFixed->setDisabled(bChecked);
  ui.comboBoxImgMoving->setDisabled(bChecked);
  ui.pushButtonRestoreOriginal->setDisabled(bChecked);
}

void DlgRegistration::AddImageToCombo(
    int comboIdx, enREGI_IMAGES option) // comboIdx 0: fixed, 1: moving
{
  switch (option) {
  case REGISTER_RAW_CBCT:
    if (m_pParent->m_cbctrecon->m_spRawReconImg != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("RAW_CBCT");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("RAW_CBCT");
      }
    }
    break;
  case REGISTER_REF_CT:
    if (m_pParent->m_cbctrecon->m_spRefCTImg != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("REF_CT");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("REF_CT");
      }
    }
    break;
  case REGISTER_MANUAL_RIGID:
    if (m_pParent->m_cbctrecon->m_spManualRigidCT != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("MANUAL_RIGID_CT");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("MANUAL_RIGID_CT");
      }
    }
    break;
  case REGISTER_AUTO_RIGID:
    if (m_pParent->m_cbctrecon->m_spAutoRigidCT != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("AUTO_RIGID_CT");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("AUTO_RIGID_CT");
      }
    }
    break;
  case REGISTER_DEFORM1:
    if (m_pParent->m_cbctrecon->m_spDeformedCT1 != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("DEFORMED_CT1");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("DEFORMED_CT1");
      }
    }
    break;
  case REGISTER_DEFORM2:
    if (m_pParent->m_cbctrecon->m_spDeformedCT2 != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("DEFORMED_CT2");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("DEFORMED_CT2");
      }
    }
    break;
  case REGISTER_DEFORM3:
    if (m_pParent->m_cbctrecon->m_spDeformedCT3 != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("DEFORMED_CT3");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("DEFORMED_CT3");
      }
    }
    break;

  case REGISTER_DEFORM_FINAL:
    if (m_pParent->m_cbctrecon->m_spDeformedCT_Final != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("DEFORMED_CT_FINAL");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("DEFORMED_CT_FINAL");
      }
    }
    break;
  case REGISTER_COR_CBCT:
    if (m_pParent->m_cbctrecon->m_spScatCorrReconImg != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("COR_CBCT");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("COR_CBCT");
      }
    }
    break;
  case REGISTER_DEFORM_SKIP_AUTORIGID:
    if (m_pParent->m_cbctrecon->m_spScatCorrReconImg != nullptr) {
      if (comboIdx == 0) {
        ui.comboBoxImgFixed->addItem("REGISTER_DEFORM_SKIP_AUTORIGID");
      } else if (comboIdx == 1) {
        ui.comboBoxImgMoving->addItem("REGISTER_DEFORM_SKIP_AUTORIGID");
      }
    }
    break;
  }
}

// externally change  combo box value
void DlgRegistration::SelectComboExternal(int idx, enREGI_IMAGES iImage) {
  QComboBox *crntCombo;

  if (idx == 0) {
    crntCombo = ui.comboBoxImgFixed;
  } else if (idx == 1) {
    crntCombo = ui.comboBoxImgMoving;
  } else {
    std::cerr << "What did you do to get here?" << std::endl;
    return;
  }

  int findIndx = -1;
  switch (iImage) {
  case REGISTER_RAW_CBCT:
    findIndx = crntCombo->findText("RAW_CBCT");
    break;
  case REGISTER_REF_CT:
    findIndx = crntCombo->findText("REF_CT");
    break;
  case REGISTER_MANUAL_RIGID:
    findIndx = crntCombo->findText("MANUAL_RIGID_CT");
    break;
  case REGISTER_AUTO_RIGID:
    findIndx = crntCombo->findText("AUTO_RIGID_CT");
    break;
  case REGISTER_DEFORM1:
    findIndx = crntCombo->findText("DEFORMED_CT1");
    break;
  case REGISTER_DEFORM2:
    findIndx = crntCombo->findText("DEFORMED_CT2");
    break;
  case REGISTER_DEFORM3:
    findIndx = crntCombo->findText("DEFORMED_CT3");
    break;
  case REGISTER_DEFORM_FINAL:
    findIndx = crntCombo->findText("DEFORMED_CT_FINAL");
    break;
  case REGISTER_COR_CBCT:
    findIndx = crntCombo->findText("COR_CBCT");
    break;
  case REGISTER_DEFORM_SKIP_AUTORIGID:
    findIndx = crntCombo->findText("REGISTER_DEFORM_SKIP_AUTORIGID");
    break;
  }
  // std::cout << "setCurrentIndx " << findIndx << std::endl;

  if (findIndx < 0) { //-1 if not found
    return;
  }

  crntCombo->setCurrentIndex(findIndx); // will call "SLT_FixedImageSelected"

  // Below are actually redendency. don't know why setCurrentIndex sometimes
  // doesn't trigger the SLT_MovingImageSelected slot func.

  QString crntStr = crntCombo->currentText();
  if (idx == 0) {
    SLT_FixedImageSelected(crntStr);
  } else if (idx == 1) {
    SLT_MovingImageSelected(crntStr);
  }
  // If it is called inside the Dlg, SLT seemed not to conneced
}

void DlgRegistration::SLT_FixedImageSelected(QString selText) {
  // QString strCrntText = ui.comboBoxImgFixed->currentText();
  LoadImgFromComboBox(
      0, selText); // here, m_spMoving and Fixed images are determined
}

void DlgRegistration::SLT_MovingImageSelected(QString selText) {
  // QString strCrntText = ui.comboBoxImgMoving->currentText();
  // std::cout << "SLT_MovingImageSelected" << std::endl;
  LoadImgFromComboBox(1, selText);
}

void DlgRegistration::UpdateVOICombobox(ctType ct_type) {
  auto struct_set =
      m_cbctregistration->m_pParent->m_structures->get_ss(ct_type);
  if (struct_set == nullptr) {
    return;
  }
  if (struct_set->slist.empty()) {
    std::cerr << "Structures not initialized yet" << std::endl;
    return;
  }
  for (auto voi : struct_set->slist) {
    ui.comboBox_VOI->addItem(QString(voi.name.c_str()));
  }
}

void DlgRegistration::SLT_RestoreMovingImg() {
  m_cbctregistration->m_pParent->RegisterImgDuplication(REGISTER_REF_CT,
                                                        REGISTER_MANUAL_RIGID);
  ui.lineEditOriginChanged->setText(QString(""));

  SelectComboExternal(1, REGISTER_MANUAL_RIGID);
}

void DlgRegistration::SLT_PreProcessCT() {
  if (!ui.checkBoxCropBkgroundCT->isChecked()) {
    std::cout << "Preprocessing is not selected." << std::endl;
    return;
  }

  int iAirThresholdShort = ui.lineEditBkDetectCT->text().toInt();

  if (m_cbctregistration->m_pParent->m_strPathPlanCTDir.length() < 3) {
    std::cout
        << "Reference CT DIR should be specified for structure based cropping"
        << std::endl;
    if ((m_spMoving == nullptr) || (m_spFixed == nullptr)) {
      return;
    }
    if (m_spFixed->GetLargestPossibleRegion().GetSize()[0] !=
            m_spMoving->GetLargestPossibleRegion().GetSize()[0] ||
        m_spFixed->GetLargestPossibleRegion().GetSize()[1] !=
            m_spMoving->GetLargestPossibleRegion().GetSize()[1] ||
        m_spFixed->GetLargestPossibleRegion().GetSize()[2] !=
            m_spMoving->GetLargestPossibleRegion().GetSize()[2]) {
      std::cout
          << "Fixed and moving image is not the same size, consider using "
             "a platimatch registration to solve this."
          << std::endl;
      return;
    }

    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(this, "No reference structures found!",
                                  "Do you wan't to attempt an auto correction "
                                  "of air and excessive circumference?",
                                  QMessageBox::Yes | QMessageBox::No);
    if (reply == QMessageBox::Yes) {
      std::cout << "Attempting automatic air filling and skin cropping..."
                << std::endl;
      m_cbctregistration->autoPreprocessCT(iAirThresholdShort, m_spFixed,
                                           m_spMoving);
    }
    return;
  }

  QString strRSName = ui.lineEditCropContourName->text();
  bool fill_bubble = ui.checkBoxFillBubbleCT->isChecked();
  int iBubbleFillingVal =
      ui.lineEditBubFillCT->text().toInt();                   // 0 = soft tissue
  int iAirFillValShort = ui.lineEditBkFillCT->text().toInt(); //-1024

  if (!m_cbctregistration->PreprocessCT(iAirThresholdShort, strRSName,
                                        fill_bubble, iBubbleFillingVal,
                                        iAirFillValShort)) {
    std::cout
        << "Error in PreprocessCT!!!scatter correction would not work out."
        << std::endl;
    m_cbctregistration->m_pParent->m_bMacroContinue = false;
  }

  show();

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);
  // if not found, just skip
  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID);

  std::cout << "FINISHED!: Pre-processing of CT image" << std::endl;

  ////Load DICOM plan
  if (m_cbctregistration->m_pParent->m_strPathPlan.isEmpty()) {
    std::cout << "No DCM plan file was found. Skipping dcm plan." << std::endl;
    return;
  }
  // QString dcmplanPath = m_pParent->m_strPathPlan;
  m_cbctregistration->LoadRTPlan(
      m_cbctregistration->m_pParent->m_strPathPlan); // fill RT_studyplan
}

void DlgRegistration::SLT_DoRegistrationDeform() {
  if ((m_spFixed == nullptr) || (m_spMoving == nullptr)) {
    return;
  }

  bool bPrepareMaskOnly = !(ui.checkBoxCropBkgroundCBCT->isChecked());

  std::cout << "0: DoRegistrationDeform: writing temporary files" << std::endl;

  // Both image type: Unsigned Short
  QString filePathFixed =
      m_cbctregistration->m_strPathPlastimatch + "/" + "fixed_deform.mha";
  QString filePathMoving =
      m_cbctregistration->m_strPathPlastimatch + "/" + "moving_deform.mha";
  QString filePathROI = m_cbctregistration->m_strPathPlastimatch + "/" +
                        "fixed_roi_DIR.mha"; // optional

  QString filePathOutput =
      m_cbctregistration->m_strPathPlastimatch + "/" + "output_deform.mha";
  QString filePathXform =
      m_cbctregistration->m_strPathPlastimatch + "/" + "xform_deform.txt";

  QString filePathOutputStage1 = m_cbctregistration->m_strPathPlastimatch +
                                 "/" + "output_deform_stage1.mha";
  QString filePathOutputStage2 = m_cbctregistration->m_strPathPlastimatch +
                                 "/" + "output_deform_stage2.mha";
  QString filePathOutputStage3 = m_cbctregistration->m_strPathPlastimatch +
                                 "/" + "output_deform_stage3.mha";

  using writerType = itk::ImageFileWriter<UShortImageType>;

  writerType::Pointer writer1 = writerType::New();
  writer1->SetFileName(filePathFixed.toLocal8Bit().constData());
  writer1->SetUseCompression(true);
  writer1->SetInput(m_spFixed);

  writerType::Pointer writer2 = writerType::New();
  writer2->SetFileName(filePathMoving.toLocal8Bit().constData());
  writer2->SetUseCompression(true);
  writer2->SetInput(m_spMoving);
  writer1->Update();
  writer2->Update();

  // Create a mask image based on the fixed sp image

  if (ui.checkBoxUseROIForDIR->isChecked()) {
    std::cout << "Creating a ROI mask for DIR.. " << std::endl;
    QString strFOVGeom = ui.lineEditFOVPos->text();

    QStringList strListFOV = strFOVGeom.split(",");
    if (strListFOV.count() == 3) {
      const auto FOV_DcmPosX = strListFOV.at(0).toFloat(); // mm
      const auto FOV_DcmPosY = strListFOV.at(1).toFloat();
      const auto FOV_Radius = strListFOV.at(2).toFloat();

      // Create Image using FixedImage sp

      // Image Pointer here
      UShortImageType::Pointer spRoiMask;
      AllocateByRef(m_spFixed, spRoiMask);
      m_cbctregistration->m_pParent->GenerateCylinderMask(
          spRoiMask, FOV_DcmPosX, FOV_DcmPosY, FOV_Radius);

      writerType::Pointer writer3 = writerType::New();
      writer3->SetFileName(filePathROI.toLocal8Bit().constData());
      writer3->SetUseCompression(true);
      writer3->SetInput(spRoiMask);
      writer3->Update();
    }
  }

  QString filePathFixed_proc = filePathFixed;

  std::cout << "1: DoRegistrationDeform: CBCT pre-processing before deformable "
               "registration"
            << std::endl;
  // std::cout << "Air region and bubble will be removed" << std::endl;

  QFileInfo info1(m_cbctregistration->m_strPathCTSkin_manRegi);
  QFileInfo info2(m_cbctregistration->m_strPathXFAutoRigid);

  if (!info1.exists() || !info2.exists()) {
    std::cout << "Fatal error! no CT skin is found or no XF auto file found. "
                 "Preprocessing will not be done. Proceeding."
              << std::endl;
  } else {
    filePathFixed_proc = m_cbctregistration->m_strPathPlastimatch + "/" +
                         "fixed_deform_proc.mha";
    // skin removal and bubble filling : output file = filePathFixed_proc
    bool bBubbleRemoval = ui.checkBoxFillBubbleCBCT->isChecked();
    double skinExp = ui.lineEditCBCTSkinCropBfDIR->text().toDouble();

    int iBubThresholdUshort = ui.lineEditBubDetectCBCT->text().toInt();
    int iBubFillUshort = ui.lineEditBubFillCBCT->text().toInt(); // 700

    m_cbctregistration->ProcessCBCT_beforeDeformRegi(
        filePathFixed, m_cbctregistration->m_strPathCTSkin_manRegi,
        filePathFixed_proc, m_cbctregistration->m_strPathXFAutoRigid,
        bBubbleRemoval, bPrepareMaskOnly, skinExp, iBubThresholdUshort,
        iBubFillUshort); // bubble filling yes

    if (bPrepareMaskOnly) {
      filePathFixed_proc = filePathFixed;
    }
  }

  std::cout << "2: DoRegistrationDeform: Creating a plastimatch command file"
            << std::endl;

  QString fnCmdRegisterRigid = "cmd_register_deform.txt";
  // QString fnCmdRegisterDeform = "cmd_register_deform.txt";
  QString pathCmdRegister =
      m_cbctregistration->m_strPathPlastimatch + "/" + fnCmdRegisterRigid;

  QString strDeformableStage1 =
      ui.lineEditArgument1->text(); // original param: 7, add output path
  QString strDeformableStage2 = ui.lineEditArgument2->text();
  QString strDeformableStage3 = ui.lineEditArgument3->text();

  strDeformableStage1.append(", ").append(filePathOutputStage1);
  strDeformableStage2.append(", ").append(filePathOutputStage2);
  strDeformableStage3.append(", ").append(filePathOutputStage3);

  bool mse = ui.radioButton_mse->isChecked();
  bool cuda = m_pParent->ui.radioButton_UseCUDA->isChecked();
  QString GradOptionStr = ui.lineEditGradOption->text();
  // For Cropped patients, FOV mask is applied.
  if (ui.checkBoxUseROIForDIR->isChecked()) {
    m_cbctregistration->GenPlastiRegisterCommandFile(
        pathCmdRegister, filePathFixed_proc, filePathMoving, filePathOutput,
        filePathXform, PLAST_BSPLINE, strDeformableStage1, strDeformableStage2,
        strDeformableStage3, filePathROI, mse, cuda, GradOptionStr);
  } else {
    m_cbctregistration->GenPlastiRegisterCommandFile(
        pathCmdRegister, filePathFixed_proc, filePathMoving, filePathOutput,
        filePathXform, PLAST_BSPLINE, strDeformableStage1, strDeformableStage2,
        strDeformableStage3, QString(""), mse, cuda, GradOptionStr);
  }
  /*void DlgRegistration::GenPlastiRegisterCommandFile(QString
  strPathCommandFile, QString strPathFixedImg, QString strPathMovingImg, QString
  strPathOutImg, QString strPathXformOut, enRegisterOption regiOption, QString
  strStageOption1, , QString strStageOption2, QString strStageOption3)*/

  std::string str_command_filepath = pathCmdRegister.toLocal8Bit().constData();

  m_cbctregistration->CallingPLMCommand(str_command_filepath);

  // Registration_parms *regp = new Registration_parms();
  // if (regp->parse_command_file (command_filepath) < 0) {
  // std::cout << "wrong use of command file" << std::endl;
  // }
  // do_registration (regp);
  // delete regp;

  std::cout << "5: DoRegistrationDeform: Registration is done" << std::endl;
  std::cout << "6: DoRegistrationDeform: Reading output image" << std::endl;

  using readerType = itk::ImageFileReader<UShortImageType>;
  QFileInfo tmpFileInfo;

  readerType::Pointer readerDefSt1 = readerType::New();

  tmpFileInfo = QFileInfo(filePathOutputStage1);
  if (tmpFileInfo.exists()) {
    readerDefSt1->SetFileName(filePathOutputStage1.toLocal8Bit().constData());
    readerDefSt1->Update();
    m_cbctregistration->m_pParent->m_spDeformedCT1 = readerDefSt1->GetOutput();
  }

  readerType::Pointer readerDefSt2 = readerType::New();
  tmpFileInfo = QFileInfo(filePathOutputStage2);
  if (tmpFileInfo.exists()) {
    readerDefSt2->SetFileName(filePathOutputStage2.toLocal8Bit().constData());
    readerDefSt2->Update();
    m_cbctregistration->m_pParent->m_spDeformedCT2 = readerDefSt2->GetOutput();
  }

  readerType::Pointer readerDefSt3 = readerType::New();
  tmpFileInfo = QFileInfo(filePathOutputStage3);
  if (tmpFileInfo.exists()) {
    readerDefSt3->SetFileName(filePathOutputStage3.toLocal8Bit().constData());
    readerDefSt3->Update();
    m_cbctregistration->m_pParent->m_spDeformedCT3 = readerDefSt3->GetOutput();
  }

  readerType::Pointer readerFinCBCT = readerType::New();
  tmpFileInfo = QFileInfo(filePathFixed_proc); // cropped image Or not
  if (tmpFileInfo.exists()) {
    readerFinCBCT->SetFileName(filePathFixed_proc.toLocal8Bit().constData());
    readerFinCBCT->Update();
    m_cbctregistration->m_pParent->m_spRawReconImg = readerFinCBCT->GetOutput();
    // std::cout << "fixed Image Path = " <<
    // filePathFixed_proc.toLocal8Bit().constData() << std::endl;
  } else {
    std::cout << "No filePathFixed_proc is available. Exit the function"
              << std::endl;
    return;
  }

  // Puncturing should be done for final Deform image

  QString strPathDeformCTFinal = filePathOutput;

  QFileInfo tmpBubFileInfo(m_cbctregistration->m_strPathMskCBCTBubble);

  if (ui.checkBoxFillBubbleCBCT->isChecked() && tmpBubFileInfo.exists()) {
    std::cout << "6B: final puncturing according to the CBCT bubble"
              << std::endl;

    /*Mask_parms parms_fill;
    strPathDeformCTFinal = m_strPathPlastimatch + "/deformCTpuncFin.mha";

    parms_fill.mask_operation = MASK_OPERATION_FILL;
    parms_fill.input_fn = filePathOutput.toLocal8Bit().constData();
    parms_fill.mask_fn = m_strPathMskCBCTBubble.toLocal8Bit().constData();
    parms_fill.output_fn = strPathDeformCTFinal.toLocal8Bit().constData();*/

    strPathDeformCTFinal =
        m_cbctregistration->m_strPathPlastimatch + "/deformCTpuncFin.mha";

    Mask_operation enMaskOp = MASK_OPERATION_FILL;
    QString input_fn = filePathOutput;
    QString mask_fn = m_cbctregistration->m_strPathMskCBCTBubble;
    QString output_fn = strPathDeformCTFinal;

    // int iBubblePunctureVal = ui.lineEditBkFillCT->text().toInt(); //0 = soft
    // tissue
    int iBubblePunctureVal =
        0; // 0 = air. deformed CT is now already a USHORT image
    int mask_value = iBubblePunctureVal;
    m_cbctregistration->plm_mask_main(enMaskOp, input_fn, mask_fn, output_fn,
                                      static_cast<float>(mask_value));
  }

  readerType::Pointer readerDeformFinal = readerType::New();
  tmpFileInfo = QFileInfo(strPathDeformCTFinal);
  if (tmpFileInfo.exists()) {
    readerDeformFinal->SetFileName(
        strPathDeformCTFinal.toLocal8Bit().constData());
    readerDeformFinal->Update();
    m_cbctregistration->m_pParent->m_spDeformedCT_Final =
        readerDeformFinal->GetOutput();
  } else {
    std::cout << "No final output is available. Exit the function" << std::endl;
    return;
  }

  std::cout << "7: DoRegistrationDeform: Reading is completed" << std::endl;

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_DEFORM_FINAL);

  std::cout << "FINISHED!: Deformable image registration. Proceed to scatter "
               "correction"
            << std::endl;
}

void DlgRegistration::SLT_PassFixedImgForAnalysis() {
  QString cur_fixed = ui.comboBoxImgFixed->currentText();
  if (m_spFixed != nullptr) {
    m_pParent->UpdateReconImage(m_spFixed, cur_fixed);
  }
}

void DlgRegistration::SLT_PassMovingImgForAnalysis() {
  QString cur_moving = ui.comboBoxImgMoving->currentText();
  if (m_spMoving != nullptr) {
    m_pParent->UpdateReconImage(m_spMoving, cur_moving);
  }
}

void DlgRegistration::SLT_DoLowerMaskIntensity() {
  if (!ui.checkBoxRemoveMaskAfterCor->isChecked()) {
    std::cout << "Error. this function is not enabled" << std::endl;
    return;
  }

  int iDiffThreshold = ui.lineEditRawCorThre->text().toInt();

  int iNoTouchThreshold = ui.lineEditiNoTouchThreshold->text().toInt();

  /*if (!m_pParent->m_spScatCorrReconImg || !m_pParent->m_spRawReconImg)
  {

      std::cout << "You need both raw and corr CBCT images" << std::endl;
      return;
  }*/
  double fInnerMargin = ui.lineEditThermoInner->text().toDouble();
  double fOuterMargin = ui.lineEditThermoOuter->text().toDouble();

  m_cbctregistration->ThermoMaskRemovingCBCT(m_spFixed, m_spMoving,
                                             iDiffThreshold, iNoTouchThreshold,
                                             fInnerMargin, fOuterMargin);

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);
  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_COR_CBCT);
}

void DlgRegistration::SLT_ExchangeRawRef() {}

void DlgRegistration::SLT_ManualMoveByDCMPlan() {
  if ((m_spFixed == nullptr) || (m_spMoving == nullptr)) {
    return;
  }

  auto final_iso_pos = m_cbctregistration->ManualMoveByDCM();

  if (final_iso_pos == nullptr) {
    std::cout << "Error!  No isocenter position was found. " << std::endl;
    return;
  }

  // ImageManualMoveOneShot(-iso_pos[0], -iso_pos[1], -iso_pos[2]);
  ImageManualMoveOneShot(final_iso_pos[0], final_iso_pos[1], final_iso_pos[2]);

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID);
}

void DlgRegistration::SLT_ManualMoveByDCMPlanOpen() {
  QString filePath = QFileDialog::getOpenFileName(
      this, "Open DCMRT Plan file",
      m_cbctregistration->m_pParent->m_strPathDirDefault, "DCMRT Plan (*.dcm)",
      nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  VEC3D planIso = m_cbctregistration->GetIsocenterDCM_FromRTPlan(filePath);

  if (planIso.x == 0.0 && planIso.y == 0.0 && planIso.z == 0.0) {
    std::cout
        << "Warning!!!! Plan iso is 0 0 0. Most likely not processed properly"
        << std::endl;
  } else {
    std::cout << "isocenter was found: " << planIso.x << ", " << planIso.y
              << ", " << planIso.z << std::endl;
  }

  if ((m_spFixed == nullptr) || (m_spMoving == nullptr)) {
    return;
  }

  if ((m_cbctregistration->m_pParent->m_spRefCTImg == nullptr) ||
      (m_cbctregistration->m_pParent->m_spManualRigidCT == nullptr)) {
    return;
  }

  // ImageManualMoveOneShot(-iso_pos[0], -iso_pos[1], -iso_pos[2]);
  ImageManualMoveOneShot(static_cast<float>(planIso.x),
                         static_cast<float>(planIso.y),
                         static_cast<float>(planIso.z));

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID);
}

void DlgRegistration::SLT_Override() {
  bool isFixed = false;
  if (ui.comboBox_imToOverride->currentText().compare(QString("Moving")) == 0) {
    isFixed = true;
  }
  if ((isFixed && (m_spFixed == nullptr)) ||
      (!isFixed && (m_spMoving == nullptr))) {
    std::cout << "The image you try to override is not loaded!" << std::endl;
    return;
  }

  int sliderPosIdxZ, sliderPosIdxY, sliderPosIdxX;

  switch (m_enViewArrange) {
  case AXIAL_FRONTAL_SAGITTAL:
    sliderPosIdxZ =
        ui.sliderPosDisp1
            ->value(); // Z corresponds to axial, Y to frontal, X to sagittal
    sliderPosIdxY = ui.sliderPosDisp2->value();
    sliderPosIdxX = ui.sliderPosDisp3->value();
    break;
  case FRONTAL_SAGITTAL_AXIAL:
    sliderPosIdxY = ui.sliderPosDisp1->value();
    sliderPosIdxX = ui.sliderPosDisp2->value();
    sliderPosIdxZ = ui.sliderPosDisp3->value();

    break;
  case SAGITTAL_AXIAL_FRONTAL:
    sliderPosIdxX = ui.sliderPosDisp1->value();
    sliderPosIdxZ = ui.sliderPosDisp2->value();
    sliderPosIdxY = ui.sliderPosDisp3->value();
    break;
  default:
    sliderPosIdxZ =
        ui.sliderPosDisp1
            ->value(); // Z corresponds to axial, Y to frontal, X to sagittal
    sliderPosIdxY = ui.sliderPosDisp2->value();
    sliderPosIdxX = ui.sliderPosDisp3->value();
    break;
  }

  // UShortImageType::SizeType img_size =
  //  m_spFixed->GetRequestedRegion().GetSize(); // 1016x1016 x z
  UShortImageType::PointType imgOrigin = m_spFixed->GetOrigin();
  UShortImageType::SpacingType imgSpacing = m_spFixed->GetSpacing();
  if (!isFixed) {
    // img_size = m_spMoving->GetRequestedRegion().GetSize(); // 1016x1016 x z
    imgOrigin = m_spMoving->GetOrigin();
    imgSpacing = m_spMoving->GetSpacing();
  }

  UShortImageType::PointType curPhysPos;
  curPhysPos[0] = imgOrigin[0] + sliderPosIdxX * imgSpacing[0]; // Z
  curPhysPos[1] = imgOrigin[1] + sliderPosIdxY * imgSpacing[1]; // Y
  curPhysPos[2] =
      imgOrigin[2] + sliderPosIdxZ * imgSpacing[2]; // Z in default setting

  UShortImageType::IndexType centerIdx{};
  if (!isFixed) {
    if (!m_spMoving->TransformPhysicalPointToIndex(curPhysPos, centerIdx)) {
      std::cerr << "Point not in fixed image!" << std::endl;
      return;
    };
  } else {
    if (!m_spFixed->TransformPhysicalPointToIndex(curPhysPos, centerIdx)) {
      std::cerr << "Point not in moving image!" << std::endl;
      return;
    };
  }

  const int radius = ui.spinBox_overrideRadius->value();
  const auto value =
      static_cast<unsigned short>(ui.spinBox_overrideValue->value() + 1024);
  size_t i = 0;
  UShortImageType::IndexType curIdx{};
  for (int curRadiusX = -radius; curRadiusX <= radius; curRadiusX++) {
    curIdx[0] = centerIdx[0] + curRadiusX;
    for (int curRadiusY = -radius; curRadiusY <= radius; curRadiusY++) {
      curIdx[1] = centerIdx[1] + curRadiusY;
      for (int curRadiusZ = -radius; curRadiusZ <= radius; curRadiusZ++) {
        curIdx[2] = centerIdx[2] + curRadiusZ;
        if (radius >= sqrt(pow(curRadiusX, 2) + pow(curRadiusY, 2) +
                           pow(curRadiusZ, 2))) {
          if (!isFixed) {
            m_spMoving->SetPixel(curIdx, value);
            i++;
          } else {
            m_spFixed->SetPixel(curIdx, value);
            i++;
          }
        }
      }
    }
  }
  std::cout << i << " pixels were overridden." << std::endl;
}

void DlgRegistration::SLT_gPMCrecalc() {
  if (m_spFixed == nullptr) {
    return;
  }

  QString plan_filepath;
  // Load dcm rtplan.
  if (ui.spinBox_NdcmPlans->value() == 1) {
    plan_filepath = QFileDialog::getOpenFileName(
        this, "Open DCMRT Plan file",
        m_cbctregistration->m_pParent->m_strPathDirDefault,
        "DCMRT Plan (*.dcm)", nullptr, nullptr);
  } else {
    plan_filepath = QFileDialog::getOpenFileName(
        this, "Open DCMRT Plan file",
        m_cbctregistration->m_pParent->m_strPathDirDefault,
        "DCMRT Plan (*.dcm)", nullptr, nullptr);
    for (auto i = 1; i < ui.spinBox_NdcmPlans->value(); i++) {
      plan_filepath =
          QString("%1,%2")
              .arg(plan_filepath)
              .arg(QFileDialog::getOpenFileName(
                  this, "Open DCMRT Plan file",
                  m_cbctregistration->m_pParent->m_strPathDirDefault,
                  "DCMRT Plan (*.dcm)", nullptr, nullptr));
    }
  }

  /* Below MUST be done externally, through command line - due to
   * incompatibilities with compilers and library versions.. Pytrip-style but
   * with somewhat better control.. (hacks to private members might be possible
   * as well:
   * http://stackoverflow.com/questions/424104/can-i-access-private-members-from-outside-the-class-without-using-friends)
   */

  // Create gPMC command line

#ifdef USE_CUDA
  enDevice gPMC_device = GPU_DEV;
#elif defined(USE_OPENCL)
  enDevice gPMC_device =
      CPU_DEV; // because I've only tested with intel graphics
#else
  enDevice gPMC_device = CPU_DEV;
#endif

  if (m_pParent->ui.radioButton_UseCPU->isChecked()) {
    gPMC_device = CPU_DEV;
  }
  const auto n_sims = ui.spinBox_Nsims->value();
  const auto n_plans = ui.spinBox_NdcmPlans->value();

  auto success = m_cbctregistration->CallingGPMCcommand(
      gPMC_device, n_sims, n_plans, plan_filepath, m_spFixed, m_spMoving,
      m_spFixedDose, m_spMovingDose);
  if (!success) {
    std::cerr << "Dose calc failed, due to the RNG in goPMC you may want to "
                 "just try again"
              << std::endl;
  }

  SLT_DrawImageWhenSliceChange();
}

void DlgRegistration::SLT_WEPLcalc() {
  // Get VOI
  auto voi_name = ui.comboBox_VOI->currentText().toStdString();

  auto gantry_angle = ui.spinBox_GantryAngle->value();
  auto couch_angle = ui.spinBox_CouchAngle->value();
  m_cbctregistration->CalculateWEPLtoVOI(voi_name, gantry_angle, couch_angle,
                                         m_spMoving);
  // Draw WEPL
  SLT_DrawImageWhenSliceChange();
}

void DlgRegistration::SLT_DoRegistrationGradient() {
  // 1) Save current image files
  if ((m_spFixed == nullptr) || (m_spMoving == nullptr)) {
    return;
  }

  if ((m_cbctregistration->m_pParent->m_spRefCTImg == nullptr) ||
      (m_cbctregistration->m_pParent->m_spManualRigidCT == nullptr)) {
    return;
  }

  std::cout << "1: writing temporary files" << std::endl;
  ui.progressBar->setValue(5);
  // Both image type: Unsigned Short
  QString filePathFixed =
      m_cbctregistration->m_strPathPlastimatch + "/" + "fixed_gradient.mha";
  QString filePathMoving =
      m_cbctregistration->m_strPathPlastimatch + "/" + "moving_gradient.mha";
  QString filePathOutput =
      m_cbctregistration->m_strPathPlastimatch + "/" + "output_gradient.mha";
  QString filePathXform =
      m_cbctregistration->m_strPathPlastimatch + "/" + "xform_gradient.txt";

  using writerType = itk::ImageFileWriter<UShortImageType>;

  writerType::Pointer writer1 = writerType::New();
  writer1->SetFileName(filePathFixed.toLocal8Bit().constData());
  writer1->SetUseCompression(true);
  writer1->SetInput(m_spFixed);

  writerType::Pointer writer2 = writerType::New();
  writer2->SetFileName(filePathMoving.toLocal8Bit().constData());
  writer2->SetUseCompression(true);
  writer2->SetInput(m_spMoving);

  writer1->Update();
  writer2->Update();

  // Preprocessing
  // 2) move CT-based skin mask on CBCT based on manual shift
  //  if (m_strPathCTSkin)

  std::cout << "2: Creating a plastimatch command file" << std::endl;
  ui.progressBar->setValue(10);

  QString fnCmdRegisterGradient = "cmd_register_gradient.txt";
  // QString fnCmdRegisterDeform = "cmd_register_deform.txt";
  QString pathCmdRegister =
      m_cbctregistration->m_strPathPlastimatch + "/" + fnCmdRegisterGradient;

  /*GenPlastiRegisterCommandFile(pathCmdRegister, filePathFixed, filePathMoving,
      filePathOutput, filePathXform, PLAST_GRADIENT, "", "", "");    */

  std::cout << "For Gradient searching only, CBCT image is a moving image, CT "
               "image is fixed image"
            << std::endl;

  bool mse = ui.radioButton_mse->isChecked();
  bool cuda = m_pParent->ui.radioButton_UseCUDA->isChecked();
  auto GradOptionStr = ui.lineEditGradOption->text();
  auto dummyStr = QString("");
  m_cbctregistration->GenPlastiRegisterCommandFile(
      pathCmdRegister, filePathMoving, filePathFixed, filePathOutput,
      filePathXform, PLAST_GRADIENT, dummyStr, dummyStr, dummyStr, dummyStr,
      mse, cuda, GradOptionStr);

  /*void DlgRegistration::GenPlastiRegisterCommandFile(QString
  strPathCommandFile, QString strPathFixedImg, QString strPathMovingImg, QString
  strPathOutImg, QString strPathXformOut, enRegisterOption regiOption, QString
  strStageOption1, , QString strStageOption2, QString strStageOption3)*/

  // const char *command_filepath = pathCmdRegister.toLocal8Bit().constData();
  std::string str_command_filepath = pathCmdRegister.toLocal8Bit().constData();

  ui.progressBar->setValue(15);
  auto trn = m_cbctregistration->CallingPLMCommandXForm(str_command_filepath);
  ui.progressBar->setValue(99); // good ol' 99%

  ImageManualMoveOneShot(static_cast<float>(-trn[0]),
                         static_cast<float>(-trn[1]),
                         static_cast<float>(-trn[2]));
  /*
std::cout << "5: Load shift values" << std::endl;

VEC3D shiftVal = GetShiftValueFromGradientXForm(filePathXform, true); //true:
inverse trans should be applied if CBCT was moving image //in mm

ImageManualMoveOneShot(shiftVal.x, shiftVal.y, shiftVal.z);
  */
  std::cout << "6: Reading is completed" << std::endl;
  ui.progressBar->setValue(100);

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID);

  ui.progressBar->setValue(0);
}

void DlgRegistration::SLT_ConfirmManualRegistration() {
  if ((m_spFixed == nullptr) || (m_spMoving == nullptr)) {
    return;
  }

  if ((m_cbctregistration->m_pParent->m_spRefCTImg == nullptr) ||
      (m_cbctregistration->m_pParent->m_spManualRigidCT == nullptr)) {
    return;
  }

  if (ui.checkBoxKeyMoving->isChecked()) {
    SLT_KeyMoving(false); // uncheck macro
  }

  // Apply post processing for raw CBCT image and generate
  std::cout << "Preprocessing for CBCT" << std::endl;

  bool bPrepareMaskOnly = !(ui.checkBoxCropBkgroundCBCT->isChecked());

  UShortImageType::PointType originBefore =
      m_cbctregistration->m_pParent->m_spRefCTImg->GetOrigin();
  UShortImageType::PointType originAfter =
      m_cbctregistration->m_pParent->m_spManualRigidCT->GetOrigin();

  double fShift[3];
  fShift[0] = (originBefore[0] - originAfter[0]); // DICOM
  fShift[1] = (originBefore[1] - originAfter[1]);
  fShift[2] = (originBefore[2] - originAfter[2]);

  std::cout << "1: writing temporary files" << std::endl;

  // Both image type: Unsigned Short
  QString filePathFixed = m_cbctregistration->m_strPathPlastimatch + "/" +
                          "fixed_rigid.mha"; // CBCT image //redundant
  QString filePathFixed_proc =
      m_cbctregistration->m_strPathPlastimatch + "/" +
      "fixed_rigid_proc.mha"; // After autoRigidbody Regi

  // writing
  using writerType = itk::ImageFileWriter<UShortImageType>;
  writerType::Pointer writer = writerType::New();
  writer->SetFileName(filePathFixed.toLocal8Bit().constData());
  writer->SetUseCompression(true);
  writer->SetInput(m_spFixed);
  writer->Update();

  std::cout << "1.A: Writing temporary files is done" << std::endl;

  QFileInfo finfoSkinFile1 = QFileInfo(m_cbctregistration->m_strPathCTSkin);
  QString strPathAlternateSkin =
      m_cbctregistration->m_strPathPlastimatch + "/" + "msk_skin_CT.mha";
  QFileInfo finfoSkinFile2 = QFileInfo(strPathAlternateSkin);

  QString strPathOriginalCTSkinMask;

  double skinExp = ui.lineEditCBCTSkinCropBfRegid->text().toDouble();
  int bkGroundValUshort = ui.lineEditBkFillCBCT->text().toInt(); // 0

  if (finfoSkinFile1.exists()) {
    strPathOriginalCTSkinMask = m_cbctregistration->m_strPathCTSkin;
    m_cbctregistration->ProcessCBCT_beforeAutoRigidRegi(
        filePathFixed, strPathOriginalCTSkinMask, filePathFixed_proc,
        &fShift[0], bPrepareMaskOnly, skinExp, bkGroundValUshort);

    if (bPrepareMaskOnly) { // currently, filePathFixed_proc == "";
      filePathFixed_proc = filePathFixed;
    }

  } else if (finfoSkinFile2.exists()) {
    std::cout << "alternative skin file will be used" << std::endl;
    strPathOriginalCTSkinMask = strPathAlternateSkin;
    m_cbctregistration->ProcessCBCT_beforeAutoRigidRegi(
        filePathFixed, strPathOriginalCTSkinMask, filePathFixed_proc,
        &fShift[0], bPrepareMaskOnly, skinExp, bkGroundValUshort);

    if (bPrepareMaskOnly) { // currently, filePathFixed_proc == "";
      filePathFixed_proc = filePathFixed;
    }
  }

  QFileInfo fInfo = QFileInfo(filePathFixed_proc);

  if (fInfo.exists() &&
      !bPrepareMaskOnly) // if fixed_rigid_proc.mha is generated successfully.
  {

    std::cout << "Trying to read file: filePathFixed_proc" << std::endl;
    // Update RawReconImg
    using readerType = itk::ImageFileReader<UShortImageType>;
    readerType::Pointer reader = readerType::New();
    reader->SetFileName(filePathFixed_proc.toLocal8Bit().constData());
    reader->Update();
    m_cbctregistration->m_pParent->m_spRawReconImg = reader->GetOutput();

    double tmpSkinMargin = ui.lineEditCBCTSkinCropBfRegid->text().toDouble();
    QString update_message =
        QString("Skin removed CBCT with margin %1 mm").arg(tmpSkinMargin);
    m_pParent->UpdateReconImage(m_cbctregistration->m_pParent->m_spRawReconImg,
                                update_message);

    std::cout << "Reading is completed" << std::endl;

    UpdateListOfComboBox(0); // combo selection signalis called
    UpdateListOfComboBox(1);

    SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
    SelectComboExternal(1, REGISTER_MANUAL_RIGID);
  }

  // Export final xform file
  QString filePathXform =
      m_cbctregistration->m_strPathPlastimatch + "/" + "xform_manual.txt";

  std::ofstream fout;
  fout.open(filePathXform.toLocal8Bit().constData());
  fout << "#Custom Transform" << std::endl;
  fout << "#Transform: Translation_only" << std::endl;
  fout << "Parameters: " << fShift[0] << " " << fShift[1] << " " << fShift[2]
       << std::endl;
  // fout << "FixedParameters: 0 0 0" << std::endl;

  fout.close();

  std::cout << "Writing manual registration transform info is done."
            << std::endl;
}

void DlgRegistration::SLT_IntensityNormCBCT() {
  float fROI_Radius = ui.lineEditNormRoiRadius->text().toFloat();

  std::cout << "Intensity is being analyzed...Please wait." << std::endl;

  float intensitySDFix = 0.0;
  float intensitySDMov = 0.0;
  float meanIntensityFix = m_cbctregistration->m_pParent->GetMeanIntensity(
      m_spFixed, fROI_Radius, &intensitySDFix);
  float meanIntensityMov = m_cbctregistration->m_pParent->GetMeanIntensity(
      m_spMoving, fROI_Radius, &intensitySDMov);

  std::cout << "Mean/SD for Fixed = " << meanIntensityFix << "/"
            << intensitySDFix << std::endl;
  std::cout << "Mean/SD for Moving = " << meanIntensityMov << "/"
            << intensitySDMov << std::endl;

  // m_pParent->ExportReconSHORT_HU(m_spMoving, QString("D:/tmpExport.mha"));

  AddConstHU(m_spFixed, static_cast<int>(meanIntensityMov - meanIntensityFix));
  // SLT_PassMovingImgForAnalysis();

  std::cout << "Intensity shifting is done! Added value = "
            << static_cast<int>(meanIntensityMov - meanIntensityFix)
            << std::endl;
  QString update_message =
      QString("Added_%1")
          .arg(static_cast<int>(meanIntensityMov - meanIntensityFix));
  m_pParent->UpdateReconImage(m_spFixed, update_message);
  SelectComboExternal(0, REGISTER_RAW_CBCT);
}
