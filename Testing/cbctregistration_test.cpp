// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include "cbctregistration_test.hpp"

#include <QString>
#include <qcombobox.h>

#include "OpenCL/ImageFilters.h"
#include "PlmWrapper.h"
#include "StructureSet.h"
#include "cbctrecon.h"
#include "cbctrecon_compute.h"
#include "cbctrecon_io.h"
#include "cbctrecon_test.hpp"
#include "cbctregistration.h"

#define FIXME_BACKGROUND_MAX (-1200)

enum enCOLOR {
  RED,
  GREEN,
};

CbctRegistrationTest::CbctRegistrationTest() {
  m_YKImgFixed = nullptr;
  m_YKImgMoving = nullptr;
  m_YKDisp = nullptr;

  m_DoseImgFixed = nullptr;
  m_DoseImgMoving = nullptr;
  m_AGDisp_Overlay = nullptr;
  ui_comboBoxImgFixed = std::make_unique<MyQComboBox>();
  ui_comboBoxImgMoving = std::make_unique<MyQComboBox>();
  ui_comboBox_VOI = std::make_unique<MyQComboBox>();
  ui_comboBox_VOItoCropBy = std::make_unique<MyQComboBox>();
}

CbctRegistrationTest::CbctRegistrationTest(CbctReconTest *parent) {
  m_pParent = parent;
  m_cbctregistration =
      std::make_unique<CbctRegistration>(parent->m_cbctrecon.get());

  m_YKImgFixed = &m_cbctregistration->m_YKImgFixed[0];
  m_YKImgMoving = &m_cbctregistration->m_YKImgMoving[0];
  m_YKDisp = &m_cbctregistration->m_YKDisp[0];

  m_DoseImgFixed = &m_cbctregistration->m_DoseImgFixed[0];
  m_DoseImgMoving = &m_cbctregistration->m_DoseImgMoving[0];
  m_AGDisp_Overlay = &m_cbctregistration->m_AGDisp_Overlay[0];
  ui_comboBoxImgFixed = std::make_unique<MyQComboBox>();
  ui_comboBoxImgMoving = std::make_unique<MyQComboBox>();
  ui_comboBox_VOI = std::make_unique<MyQComboBox>();
  ui_comboBox_VOItoCropBy = std::make_unique<MyQComboBox>();
}

void CbctRegistrationTest::initCbctRegistrationTest(QString &strDCMUID) {
  m_cbctregistration->SetPlmOutputDir(strDCMUID);

  const UShortImageType::Pointer spNull;
  // unlink all of the pointers
  // m_pParent->m_spReconImg->Delete(); //fixed image // ID: RawCBCT
  auto p_parent = m_cbctregistration->m_pParent;
  p_parent->m_spRefCTImg = spNull;
  p_parent->m_spManualRigidCT =
      spNull; // copied from RefCTImg; ID: RefCT --> Moving Img, cloned
  p_parent->m_spAutoRigidCT = spNull; // ID: AutoRigidCT
  p_parent->m_spDeformedCT1 = spNull; // Deformmation will be carried out based
                                      // on Moving IMage of GUI //AutoDeformCT1
  p_parent->m_spDeformedCT2 = spNull; // AutoDeformCT2
  p_parent->m_spDeformedCT3 = spNull; // AutoDeformCT3
  p_parent->m_spDeformedCT_Final = spNull; // AutoDeformCT3

  this->ui_checkBoxKeyMoving = false;

  // show();

  UpdateListOfComboBox(0);
  UpdateListOfComboBox(1);
  // if not found, just skip
  SelectComboExternal(0, REGISTER_RAW_CBCT);     // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID); // WILL BE IGNORED
}

void CbctRegistrationTest::LoadImgFromComboBox(
    const int idx,
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
    m_spFixed = spTmpImg.GetPointer();

  } else if (idx == 1) {
    m_spMoving = spTmpImg.GetPointer();
    // std::cout << "idx: " << idx << "m_spMoving"  << m_spMoving << std::endl;
  }
}

void CbctRegistrationTest::LoadVOIFromComboBox(
    int /*idx*/, QString &strSelectedComboTxt) const {

  auto ct_type = PLAN_CT;
  const auto ct = this->ui_comboBoxImgMoving->currentText().toStdString();
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
}

// search  for the  main data, if there  is, add  the predefined name to the
// combobox
void CbctRegistrationTest::UpdateListOfComboBox(const int idx) const {
  MyQComboBox *crntCombo;

  if (idx == 0) {
    crntCombo = this->ui_comboBoxImgFixed.get();
  } else {
    crntCombo = this->ui_comboBoxImgMoving.get();
  }

  // remove all the list
  crntCombo->clear();
  const auto p_parent = m_cbctregistration->m_pParent;

  if (p_parent->m_spRawReconImg != nullptr) {
    crntCombo->addItem("RAW_CBCT");
  }

  if (p_parent->m_spRefCTImg != nullptr) {
    crntCombo->addItem("REF_CT");
  }

  if (p_parent->m_spManualRigidCT != nullptr) {
    crntCombo->addItem("MANUAL_RIGID_CT");
  }

  if (p_parent->m_spAutoRigidCT != nullptr) {
    crntCombo->addItem("AUTO_RIGID_CT");
  }

  if (p_parent->m_spDeformedCT1 != nullptr) {
    crntCombo->addItem("DEFORMED_CT1");
  }

  if (p_parent->m_spDeformedCT2 != nullptr) {
    crntCombo->addItem("DEFORMED_CT2");
  }

  if (p_parent->m_spDeformedCT3 != nullptr) {
    crntCombo->addItem("DEFORMED_CT3");
  }

  if (p_parent->m_spDeformedCT_Final != nullptr) {
    crntCombo->addItem("DEFORMED_CT_FINAL");
  }

  if (p_parent->m_spScatCorrReconImg != nullptr) {
    crntCombo->addItem("COR_CBCT");
  }
}

void CbctRegistrationTest::SLT_RestoreImageSingle() const {
  const auto mainWndIdx = 0;
  m_YKDisp[mainWndIdx].SetZoom(1.0);
  m_YKDisp[mainWndIdx].SetOffset(0, 0);
}

void CbctRegistrationTest::SLT_RestoreImageAll() const {
  for (auto &i : m_cbctregistration->m_YKDisp) {
    i.SetZoom(1.0);
    i.SetOffset(0, 0);
  }
}
// Auto registration
void CbctRegistrationTest::SLT_DoRegistrationRigid() // plastimatch auto
                                                     // registration
{
  // 1) Save current image files
  if (m_spFixed == nullptr || m_spMoving == nullptr) {
    return;
  }

  // Before the registration, remove background and fill the air-bubble of the
  // cbct  1) calculate manual shift value  PreProcess_CBCT(); // remove skin
  // and fill bubble before registration

  if (m_pParent->m_cbctrecon->m_spRefCTImg == nullptr ||
      m_pParent->m_cbctrecon->m_spManualRigidCT == nullptr) {
    return;
  }

  if (this->ui_checkBoxKeyMoving) {
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

  std::cout << "1: writing temporary files" << std::endl;

  // Both image type: Unsigned Short
  auto filePathFixed =
      m_cbctregistration->m_strPathPlastimatch + "/" + "fixed_rigid.mha";
  auto filePathMoving =
      m_cbctregistration->m_strPathPlastimatch + "/" + "moving_rigid.mha";
  auto filePathOutput =
      m_cbctregistration->m_strPathPlastimatch + "/" + "output_rigid.mha";
  const auto filePathXform =
      m_cbctregistration->m_strPathPlastimatch + "/" + "xform_rigid.txt";
  auto filePathROI = m_cbctregistration->m_strPathPlastimatch + "/" +
                     "fixed_roi_rigid.mha"; // optional

  using writerType = itk::ImageFileWriter<UShortImageType>;

  auto writer1 = writerType::New();
  writer1->SetFileName(filePathFixed.toLocal8Bit().constData());
  writer1->SetUseCompression(true);
  writer1->SetInput(m_spFixed);

  auto writer2 = writerType::New();
  writer2->SetFileName(filePathMoving.toLocal8Bit().constData());
  writer2->SetUseCompression(true);
  writer2->SetInput(m_spMoving);

  writer1->Update();
  writer2->Update();

  if (this->ui_checkBoxUseROIForRigid) {
    std::cout << "Creating a ROI mask for Rigid registration " << std::endl;
    auto strFOVGeom = this->ui_lineEditFOVPos;

    auto strListFOV = strFOVGeom.split(",");
    if (strListFOV.count() == 3) {
      const auto FOV_DcmPosX = strListFOV.at(0).toFloat(); // mm
      const auto FOV_DcmPosY = strListFOV.at(1).toFloat();
      const auto FOV_Radius = strListFOV.at(2).toFloat();

      // Create Image using FixedImage sp

      // Image Pointer here
      UShortImageType::Pointer spRoiMask;
      AllocateByRef<UShortImageType, UShortImageType>(m_spFixed, spRoiMask);
      m_pParent->m_cbctrecon->GenerateCylinderMask(spRoiMask, FOV_DcmPosX,
                                                   FOV_DcmPosY, FOV_Radius);

      auto writer3 = writerType::New();
      writer3->SetFileName(filePathROI.toLocal8Bit().constData());
      writer3->SetUseCompression(true);
      writer3->SetInput(spRoiMask);
      writer3->Update();
    }
    std::cout << "2.A:  ROI-based Rigid body registration will be done"
              << std::endl;
  } else {
    filePathROI = QString("");
  }

  // Preprocessing
  // 2) move CT-based skin mask on CBCT based on manual shift
  //  if (m_strPathCTSkin)

  // m_strPathCTSkin is prepared after PreProcessCT()

  std::cout << "2: Creating a plastimatch command file" << std::endl;

  const auto fnCmdRegisterRigid = QString("cmd_register_rigid.txt");
  // QString fnCmdRegisterDeform = "cmd_register_deform.txt";
  auto pathCmdRegister =
      m_cbctregistration->m_strPathPlastimatch + "/" + fnCmdRegisterRigid;

  const auto strDummy = QString("");

  const auto mse = this->ui_radioButton_mse;
#ifdef USE_CUDA
  const auto cuda = true;
#else
  const auto cuda = false;
#endif
  auto GradOptionStr = this->ui_lineEditGradOption;
  // For Cropped patients, FOV mask is applied.

  m_cbctregistration->GenPlastiRegisterCommandFile(
      pathCmdRegister, filePathFixed, filePathMoving, filePathOutput,
      filePathXform, PLAST_RIGID, strDummy, strDummy, strDummy, filePathROI,
      mse, cuda, GradOptionStr);

  std::string str_command_filepath = pathCmdRegister.toLocal8Bit().constData();

  m_cbctregistration->CallingPLMCommand(str_command_filepath);

  std::cout << "5: Reading output image-CT" << std::endl;

  auto reader = itk::ImageFileReader<UShortImageType>::New();
  reader->SetFileName(filePathOutput.toLocal8Bit().constData());
  reader->Update();
  m_pParent->m_cbctrecon->m_spAutoRigidCT = reader->GetOutput();

  std::cout << "6: Reading is completed" << std::endl;

  const QFile xform_file(filePathXform);

  // If rigid structure already exists copy it to plan-ct structure before
  // applying rigid again.
  if (!m_cbctregistration->m_pParent->m_structures->is_ss_null<RIGID_CT>()) {
    auto rigid_ss = std::make_unique<Rtss_modern>(
        *(m_cbctregistration->m_pParent->m_structures->get_ss(RIGID_CT)));
    m_cbctregistration->m_pParent->m_structures->set_planCT_ss(
        std::move(rigid_ss));
  }

  const auto transform_success =
      m_cbctregistration->m_pParent->m_structures->ApplyTransformTo<PLAN_CT>(
          xform_file);
  if (transform_success) {
    std::cout << "7: Contours registered" << std::endl;
  }

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_AUTO_RIGID);

  m_cbctregistration->m_strPathXFAutoRigid = filePathXform; // for further use
}

void CbctRegistrationTest::ImageManualMove(const int direction,
                                           const double resol) const {
  if (m_spMoving == nullptr) {
    return;
  }

  // USHORT_ImageType::SizeType imgSize =
  // m_spMoving->GetRequestedRegion().GetSize(); //1016x1016 x z
  auto imgOrigin = m_spMoving->GetOrigin();
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

  // Display relative movement
  // Starting point? RefCT image
  // Only Valid when Moving image is the ManualMove
  auto imgOriginRef = m_pParent->m_cbctrecon->m_spRefCTImg->GetOrigin();

  QString strDelta;
  strDelta.sprintf(
      "delta(mm): %3.1f, %3.1f, %3.1f", imgOrigin[0] - imgOriginRef[0],
      imgOrigin[1] - imgOriginRef[1], imgOrigin[2] - imgOriginRef[2]);
}

// change origin of moving image by shift value acqui_ed from gradient
// registration
void CbctRegistrationTest::ImageManualMoveOneShot(const float shiftX,
                                                  const float shiftY,
                                                  const float shiftZ) const
// DICOM coordinate
{
  if (m_spMoving == nullptr) {
    return;
  }

  using USPointType = UShortImageType::PointType;
  auto imgOrigin = m_spMoving->GetOrigin();
  imgOrigin[0] = imgOrigin[0] - static_cast<USPointType::ValueType>(shiftX);
  imgOrigin[1] = imgOrigin[1] - static_cast<USPointType::ValueType>(shiftY);
  imgOrigin[2] = imgOrigin[2] - static_cast<USPointType::ValueType>(shiftZ);

  m_spMoving->SetOrigin(imgOrigin);

  // Display relative movement
  // Starting point? RefCT image
  // Only Valid when Moving image is the ManualMove
  auto imgOriginRef = m_pParent->m_cbctrecon->m_spRefCTImg->GetOrigin();

  QString strDelta;
  strDelta.sprintf(
      "delta(mm): %3.1f, %3.1f, %3.1f", imgOrigin[0] - imgOriginRef[0],
      imgOrigin[1] - imgOriginRef[1], imgOrigin[2] - imgOriginRef[2]);
  std::cerr << strDelta.toStdString() << "\n";
}

void CbctRegistrationTest::SLT_KeyMoving(
    const bool bChecked) // Key Moving check box
{
  if (bChecked) {
    SelectComboExternal(1,
                        REGISTER_MANUAL_RIGID); // should be
  }
}

void CbctRegistrationTest::AddImageToCombo(const int comboIdx,
                                           const enREGI_IMAGES option) const
// comboIdx 0: fixed, 1: moving
{
  switch (option) {
  case REGISTER_RAW_CBCT:
    if (m_pParent->m_cbctrecon->m_spRawReconImg != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("RAW_CBCT");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("RAW_CBCT");
      }
    }
    break;
  case REGISTER_REF_CT:
    if (m_pParent->m_cbctrecon->m_spRefCTImg != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("REF_CT");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("REF_CT");
      }
    }
    break;
  case REGISTER_MANUAL_RIGID:
    if (m_pParent->m_cbctrecon->m_spManualRigidCT != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("MANUAL_RIGID_CT");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("MANUAL_RIGID_CT");
      }
    }
    break;
  case REGISTER_AUTO_RIGID:
    if (m_pParent->m_cbctrecon->m_spAutoRigidCT != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("AUTO_RIGID_CT");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("AUTO_RIGID_CT");
      }
    }
    break;
  case REGISTER_DEFORM1:
    if (m_pParent->m_cbctrecon->m_spDeformedCT1 != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("DEFORMED_CT1");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("DEFORMED_CT1");
      }
    }
    break;
  case REGISTER_DEFORM2:
    if (m_pParent->m_cbctrecon->m_spDeformedCT2 != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("DEFORMED_CT2");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("DEFORMED_CT2");
      }
    }
    break;
  case REGISTER_DEFORM3:
    if (m_pParent->m_cbctrecon->m_spDeformedCT3 != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("DEFORMED_CT3");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("DEFORMED_CT3");
      }
    }
    break;

  case REGISTER_DEFORM_FINAL:
    if (m_pParent->m_cbctrecon->m_spDeformedCT_Final != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("DEFORMED_CT_FINAL");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("DEFORMED_CT_FINAL");
      }
    }
    break;
  case REGISTER_COR_CBCT:
    if (m_pParent->m_cbctrecon->m_spScatCorrReconImg != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("COR_CBCT");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("COR_CBCT");
      }
    }
    break;
  case REGISTER_DEFORM_SKIP_AUTORIGID:
    if (m_pParent->m_cbctrecon->m_spScatCorrReconImg != nullptr) {
      if (comboIdx == 0) {
        this->ui_comboBoxImgFixed->addItem("REGISTER_DEFORM_SKIP_AUTORIGID");
      } else if (comboIdx == 1) {
        this->ui_comboBoxImgMoving->addItem("REGISTER_DEFORM_SKIP_AUTORIGID");
      }
    }
    break;
  }
}

// externally change  combo box value
void CbctRegistrationTest::SelectComboExternal(const int idx,
                                               const enREGI_IMAGES iImage) {
  MyQComboBox *crntCombo;

  if (idx == 0) {
    crntCombo = this->ui_comboBoxImgFixed.get();
  } else if (idx == 1) {
    crntCombo = this->ui_comboBoxImgMoving.get();
  } else {
    std::cerr << "What did you do to get here?" << std::endl;
    return;
  }

  auto findIndx = -1;
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

  const auto crntStr = crntCombo->currentText();
  if (idx == 0) {
    SLT_FixedImageSelected(crntStr);
  } else if (idx == 1) {
    SLT_MovingImageSelected(crntStr);
  }
  // If it is called inside the Dlg, SLT seemed not to conneced
}

ctType get_ctType(const QString &selText) {
  if (selText.compare("REF_CT") == 0) {
    return PLAN_CT;
  }
  if (selText.compare("MANUAL_RIGID_CT") == 0 ||
      selText.compare("AUTO_RIGID_CT") == 0) {
    return RIGID_CT;
  }
  if (selText.compare("DEFORMED_CT1") == 0 ||
      selText.compare("DEFORMED_CT2") == 0 ||
      selText.compare("DEFORMED_CT3") == 0 ||
      selText.compare("DEFORMED_CT_FINAL") == 0) {
    return DEFORM_CT;
  }
  return PLAN_CT;
}

void CbctRegistrationTest::SLT_FixedImageSelected(QString selText) {
  // QString strCrntText = this->ui_comboBoxImgFixed->currentText();
  LoadImgFromComboBox(
      0, selText); // here, m_spMoving and Fixed images are determined
}

void CbctRegistrationTest::SLT_MovingImageSelected(QString selText) {
  // QString strCrntText = this->ui_comboBoxImgMoving->currentText();
  // std::cout << "SLT_MovingImageSelected" << std::endl;
  LoadImgFromComboBox(1, selText);
  const auto cur_ct = get_ctType(selText);
  UpdateVOICombobox(cur_ct);
}

void CbctRegistrationTest::UpdateVOICombobox(const ctType ct_type) const {
  auto struct_set =
      m_cbctregistration->m_pParent->m_structures->get_ss(ct_type);
  if (struct_set == nullptr) {
    return;
  }
  if (struct_set->slist.empty()) {
    std::cerr << "Structures not initialized yet" << std::endl;
    return;
  }
  this->ui_comboBox_VOI->clear();
  for (const auto &voi : struct_set->slist) {
    this->ui_comboBox_VOI->addItem(QString(voi.name.c_str()));
    this->ui_comboBox_VOItoCropBy->addItem(QString(voi.name.c_str()));
  }
}

void CbctRegistrationTest::SLT_RestoreMovingImg() {
  m_cbctregistration->m_pParent->RegisterImgDuplication(REGISTER_REF_CT,
                                                        REGISTER_MANUAL_RIGID);
  SelectComboExternal(1, REGISTER_MANUAL_RIGID);
}

void CbctRegistrationTest::SLT_PreProcessCT() {
  if (!this->ui_checkBoxCropBkgroundCT) {
    std::cout << "Preprocessing is not selected." << std::endl;
    return;
  }

  const auto iAirThresholdShort = this->ui_spinBoxBkDetectCT;

  if (this->ui_comboBox_VOItoCropBy->count() < 1) {
    std::cout
        << "Reference CT DIR should be specified for structure based cropping"
        << std::endl;
    if (m_spMoving == nullptr || m_spFixed == nullptr) {
      return;
    }
    const auto fixed_size = m_spFixed->GetLargestPossibleRegion().GetSize();
    const auto moving_size = m_spMoving->GetLargestPossibleRegion().GetSize();
    if (fixed_size[0] != moving_size[0] || fixed_size[1] != moving_size[1] ||
        fixed_size[2] != moving_size[2]) {
      std::cout
          << "Fixed and moving image is not the same size, consider using "
             "a platimatch registration to solve this."
          << std::endl;
      return;
    }

    const auto reply = true;
    /*QMessageBox::question(this, "No reference structures found!",
                          "Do you wan't to attempt an auto correction "
                          "of air and excessive circumference?",
                          QMessageBox::Yes | QMessageBox::No);*/
    if (reply) {
      std::cout << "Attempting automatic air filling and skin cropping..."
                << std::endl;
      m_cbctregistration->autoPreprocessCT(iAirThresholdShort, m_spFixed,
                                           m_spMoving);
    }
    return;
  }

  const auto strRSName = this->ui_comboBox_VOItoCropBy->currentText();
  const auto fill_bubble = this->ui_checkBoxFillBubbleCT;
  const auto iBubbleFillingVal =
      this->ui_spinBoxBubFillCT;                          // 1000 = soft tissue
  const auto iAirFillValShort = this->ui_spinBoxBkFillCT; // 0 = air

  const auto cur_ct_text = ui_comboBoxImgMoving->currentText();
  const auto cur_ct = get_ctType(cur_ct_text);
  const auto &rt_structs =
      m_cbctregistration->m_pParent->m_structures->get_ss(cur_ct);

  QString image_str;
  if (ui_comboBoxImToCropFill.compare("Moving") == 0) {
    image_str = ui_comboBoxImgMoving->currentText();
  } else {
    image_str = ui_comboBoxImgFixed->currentText();
  }

  auto &image = m_cbctregistration->get_image_from_combotext(image_str);

  if (!m_cbctregistration->PreprocessCT(image, iAirThresholdShort, rt_structs,
                                        strRSName, fill_bubble,
                                        iBubbleFillingVal, iAirFillValShort)) {
    std::cout
        << "Error in PreprocessCT!!!scatter correction would not work out."
        << std::endl;
    m_cbctregistration->m_pParent->m_bMacroContinue = false;
  }

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

void CbctRegistrationTest::SLT_DoRegistrationDeform() {
  if (m_spFixed == nullptr || m_spMoving == nullptr) {
    return;
  }

  const auto bPrepareMaskOnly = !this->ui_checkBoxCropBkgroundCT;

  std::cout << "0: DoRegistrationDeform: writing temporary files" << std::endl;

  // Both image type: Unsigned Short
  auto filePathFixed =
      m_cbctregistration->m_strPathPlastimatch + "/" + "fixed_deform.mha";
  auto filePathMoving =
      m_cbctregistration->m_strPathPlastimatch + "/" + "moving_deform.mha";
  auto filePathROI = m_cbctregistration->m_strPathPlastimatch + "/" +
                     "fixed_roi_DIR.mha"; // optional

  const auto filePathOutput =
      m_cbctregistration->m_strPathPlastimatch + "/" + "output_deform.mha";
  const auto filePathXform =
      m_cbctregistration->m_strPathPlastimatch + "/" + "xform_deform.txt";

  auto filePathOutputStage1 = m_cbctregistration->m_strPathPlastimatch + "/" +
                              "output_deform_stage1.mha";
  auto filePathOutputStage2 = m_cbctregistration->m_strPathPlastimatch + "/" +
                              "output_deform_stage2.mha";
  auto filePathOutputStage3 = m_cbctregistration->m_strPathPlastimatch + "/" +
                              "output_deform_stage3.mha";

  using writerType = itk::ImageFileWriter<UShortImageType>;

  auto writer1 = writerType::New();
  writer1->SetFileName(filePathFixed.toLocal8Bit().constData());
  writer1->SetUseCompression(true);
  writer1->SetInput(m_spFixed);

  auto writer2 = writerType::New();
  writer2->SetFileName(filePathMoving.toLocal8Bit().constData());
  writer2->SetUseCompression(true);
  writer2->SetInput(m_spMoving);
  writer1->Update();
  writer2->Update();

  // Create a mask image based on the fixed sp image

  if (this->ui_checkBoxUseROIForDIR) {
    std::cout << "Creating a ROI mask for DIR.. " << std::endl;
    auto strFOVGeom = this->ui_lineEditFOVPos;

    auto strListFOV = strFOVGeom.split(",");
    if (strListFOV.count() == 3) {
      const auto FOV_DcmPosX = strListFOV.at(0).toFloat(); // mm
      const auto FOV_DcmPosY = strListFOV.at(1).toFloat();
      const auto FOV_Radius = strListFOV.at(2).toFloat();

      // Create Image using FixedImage sp

      // Image Pointer here
      UShortImageType::Pointer spRoiMask;
      AllocateByRef<UShortImageType, UShortImageType>(m_spFixed, spRoiMask);
      m_cbctregistration->m_pParent->GenerateCylinderMask(
          spRoiMask, FOV_DcmPosX, FOV_DcmPosY, FOV_Radius);

      auto writer3 = writerType::New();
      writer3->SetFileName(filePathROI.toLocal8Bit().constData());
      writer3->SetUseCompression(true);
      writer3->SetInput(spRoiMask);
      writer3->Update();
    }
  }

  auto filePathFixed_proc = filePathFixed;

  std::cout << "1: DoRegistrationDeform: CBCT pre-processing before deformable "
               "registration"
            << std::endl;
  // std::cout << "Air region and bubble will be removed" << std::endl;

  QFileInfo info1(m_cbctregistration->m_strPathCTSkin_autoRegi);
  QFileInfo info2(m_cbctregistration->m_strPathXFAutoRigid);

  if (!info1.exists() || !info2.exists()) {
    std::cout << "Fatal error! no CT skin is found or no XF auto file found. "
                 "Preprocessing will not be done. Proceeding."
              << std::endl;
  } else {
    filePathFixed_proc = m_cbctregistration->m_strPathPlastimatch + "/" +
                         "fixed_deform_proc.mha";
    // skin removal and bubble filling : output file = filePathFixed_proc
    const auto bBubbleRemoval = this->ui_checkBoxFillBubbleCBCT;
    const auto skinExp = this->ui_lineEditCBCTSkinCropBfDIR.toDouble();

    const auto iBubThresholdUshort = this->ui_spinBoxBkDetectCT; // 500
    const auto iBubFillUshort = this->ui_spinBoxBubFillCT;       // 1000

    m_cbctregistration->ProcessCBCT_beforeDeformRegi(
        filePathFixed, m_cbctregistration->m_strPathCTSkin_autoRegi,
        filePathFixed_proc, m_cbctregistration->m_strPathXFAutoRigid,
        bBubbleRemoval, bPrepareMaskOnly, skinExp, iBubThresholdUshort,
        iBubFillUshort); // bubble filling yes

    if (bPrepareMaskOnly) {
      filePathFixed_proc = filePathFixed;
    }
  }

  std::cout << "2: DoRegistrationDeform: Creating a plastimatch command file"
            << std::endl;

  const auto fnCmdRegisterRigid = QString("cmd_register_deform.txt");
  // QString fnCmdRegisterDeform = "cmd_register_deform.txt";
  auto pathCmdRegister =
      m_cbctregistration->m_strPathPlastimatch + "/" + fnCmdRegisterRigid;

  auto strDeformableStage1 =
      this->ui_lineEditArgument1; // original param: 7, add output path
  auto strDeformableStage2 = this->ui_lineEditArgument2;
  auto strDeformableStage3 = this->ui_lineEditArgument3;

  strDeformableStage1.append(", ").append(filePathOutputStage1);
  strDeformableStage2.append(", ").append(filePathOutputStage2);
  strDeformableStage3.append(", ").append(filePathOutputStage3);

  const auto mse = this->ui_radioButton_mse;
  const auto cuda = this->ui_radioButton_UseCUDA;
  auto GradOptionStr = this->ui_lineEditGradOption;
  // For Cropped patients, FOV mask is applied.
  if (this->ui_checkBoxUseROIForDIR) {
    m_cbctregistration->GenPlastiRegisterCommandFile(
        pathCmdRegister, filePathFixed_proc, filePathMoving, filePathOutput,
        filePathXform, PLAST_BSPLINE, strDeformableStage1, strDeformableStage2,
        strDeformableStage3, filePathROI, mse, cuda, GradOptionStr);
  } else {
    m_cbctregistration->GenPlastiRegisterCommandFile(
        pathCmdRegister, filePathFixed_proc, filePathMoving, filePathOutput,
        filePathXform, PLAST_BSPLINE, strDeformableStage1, strDeformableStage2,
        strDeformableStage3, QString(), mse, cuda, GradOptionStr);
  }
  /*void CbctRegistrationTest::GenPlastiRegisterCommandFile(QString
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

  auto readerDefSt1 = readerType::New();

  auto tmpFileInfo = QFileInfo(filePathOutputStage1);
  if (tmpFileInfo.exists()) {
    readerDefSt1->SetFileName(filePathOutputStage1.toLocal8Bit().constData());
    readerDefSt1->Update();
    m_cbctregistration->m_pParent->m_spDeformedCT1 = readerDefSt1->GetOutput();
  }

  auto readerDefSt2 = readerType::New();
  tmpFileInfo = QFileInfo(filePathOutputStage2);
  if (tmpFileInfo.exists()) {
    readerDefSt2->SetFileName(filePathOutputStage2.toLocal8Bit().constData());
    readerDefSt2->Update();
    m_cbctregistration->m_pParent->m_spDeformedCT2 = readerDefSt2->GetOutput();
  }

  auto readerDefSt3 = readerType::New();
  tmpFileInfo = QFileInfo(filePathOutputStage3);
  if (tmpFileInfo.exists()) {
    readerDefSt3->SetFileName(filePathOutputStage3.toLocal8Bit().constData());
    readerDefSt3->Update();
    m_cbctregistration->m_pParent->m_spDeformedCT3 = readerDefSt3->GetOutput();
  }

  auto readerFinCBCT = readerType::New();
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

  auto strPathDeformCTFinal = filePathOutput;

  QFileInfo tmpBubFileInfo(m_cbctregistration->m_strPathMskCBCTBubble);

  if (this->ui_checkBoxFillBubbleCBCT && tmpBubFileInfo.exists()) {
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

    const auto enMaskOp = MASK_OPERATION_FILL;
    auto input_fn = filePathOutput;
    auto mask_fn = m_cbctregistration->m_strPathMskCBCTBubble;
    auto output_fn = strPathDeformCTFinal;

    // int iBubblePunctureVal = this->ui_lineEditBkFillCT->text().toInt(); //0 =
    // soft tissue
    const auto iBubblePunctureVal =
        0; // 0 = air. deformed CT is now already a USHORT image
    const auto mask_value = iBubblePunctureVal;
    m_cbctregistration->plm_mask_main(enMaskOp, input_fn, mask_fn, output_fn,
                                      static_cast<float>(mask_value));
  }

  auto readerDeformFinal = readerType::New();
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

  const QFile xform_file(filePathXform);
  m_cbctregistration->m_pParent->m_structures->ApplyTransformTo<RIGID_CT>(
      xform_file);

  std::cout << "8: Contours deformed" << std::endl;

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_DEFORM_FINAL);

  std::cout << "FINISHED!: Deformable image registration. Proceed to scatter "
               "correction"
            << std::endl;
}

void CbctRegistrationTest::SLT_DoLowerMaskIntensity() {
  if (!this->ui_checkBoxRemoveMaskAfterCor) {
    std::cout << "Error. this function is not enabled" << std::endl;
    return;
  }

  const auto iDiffThreshold = this->ui_lineEditRawCorThre.toInt();

  const auto iNoTouchThreshold = this->ui_lineEditiNoTouchThreshold.toInt();

  /*if (!m_pParent->m_spScatCorrReconImg || !m_pParent->m_spRawReconImg)
  {

      std::cout << "You need both raw and corr CBCT images" << std::endl;
      return;
  }*/
  const auto fInnerMargin = this->ui_lineEditThermoInner.toDouble();
  const auto fOuterMargin = this->ui_lineEditThermoOuter.toDouble();

  m_cbctregistration->ThermoMaskRemovingCBCT(m_spFixed, m_spMoving,
                                             iDiffThreshold, iNoTouchThreshold,
                                             fInnerMargin, fOuterMargin);

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);
  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_COR_CBCT);
}

void CbctRegistrationTest::SLT_ManualMoveByDCMPlan() {
  if (m_spFixed == nullptr || m_spMoving == nullptr) {
    return;
  }

  const auto final_iso_pos = m_cbctregistration->ManualMoveByDCM();

  if (final_iso_pos == nullptr) {
    std::cout << "Error!  No isocenter position was found. " << std::endl;
    return;
  }

  ImageManualMoveOneShot(final_iso_pos[0], final_iso_pos[1], final_iso_pos[2]);

  const auto trn_vec = FloatVector{static_cast<float>(final_iso_pos[0]),
                                   static_cast<float>(final_iso_pos[1]),
                                   static_cast<float>(final_iso_pos[2])};
  auto &structs = m_cbctregistration->m_pParent->m_structures;
  structs->ApplyVectorTransform_InPlace<PLAN_CT>(trn_vec);

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID);
}

void CbctRegistrationTest::SLT_ManualMoveByDCMPlanOpen(QString &filePath) {
  const auto p_parent = m_cbctregistration->m_pParent;

  if (filePath.length() < 1) {
    return;
  }

  const auto planIso = m_cbctregistration->GetIsocenterDCM_FromRTPlan(filePath);

  if (fabs(planIso.x) < 0.001 && fabs(planIso.y) < 0.001 &&
      fabs(planIso.z) < 0.001) {
    std::cout
        << "Warning!!!! Plan iso is 0 0 0. Most likely not processed properly"
        << std::endl;
  } else {
    std::cout << "isocenter was found: " << planIso.x << ", " << planIso.y
              << ", " << planIso.z << std::endl;
  }

  if (m_spFixed == nullptr || m_spMoving == nullptr) {
    return;
  }

  if (p_parent->m_spRefCTImg == nullptr ||
      p_parent->m_spManualRigidCT == nullptr) {
    return;
  }

  ImageManualMoveOneShot(static_cast<float>(planIso.x),
                         static_cast<float>(planIso.y),
                         static_cast<float>(planIso.z));

  const auto trn_vec =
      FloatVector{static_cast<float>(planIso.x), static_cast<float>(planIso.y),
                  static_cast<float>(planIso.z)};
  auto &structs = m_cbctregistration->m_pParent->m_structures;
  structs->ApplyVectorTransform_InPlace<PLAN_CT>(trn_vec);

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID);
}

void CbctRegistrationTest::SLT_Override(const int sliderPosIdxX,
                                        const int sliderPosIdxY,
                                        const int sliderPosIdxZ,
                                        const size_t radius,
                                        const int new_value,
                                        const QString &img_to_override) const {
  auto isFixed = false;
  if (img_to_override.compare(QString("Moving")) == 0) {
    isFixed = true;
  }
  if ((isFixed && m_spFixed == nullptr) ||
      (!isFixed && m_spMoving == nullptr)) {
    std::cout << "The image you try to override is not loaded!\n";
    return;
  }

  // UShortImageType::SizeType img_size =
  //  m_spFixed->GetRequestedRegion().GetSize(); // 1016x1016 x z
  auto imgOrigin = m_spFixed->GetOrigin();
  auto imgSpacing = m_spFixed->GetSpacing();
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
    }
  } else {
    if (!m_spFixed->TransformPhysicalPointToIndex(curPhysPos, centerIdx)) {
      std::cerr << "Point not in moving image!" << std::endl;
      return;
    }
  }

  const auto value = static_cast<unsigned short>(new_value + 1024);
  size_t i = 0;
  UShortImageType::IndexType curIdx{};
  const auto s_radius = static_cast<int>(radius);
  for (auto curRadiusX = -s_radius; curRadiusX <= s_radius; curRadiusX++) {
    curIdx[0] = centerIdx[0] + curRadiusX;
    for (auto curRadiusY = -s_radius; curRadiusY <= s_radius; curRadiusY++) {
      curIdx[1] = centerIdx[1] + curRadiusY;
      for (auto curRadiusZ = -1; curRadiusZ <= 1; curRadiusZ++) {
        curIdx[2] = centerIdx[2] + curRadiusZ;
        if (static_cast<double>(radius) >=
            sqrt(pow(curRadiusX, 2) + pow(curRadiusY, 2) +
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

void CbctRegistrationTest::SLT_DoEclRegistration(
    const DoubleVector &translation_vec, const DoubleVector &rotation_vec) {
  auto &ct_img = this->m_spFixed;

  saveImageAsMHA<UShortImageType>(ct_img, "fixed_bkp.mha");
  this->m_spFixed = CbctRegistration::MoveByEclRegistration(
      translation_vec, rotation_vec, ct_img);
}

void CbctRegistrationTest::SLT_ResetEclRegistration() {
  const auto filename = std::string("fixed_bkp.mha");
  this->m_spFixed = loadMHAImageAs<UShortImageType>(filename);
}

void CbctRegistrationTest::SLT_gPMCrecalc(std::vector<QString> &dcm_plans,
                                          const size_t n_sims) {
  if (m_spFixed == nullptr) {
    return;
  }

  QString plan_filepath;
  // Load dcm rtplan.
  if (dcm_plans.size() == 1) {
    plan_filepath = dcm_plans.at(0);
  } else {
    plan_filepath = dcm_plans.at(0);
    for (auto i = 1U; i < dcm_plans.size(); ++i) {
      plan_filepath = QString("%1,%2").arg(plan_filepath, dcm_plans.at(i));
    }
  }

  // Create gPMC command line

#ifdef USE_CUDA
  auto gPMC_device = GPU_DEV;
#elif defined(USE_OPENCL)
  enDevice gPMC_device =
      CPU_DEV; // because I've only tested with intel graphics
#else
  auto gPMC_device = CPU_DEV;
#endif

  if (this->ui_radioButton_UseCPU) {
    gPMC_device = CPU_DEV;
  }
  const auto n_plans = dcm_plans.size();

  const auto success = m_cbctregistration->CallingGPMCcommand(
      gPMC_device, n_sims, n_plans, plan_filepath, m_spFixed, m_spMoving,
      m_spFixedDose, m_spMovingDose);
  if (!success) {
    std::cerr << "Dose calc failed, due to the RNG in goPMC you may want to "
                 "just try again"
              << std::endl;
  }
}

void CbctRegistrationTest::SLT_WEPLcalc(const int gantry_angle,
                                        const int couch_angle) const {
  // Get VOI
  const auto voi_name = this->ui_comboBox_VOI->currentText().toStdString();

  const auto ct_type = get_ctType(ui_comboBoxImgMoving->currentText());
  const auto ss = m_pParent->m_cbctrecon->m_structures->get_ss(ct_type);
  m_cbctregistration->cur_voi = ss->get_roi_by_name(voi_name);

  const auto wepl_voi =
      CalculateWEPLtoVOI(m_cbctregistration->cur_voi.get(), gantry_angle,
                         couch_angle, m_spMoving, m_spFixed);
  m_cbctregistration->WEPL_voi = std::make_unique<Rtss_roi_modern>(*wepl_voi);
  // Draw WEPL
}

void CbctRegistrationTest::SLT_DoRegistrationGradient() {
  // 1) Save current image files
  if (m_spFixed == nullptr || m_spMoving == nullptr) {
    return;
  }

  if (m_cbctregistration->m_pParent->m_spRefCTImg == nullptr ||
      m_cbctregistration->m_pParent->m_spManualRigidCT == nullptr) {
    return;
  }

  std::cout << "1: writing temporary files" << std::endl;
  // Both image type: Unsigned Short
  auto filePathFixed =
      m_cbctregistration->m_strPathPlastimatch + "/" + "fixed_gradient.mha";
  auto filePathMoving =
      m_cbctregistration->m_strPathPlastimatch + "/" + "moving_gradient.mha";
  const auto filePathOutput =
      m_cbctregistration->m_strPathPlastimatch + "/" + "output_gradient.mha";
  const auto filePathXform =
      m_cbctregistration->m_strPathPlastimatch + "/" + "xform_gradient.txt";

  using writerType = itk::ImageFileWriter<UShortImageType>;

  auto writer1 = writerType::New();
  writer1->SetFileName(filePathFixed.toLocal8Bit().constData());
  writer1->SetUseCompression(true);
  writer1->SetInput(m_spFixed);

  auto writer2 = writerType::New();
  writer2->SetFileName(filePathMoving.toLocal8Bit().constData());
  writer2->SetUseCompression(true);
  writer2->SetInput(m_spMoving);

  writer1->Update();
  writer2->Update();

  // Preprocessing
  // 2) move CT-based skin mask on CBCT based on manual shift
  //  if (m_strPathCTSkin)

  std::cout << "2: Creating a plastimatch command file" << std::endl;

  const auto fnCmdRegisterGradient = QString("cmd_register_gradient.txt");
  // QString fnCmdRegisterDeform = "cmd_register_deform.txt";
  auto pathCmdRegister =
      m_cbctregistration->m_strPathPlastimatch + "/" + fnCmdRegisterGradient;

  /*GenPlastiRegisterCommandFile(pathCmdRegister, filePathFixed, filePathMoving,
      filePathOutput, filePathXform, PLAST_GRADIENT, "", "", "");    */

  std::cout << "For Gradient searching only, CBCT image is a moving image, CT "
               "image is fixed image"
            << std::endl;

  const auto mse = this->ui_radioButton_mse;
  const auto cuda = this->ui_radioButton_UseCUDA;
  auto GradOptionStr = this->ui_lineEditGradOption;
  const auto dummyStr = QString("");
  m_cbctregistration->GenPlastiRegisterCommandFile(
      pathCmdRegister, filePathMoving, filePathFixed, filePathOutput,
      filePathXform, PLAST_GRADIENT, dummyStr, dummyStr, dummyStr, dummyStr,
      mse, cuda, GradOptionStr);

  // const char *command_filepath = pathCmdRegister.toLocal8Bit().constData();
  std::string str_command_filepath = pathCmdRegister.toLocal8Bit().constData();

  auto trn = m_cbctregistration->CallingPLMCommandXForm(str_command_filepath);

  const auto trn_vec =
      FloatVector{static_cast<float>(-trn[0]), static_cast<float>(-trn[1]),
                  static_cast<float>(-trn[-2])};

  auto &structs = m_cbctregistration->m_pParent->m_structures;
  structs->ApplyVectorTransform_InPlace<PLAN_CT>(trn_vec);

  ImageManualMoveOneShot(static_cast<float>(-trn[0]),
                         static_cast<float>(-trn[1]),
                         static_cast<float>(-trn[2]));

  std::cout << "6: Reading is completed" << std::endl;

  UpdateListOfComboBox(0); // combo selection signalis called
  UpdateListOfComboBox(1);

  SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  SelectComboExternal(1, REGISTER_MANUAL_RIGID);
}

void CbctRegistrationTest::SLT_ConfirmManualRegistration() {
  if (m_spFixed == nullptr || m_spMoving == nullptr) {
    return;
  }

  auto p_parent = m_cbctregistration->m_pParent;
  if (p_parent->m_spRefCTImg == nullptr ||
      p_parent->m_spManualRigidCT == nullptr) {
    return;
  }

  if (this->ui_checkBoxKeyMoving) {
    SLT_KeyMoving(false); // uncheck macro
  }

  // Apply post processing for raw CBCT image and generate
  std::cout << "Preprocessing for CBCT" << std::endl;

  auto originBefore = p_parent->m_spRefCTImg->GetOrigin();
  auto originAfter = p_parent->m_spManualRigidCT->GetOrigin();

  double fShift[3];
  fShift[0] = originBefore[0] - originAfter[0]; // DICOM
  fShift[1] = originBefore[1] - originAfter[1];
  fShift[2] = originBefore[2] - originAfter[2];

  const auto trn_vec = FloatVector{static_cast<float>(-fShift[0]),
                                   static_cast<float>(-fShift[1]),
                                   static_cast<float>(-fShift[2])};
  m_cbctregistration->m_pParent->m_structures
      ->ApplyVectorTransform_InPlace<PLAN_CT>(trn_vec);

  if (this->ui_checkBoxCropBkgroundCT) {

    // RIGID_CT structs should be created above
    auto structures =
        m_cbctregistration->m_pParent->m_structures->get_ss(RIGID_CT);
    const auto voi_name = this->ui_comboBox_VOItoCropBy->currentText();
    if (structures != nullptr) {
      const auto voi = structures->get_roi_ref_by_name(voi_name.toStdString());
      OpenCL_crop_by_struct_InPlace(p_parent->m_spRawReconImg, voi);
    }

    /*auto update_message = QString("CBCT cropped outside of ") + voi_name;
    m_pParent->UpdateReconImage(p_parent->m_spRawReconImg, update_message);*/

    UpdateListOfComboBox(0); // combo selection signalis called
    UpdateListOfComboBox(1);

    SelectComboExternal(0, REGISTER_RAW_CBCT); // will call fixedImageSelected
    SelectComboExternal(1, REGISTER_MANUAL_RIGID);
  }

  // Export final xform file
  auto filePathXform =
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

  // this->ui.pushButtonConfirmManualRegi->setDisabled(true);
}

void CbctRegistrationTest::SLT_IntensityNormCBCT(const float fROI_Radius) {

  std::cout << "Intensity is being analyzed...Please wait." << std::endl;

  float intensitySDFix = 0.0;
  float intensitySDMov = 0.0;
  const auto meanIntensityFix = m_cbctregistration->m_pParent->GetMeanIntensity(
      m_spFixed, fROI_Radius, &intensitySDFix);
  const auto meanIntensityMov = m_cbctregistration->m_pParent->GetMeanIntensity(
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
  auto update_message =
      QString("Added_%1")
          .arg(static_cast<int>(meanIntensityMov - meanIntensityFix));
  SelectComboExternal(0, REGISTER_RAW_CBCT);
}
