#include "cbctrecon_mainwidget.h"

// Std
#include <iostream>
#include <thread>
#include <vector>

// Qt
#include <QtWidgets/QMainWindow>
#include <qclipboard.h>
#include <qfiledialog.h>
#include <qinputdialog.h>
#include <qmessagebox.h>
#include <qstandarditemmodel.h>
#include <qstring.h>
#include <qtimer.h>

// ITK
// #include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"
#include "itkTimeProbe.h"

// PLM
#include "itk_mask.h"
#include "mha_io.h"
#include "nki_io.h"

// Local
#include "DlgExternalCommand.h"
#include "DlgRegistration.h"
#include "cbctrecon.h"
#include "cbctrecon_compute.h"
#include "cbctrecon_io.h"
#include "cbctregistration.h"
#include "qcustomplot.h"

#pragma GCC poison new

CbctReconWidget::CbctReconWidget(QWidget *parent, const Qt::WindowFlags flags)
    : QMainWindow(parent, flags) {

  this->ui.setupUi(this);

  // Disable cuda & opencl as defaults
  this->ui.radioButton_UseCPU->setChecked(true);
  this->ui.radioButton_UseCUDA->setDisabled(true);
  this->ui.radioButton_UseOpenCL->setDisabled(true);

  this->ui.pushButton_DoRecon->setDisabled(true);

  connect(this->ui.labelImageRaw, SIGNAL(Mouse_Move()), this,
          SLOT(SLT_DataProbeProj()));
  connect(this->ui.labelImageRaw, SIGNAL(Mouse_Pressed_Left()), this,
          SLOT(SLT_CalculateROI_Proj())); // added

  connect(this->ui.labelReconImage, SIGNAL(Mouse_Move()), this,
          SLOT(SLT_DataProbeRecon()));
  connect(this->ui.labelReconImage, SIGNAL(Mouse_Pressed_Left()), this,
          SLOT(SLT_CalculateROI_Recon()));

  // Mouse_Left

  // connect(this->ui.MyLabel, SIGNAL(Mouse_Pos()), this, SLOT(MouseCurrPos()));

#if USE_CUDA
  this->ui.radioButton_UseCUDA->setDisabled(false);
  this->ui.radioButton_UseCUDA->setChecked(true);
#endif

#if USE_OPENCL_PLM
  this->ui.radioButton_UseOpenCL->setDisabled(false);
#endif

  this->m_cbctrecon = std::make_unique<CbctRecon>();
  m_dlgRegistration = std::make_unique<DlgRegistration>(this);
  m_cbctregistration = m_dlgRegistration->m_cbctregistration.get();

  QString tmp_folder("tmp");
  init_DlgRegistration(tmp_folder); // to Setup plastimatch folder. this is
                                    // useful if registration will be only done
  // m_pDlgHistogram = new DlgHistogram(this);
  this->m_dlgExternalCommand = std::make_unique<DlgExternalCommand>(this);

  // 20141017 QTIMER for sync
  m_Timer = std::make_unique<QTimer>(this);
  connect(m_Timer.get(), SIGNAL(timeout()), this, SLOT(SLT_TimerEvent()));
  m_busyTimer = false;

  this->m_cbctrecon->m_strPathDirDefault =
      R"(D:\Program_data\01_20140827_CBCT_All\04_patient_phan_pelvis_M\IMAGES\img_1.3.46.423632.135786.1409186054.9_M20mAs6440)";

  auto strPathCurAppDir = QDir::currentPath(); // should be same as .exe file
  std::cout << "Current app path= "
            << strPathCurAppDir.toLocal8Bit().constData() << std::endl;
  this->m_cbctrecon->m_strPathDefaultConfigFile =
      strPathCurAppDir + "/" + "DefaultConfig.cfg";

  if (!LoadCurrentSetting(
          this->m_cbctrecon->m_strPathDefaultConfigFile)) // Update GUI
  {
    std::cout << "DefaultConfig.cfg is not found in the application folder. A "
                 "new one will be created"
              << std::endl;
    if (!SaveCurrentSetting(this->m_cbctrecon->m_strPathDefaultConfigFile)) {
      std::cout << "Error in SaveCurrentSetting" << std::endl;
    }
  }
}

void CbctReconWidget::init_DlgRegistration(QString &str_dcm_uid) const
// init dlgRegistrations
{
  this->m_dlgRegistration->initDlgRegistration(
      str_dcm_uid); // NULLing all temporary spImage
}

void CbctReconWidget::SLT_LoadRawImages() { LoadRawHisImages(); }

void CbctReconWidget::SLT_DrawRawImages() const {
  const auto crntIdx = this->ui.spinBoxImgIdx->value();

  if (crntIdx >= this->m_cbctrecon->m_iImgCnt) {
    return;
  }

  const auto windowMin = this->ui.sliderRawMin->value();
  const auto windowMax = this->ui.sliderRawMax->value();

  auto tmpInfo =
      QFileInfo(this->m_cbctrecon->m_arrYKImage[crntIdx].m_strFilePath);
  this->ui.lineEditFileName->setText(tmpInfo.fileName());

  const auto width = this->m_cbctrecon->m_arrYKImage[crntIdx].m_iWidth;
  const auto height = this->m_cbctrecon->m_arrYKImage[crntIdx].m_iHeight;
  this->m_cbctrecon->m_dspYKImgProj->CreateImage(width, height, 0);
  this->m_cbctrecon->m_dspYKImgProj->CopyFromBuffer(
      this->m_cbctrecon->m_arrYKImage[crntIdx].m_pData, width, height);

  this->m_cbctrecon->m_dspYKImgProj->FillPixMapMinMax(windowMin, windowMax);
  this->ui.labelImageRaw->SetBaseImage(this->m_cbctrecon->m_dspYKImgProj.get());
  this->ui.labelImageRaw->update();
}

void CbctReconWidget::SLT_DrawProjImages() {
  if (this->m_cbctrecon->m_dspYKImgProj == nullptr) {
    return;
  }

  if (this->m_cbctrecon->m_iImgCnt > 0) {
    SLT_DrawRawImages();
    //		SLT_DrawGraph();
    SLT_UpdateTable();
    return;
  }

  // int iReqSlice = this->ui.spinBoxImgIdx->value();

  if (!this->m_cbctrecon->FillProjForDisplay(this->ui.spinBoxImgIdx->value())) {
    return;
  }

  this->m_cbctrecon->m_dspYKImgProj->FillPixMapMinMax(
      this->ui.sliderRawMin->value(), this->ui.sliderRawMax->value());

  this->ui.labelImageRaw->SetBaseImage(this->m_cbctrecon->m_dspYKImgProj.get());
  this->ui.labelImageRaw->update();

  SLT_UpdateTable();
}

void CbctReconWidget::SLT_FileNameHex2Dec() {
  auto files = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open",
      this->m_cbctrecon->m_strPathDirDefault, "projection images (*.his)");

  const auto cnt = files.size();
  if (cnt <= 0) {
    return;
  }

  const auto strMsg = "Original file names will be gone. Backup is strongly "
                      "recommended. Continue?";

  QMessageBox msgBox;
  msgBox.setText(strMsg);
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);

  const auto res = msgBox.exec();

  if (res == QMessageBox::Yes) {
    this->m_cbctrecon->RenameFromHexToDecimal(files);
  }
}

void CbctReconWidget::SLT_MakeElektaXML() {
  // Define IMAGE.DBF path
  auto filePath_ImageDBF = QFileDialog::getOpenFileName(
      this, "SelectIMAGE.DBF file", this->m_cbctrecon->m_strPathDirDefault,
      "Elekta DB file (*.dbf)", nullptr, nullptr);

  if (filePath_ImageDBF.length() < 2) {
    return;
  }

  auto filePath_FrameDBF = QFileDialog::getOpenFileName(
      this, "Select FRAME.DBF file", this->m_cbctrecon->m_strPathDirDefault,
      "Elekta DB file (*.dbf)", nullptr, nullptr);

  if (filePath_FrameDBF.length() < 2) {
    return;
  }

  QString DICOM_UID;
  QInputDialog inputDlg;

  bool ok;
  auto text = QInputDialog::getText(
      this, "Input Dialog", "DICOM_UID:", QLineEdit::Normal, "DICOM_UID", &ok);

  if (ok && !text.isEmpty()) {
    DICOM_UID = text;
  }

  if (DICOM_UID.length() < 2) {
    return;
  }

  auto genFilePath =
      MakeElektaXML(filePath_ImageDBF, filePath_FrameDBF, DICOM_UID);
  std::cout << "Generated ElektaXML path: "
            << genFilePath.toLocal8Bit().constData() << std::endl;
}

void CbctReconWidget::SLT_OpenOffsetFile() {
  // QString strPath = QFileDialog::getOpenFileNames(this,"Select one or more
  // files to open","/home","Images (*.raw)");
  auto strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open",
      this->m_cbctrecon->m_strPathDirDefault, "raw image (*.raw)");

  if (strPath.length() <= 1) {
    return;
  }
  auto stdstr_path = strPath.toStdString();
  this->m_cbctrecon->LoadCalibData(stdstr_path, OFFSET_CALIB);

  this->ui.lineEdit_offsetPath->setText(strPath);
}

void CbctReconWidget::SLT_OpenGainFile() {
  auto strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open",
      this->m_cbctrecon->m_strPathDirDefault, "raw image (*.raw)");

  if (strPath.length() <= 1) {
    return;
  }

  this->ui.lineEdit_gainPath->setText(strPath);
  auto stdstr_path = strPath.toStdString();
  this->m_cbctrecon->LoadCalibData(stdstr_path, GAIN_CALIB);
}

void CbctReconWidget::SLT_OpenBadpixelFile() {
  auto strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open",
      this->m_cbctrecon->m_strPathDirDefault, "bad pixel map file (*.pmf)");

  if (strPath.length() <= 1) {
    return;
  }

  this->ui.lineEdit_badpixelPath->setText(strPath);
  auto stdstr_path = strPath.toStdString();
  this->m_cbctrecon->LoadCalibData(stdstr_path, BADPIXEL_CALIB);
  // m_pImgGain->LoadRawImage(strPath.toLocal8Bit(),IMG_WIDTH, IMG_HEIGHT);
}

void CbctReconWidget::SLT_ApplyCalibration() const {
  if (this->m_cbctrecon->m_iImgCnt < 1) {
    return;
  }

  const auto bDarkCorrApply = this->ui.checkBox_offsetOn->isChecked();
  const auto bGainCorrApply = this->ui.checkBox_gainOn->isChecked();
  const auto bDefectMapApply = this->ui.checkBox_badpixelOn->isChecked();
  for (auto i = 0; i < this->m_cbctrecon->m_iImgCnt; i++) {
    this->m_cbctrecon->CorrectSingleFile(
        &this->m_cbctrecon->m_arrYKImage[i], bDarkCorrApply, bGainCorrApply,
        bDefectMapApply); // pixel value will be changed
  }
  SLT_DrawRawImages();
}

void CbctReconWidget::SLT_DrawReconImage() {
  if (this->m_cbctrecon->m_dspYKReconImage == nullptr) {
    return;
  }

  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    std::cout << "no recon image to be displayed" << std::endl;
    return;
  }

  //  The ExtractImageFilter type is instantiated using the input and
  //  output image types. A filter object is created with the New()
  //  method and assigned to a SmartPointer.

  using ExtractFilterType =
      itk::ExtractImageFilter<UShortImageType, UShortImage2DType>;
  auto extractFilter = ExtractFilterType::New();

  using DuplicatorType = itk::ImageDuplicator<UShortImageType>;
  auto duplicator = DuplicatorType::New();
  duplicator->SetInputImage(this->m_cbctrecon->m_spCrntReconImg);
  duplicator->Update();
  const auto clonedImage = duplicator->GetOutput();

  extractFilter->SetDirectionCollapseToSubmatrix();

  auto crnt_region_3d = clonedImage->GetBufferedRegion();

  //  We take the size from the region and collapse the size in the $Z$
  //  component by setting its value to $1$.

  // Get Image Size and Extraction Index info.
  auto size = crnt_region_3d.GetSize();
  size[2] = 0; // z size number = 0 --> should not be 1

  auto start = crnt_region_3d.GetIndex();
  const auto iSliceNumber = this->ui.spinBoxReconImgSliceNo->value();
  start[2] = iSliceNumber; // 60

  const auto originZ = this->m_cbctrecon->m_spCrntReconImg->GetOrigin()[2];
  const auto spacingZ = this->m_cbctrecon->m_spCrntReconImg->GetSpacing()[2];
  const auto posZ = originZ + iSliceNumber * spacingZ;

  const auto strPosZ = QString("%1").arg(posZ, 0, 'f', 2);
  // strPosZ.sprintf("%4.2f", posZ);
  this->ui.lineEdit_CurrentPosZ->setText(strPosZ);

  // Define a region to generate
  UShortImageType::RegionType desiredRegion;
  desiredRegion.SetSize(size);   // 410 410 0
  desiredRegion.SetIndex(start); // 0 0 60

  // Error occurred here --> sloved by crntRegion3D =
  // m_spReconImg->GetBufferedRegion();
  extractFilter->SetExtractionRegion(desiredRegion); // error

  //  Below we connect the reader, filter and writer to form the data
  //  processing pipeline.
  extractFilter->SetInput(clonedImage);

  extractFilter->Update();

  UShortImage2DType::Pointer pCrnt2D = extractFilter->GetOutput();
  this->m_cbctrecon->m_dspYKReconImage = YK16GrayImage::CopyItkImage2YKImage(
      pCrnt2D,
      std::move(this->m_cbctrecon->m_dspYKReconImage)); // dimension should be
                                                        // same automatically.

  // m_dspYKReconImage->SaveDataAsRaw("D:\\RawFile.raw"); //410 410 OK

  const auto physPosX = this->ui.lineEdit_PostFOV_X->text().toFloat();
  const auto physPosY = this->ui.lineEdit_PostFOV_Y->text().toFloat();
  const auto physRadius = this->ui.lineEdit_PostFOV_R->text().toFloat();
  const auto physTablePosY = this->ui.lineEdit_PostTablePosY->text().toFloat();
  this->m_cbctrecon->PostApplyFOVDispParam(physPosX, physPosY, physRadius,
                                           physTablePosY);
  // SLT_UpdatePostProcDispObj();

  if (this->ui.checkBox_PostDispObjOn->isChecked()) {
    this->m_cbctrecon->m_dspYKReconImage->m_bDrawFOVCircle = true;
    this->m_cbctrecon->m_dspYKReconImage->m_bDrawTableLine = true;
  }

  else {
    this->m_cbctrecon->m_dspYKReconImage->m_bDrawFOVCircle = false;
    this->m_cbctrecon->m_dspYKReconImage->m_bDrawTableLine = false;
  }

  this->m_cbctrecon->m_dspYKReconImage->FillPixMapMinMax(
      this->ui.sliderReconImgMin->value(), this->ui.sliderReconImgMax->value());
  this->ui.labelReconImage->SetBaseImage(
      this->m_cbctrecon->m_dspYKReconImage.get());
  this->ui.labelReconImage->update();

  // SLT_DrawGraph();
  SLT_UpdateTable();
}

void CbctReconWidget::SLT_OpenElektaGeomFile() {
  auto strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open",
      this->m_cbctrecon->m_strPathDirDefault, "Geometry file (*.xml)");

  if (strPath.length() <= 1) {
    return;
  }

  this->ui.lineEdit_ElektaGeomPath->setText(strPath);
}

void CbctReconWidget::SLT_SetOutputPath() {
  auto strPath = QFileDialog::getSaveFileName(
      this, "File path to save", "D:\\", "meta 3D image data (*.mha)", nullptr,
      nullptr); // Filename don't need to exist

  if (strPath.length() <= 1) {
    return;
  }

  this->ui.lineEdit_OutputFilePath->setText(strPath);
}

FDK_options CbctReconWidget::getFDKoptions() const {
  FDK_options fdk_options;

  fdk_options.TruncCorFactor =
      this->ui.lineEdit_Ramp_TruncationCorrection->text().toDouble();
  fdk_options.HannCutX = this->ui.lineEdit_Ramp_HannCut->text().toDouble();
  fdk_options.HannCutY = this->ui.lineEdit_Ramp_HannCutY->text().toDouble();
  fdk_options.CosCut = this->ui.lineEdit_Ramp_CosineCut->text().toDouble();
  fdk_options.HammCut = this->ui.lineEdit_Ramp_Hamming->text().toDouble();

  fdk_options.displacedDetectorFilter = this->ui.checkBox_UseDDF->isChecked();
  fdk_options.updateAfterDDF =
      this->ui.checkBox_UpdateAfterFiltering->isChecked();
  fdk_options.ParkerShortScan = this->ui.checkBox_UsePSSF->isChecked();

  fdk_options.ct_spacing[0] = this->ui.lineEdit_outImgSp_AP->text().toDouble();
  fdk_options.ct_spacing[1] = this->ui.lineEdit_outImgSp_SI->text().toDouble();
  fdk_options.ct_spacing[2] = this->ui.lineEdit_outImgSp_LR->text().toDouble();

  fdk_options.ct_size[0] = this->ui.lineEdit_outImgDim_AP->text().toInt();
  fdk_options.ct_size[1] = this->ui.lineEdit_outImgDim_SI->text().toInt();
  fdk_options.ct_size[2] = this->ui.lineEdit_outImgDim_LR->text().toInt();

  fdk_options.medianRadius[0] =
      this->ui.lineEdit_PostMedSizeX->text().toInt(); // radius along x
  fdk_options.medianRadius[1] =
      this->ui.lineEdit_PostMedSizeY->text().toInt(); // radius along y
  fdk_options.medianRadius[2] =
      this->ui.lineEdit_PostMedSizeZ->text().toInt(); // radius along z
  fdk_options.medianFilter = this->ui.checkBox_PostMedianOn->isChecked();

  fdk_options.outputFilePath = this->ui.lineEdit_OutputFilePath->text();

  return fdk_options;
}

void CbctReconWidget::SLT_DoReconstruction() {
  const auto fdk_options = getFDKoptions();

  itk::TimeProbe reconTimeProbe;
  reconTimeProbe.Start();

  if (this->ui.radioButton_UseCUDA->isChecked()) {
    this->m_cbctrecon->DoReconstructionFDK<CUDA_DEVT>(REGISTER_RAW_CBCT,
                                                      fdk_options);
  } else if (this->ui.radioButton_UseOpenCL->isChecked()) {
    this->m_cbctrecon->DoReconstructionFDK<OPENCL_DEVT>(REGISTER_RAW_CBCT,
                                                        fdk_options);
  } else {
    this->m_cbctrecon->DoReconstructionFDK<CPU_DEVT>(REGISTER_RAW_CBCT,
                                                     fdk_options);
  }

  reconTimeProbe.Stop();
  std::cout << "It took " << reconTimeProbe.GetMean() << ' '
            << reconTimeProbe.GetUnit() << std::endl;
  this->ui.lineEdit_ReconstructionTime->setText(
      QString("%1").arg(reconTimeProbe.GetMean()));

  this->ui.spinBoxReconImgSliceNo->setMinimum(0);
  this->ui.spinBoxReconImgSliceNo->setMaximum(fdk_options.ct_size[1] - 1);
  this->ui.spinBoxReconImgSliceNo->setValue(qRound(
      fdk_options.ct_size[1] / 2.0)); // DrawReconImage is called automatically

  SLT_DrawProjImages();

  QString update_text("RAW_CBCT");
  UpdateReconImage(this->m_cbctrecon->m_spCrntReconImg, update_text);

  m_dlgRegistration->UpdateListOfComboBox(0); // combo selection
                                              // signalis called
  m_dlgRegistration->UpdateListOfComboBox(1);
  // m_pDlgRegistration->SelectComboExternal(0, REGISTER_RAW_CBCT); // will call
  // fixedImageSelected  m_pDlgRegistration->SelectComboExternal(1,
  // REGISTER_RAW_CBCT );

  // After first reconstruction, set Median size to 0 0 1 for scatter corrected
  // solution
  /* this->ui.lineEdit_PostMedSizeX->setText(QString("%1").arg(0.0));
  this->ui.lineEdit_PostMedSizeY->setText(QString("%1").arg(0.0));
  this->ui.lineEdit_PostMedSizeZ->setText(QString("%1").arg(1.0));*/
}

void CbctReconWidget::SLT_InitializeGraphLim() const {
  // Set Max Min at graph
  if (this->ui.radioButton_graph_proj->isChecked()) {
    if (this->m_cbctrecon->m_iImgCnt > 0) // if indep raw his images are loaded
    {
      const auto horLen = this->m_cbctrecon->m_dspYKImgProj->m_iWidth;
      // int verLen = m_dspYKImgProj->m_iHeight;

      // set edit maxium min
      const auto strXMin = QString("%1").arg(horLen);
      this->ui.lineEditXMin->setText("0");
      this->ui.lineEditXMax->setText(strXMin);

      const auto strYMin =
          QString("%1").arg(this->m_cbctrecon->m_fProjImgValueMin, 0, 'f', 1);
      const auto strYMax =
          QString("%1").arg(this->m_cbctrecon->m_fProjImgValueMax, 0, 'f', 1);

      this->ui.lineEditYMin->setText(strYMin);
      this->ui.lineEditYMax->setText(strYMax);
    }

    if (this->m_cbctrecon->m_spProjImg3DFloat == nullptr) {
      return;
    }

    const auto horLen =
        this->m_cbctrecon->m_spProjImg3DFloat->GetBufferedRegion().GetSize()[0];
    // int verLen = m_spProjImg3DFloat->GetBufferedRegion().GetSize()[1];

    // set edit maxium min
    const auto strXMin = QString("%1").arg(horLen);
    this->ui.lineEditXMin->setText("0");
    this->ui.lineEditXMax->setText(strXMin);

    const auto strYMin =
        QString("%1").arg(this->m_cbctrecon->m_fProjImgValueMin, 0, 'f', 1);
    const auto strYMax =
        QString("%1").arg(this->m_cbctrecon->m_fProjImgValueMax, 0, 'f', 1);

    this->ui.lineEditYMin->setText(strYMin);
    this->ui.lineEditYMax->setText(strYMax);
  } else if (this->ui.radioButton_graph_recon->isChecked()) {
    if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
      return;
    }

    const auto horLen =
        this->m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize()[0];
    // int verLen = m_spCrntReconImg->GetBufferedRegion().GetSize()[1];

    // set edit maxium min

    const auto strXMax = QString("%1").arg(horLen);
    this->ui.lineEditXMin->setText("0");
    this->ui.lineEditXMax->setText(strXMax);

    const auto strYMin = QString("%1").arg(0.0, 0, 'f', 1);
    const auto strYMax = QString("%1").arg(2000.0, 0, 'f', 1);

    this->ui.lineEditYMin->setText(strYMin);
    this->ui.lineEditYMax->setText(strYMax);
  }
}

void CbctReconWidget::SLT_CopyTableToClipBoard() const {
  // qApp->clipboard()->clear();
  auto clipboard = QApplication::clipboard();
  clipboard->clear();

  QStringList list;

  const auto rowCnt = m_pTableModel->rowCount();
  const auto columnCnt = m_pTableModel->columnCount();

  list << "\n";
  // for (int i = 0 ; i < columnCnt ; i++)
  //{
  const auto tmpInfo = QFileInfo(this->ui.lineEdit_Cur3DFileName->text());
  // list << "Index";
  list << tmpInfo.baseName();
  list << "\n";

  list << "Pos(mm)";
  list << "Intensity";
  list << "\n";

  for (auto j = 0; j < rowCnt; j++) {
    for (auto i = 0; i < columnCnt; i++) {
      auto item = m_pTableModel->item(j, i);
      list << item->text();
    }
    list << "\n";
  }

  clipboard->setText(list.join("\t"));
  // qApp->clipboard()->setText(list.join("\t"));
}

void CbctReconWidget::SLT_SetHisDir() // Initialize all image buffer
{
  // Initializing..

  // Set folder --> then use RTK HIS Reader
  auto dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), this->m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.length() <= 1) {
    return;
  }

  this->ui.lineEdit_HisDirPath->setText(dirPath);

  this->m_cbctrecon->SetProjDir(dirPath);
  init_DlgRegistration(this->m_cbctrecon->m_strDCMUID);

  this->ui.lineEdit_ElektaGeomPath->setText(
      this->m_cbctrecon->m_strPathGeomXML);
  this->ui.lineEdit_PathCBCTSkinPath->setText(
      this->m_cbctrecon->m_strPathRS_CBCT);

  float kVp = 0.0;
  float mA = 0.0;
  float ms = 0.0;
  GetXrayParamFromINI(this->m_cbctrecon->m_strPathElektaINI, kVp, mA, ms);

  if (kVp * mA * ms != 0) {
    // update GUI
    std::cout << "Updating current mAs setting from INI file: "
              << "kVp= " << kVp << ", mA= " << mA << ", ms= " << ms
              << std::endl;
  }
  this->ui.lineEdit_CurmAs->setText(QString("%1, %2").arg(mA).arg(ms));

  VEC3D couch_trans = {-999, -999,
                       -999}; // mm. In the text file, these values are in cm.
  VEC3D couch_rot = {-999, -999,
                     -999}; // mm. In the text file, these values are in cm.

  const auto res = GetCouchShiftFromINIXVI(
      this->m_cbctrecon->m_strPathElektaINIXVI2, &couch_trans, &couch_rot);

  if (res) {
    const auto strTransX = QString::number(couch_trans.x, 'f', 1);
    const auto strTransY = QString::number(couch_trans.y, 'f', 1);
    const auto strTransZ = QString::number(couch_trans.z, 'f', 1);
    const auto strTransAll = strTransX + "," + strTransY + "," + strTransZ;

    const auto strRotX = QString::number(couch_rot.x, 'f', 1);
    const auto strRotY = QString::number(couch_rot.y, 'f', 1);
    const auto strRotZ = QString::number(couch_rot.z, 'f', 1);

    const auto strRotAll = strRotX + "," + strRotY + "," + strRotZ;

    this->ui.lineEdit_CouchTrans->setText(strTransAll);
    this->ui.lineEdit_CouchRot->setText(strRotAll);
  } else {
    this->ui.lineEdit_CouchTrans->setText("Not available");
    this->ui.lineEdit_CouchRot->setText("Not available");
  }

  this->m_cbctrecon->m_vSelectedFileNames.clear();

  std::cout << "Push Load button to load projection images" << std::endl;
}

QString getBowtiePath(QWidget *parent, const QDir &calDir) {
  return QFileDialog::getOpenFileName(
      parent, "Find air(+bowtie) filter image for subtraction",
      calDir.absolutePath(), "Projection (*.xim)", nullptr, nullptr);
}

std::tuple<bool, bool> CbctReconWidget::probeUser(const QString &guessDir) {

  auto dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open CT DICOM Directory"), guessDir,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  auto dcm_success = false;
  if (!(dirPath.length() <= 1)) {

    if (this->m_cbctrecon->ReadDicomDir(dirPath)) {

      m_dlgRegistration->UpdateVOICombobox(PLAN_CT);
      // UpdateReconImage(m_spRefCTImg, QString("DICOM reference image"));

      this->m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT,
                                                REGISTER_MANUAL_RIGID);
      dcm_success = true;
    }
  }
  auto instaRecon = false;
  const auto reply =
      QMessageBox::question(this, "Instant Recon?",
                            "Do you want to reconstruct projections as soon as "
                            "they are loaded?\n(Using the current settings)",
                            QMessageBox::Yes | QMessageBox::No);
  if (reply == QMessageBox::Yes) {
    instaRecon = true;
  }

  return std::make_tuple(instaRecon, dcm_success);
}

FilterReaderType::Pointer
CbctReconWidget::ReadBowtieFileWhileProbing(const QString &proj_path,
                                            std::tuple<bool, bool> &answers) {

  auto bowtiereader =
      FilterReaderType::New(); // we use is because we need the projections to
                               // be in the same unit (order of magnitude)

  QDir guessDir(proj_path + QString("/../"));

  const auto calDir(proj_path + QString("/Calibrations/"));

  QString bowtiePath;

  switch (this->m_cbctrecon->m_projFormat) {
  case XIM_FORMAT:
    bowtiePath = getBowtiePath(this, calDir);
    if (bowtiePath.length() > 1) {
      std::cout << "loading bowtie-filter..." << std::endl;
      std::vector<std::string> filepath;
      filepath.push_back(bowtiePath.toStdString());
      bowtiereader->SetFileNames(filepath);
      // std::thread calc_thread_bowtie(read_bowtie_projection, bowtiereader);
      std::thread calc_thread_bowtie(
          [&bowtiereader] { bowtiereader->Update(); });
      answers = probeUser(guessDir.absolutePath());
      calc_thread_bowtie.join();
    } else {
      answers = probeUser(
          guessDir.absolutePath()); // looks ugly, but allows threading
    }
    break;
  default:
    answers = probeUser(guessDir.absolutePath()); // ^^^
    break;
  }
  if (bowtiePath.length() > 1) {
    return bowtiereader;
  }
  return nullptr;
}

void CbctReconWidget::SLT_LoadSelectedProjFiles() // main loading fuction for
                                                  // projection images
{
  this->ui.pushButton_DoRecon->setDisabled(true);
  // 1) Get all projection file names
  auto dirPath = this->ui.lineEdit_HisDirPath->text();
  //.toLocal8Bit().constData();

  if (!QFile::exists(dirPath)) {
    std::cout << "Projection file directory was not found. Retry." << std::endl;
    return;
  }

  auto names = this->m_cbctrecon->GetProjFileNames(dirPath);

  if (!this->m_cbctrecon->IsFileNameOrderCorrect(names) &&
      this->m_cbctrecon->m_projFormat != XIM_FORMAT) {
    std::cout << "Check the file name order" << std::endl;
    QMessageBox::warning(this, "warning", "Error on File Name Sorting!");
    return;
  }

  std::cout << "File name order was cross-checked and found to be OK!"
            << std::endl;

  const auto fullCnt = names.size();
  if (fullCnt <= 0) {
    std::cout << "No projection file was found. Retry." << std::endl;
    return;
  }

  std::cout << fullCnt << "  projection files were found." << std::endl;

  // 2) Elekta Geometry file
  const auto geomPath = this->ui.lineEdit_ElektaGeomPath->text();
  QFileInfo geomFileInfo(geomPath);
  //! QFile::exists(geomPath)

  if (!this->m_cbctrecon->LoadGeometry(geomFileInfo, names)) {
    if (!this->m_cbctrecon->m_strError.isEmpty()) {
      QMessageBox::critical(this, "LoadXVIGeometryFile",
                            this->m_cbctrecon->m_strError, QMessageBox::Ok);
    }
  }

  const auto iFullGeoDataSize =
      this->m_cbctrecon->m_spFullGeometry->GetGantryAngles().size();
  if (iFullGeoDataSize < 1) {
    std::cout << "Not enough projection image (should be > 0)" << std::endl;
    return;
  }

  if (iFullGeoDataSize != fullCnt) {
    if (this->m_cbctrecon->m_projFormat != XIM_FORMAT) {
      std::cout << "Size of geometry data and file numbers are not same! Check "
                   "and retry"
                << std::endl;
      return;
    }

    const auto reply = QMessageBox::question(
        this, "Mismatch in number of files and Geometry information!",
        "Mismatch in number of files and Geometry information!\nHowever, Xim "
        "detected, so it may be safe to continue anyway?",
        QMessageBox::Yes | QMessageBox::No);
    if (reply == QMessageBox::Yes) {
      std::cout << "continuing despite warning..." << std::endl;
    } else {
      return;
    }
  }

  auto angle_gaps = this->m_cbctrecon->m_spFullGeometry->GetAngularGaps(
      this->m_cbctrecon->m_spFullGeometry->GetSourceAngles());

  auto sum_gap =
      std::accumulate(std::begin(angle_gaps), std::end(angle_gaps), 0.0);
  sum_gap /= itk::Math::pi * 180.0;
  const auto mean_gap = sum_gap / angle_gaps.size();

  std::cout << "AngularGaps Sum (deg):" << sum_gap
            << ", Mean (deg): " << mean_gap << std::endl;

  const auto gantryAngleInterval =
      this->ui.lineEdit_ManualProjAngleGap->text().toDouble();

  // if (this->ui.Radio_KeepOriginalAngles->isChecked())
  if (this->ui.Radio_ManualProjAngleGap->isChecked()) {
    // bManualGap = true;
    // std::cout << "Input angle gap in deg: " ;
    // cin >> gantryAngleInterval;

    if (gantryAngleInterval < mean_gap) {
      std::cout << "Angle gap size is too small. Terminating the app"
                << std::endl;
      return;
      // bManualGap = false;
    }
  }

  const auto exclude_ids = this->m_cbctrecon->GetExcludeProjFiles(
      this->ui.Radio_ManualProjAngleGap->isChecked(), gantryAngleInterval);

  this->m_cbctrecon->LoadSelectedProj(exclude_ids, names);

  // Reads the cone beam projections
  using ReaderType = rtk::ProjectionsReader<FloatImageType>;
  auto reader = ReaderType::New();
  reader->SetFileNames(this->m_cbctrecon->m_vSelectedFileNames);
  // TRY_AND_EXIT_ON_ITK_EXCEPTION(
  // std::thread calc_thread(read_projections, reader);
  std::thread calc_thread([&reader]() { reader->Update(); });
  // calc_thread.detach();

  std::cout << "Reader detached from main thread" << std::endl;

  // After reading the whole file,
  // HIS header should be saved
  this->m_cbctrecon->saveHisHeader();

  //  Insta Recon, Dcm read
  const auto geopath = geomFileInfo.absolutePath();
  std::tuple<bool, bool> answers;
  auto bowtie_reader = ReadBowtieFileWhileProbing(geopath, answers);

  calc_thread.join();
  std::cout << "Reader re-attached to main thread" << std::endl;

  if (bowtie_reader != nullptr) {
    ApplyBowtie(reader, bowtie_reader);
  }
  if (this->m_cbctrecon->m_projFormat == HND_FORMAT) {
    std::cout << "Fitted bowtie-filter correction ongoing..." << std::endl;
    SLT_DoBowtieCorrection();
  }

  saveImageAsMHA<FloatImageType>(this->m_cbctrecon->m_spProjImg3DFloat);
  auto res_factor = this->ui.lineEdit_DownResolFactor->text().toDouble();
  if (!this->m_cbctrecon->ResampleProjections(res_factor)) { // 0.5
    // reset factor if image was not resampled
    this->ui.lineEdit_DownResolFactor->setText("1.0");
  }

  this->m_cbctrecon->ConvertLineInt2Intensity(
      this->m_cbctrecon->m_spProjImg3DFloat,
      this->m_cbctrecon->m_spProjImgRaw3D,
      65535); // if X not 1024 == input size: out_offset =
              // in_offset + (1024*res_f -
              // X*res_f)*out_spacing     <- will still
              // break down at fw_projection

  this->ui.pushButton_DoRecon->setEnabled(true);

  this->ui.spinBoxImgIdx->setMinimum(0);
  this->ui.spinBoxImgIdx->setMaximum(
      this->m_cbctrecon->m_vSelectedFileNames.size() - 1);
  this->ui.spinBoxImgIdx->setValue(0); // it doesn't call Draw Event .. don't
                                       // know why.

  this->m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image

  SLT_InitializeGraphLim();

  this->SLT_DrawProjImages(); // Update Table is called

  if (!std::get<0>(answers)) { // instaRecon
    std::cout
        << "FINISHED!: Loading projection files. Proceed to reconstruction"
        << std::endl;
  } else {
    SLT_DoReconstruction();
  }

  if (std::get<1>(answers)) { // CT DCM dir was found
    SLT_ViewRegistration();
  }
}

void CbctReconWidget::SLT_DataProbeProj() const {
  const auto dspWidth = this->ui.labelImageRaw->width();
  const auto dspHeight = this->ui.labelImageRaw->height();

  if (this->m_cbctrecon->m_iImgCnt >
      0) // there is indep loaded projection files
  {
    const auto dataWidth = this->m_cbctrecon->m_dspYKImgProj->m_iWidth;
    const auto dataHeight = this->m_cbctrecon->m_dspYKImgProj->m_iHeight;

    const auto dataX = qRound(this->ui.labelImageRaw->x /
                              static_cast<double>(dspWidth * dataWidth));
    const auto dataY = qRound(this->ui.labelImageRaw->y /
                              static_cast<double>(dspHeight * dataHeight));
    const auto dataZ = this->ui.spinBoxImgIdx->value();
    const auto fProbeValue = static_cast<double>(
        this->m_cbctrecon->m_dspYKImgProj->m_pData[dataWidth * dataY + dataX]);
    const auto dspText = QString("(%1, %2, %3): %4")
                             .arg(dataX)
                             .arg(dataY)
                             .arg(dataZ)
                             .arg(fProbeValue, 0, 'f', 2);
    this->ui.lineEdit_DataProbe_Proj->setText(dspText);
  } else {
    if (this->m_cbctrecon->m_spProjImg3DFloat == nullptr) {
      return;
    }

    const auto dataWidth = static_cast<int>(
        this->m_cbctrecon->m_spProjImg3DFloat->GetBufferedRegion()
            .GetSize()[0]);
    const auto dataHeight = static_cast<int>(
        this->m_cbctrecon->m_spProjImg3DFloat->GetBufferedRegion()
            .GetSize()[1]);

    // int crntIdx = this->ui.spinBoxImgIdx->value();
    // These are displayed data (just index data)
    const auto dataX = qRound(this->ui.labelImageRaw->x /
                              static_cast<double>(dspWidth * dataWidth));
    const auto dataY = qRound(this->ui.labelImageRaw->y /
                              static_cast<double>(dspHeight * dataHeight));
    const auto dataZ = this->ui.spinBoxImgIdx->value();

    // fProbeValue = m_dspYKImgProj->m_pData[dataWidth*dataY +
    // dataX]/m_multiplyFactor;
    const auto fProbeValue =
        static_cast<double>(this->m_cbctrecon->m_dspYKImgProj
                                ->m_pData[dataWidth * dataY + dataX]) /
            this->m_cbctrecon->m_multiplyFactor +
        this->m_cbctrecon->m_fProjImgValueMin;
    const auto dspText = QString("(%1, %2, %3): %4")
                             .arg(dataX)
                             .arg(dataY)
                             .arg(dataZ)
                             .arg(fProbeValue, 0, 'f', 2);
    this->ui.lineEdit_DataProbe_Proj->setText(dspText);
  }
}

void CbctReconWidget::SLT_DataProbeRecon() const {
  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  const auto dspWidth = this->ui.labelReconImage->width();
  const auto dspHeight = this->ui.labelReconImage->height();

  const auto dataWidth =
      this->m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize()[0];
  const auto dataHeight =
      this->m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize()[1];

  // int crntIdx = this->ui.spinBoxImgIdx->value();
  // These are displayed data (just index data)

  const auto dataX =
      std::lround(this->ui.labelReconImage->x / dspWidth * dataWidth);
  const auto dataY =
      std::lround(this->ui.labelReconImage->y / dspHeight * dataHeight);
  const auto dataZ = this->ui.spinBoxReconImgSliceNo->value();

  const auto iProbeValue =
      this->m_cbctrecon->m_dspYKReconImage->m_pData[dataWidth * dataY + dataX];
  // unsigned short iProbeValue = GetValueFrom3DImageUshort(dataX, dataY, dataZ,
  // m_spReconImg);

  const auto zero = QChar('0');
  const auto dspText = QString("(%1, %2, %3): %4")
                           .arg(dataX, 3, 10, zero)
                           .arg(dataY, 3, 10, zero)
                           .arg(dataZ, 3, 10, zero)
                           .arg(iProbeValue);
  // dspText.sprintf("(%03d, %03d, %03d): %d", dataX, dataY, dataZ,
  // iProbeValue);
  this->ui.lineEdit_DataProbe_Recon->setText(dspText);
}

void CbctReconWidget::SLT_DrawGraph() const
// based on profile
{
  if (m_pTableModel == nullptr) {
    return;
  }

  // Draw only horizontal, center

  QVector<double> vAxisX; // can be rows or columns
  QVector<double> vAxisY;

  // QStandardItemModel 	m_pTableModel.item()
  const auto dataLen = m_pTableModel->rowCount();

  if (dataLen < 1) {
    return;
  }

  // std::cout << "check graph 1" << std::endl;
  this->ui.customPlot->clearGraphs();

  auto minX = 9999.0;
  auto maxX = -1.0;

  for (auto i = 0; i < dataLen; i++) {
    auto tableItem1 = m_pTableModel->item(i, 0);
    auto tableItem2 = m_pTableModel->item(i, 1);
    auto tableVal1 = tableItem1->text().toDouble();
    auto tableVal2 = tableItem2->text().toDouble();

    if (minX > tableVal1) {
      minX = tableVal1;
    }
    if (maxX < tableVal1) {
      maxX = tableVal1;
    }

    vAxisX.push_back(tableVal1);
    vAxisY.push_back(tableVal2);
  }

  // std::cout << "check graph 2" << std::endl;

  this->ui.customPlot->addGraph();
  this->ui.customPlot->graph(0)->setData(vAxisX, vAxisY);
  this->ui.customPlot->graph(0)->setPen(QPen(Qt::blue));
  this->ui.customPlot->graph(0)->setName("Image profile");

  this->ui.lineEditXMin->setText(QString("%1").arg(minX));
  this->ui.lineEditXMax->setText(QString("%1").arg(maxX));

  const auto tmpXMin = this->ui.lineEditXMin->text().toDouble();
  const auto tmpXMax = this->ui.lineEditXMax->text().toDouble();
  const auto tmpYMin = this->ui.lineEditYMin->text().toDouble();
  const auto tmpYMax = this->ui.lineEditYMax->text().toDouble();

  // std::cout << "check graph 3" << std::endl;

  this->ui.customPlot->xAxis->setRange(tmpXMin, tmpXMax);
  this->ui.customPlot->yAxis->setRange(tmpYMin, tmpYMax);

  this->ui.customPlot->xAxis->setLabel("mm");
  this->ui.customPlot->yAxis->setLabel("Intensity");
  this->ui.customPlot->setTitle("Image Profile");
  auto titleFont = font();
  titleFont.setPointSize(10);
  this->ui.customPlot->setTitleFont(titleFont);

  // std::cout << "check graph 4" << std::endl;

  this->ui.customPlot->legend->setVisible(false);
  auto legendFont = font();   // start out with MainWindow's font..
  legendFont.setPointSize(9); // and make a bit smaller for legend
  this->ui.customPlot->legend->setFont(legendFont);
  this->ui.customPlot->legend->setPositionStyle(QCPLegend::psTopRight);
  this->ui.customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 200)));

  // std::cout << "check graph 5" << std::endl;
  this->ui.customPlot->replot();

  // SLT_UpdateTable();
}

void CbctReconWidget::SLT_UpdateTable() {
  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    this->ui.radioButton_graph_proj->setChecked(true);
  }

  // std::cout << "check 1" << std::endl;
  YK16GrayImage *pYKImg;
  auto fMultiPlyFactor = 1.0;
  auto fMinValue = 0.0;

  if (this->ui.radioButton_graph_proj->isChecked()) {
    pYKImg = this->m_cbctrecon->m_dspYKImgProj
                 .get(); // you may look, but no touching!

    if (this->m_cbctrecon->m_iImgCnt > 0) { // if indep image
      fMultiPlyFactor = 1.0;
    } else {
      fMultiPlyFactor = this->m_cbctrecon->m_multiplyFactor;
      fMinValue = this->m_cbctrecon->m_fProjImgValueMin;
    }
  } else {
    pYKImg = this->m_cbctrecon->m_dspYKReconImage.get();
    fMinValue = 0.0;
  }
  if (pYKImg == nullptr) {
    return;
  }

  // std::cout << "check 2" << std::endl;

  // std::cout << "check 3" << std::endl;
  auto columnSize = 2;
  auto rowSize = pYKImg->m_iHeight;

  /// int rowSize = pYKImg->m_iWidth;

  if (this->ui.radioButton_Profile_Hor->isChecked()) {
    // columnSize = 2;
    rowSize = pYKImg->m_iWidth;
  }

  // std::cout << "check 4" << std::endl;
  m_pTableModel.reset();
  m_pTableModel = std::make_unique<QStandardItemModel>(
      rowSize, columnSize, this); // 2 Rows and 3 Columns

  // for (int i = 0 ; i<columnSize ; i++)
  //{
  // QFileInfo tmpInfo = QFileInfo(m_arrYKImage[i].m_strFilePath);
  // m_pTableModel->setHorizontalHeaderItem(0, new
  // QStandardItem(QString("Index"))); m_pTableModel->setHorizontalHeaderItem(0,
  // new QStandardItem(QString("Profile")));
  auto pos_item = std::make_unique<QStandardItem>(QString("Position(mm)"));
  auto val_item = std::make_unique<QStandardItem>(QString("Value"));

  m_pTableModel->setHorizontalHeaderItem(0, pos_item.get());
  m_pTableModel->setHorizontalHeaderItem(1, val_item.get());
  //}

  // std::cout << "check 5" << std::endl;
  // int width = pYKImg->m_iWidth;
  // int height = pYKImg->m_iHeight;
  // int fixedY = qRound(height / 2.0);

  auto originX = 0.0;
  auto originY = 0.0;
  auto spacingX = 1.0;
  auto spacingY = 1.0;

  if (!this->ui.radioButton_graph_proj->isChecked()) {
    if (this->m_cbctrecon->m_spCrntReconImg != nullptr) {
      auto tmpOrigin = this->m_cbctrecon->m_spCrntReconImg->GetOrigin();
      auto tmpSpacing = this->m_cbctrecon->m_spCrntReconImg->GetSpacing();
      originX = tmpOrigin[0];
      originY = tmpOrigin[1];
      spacingX = tmpSpacing[0];
      spacingY = tmpSpacing[1];
    }
  }

  // std::cout << "check 6" << std::endl;

  QVector<qreal> vPos;
  if (this->ui.radioButton_Profile_Hor->isChecked()) {
    for (auto i = 0; i < rowSize; i++) {
      vPos.push_back(originX + i * spacingX);
    }
  } else {
    for (auto i = 0; i < rowSize; i++) {
      vPos.push_back(originY + i * spacingY);
    }
  }

  QVector<qreal> vProfile;
  if (this->ui.radioButton_Profile_Hor->isChecked()) {
    pYKImg->GetProfileData(vProfile, DIRECTION_HOR);
  } else {
    pYKImg->GetProfileData(vProfile, DIRECTION_VER);
  }

  // int i = fixedY;
  for (auto i = 0; i < rowSize; i++) {
    const auto tmpVal1 = vPos[i];
    auto xpos_item =
        std::make_unique<QStandardItem>(QString("%1").arg(tmpVal1));
    m_pTableModel->setItem(i, 0, xpos_item.get());

    const auto tmpVal2 = vProfile[i] / fMultiPlyFactor + fMinValue;
    auto profval_item =
        std::make_unique<QStandardItem>(QString("%1").arg(tmpVal2));
    m_pTableModel->setItem(i, 1, profval_item.get());
  }

  this->ui.tableViewReconImgProfile->setModel(
      m_pTableModel.get()); // also for proj

  // std::cout << "check 7" << std::endl;
  SLT_DrawGraph();
}

// Mouse Left Click
void CbctReconWidget::SLT_CalculateROI_Recon() {
  if (this->m_cbctrecon->m_dspYKReconImage == nullptr) {
    return;
  }

  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }
  auto &CrntReconImg = this->m_cbctrecon->m_spCrntReconImg;
  auto &dspReconImg = this->m_cbctrecon->m_dspYKReconImage;

  const auto dspWidth = this->ui.labelReconImage->width();
  const auto dspHeight = this->ui.labelReconImage->height();

  const auto dataWidth = dspReconImg->m_iWidth;
  const auto dataHeight = dspReconImg->m_iHeight;
  if (dataWidth * dataHeight == 0) {
    return;
  }

  // int crntIdx = this->ui.spinBoxImgIdx->value();
  // These are displayed data (just index data)

  const auto dataX = qRound(this->ui.labelReconImage->x /
                            static_cast<double>(dspWidth * dataWidth));
  const auto dataY = qRound(this->ui.labelReconImage->y /
                            static_cast<double>(dspHeight * dataHeight));
  const auto dataZ = this->ui.spinBoxReconImgSliceNo->value();

  const auto originX = static_cast<double>(CrntReconImg->GetOrigin()[0]);
  const auto originY = static_cast<double>(CrntReconImg->GetOrigin()[1]);
  const auto originZ = static_cast<double>(CrntReconImg->GetOrigin()[2]);

  const auto spacingX = static_cast<double>(CrntReconImg->GetSpacing()[0]);
  const auto spacingY = static_cast<double>(CrntReconImg->GetSpacing()[1]);
  const auto spacingZ = static_cast<double>(CrntReconImg->GetSpacing()[2]);

  const auto posX = originX + dataX * spacingX;
  const auto posY = originY + dataY * spacingY;
  const auto posZ = originZ + dataZ * spacingZ;

  const auto tmpStr1 = QString("%1").arg(posX, 0, 'f', 2);
  const auto tmpStr2 = QString("%1").arg(posY, 0, 'f', 2);
  const auto tmpStr3 = QString("%1").arg(posZ, 0, 'f', 2);
  this->ui.lineEdit_ForcedProbePosX->setText(tmpStr1);
  this->ui.lineEdit_ForcedProbePosY->setText(tmpStr2);
  this->ui.lineEdit_ForcedProbePosZ->setText(tmpStr3);

  dspReconImg->SetProfileProbePos(dataX, dataY);
  if (this->ui.radioButton_Profile_Hor->isChecked()) {
    dspReconImg->m_bDrawProfileX = true;
    dspReconImg->m_bDrawProfileY = false;
  } else {
    dspReconImg->m_bDrawProfileX = false;
    dspReconImg->m_bDrawProfileY = true;
  }

  // dspReconImg value itself
  const auto ROI_size = this->ui.lineEdit_ROI_size->text().toInt();
  if (ROI_size < 0) {
    return;
  }

  if (ROI_size > 0) {
    dspReconImg->setROI(
        qRound(dataX - ROI_size / 2.0), qRound(dataY - ROI_size / 2.0),
        qRound(dataX + ROI_size / 2.0), qRound(dataY + ROI_size / 2.0));
    dspReconImg->CalcImageInfo_ROI();
    dspReconImg->DrawROIOn(true);

    // m_dspYKImgProj->m_pData[iNumWidth + width*iNumHeight] = (unsigned
    // short)((tmpVal- m_fProjImgValueMin)*m_multiplyFactor);

    const auto strMean =
        QString("%1").arg(dspReconImg->m_fPixelMean_ROI, 0, 'f', 2);
    const auto strSD =
        QString("%1").arg(dspReconImg->m_fPixelSD_ROI, 0, 'f', 2);
    this->ui.lineEdit_ROI_mean->setText(strMean);
    this->ui.lineEdit_ROI_SD->setText(strSD);
  } else {
    dspReconImg->DrawROIOn(false);
  }

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_CalculateROI_Proj() {
  if (this->m_cbctrecon->m_dspYKImgProj == nullptr) {
    return;
  }
  auto &dspProjImg = this->m_cbctrecon->m_dspYKImgProj;

  const auto dspWidth = this->ui.labelImageRaw->width();
  const auto dspHeight = this->ui.labelImageRaw->height();

  const auto dataWidth = dspProjImg->m_iWidth;
  const auto dataHeight = dspProjImg->m_iHeight;
  if (dataWidth * dataHeight == 0) {
    return;
  }

  // int crntIdx = this->ui.spinBoxImgIdx->value();
  // These are displayed data (just index data)
  const auto dataX = qRound(this->ui.labelImageRaw->x /
                            static_cast<double>(dspWidth * dataWidth));
  const auto dataY = qRound(this->ui.labelImageRaw->y /
                            static_cast<double>(dspHeight * dataHeight));
  const auto dataZ = this->ui.spinBoxImgIdx->value();

  const auto originX = 0.0; // 0
  const auto originY = 0.0; // 0
  const auto originZ = 0.0; // 0

  const auto spacingX = 1.0; // 1
  const auto spacingY = 1.0; // 1
  const auto spacingZ = 1.0; // 1

  const auto posX = originX + dataX * spacingX;
  const auto posY = originY + dataY * spacingY;
  const auto posZ = originZ + dataZ * spacingZ;

  this->ui.lineEdit_ForcedProbePosX->setText(QString("%1").arg(posX));
  this->ui.lineEdit_ForcedProbePosY->setText(QString("%1").arg(posY));
  this->ui.lineEdit_ForcedProbePosZ->setText(QString("%1").arg(posZ));

  dspProjImg->SetProfileProbePos(dataX, dataY);

  if (this->ui.radioButton_Profile_Hor->isChecked()) {
    dspProjImg->m_bDrawProfileX = true;
    dspProjImg->m_bDrawProfileY = false;
  } else {
    dspProjImg->m_bDrawProfileX = false;
    dspProjImg->m_bDrawProfileY = true;
  }

  // m_dspYKReconImage value itself
  const auto ROI_size = this->ui.lineEdit_ROI_size->text().toInt();
  if (ROI_size < 0) {
    return;
  }

  if (ROI_size > 0) {
    dspProjImg->setROI(
        qRound(dataX - ROI_size / 2.0), qRound(dataY - ROI_size / 2.0),
        qRound(dataX + ROI_size / 2.0), qRound(dataY + ROI_size / 2.0));
    dspProjImg->CalcImageInfo_ROI();
    dspProjImg->DrawROIOn(true);
    // strMean.sprintf("%5.1f", dspProjImg->m_fPixelMean_ROI);
    const auto strMean = QString("%1").arg(
        dspProjImg->m_fPixelMean_ROI / this->m_cbctrecon->m_multiplyFactor +
            this->m_cbctrecon->m_fProjImgValueMin,
        0, 'f', 2);

    const auto strSD = QString("%1").arg(
        dspProjImg->m_fPixelSD_ROI / this->m_cbctrecon->m_multiplyFactor, 0,
        'f', 2);
    this->ui.lineEdit_ROI_mean->setText(strMean);
    this->ui.lineEdit_ROI_SD->setText(strSD);
  } else {
    dspProjImg->DrawROIOn(false);
  }

  SLT_DrawProjImages();
}

void CbctReconWidget::SLT_GoForcedProbePos() // when forced probe button was
                                             // clicked
{
  const auto fForcedProbePosX = this->ui.lineEdit_ForcedProbePosX->text()
                                    .toDouble(); // data is the reference
  const auto fForcedProbePosY =
      this->ui.lineEdit_ForcedProbePosY->text().toDouble();
  const auto fForcedProbePosZ =
      this->ui.lineEdit_ForcedProbePosZ->text().toDouble();

  // First change the scene acc to Z value
  double originX, originY, originZ;
  double spacingX, spacingY, spacingZ;

  int dataX, dataY;

  if (this->ui.radioButton_graph_proj->isChecked()) {
    if (this->m_cbctrecon->m_spProjImg3DFloat == nullptr) {
      return;
    }
    auto &ProjImg3D = this->m_cbctrecon->m_spProjImg3DFloat;

    originX = ProjImg3D->GetOrigin()[0];
    originY = ProjImg3D->GetOrigin()[1];
    originZ = ProjImg3D->GetOrigin()[2];

    spacingX = ProjImg3D->GetSpacing()[0];
    spacingY = ProjImg3D->GetSpacing()[1];
    spacingZ = ProjImg3D->GetSpacing()[2];

    const auto sliceIdx = qRound((fForcedProbePosZ - originZ) / spacingZ);

    if (sliceIdx < 0 || sliceIdx >= this->m_cbctrecon->m_iImgCnt) {
      return;
    }

    this->ui.spinBoxImgIdx->setValue(sliceIdx); // Draw function is called

    const auto dspWidth = this->ui.labelImageRaw->width();
    const auto dspHeight = this->ui.labelImageRaw->height();

    const auto dataWidth = this->m_cbctrecon->m_dspYKImgProj->m_iWidth;
    const auto dataHeight = this->m_cbctrecon->m_dspYKImgProj->m_iHeight;

    dataX = qRound((fForcedProbePosX - originX) / spacingX);
    dataY = qRound((fForcedProbePosY - originY) / spacingY);

    if (dataX < 0 || dataX >= dataWidth || dataY < 0 || dataY >= dataHeight) {
      return;
    }

    this->ui.labelImageRaw->x =
        qRound(dataX / static_cast<double>(dataWidth) * dspWidth);
    this->ui.labelImageRaw->y =
        qRound(dataY / static_cast<double>(dataHeight) * dspHeight);

    SLT_CalculateROI_Proj();
  } else if (this->ui.radioButton_graph_recon->isChecked()) {
    if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
      return;
    }
    auto &CrntReconImg = this->m_cbctrecon->m_spCrntReconImg;
    originX = CrntReconImg->GetOrigin()[0];
    originY = CrntReconImg->GetOrigin()[1];
    originZ = CrntReconImg->GetOrigin()[2];

    spacingX = CrntReconImg->GetSpacing()[0];
    spacingY = CrntReconImg->GetSpacing()[1];
    spacingZ = CrntReconImg->GetSpacing()[2];

    const auto sliceIdx = qRound((fForcedProbePosZ - originZ) / spacingZ);

    if (sliceIdx < 0 ||
        sliceIdx >=
            static_cast<int>(CrntReconImg->GetBufferedRegion().GetSize()[2])) {
      return;
    }

    this->ui.spinBoxReconImgSliceNo->setValue(
        sliceIdx); // Draw function is called

    const auto dspWidth = this->ui.labelReconImage->width();
    const auto dspHeight = this->ui.labelReconImage->height();

    const auto dataWidth = this->m_cbctrecon->m_dspYKReconImage->m_iWidth;
    const auto dataHeight = this->m_cbctrecon->m_dspYKReconImage->m_iHeight;

    dataX = qRound((fForcedProbePosX - originX) / spacingX);
    dataY = qRound((fForcedProbePosY - originY) / spacingY);

    if (dataX < 0 || dataX >= dataWidth || dataY < 0 || dataY >= dataHeight) {
      return;
    }

    this->ui.labelReconImage->x =
        qRound(dataX / static_cast<double>(dataWidth) * dspWidth);
    this->ui.labelReconImage->y =
        qRound(dataY / static_cast<double>(dataHeight) * dspHeight);

    SLT_CalculateROI_Recon();
  }
}

void CbctReconWidget::SLT_PostApplyFOVDispParam() {
  // this->m_cbctrecon->PostApplyFOVDispParam();
  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_DoPostProcessing() {
  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position

  const auto physPosX = this->ui.lineEdit_PostFOV_X->text().toFloat();
  const auto physPosY = this->ui.lineEdit_PostFOV_Y->text().toFloat();

  const auto physRadius = this->ui.lineEdit_PostFOV_R->text().toFloat();
  const auto physTablePosY = this->ui.lineEdit_PostTablePosY->text().toFloat();

  std::cout << "YKDEBUG " << physPosX << "," << physPosY << "," << physRadius
            << "," << physTablePosY << std::endl;

  this->m_cbctrecon->CropFOV3D(this->m_cbctrecon->m_spCrntReconImg, physPosX,
                               physPosY, physRadius, physTablePosY);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_PostProcCropInv() {
  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position

  const auto physPosX = this->ui.lineEdit_PostFOV_X->text().toDouble();
  const auto physPosY = this->ui.lineEdit_PostFOV_Y->text().toDouble();

  const auto physRadius = this->ui.lineEdit_PostFOV_R->text().toDouble();
  // double physTablePosY = this->ui.lineEdit_PostTablePosY->text().toDouble();

  auto origin = this->m_cbctrecon->m_spCrntReconImg->GetOrigin();
  auto spacing = this->m_cbctrecon->m_spCrntReconImg->GetSpacing();
  // UShortImageType::SizeType size =
  // m_spCrntReconImg->GetBufferedRegion().GetSize();

  // itk::ImageSliceConstIteratorWithIndex<FloatImageType> it (m_spReconImg,
  // m_spReconImg->GetRequestedRegion());
  itk::ImageSliceIteratorWithIndex<UShortImageType> it(
      this->m_cbctrecon->m_spCrntReconImg,
      this->m_cbctrecon->m_spCrntReconImg->GetRequestedRegion());

  // ImageSliceConstIteratorWithIndex<ImageType> it( image,
  // image->GetRequestedRegion() );
  // UShortImageType::SizeType imgSize =
  // m_spCrntReconImg->GetRequestedRegion().GetSize(); //1016x1016 x z

  // int width = imgSize[0];
  // int height = imgSize[1];

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  auto iNumSlice = 0;

  while (!it.IsAtEnd()) {
    auto iPosY = 0;
    while (!it.IsAtEndOfSlice()) {
      auto iPosX = 0;
      while (!it.IsAtEndOfLine()) {
        // Calculate physical position

        const auto crntPhysX = iPosX * static_cast<double>(spacing[0]) +
                               static_cast<double>(origin[0]);
        const auto crntPhysY = iPosY * static_cast<double>(spacing[1]) +
                               static_cast<double>(origin[1]);

        // crop inside of FOV
        if (pow(crntPhysX - physPosX, 2.0) + pow(crntPhysY - physPosY, 2.0) <
            pow(physRadius, 2.0)) {
          it.Set(0);
        }
        // if (crntPhysY >= physTablePosY) //table cropping = same
        //{
        //    it.Set(0);
        //}
        ++it;
        iPosX++;
      }
      it.NextLine();
      iPosY++;
    }
    it.NextSlice();
    iNumSlice++;
  }

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_ExportALL_DCM_and_SHORT_HU_and_calc_WEPL() {

  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  auto dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), this->m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.isEmpty()) {
    return;
  }

  // Get current folder
  QDir crntDir(dirPath);

  QInputDialog inputDlg;

  bool ok;
  auto textInput = QInputDialog::getText(
      this, "Input Dialog", "Set Patient ID and Name", QLineEdit::Normal,
      "PatientID_LastName_FirstName", &ok);

  // QString strEndFix = "YKP";
  QString strPatientID;
  QString strLastName;
  QString strFirstName;

  if (ok && !textInput.isEmpty()) {
    auto strListPtInfo = textInput.split("_");

    if (strListPtInfo.count() >= 3) {
      strPatientID = strListPtInfo.at(0);
      strLastName = strListPtInfo.at(1);
      strFirstName = strListPtInfo.at(2);
    } else if (strListPtInfo.count() == 2) {
      strPatientID = strListPtInfo.at(0);
      strLastName = strListPtInfo.at(1);
    } else if (strListPtInfo.count() == 1) {
      strPatientID = strListPtInfo.at(0);
    } else {
      strPatientID = this->m_cbctrecon->m_strDCMUID;
    }
    // strPatientID = m_strDCMUID + "_" + strEndFix;
  } else {
    strPatientID = this->m_cbctrecon->m_strDCMUID;
  }

  if (strPatientID.isEmpty()) {
    return;
  }

  for (auto i = 0; i < this->m_dlgRegistration->ui.comboBoxImgFixed->count();
       i++) {
    this->m_dlgRegistration->ui.comboBoxImgFixed->setCurrentIndex(i);
    auto strDirName =
        this->m_dlgRegistration->ui.comboBoxImgFixed->currentText();
    const auto tmpResult =
        crntDir.mkdir(strDirName); // what if the directory exists?
    if (!tmpResult) {
      std::cout
          << "DICOM dir seems to exist already. Files will be overwritten."
          << std::endl;
    }

    auto strSavingFolder = dirPath + "/" + strDirName;
    auto strFullName = strLastName + ", " + strFirstName;
    m_dlgRegistration->LoadImgFromComboBox(0, strDirName);

    SaveUSHORTAsSHORT_DICOM(m_dlgRegistration->m_spFixed, strPatientID,
                            strFullName, strSavingFolder);
    auto mhaFileName = strSavingFolder + "/" + strDirName + ".mha";
    ExportReconSHORT_HU(m_dlgRegistration->m_spFixed, mhaFileName);
  }
  SLT_GeneratePOIData();
  const auto angle_end_one = QString("1");
  this->ui.lineEdit_AngEnd->setText(angle_end_one);
  SLT_ExportAngularWEPL_byFile();
}

void CbctReconWidget::SLT_ExportReconSHORT_HU() {
  auto strPath = QFileDialog::getSaveFileName(this, "Save Image", "",
                                              "signed short meta image (*.mha)",
                                              nullptr, nullptr);
  if (strPath.length() <= 1) {
    return;
  }
  ExportReconSHORT_HU(this->m_cbctrecon->m_spCrntReconImg, strPath);
}

void CbctReconWidget::SLT_DoBHC() {

  std::cout << "Beam hardening correction is under progress.." << std::endl;
  this->m_cbctrecon->DoBeamHardeningCorrection(); // only for m_spProjImg3D
  this->m_cbctrecon->SetMaxAndMinValueOfProjectionImage();

  SLT_DrawProjImages();
}

void CbctReconWidget::SLT_DoBowtieCorrection() {
  if (this->m_cbctrecon->m_spProjImg3DFloat == nullptr) {
    return;
  }

  if (this->m_cbctrecon->m_projFormat != HND_FORMAT) {
    std::cout
        << "Bow tie filtering should not be used for His data or Xim data!!"
        << std::endl;
    return;
  }

  const auto strList = this->ui.comboBox_fBTcor->currentText().split(';');

  this->m_cbctrecon->BowtieByFit(this->ui.checkBox_Fullfan->isChecked(),
                                 strList);

  this->m_cbctrecon->SetMaxAndMinValueOfProjectionImage();
  SLT_DrawProjImages();
  std::cout << "Bow-tie correction done." << std::endl;
}

void CbctReconWidget::SLT_ViewRegistration() const
// default showing function
{
  m_dlgRegistration->UpdateListOfComboBox(0); // combo selection
                                              // signalis called
  m_dlgRegistration->UpdateListOfComboBox(1);
  m_dlgRegistration->show();
}

void CbctReconWidget::SLT_ViewHistogram() // default showing function
{
  /*
  m_pDlgHistogram->show();

  if (m_pDlgRegistration->m_spMoving) {
          ForwardProjection(m_pDlgRegistration->m_spMoving, m_spCustomGeometry,
  m_spProjImgCT3D, false); //final moving image
  }

  m_pDlgHistogram->SLT_DrawGraph();
  */
}

// output spProjCT3D => intensity value, not line integral
void CbctReconWidget::ForwardProjection(UShortImageType::Pointer &spVolImg3D,
                                        GeometryType::Pointer &spGeometry,
                                        UShortImageType::Pointer &spProjCT3D,
                                        const bool bSave, const bool use_cuda) {
  if (spVolImg3D == nullptr) {
    std::cout << "ERROR! No 3D-CT file. Load 3D CT file first" << std::endl;
    return;
  }

  if (this->m_cbctrecon->m_iCntSelectedProj < 1 && bSave) {
    std::cout << "Error! No projection image is loaded" << std::endl;
    return;
  }

  if (spGeometry->GetGantryAngles().empty()) {
    std::cout << "No geometry!" << std::endl;
    return;
  }

#ifndef USE_CUDA
  this->m_cbctrecon->CPU_ForwardProjection(spVolImg3D, spGeometry,
                                           spProjCT3D); // final moving image
#else
  if (use_cuda) {
    this->m_cbctrecon->CUDA_ForwardProjection(spVolImg3D, spGeometry,
                                              spProjCT3D); // final moving image
  } else {
    this->m_cbctrecon->CPU_ForwardProjection(spVolImg3D, spGeometry,
                                             spProjCT3D); // final moving image
  }
#endif // !USE_CUDA
  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout
        << " Patient DIR Path: "
        << this->m_cbctrecon->m_strPathPatientDir.toLocal8Bit().constData()
        << std::endl;

    auto manuallySelectedDir = false; // <- just to make sure I don't break
                                      // usecases of the older version.
    if (this->m_cbctrecon->m_strPathPatientDir.isEmpty()) {
      std::cout << "File save error!: No patient DIR name" << std::endl;

      this->m_cbctrecon->m_strPathPatientDir =
          QFileDialog::getExistingDirectory(
              this, tr("Open Directory"), ".",
              QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
      if (this->m_cbctrecon->m_strPathPatientDir.length() <= 1) {
        return;
      }
      manuallySelectedDir = true;
    }

    // Get current folder
    const auto subdir_images("IMAGES");
    const auto strCrntDir = this->m_cbctrecon->m_strPathPatientDir + "/" +
                            subdir_images; // current Proj folder

    // Make a sub directory
    QDir crntDir(strCrntDir);

    if (!crntDir.exists()) {
      if (manuallySelectedDir) {
        QDir current_dir(this->m_cbctrecon->m_strPathPatientDir);
        const auto success = current_dir.mkdir(subdir_images);
        if (!success) {
          std::cerr << "Could not create subfolder IMAGES in given directory"
                    << std::endl;
          return;
        }
      } else {
        std::cout << "File save error: The specified folder does not exist."
                  << std::endl;
        return;
      }
    }

    const auto fwdDirName = "fwd_" + this->m_cbctrecon->m_strDCMUID;

    const auto tmpResult =
        crntDir.mkdir(fwdDirName); // what if the directory exists?

    if (!tmpResult) {
      std::cout << "FwdProj directory seems to exist already. Files will be "
                   "overwritten."
                << std::endl;
    }

    auto strSavingFolder = strCrntDir + "/" + fwdDirName;
    this->m_cbctrecon->SaveProjImageAsHIS(
        spProjCT3D, this->m_cbctrecon->m_arrYKBufProj, strSavingFolder,
        this->m_cbctrecon->m_fResampleF);
  }
}

void CbctReconWidget::SLT_DoScatterCorrection_APRIORI() {

  const auto bExportProj_Fwd = this->ui.checkBox_ExportFwd->isChecked();
  const auto bExportProj_Scat = this->ui.checkBox_ExportScat->isChecked();
  const auto bExportProj_Cor = this->ui.checkBox_ExportCor->isChecked();

  // ForwardProjection(m_spRefCTImg, m_spCustomGeometry, m_spProjImgCT3D,
  // false); //final moving image
  if (m_dlgRegistration->m_spMoving != nullptr) {
    ForwardProjection(
        m_dlgRegistration->m_spMoving, this->m_cbctrecon->m_spCustomGeometry,
        this->m_cbctrecon->m_spProjImgCT3D, bExportProj_Fwd,
        this->ui.radioButton_UseCUDA->isChecked()); // final moving image
  } else if (this->m_cbctrecon->m_spRefCTImg != nullptr) {
    std::cout << "No Moving image in Registration is found. Ref CT image will "
                 "be used instead"
              << std::endl;
    ForwardProjection(
        this->m_cbctrecon->m_spRefCTImg, this->m_cbctrecon->m_spCustomGeometry,
        this->m_cbctrecon->m_spProjImgCT3D, bExportProj_Fwd,
        this->ui.radioButton_UseCUDA->isChecked()); // final moving image
  } else {
    std::cout << "Error!: No ref image for forward projection is found."
              << std::endl;
    return;
  }

  // YKTEMP
  std::cout
      << "ProjImgCT Size = "
      << this->m_cbctrecon->m_spProjImgCT3D->GetBufferedRegion().GetSize()[0]
      << ", "
      << this->m_cbctrecon->m_spProjImgCT3D->GetBufferedRegion().GetSize()[1]
      << ", "
      << this->m_cbctrecon->m_spProjImgCT3D->GetBufferedRegion().GetSize()[2]
      << std::endl;
  std::cout << "ProjImgCT origin = "
            << this->m_cbctrecon->m_spProjImgCT3D->GetOrigin()[0] << ", "
            << this->m_cbctrecon->m_spProjImgCT3D->GetOrigin()[1] << ", "
            << this->m_cbctrecon->m_spProjImgCT3D->GetOrigin()[2] << std::endl;
  std::cout << "ProjImgCT spacing = "
            << this->m_cbctrecon->m_spProjImgCT3D->GetSpacing()[0] << ", "
            << this->m_cbctrecon->m_spProjImgCT3D->GetSpacing()[1] << ", "
            << this->m_cbctrecon->m_spProjImgCT3D->GetSpacing()[2] << std::endl;

  // double scaResam = this->ui.lineEdit_scaResam->text().toDouble();
  const auto scaMedian = this->ui.lineEdit_scaMedian->text().toDouble();
  const auto scaGaussian = this->ui.lineEdit_scaGaussian->text().toDouble();

  std::cout << "Generating scatter map is ongoing..." << std::endl;

  this->m_cbctrecon->GenScatterMap_PriorCT(
      this->m_cbctrecon->m_spProjImgRaw3D, this->m_cbctrecon->m_spProjImgCT3D,
      this->m_cbctrecon->m_spProjImgScat3D, scaMedian, scaGaussian,
      this->m_cbctrecon->m_iFixedOffset_ScatterMap,
      bExportProj_Scat); // void GenScatterMap2D_PriorCT()

  std::cout << "To account for the mAs values, the intensity scale factor of "
            << GetRawIntensityScaleFactor(this->m_cbctrecon->m_strRef_mAs,
                                          this->m_cbctrecon->m_strCur_mAs)
            << "was multiplied during scatter correction to avoid negative "
               "scatter"
            << std::endl;

  this->ui.lineEdit_CurmAs->setText(this->m_cbctrecon->m_strCur_mAs);
  this->ui.lineEdit_RefmAs->setText(this->m_cbctrecon->m_strRef_mAs);

  this->m_cbctrecon->m_spProjImgCT3D->Initialize(); // memory saving

  std::cout << "Scatter correction is in progress..." << std::endl;

  const auto postScatMedianSize =
      this->ui.lineEdit_scaPostMedian->text().toInt();
  this->m_cbctrecon->ScatterCorr_PrioriCT(
      this->m_cbctrecon->m_spProjImgRaw3D, this->m_cbctrecon->m_spProjImgScat3D,
      this->m_cbctrecon->m_spProjImgCorr3D,
      this->m_cbctrecon->m_iFixedOffset_ScatterMap, postScatMedianSize,
      bExportProj_Cor);
  this->m_cbctrecon->m_spProjImgScat3D->Initialize(); // memory saving

  std::cout << "AfterCorrectionMacro is ongoing..." << std::endl;

  // Update UI
  this->ui.pushButton_DoRecon->setEnabled(true);
  this->ui.spinBoxImgIdx->setMinimum(0);
  const auto iSizeZ =
      this->m_cbctrecon->m_spProjImg3DFloat->GetRequestedRegion().GetSize()[2];
  this->ui.spinBoxImgIdx->setMaximum(iSizeZ - 1);
  this->ui.spinBoxImgIdx->setValue(0);
  this->m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();
  SLT_DrawProjImages(); // Update Table is called

  auto fdk_options = getFDKoptions();

  this->m_cbctrecon->AfterScatCorrectionMacro(
      this->ui.radioButton_UseCUDA->isChecked(),
      this->ui.radioButton_UseOpenCL->isChecked(),
      this->ui.checkBox_ExportVolDICOM->isChecked(), fdk_options);

  // Skin removal (using CT contour w/ big margin)
  std::cout
      << "Post  FDK reconstruction is done. Moving on to post skin removal"
      << std::endl;

  m_cbctregistration->PostSkinRemovingCBCT(this->m_cbctrecon->m_spRawReconImg);
  m_cbctregistration->PostSkinRemovingCBCT(
      this->m_cbctrecon->m_spScatCorrReconImg);

  // 20151208 Removal of high intensity skin mask
  // Main issue: raw CBCT projection includes mask, deformed CT doesn't include
  // mask. In case of weight loss, mask signal is independent from skin contour,
  // but deformed CT cannot have that signal.  Therefore, after the subtraction
  // (CBCTcor projections), there is always a big peak. DIR quality doesn't
  // matter because it cannot 'create' mask signal anyway.  Assumption: near the
  // skin contour, this kind of discrepancy is not expected.
  // m_pDlgRegistration->ThermoMaskRemovingCBCT(m_spRawReconImg,
  // m_spScatCorrReconImg, threshold_HU);

  m_dlgRegistration->UpdateListOfComboBox(0); // combo selection signalis
                                              // called
  m_dlgRegistration->UpdateListOfComboBox(1);
  m_dlgRegistration->SelectComboExternal(
      0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  m_dlgRegistration->SelectComboExternal(1, REGISTER_COR_CBCT);

  m_dlgRegistration
      ->SLT_DoLowerMaskIntensity(); // it will check the check button.

  SLT_DrawProjImages();

  std::cout << "Updating ReconImage..";
  auto updated_text = QString("Scatter corrected CBCT");
  UpdateReconImage(this->m_cbctrecon->m_spScatCorrReconImg,
                   updated_text); // main GUI update

  std::cout << "FINISHED!Scatter correction: CBCT DICOM files are saved"
            << std::endl;
}

// called whenver recon 3D image for display changes.
void CbctReconWidget::UpdateReconImage(UShortImageType::Pointer &spNewImg,
                                       QString &fileName) {
  this->m_cbctrecon->m_spCrntReconImg = spNewImg;

  const auto origin_new = this->m_cbctrecon->m_spCrntReconImg->GetOrigin();
  const auto spacing_new = this->m_cbctrecon->m_spCrntReconImg->GetSpacing();
  const auto size_new =
      this->m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize();

  std::cout << "New Origin" << origin_new << std::endl;
  std::cout << "New spacing" << spacing_new << std::endl;
  std::cout << "New size" << size_new << std::endl;

  this->ui.lineEdit_Cur3DFileName->setText(fileName);

  auto size =
      this->m_cbctrecon->m_spCrntReconImg->GetRequestedRegion().GetSize();

  this->m_cbctrecon->m_dspYKReconImage->CreateImage(size[0], size[1], 0);

  disconnect(this->ui.spinBoxReconImgSliceNo, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawReconImage()));

  this->ui.spinBoxReconImgSliceNo->setMinimum(0);
  this->ui.spinBoxReconImgSliceNo->setMaximum(size[2] - 1);

  const auto initVal = qRound((size[2] - 1) / 2.0);
  // SLT_DrawReconImage(); //Update Table, Update Graph

  // m_dspYKReconImage->CreateImage(size_trans[0], size_trans[1],0);
  SLT_InitializeGraphLim();

  this->ui.spinBoxReconImgSliceNo->setValue(initVal);
  this->ui.radioButton_graph_recon->setChecked(true);

  connect(this->ui.spinBoxReconImgSliceNo, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawReconImage()));

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_TempAudit() const {
  if (this->m_cbctrecon->m_spRawReconImg != nullptr) {
    std::cout << "m_spRawReconImg " << this->m_cbctrecon->m_spRawReconImg
              << std::endl;
  }

  if (this->m_cbctrecon->m_spRefCTImg != nullptr) {
    std::cout << "m_spRefCTImg " << this->m_cbctrecon->m_spRefCTImg
              << std::endl;
  }

  if (this->m_cbctrecon->m_spCrntReconImg != nullptr) {
    std::cout << "m_spCrntReconImg " << this->m_cbctrecon->m_spCrntReconImg
              << std::endl;
  }
}

void CbctReconWidget::SLT_LoadPlanCT_USHORT() {
  // typedef itk::ImageFileWriter<FloatImageType> WriterType;
  using ReaderType = itk::ImageFileReader<UShortImageType>;
  auto reader = ReaderType::New();

  auto fileName = QFileDialog::getOpenFileName(
      this, "Open Image", this->m_cbctrecon->m_strPathDirDefault,
      "Projection file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  reader->SetFileName(fileName.toLocal8Bit().constData());
  reader->Update();

  this->m_cbctrecon->m_spRefCTImg = reader->GetOutput();
  auto ref_ct = QString("RefCT");
  UpdateReconImage(this->m_cbctrecon->m_spRefCTImg, ref_ct);

  this->m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT,
                                            REGISTER_MANUAL_RIGID);
}

void CbctReconWidget::SLT_CalcAndSaveAngularWEPL() // single point
{
  std::vector<WEPLData> vOutputWEPL;

  const auto fAngleGap = this->ui.lineEdit_WEPL_AngRes->text().toDouble();
  const auto fAngleStart = this->ui.lineEdit_AngStart->text().toDouble();
  const auto fAngleEnd = this->ui.lineEdit_AngEnd->text().toDouble();

  const auto cur_poi =
      VEC3D{this->ui.lineEdit_ForcedProbePosX->text().toDouble(), // in mm
            this->ui.lineEdit_ForcedProbePosY->text().toDouble(),
            this->ui.lineEdit_ForcedProbePosZ->text().toDouble()};

  this->m_cbctrecon->GetAngularWEPL_SinglePoint(
      this->m_cbctrecon->m_spCrntReconImg, fAngleGap, fAngleStart, fAngleEnd,
      cur_poi, 0, vOutputWEPL, true);
  std::cout << "Computed WEPL points: " << vOutputWEPL.size() << std::endl;

  // export arrWEPL
  auto filePath = QFileDialog::getSaveFileName(
      this, "Save data", "", "txt image file (*.txt)", nullptr,
      nullptr); // Filename don't need to exist

  if (filePath.length() < 1) {
    return;
  }

  std::ofstream fout;
  fout.open(filePath.toLocal8Bit().constData());

  fout << "Angle"
       << "	"
       << "WEPL(mm)" << std::endl;

  for (auto &i : vOutputWEPL) {
    fout << i.fGanAngle << "	" << i.fWEPL << std::endl;
  }

  fout.close();
  std::cout << "Saving angular WEPL is completed" << std::endl;
}

// This function deals with the current projection image converted by
// ElektaProjReader (intensity to lineintegral is already done. now the type is
// float, 3D)
void CbctReconWidget::SLT_DoScatterCorrectionUniform() {
  if (this->m_cbctrecon->m_spProjImg3DFloat == nullptr) {
    return;
  }

  UShortImageType::Pointer spIntensityRaw;
  this->m_cbctrecon->ConvertLineInt2Intensity(
      this->m_cbctrecon->m_spProjImg3DFloat, spIntensityRaw, 65535);

  using ScatterFilterType =
      rtk::BoellaardScatterCorrectionImageFilter<UShortImageType,
                                                 UShortImageType>;

  auto spScatFilter = ScatterFilterType::New();

  const auto airThre = this->ui.lineEdit_uniAirThre->text().toDouble();
  const auto scat2PrimRatio = this->ui.lineEdit_uniSPR->text().toDouble();
  const auto nonNagativity = this->ui.lineEdit_uniNegativity->text().toDouble();

  std::cout << "Boallaard uniform scatter correction is being applied"
            << std::endl;
  std::cout << "Air threshold: " << airThre << std::endl;
  std::cout << "Scatter to Primary ratio(will be automatically adjusted by the "
               "minimum intensity): "
            << scat2PrimRatio << std::endl;
  std::cout << "NonNegativityConstraintThreshold: " << nonNagativity
            << std::endl;

  spScatFilter->SetInput(spIntensityRaw);
  spScatFilter->SetAirThreshold(airThre);
  spScatFilter->SetScatterToPrimaryRatio(scat2PrimRatio);
  spScatFilter->SetNonNegativityConstraintThreshold(nonNagativity);

  spScatFilter->Update(); // errir occured

  UShortImageType::Pointer spIntensityUniformCorr = spScatFilter->GetOutput();

  this->m_cbctrecon->ConvertIntensity2LineInt(
      spIntensityUniformCorr, this->m_cbctrecon->m_spProjImg3DFloat, 65535);

  // ConvertLineInt2Intensity(m_spProjImg3DFloat, m_spProjImgRaw3D, 65535);

  this->m_cbctrecon->m_spProjImgRaw3D = spIntensityUniformCorr;

  this->ui.spinBoxImgIdx->setValue(0);
  this->m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();

  this->SLT_DrawProjImages(); // Update Table is called

  std::cout << "FINISHED!: Uniform Scatter correction for raw projection "
               "images (Boallaard method) is completed. Proceed to "
               "reconstruction"
            << std::endl;
}

void CbctReconWidget::SLT_FileExportShortDICOM_CurrentImg() {

  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  auto dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), this->m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.isEmpty()) {
    return;
  }

  // Get current folder
  QDir crntDir(dirPath);

  QInputDialog inputDlg;

  bool ok;
  auto textInput = QInputDialog::getText(
      this, "Input Dialog", "Set Patient ID and Name", QLineEdit::Normal,
      "PatientID_LastName_FirstName", &ok);

  // QString strEndFix = "YKP";
  QString strPatientID;
  QString strLastName;
  QString strFirstName;

  if (ok && !textInput.isEmpty()) {
    auto strListPtInfo = textInput.split("_");

    if (strListPtInfo.count() >= 3) {
      strPatientID = strListPtInfo.at(0);
      strLastName = strListPtInfo.at(1);
      strFirstName = strListPtInfo.at(2);
    } else if (strListPtInfo.count() == 2) {
      strPatientID = strListPtInfo.at(0);
      strLastName = strListPtInfo.at(1);
    } else if (strListPtInfo.count() == 1) {
      strPatientID = strListPtInfo.at(0);
    } else {
      strPatientID = this->m_cbctrecon->m_strDCMUID;
    }
    // strPatientID = m_strDCMUID + "_" + strEndFix;
  } else {
    strPatientID = this->m_cbctrecon->m_strDCMUID;
  }

  if (strPatientID.isEmpty()) {
    return;
  }

  const auto strDirName = strPatientID + "_DCM";
  const auto tmpResult =
      crntDir.mkdir(strDirName); // what if the directory exists?
  if (!tmpResult) {
    std::cout << "DICOM dir seems to exist already. Files will be overwritten."
              << std::endl;
  }

  auto strSavingFolder = dirPath + "/" + strDirName;
  auto strFullName = strLastName + ", " + strFirstName;
  SaveUSHORTAsSHORT_DICOM(this->m_cbctrecon->m_spCrntReconImg, strPatientID,
                          strFullName, strSavingFolder);
}

void CbctReconWidget::SLT_AddConstHUToCurImg() {
  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }
  const auto addingVal = this->ui.lineEdit_AddConstHU->text().toInt();
  AddConstHU(this->m_cbctrecon->m_spCrntReconImg, addingVal);
  auto updated_text = QString("Added%1").arg(addingVal);
  UpdateReconImage(this->m_cbctrecon->m_spCrntReconImg, updated_text);
}

void CbctReconWidget::SLT_SetCBCTSkinRSPath() {
  auto strPath = QFileDialog::getOpenFileName(
      this, "Open RS file", this->m_cbctrecon->m_strPathDirDefault,
      "DICOM RS (*.dcm)", nullptr, nullptr);

  if (strPath.length() <= 1) {
    return;
  }

  this->ui.lineEdit_PathCBCTSkinPath->setText(strPath);
}

void CbctReconWidget::SLT_CropSkinUsingThreshold() {

  auto update_text = QString("Thresh-based skin cropped image");
  const auto thresh = this->ui.lineEdit_Threshold->text().toInt();
  const auto erode = this->ui.lineEdit_ErodeRadius->text().toInt();
  const auto dilate = this->ui.lineEdit_DilateRadius->text().toInt();

  const auto imgType =
      this->m_cbctrecon->CropSkinUsingThreshold(thresh, erode, dilate);

  switch (imgType) {
  case 1:
    UpdateReconImage(this->m_cbctrecon->m_spRawReconImg, update_text);
    break;
  case 2:
    UpdateReconImage(this->m_cbctrecon->m_spRefCTImg, update_text);
    break;
  case 3:
    UpdateReconImage(this->m_cbctrecon->m_spScatCorrReconImg, update_text);
    break;
  default:
    std::cerr << "WTF!?" << std::endl;
    break;
  }
}

void CbctReconWidget::SLT_CropSkinUsingRS() {
  auto strPathRS = this->ui.lineEdit_PathCBCTSkinPath->text();
  if (strPathRS.length() < 1) {
    return;
  }

  const auto croppingMargin = this->ui.lineEdit_SkinMargin->text().toDouble();
  auto update_text = QString("RS-based skin cropped image");
  if (this->m_cbctrecon->m_spCrntReconImg ==
      this->m_cbctrecon->m_spRawReconImg) {
    m_cbctregistration->CropSkinUsingRS(this->m_cbctrecon->m_spRawReconImg,
                                        strPathRS, croppingMargin);
    UpdateReconImage(this->m_cbctrecon->m_spRawReconImg, update_text);
  } else if (this->m_cbctrecon->m_spCrntReconImg ==
             this->m_cbctrecon->m_spRefCTImg) {
    m_cbctregistration->CropSkinUsingRS(this->m_cbctrecon->m_spRefCTImg,
                                        strPathRS, croppingMargin);
    UpdateReconImage(this->m_cbctrecon->m_spRefCTImg, update_text);
  } else if (this->m_cbctrecon->m_spCrntReconImg ==
             this->m_cbctrecon->m_spScatCorrReconImg) {
    m_cbctregistration->CropSkinUsingRS(this->m_cbctrecon->m_spScatCorrReconImg,
                                        strPathRS, croppingMargin);
    UpdateReconImage(this->m_cbctrecon->m_spScatCorrReconImg, update_text);
  }
}

void CbctReconWidget::SLT_ExportAngularWEPL_byFile() {
  // export arrWEPL
  auto filePath = QFileDialog::getSaveFileName(
      this, "Save data", this->m_cbctrecon->m_strPathDirDefault,
      "txt image file (*.txt)", nullptr,
      nullptr); // Filename don't need to exist

  if (filePath.length() < 1) {
    return;
  }

  const auto fAngleStart = this->ui.lineEdit_AngStart->text().toDouble();
  const auto fAngleEnd = this->ui.lineEdit_AngEnd->text().toDouble();
  const auto fAngleGap = this->ui.lineEdit_WEPL_AngRes->text().toDouble();

  this->m_cbctrecon->ExportAngularWEPL_byFile(filePath, fAngleStart, fAngleEnd,
                                              fAngleGap);
}

void CbctReconWidget::SLT_GeneratePOIData() const
// it fills m_vPOI_DCM
{
  this->m_cbctrecon->GeneratePOIData(
      this->ui.checkBox_AP->isChecked(),
      this->ui.lineEdit_PostTablePosY->text().toDouble());
}

void CbctReconWidget::SLT_LoadPOIData() // it fills m_vPOI_DCM
{
  if (!this->m_cbctrecon->m_vPOI_DCM.empty()) {
    this->m_cbctrecon->m_vPOI_DCM.clear();
  }

  auto filePath = QFileDialog::getOpenFileName(
      this, "POI data file", this->m_cbctrecon->m_strPathDirDefault,
      "POI data file (*.txt)", nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  std::ifstream fin;
  fin.open(filePath.toLocal8Bit().constData(), std::ios::in);
  if (fin.fail()) {
    return;
  }

  char str[MAX_LINE_LENGTH];
  // File format:
  // no header

  // x1	y1	z1, in mm
  // x2	y2	z2,
  //..

  while (!fin.eof()) {
    VEC3D fPOI{};
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString strLine(&str[0]);

    // QStringList strList = strLine.split(' ');//tab
    auto strList = strLine.split('\t'); // tab
    // third one is the organ name
    if (strList.length() < 3) {
      std::cout << "abnormal file expression." << std::endl;
      break;
    }
    fPOI.x = strList.at(0).toDouble();
    fPOI.y = strList.at(1).toDouble();
    fPOI.z = strList.at(2).toDouble();

    this->m_cbctrecon->m_vPOI_DCM.push_back(fPOI);
  }
  for (auto i = 0; i < static_cast<int>(this->m_cbctrecon->m_vPOI_DCM.size());
       i++) {
    std::cout << "Data " << i << "	"
              << this->m_cbctrecon->m_vPOI_DCM.at(i).x << ", "
              << this->m_cbctrecon->m_vPOI_DCM.at(i).y << ", "
              << this->m_cbctrecon->m_vPOI_DCM.at(i).z << std::endl;
  }
  std::cout << "POI data has been loaded. "
            << this->m_cbctrecon->m_vPOI_DCM.size() << " data points are read"
            << std::endl;
  fin.close();
}

void CbctReconWidget::SLT_StartSyncFromSharedMem() {
  // int msInterval = this->ui.lineEditTimerInterval->text().toInt();

  // if (msInterval > 0 && msInterval < 5000)
  //{
  // if (m_arrYKImage != NULL)
  // {
  // delete[] m_arrYKImage;
  // m_arrYKImage = NULL;
  // m_iImgCnt = 0;
  // }

  // m_iImgCnt = 1;
  // m_arrYKImage = new YK16GrayImage [m_iImgCnt];
  // m_arrYKImage[0].CreateImage(1024,1024,0);
  //
  // m_busyTimer =false;
  // m_Timer->start(msInterval);
  //}
}

void CbctReconWidget::SLT_StopSyncFromSharedMem() {
  // m_Timer->stop();
#ifndef _WIN32
  std::cerr << "Function not implemented for non-windows systems!!"
            << std::endl;
  return;
#else
  // HANDLE hSemaphore = OpenSemaphore(SYNCHRONIZE ,FALSE, "YKSemaphore");
  // Option SYNCHRONIZE doesn't work! you cannot release Semaphore due to the
  // access is denied (GetLastError 5)
  const auto hSemaphore =
      OpenSemaphore(SEMAPHORE_ALL_ACCESS, FALSE, "YKSemaphore");
  // increase counter
  //  LONG prev_counter;
  // ReleaseSemaphore(hSemaphore, 1, &prev_counter);
  // decrease counter
  // ReleaseSemaphore(hSemaphore, 1, &prev_counter);

  std::ofstream fout;
  fout.open("E:\\SemphoreLogC++.txt");

  const auto max_cnt = 200;

  auto cnt = 0;
  while (cnt < max_cnt) {
    cnt++;

    const auto dwWaitResult = WaitForSingleObject(hSemaphore, 5);

    switch (dwWaitResult) {
      // The semaphore object was signaled.
    case WAIT_OBJECT_0: // 0 -->This got the semaphore token
      // PERFORM TASK
      fout << "Thread " << GetCurrentThreadId() << ": wait succeeded "
           << std::endl;
      Sleep(10);

      if (ReleaseSemaphore(hSemaphore,    // handle to semaphore
                           1,             // increase count by one
                           nullptr) == 0) // not interested in previous count
      {
        fout << "ReleaseSemaphore error: " << GetLastError()
             << "\n"; // Access is denied.
      }
      break;

    case WAIT_TIMEOUT: // no need of release //258
      // Sleep(50);
      fout << "Thread " << GetCurrentThreadId() << ": wait timed out "
           << std::endl;
      break;
    default:
      break;
    }
    //	ReleaseSemaphore(hSemaphore, 1, &prev_counter);
  }

  fout.close();
#endif
}

void CbctReconWidget::SLT_TimerEvent() {

#ifndef _WIN32
  std::cerr << "Function not implemented for non-windows systems!!"
            << std::endl;
  return;
#else
  if (m_busyTimer) {
    return;
  }

  if (this->m_cbctrecon->m_arrYKImage.empty() ||
      this->m_cbctrecon->m_iImgCnt != 1) {
    return;
  }

  m_busyTimer = true;

  // Look into the shared mem
  TCHAR szName[] = TEXT("YKSharedMemory");
  void *const handle = OpenFileMapping(FILE_MAP_READ, FALSE, szName);

  if (handle == nullptr) {
    std::cout << "Cannot open Mapped file" << std::endl;
    SLT_StopSyncFromSharedMem();
    return;
  }

  const auto size = 1024 * 1024 * 2;
  const auto pix_size = static_cast<int>(size / 2.0);
  const auto char_buf = static_cast<unsigned char *>(
      MapViewOfFile(handle, FILE_MAP_READ, 0, 0, size));

  if (char_buf == nullptr) {
    std::cout << "Shared memory was not read. Timer will be stopped"
              << std::endl;
    SLT_StopSyncFromSharedMem();
    CloseHandle(handle);
    delete[] char_buf;
    return;
  }

  // byte array to unsigned short
  // assuming little endian

  // unsigned short* imgBuf = new unsigned short [pix_size];

  for (auto i = 0; i < pix_size; i++) {
    const auto idxA = i * 2 + 1;
    const auto idxB = i * 2;
    // 0: 1,0  1: 3,2 ...
    this->m_cbctrecon->m_arrYKImage.at(0).m_pData[i] =
        (char_buf[idxA] << 8) | char_buf[idxB]; // little endian
  }

  this->ui.spinBoxImgIdx->setValue(0);
  SLT_DrawRawImages();

  CloseHandle(handle);

  m_busyTimer = false;
#endif
}

void CbctReconWidget::SLTM_ViewExternalCommand() const {
  this->m_dlgExternalCommand->show();
}

void CbctReconWidget::SLTM_LoadDICOMdir() {
  auto dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), this->m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.length() <= 1) {
    return;
  }

  if (this->m_cbctrecon->ReadDicomDir(dirPath)) {

    m_dlgRegistration->UpdateVOICombobox(PLAN_CT);
    auto update_text = QString("DICOM reference image");
    UpdateReconImage(this->m_cbctrecon->m_spRefCTImg, update_text);

    this->m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT,
                                              REGISTER_MANUAL_RIGID);
  }
}

void CbctReconWidget::SLTM_LoadRTKoutput() {
  auto filePath = QFileDialog::getOpenFileName(
      this, "Open Image", this->m_cbctrecon->m_strPathDirDefault,
      "rtk output float image (*.mha)", nullptr, nullptr);
  this->m_cbctrecon->LoadExternalFloatImage(filePath, true);
  QFileInfo outFileInfo(filePath);
  auto strCrntFileName = outFileInfo.fileName();
  UpdateReconImage(this->m_cbctrecon->m_spRawReconImg, strCrntFileName);
}

// Only can be used for m_spRawRecon // NOT USED AT ALL?
void CbctReconWidget::FileExportByGUI() const
// USHORT
{
  auto outputFilePath = this->ui.lineEdit_OutputFilePath->text();
  QFileInfo outFileInfo(outputFilePath);
  auto outFileDir = outFileInfo.absoluteDir();

  // bool b = outFileDir.exists();
  // QString tmpPath = outFileDir.absolutePath();

  if (outputFilePath.length() < 2 || !outFileDir.exists()) {
    std::cout << "No available output path. Should be exported later"
              << std::endl;
  } else {
    using WriterType = itk::ImageFileWriter<UShortImageType>;
    auto writer = WriterType::New();
    writer->SetFileName(outputFilePath.toLocal8Bit().constData());
    writer->SetUseCompression(true); // not exist in original code (rtkfdk)
    writer->SetInput(this->m_cbctrecon->m_spRawReconImg);

    std::cout << "Writing the image to: "
              << outputFilePath.toLocal8Bit().constData() << std::endl;

    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());

    std::cout << std::endl;
    std::cout << "The output image was successfully saved" << std::endl;
  }
}

void CbctReconWidget::SLT_OutPathEdited() const {
  if (!this->ui.lineEdit_OutputFilePath->text().isEmpty()) {
    this->ui.lineEdit_outImgDim_LR->setEnabled(true);
    this->ui.lineEdit_outImgDim_AP->setEnabled(true);
    this->ui.lineEdit_outImgDim_SI->setEnabled(true);
    this->ui.lineEdit_outImgSp_LR->setEnabled(true);
    this->ui.lineEdit_outImgSp_AP->setEnabled(true);
    this->ui.lineEdit_outImgSp_SI->setEnabled(true);
  } else {
    this->ui.lineEdit_outImgDim_LR->setEnabled(false);
    this->ui.lineEdit_outImgDim_AP->setEnabled(false);
    this->ui.lineEdit_outImgDim_SI->setEnabled(false);
    this->ui.lineEdit_outImgSp_LR->setEnabled(false);
    this->ui.lineEdit_outImgSp_AP->setEnabled(false);
    this->ui.lineEdit_outImgSp_SI->setEnabled(false);
  }
}

void CbctReconWidget::SLT_MedianFilterDoNow() {
  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  /*QString strCrntFileName;
  QFileInfo outFileInfo(strPath);
  strCrntFileName = outFileInfo.fileName();*/

  UShortImageType::SizeType indexRadius{};
  indexRadius[0] =
      this->ui.lineEdit_PostMedSizeX->text().toInt(); // radius along x
  indexRadius[1] =
      this->ui.lineEdit_PostMedSizeY->text().toInt(); // radius along y
  indexRadius[2] =
      this->ui.lineEdit_PostMedSizeZ->text().toInt(); // radius along y

  if (indexRadius[0] != 0 || indexRadius[1] != 0 || indexRadius[2] != 0) {
    this->m_cbctrecon->MedianFilterByGUI(indexRadius);

    auto prevFileName = this->ui.lineEdit_Cur3DFileName->text();
    UpdateReconImage(this->m_cbctrecon->m_spCrntReconImg,
                     prevFileName.append("_med"));
  } else {
    std::cout << "Not valid median window" << std::endl;
  }
}

void CbctReconWidget::SLT_Export2DDose_TIF() // 2D dose from current displayed
                                             // image of reconstruction
{
  if (this->m_cbctrecon->m_dspYKReconImage == nullptr) {
    return;
  }

  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  auto strPath = QFileDialog::getSaveFileName(this, "Save Image", "",
                                              "signed short meta image (*.tif)",
                                              nullptr, nullptr);
  if (strPath.length() <= 1) {
    return;
  }

  const auto originLeft =
      static_cast<double>(this->m_cbctrecon->m_spCrntReconImg->GetOrigin()[0]);
  const auto originTop = static_cast<double>(
      this->m_cbctrecon->m_spCrntReconImg->GetOrigin()[1]); // not sure...

  const auto spacingX =
      static_cast<double>(this->m_cbctrecon->m_spCrntReconImg->GetSpacing()[0]);
  const auto spacingY = static_cast<double>(
      this->m_cbctrecon->m_spCrntReconImg->GetSpacing()[1]); // not sure...

  if (!SaveDoseGrayImage(strPath.toLocal8Bit().constData(),
                         this->m_cbctrecon->m_dspYKReconImage->m_iWidth,
                         this->m_cbctrecon->m_dspYKReconImage->m_iHeight,
                         spacingX, spacingY, originLeft, originTop,
                         this->m_cbctrecon->m_dspYKReconImage->m_pData)) {
    std::cout << "Failed in save gray dose file" << std::endl;
  } else {
    std::cout << "image exported successfully." << std::endl;
  }
}
void CbctReconWidget::SLTM_Export2DDoseMapAsMHA() {

  auto strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "itk compatible meta image (*.mha)", nullptr,
      nullptr);

  this->m_cbctrecon->Export2DDoseMapAsMHA(strPath);
}

void CbctReconWidget::SLTM_ExportProjGeometryTXT() {

  auto strPath = QFileDialog::getSaveFileName(this, "Save text file", "",
                                              "text (*.txt)", nullptr, nullptr);

  this->m_cbctrecon->ExportProjGeometryTXT(strPath);
}

void CbctReconWidget::SLTM_ForwardProjection() {
  if (this->m_cbctrecon->m_spRawReconImg == nullptr) {
    return;
  }

  auto crntGeometry = GeometryType::New();

  if (this->m_cbctrecon->m_spCustomGeometry == nullptr) {
    std::cout << "No geometry is ready. moving on to 360 projection"
              << std::endl;

    const auto curSID = 1000.0;
    const auto curSDD = 1536.0;
    const auto curProjOffsetX = 0.0;
    const auto curProjOffsetY = 0.0;
    const auto curOutOfPlaneAngles = 0.0;
    const auto curInPlaneAngles = 0.0;
    const auto curSrcOffsetX = 0.0;
    const auto curSrcOffsetY = 0.0;

    // double startAngle = 180.0; //kV = 270.0, CW
    const auto startAngle = 270; // kV = 360.0, CW
    // int NumOfProj = 360;
    const auto NumOfProj = 1;

    for (auto i = 0; i < NumOfProj; i++) {
      auto cur_mv_gantry_angle = startAngle + i;
      if (cur_mv_gantry_angle > 360.0) {
        cur_mv_gantry_angle = cur_mv_gantry_angle - 360.0;
      }
      // AddProjection: current CBCT software version only requires MV gantry
      // angle!!!
      crntGeometry->AddProjection(
          curSID, curSDD, cur_mv_gantry_angle, curProjOffsetX,
          curProjOffsetY,                        // Flexmap
          curOutOfPlaneAngles, curInPlaneAngles, // In elekta, these are 0
          curSrcOffsetX, curSrcOffsetY);         // In elekta, these are 0
    }

    ForwardProjection(this->m_cbctrecon->m_spRawReconImg, crntGeometry,
                      this->m_cbctrecon->m_spProjImgRaw3D, false,
                      this->ui.radioButton_UseCUDA->isChecked());
    // Save proj3D;

    // QString outputPath = "D:/ProjTemplate.mha";
    // QString outputPath = "D:/2D3DRegi/FwdProj_0.mha";
    auto outputPath = QFileDialog::getSaveFileName(
        this, "File path to save", this->m_cbctrecon->m_strPathDirDefault,
        "Projection stack (*.mha)", nullptr,
        nullptr); // Filename don't need to exist
    if (outputPath.length() <= 1) {
      return;
    }

    using WriterType = itk::ImageFileWriter<UShortImageType>;
    auto writer = WriterType::New();
    writer->SetFileName(outputPath.toLocal8Bit().constData());
    // writer->SetUseCompression(true);
    writer->SetUseCompression(true); // for plastimatch
    writer->SetInput(this->m_cbctrecon->m_spProjImgRaw3D);
    writer->Update();

    return;
  }
  // if there is a geometry

  const auto cntProj =
      this->m_cbctrecon->m_spCustomGeometry->GetGantryAngles().size();

  if (cntProj < 1) {
    std::cout << "ERROR: geometry is not ready" << std::endl;
    return;
  }

  /*  QMessageBox msgBox;
    msgBox.setText("Do you want to override panel shifts with 0 before forward
    projection?"); msgBox.setStandardButtons(QMessageBox::Yes |
    QMessageBox::No);

    int bOverridePanelShift = msgBox.exec();*/

  // if (res == QMessageBox::Yes)

  // Regenerate geometry object

  for (auto i = 0U; i < cntProj; i++) {
    const auto curSID =
        this->m_cbctrecon->m_spCustomGeometry->GetSourceToIsocenterDistances()
            .at(i);
    const auto curSDD =
        this->m_cbctrecon->m_spCustomGeometry->GetSourceToDetectorDistances()
            .at(i);
    auto curGantryAngle =
        this->m_cbctrecon->m_spCustomGeometry->GetGantryAngles().at(i);
    const auto kVAng = curGantryAngle * 360. / (2. * itk::Math::pi);
    auto MVAng =
        kVAng - (this->m_cbctrecon->m_projFormat == HIS_FORMAT ? 0.0 : 90.0);
    if (MVAng < 0.0) {
      MVAng = MVAng + 360.0;
    }
    curGantryAngle = MVAng;

    const auto curProjOffsetX =
        this->m_cbctrecon->m_spCustomGeometry->GetProjectionOffsetsX().at(i);
    const auto curProjOffsetY =
        this->m_cbctrecon->m_spCustomGeometry->GetProjectionOffsetsY().at(i);

    const auto curOutOfPlaneAngles =
        this->m_cbctrecon->m_spCustomGeometry->GetOutOfPlaneAngles().at(i);
    const auto curInPlaneAngles =
        this->m_cbctrecon->m_spCustomGeometry->GetInPlaneAngles().at(i);

    const auto curSrcOffsetX =
        this->m_cbctrecon->m_spCustomGeometry->GetSourceOffsetsX().at(i);
    const auto curSrcOffsetY =
        this->m_cbctrecon->m_spCustomGeometry->GetSourceOffsetsY().at(i);

    // if (bOverridePanelShift)
    //{
    //    /*curProjOffsetX = 0.0;
    //    curProjOffsetY = 0.0;*/
    //    curProjOffsetX = 0.0;//half fan, shifted toward patient right when kV
    //    source is 0 curProjOffsetY = 0.0;//shifted toward superior
    //}

    crntGeometry->AddProjection(
        curSID, curSDD, curGantryAngle, curProjOffsetX,
        curProjOffsetY,                        // Flexmap
        curOutOfPlaneAngles, curInPlaneAngles, // In elekta, these are 0
        curSrcOffsetX, curSrcOffsetY);         // In elekta, these are 0
  }

  ForwardProjection(this->m_cbctrecon->m_spRawReconImg, crntGeometry,
                    this->m_cbctrecon->m_spProjImgRaw3D, true,
                    this->ui.radioButton_UseCUDA->isChecked());

  // Export geometry txt
  /* QString strPath = QFileDialog::getSaveFileName(this, "Save geometry file
   for forward projection", "", "text (*.txt)", 0, 0);

   if (strPath.length() <= 1)
   return;

   std::vector<double>::const_iterator itAng, itShiftX, itShiftY;

   int cntAngle = crntGeometry->GetGantryAngles().size();
   int cntShiftX = crntGeometry->GetProjectionOffsetsX().size();
   int cntShiftY = crntGeometry->GetProjectionOffsetsY().size();

   if (cntAngle <= 0)
   {
   std::cout << "Error! no angle std::vector is found" << std::endl;
   return;
   }


   if (cntAngle != cntShiftX || cntAngle != cntShiftY)
   {
   std::cout << "Error! Angle number and shift number are not matching." <<
   std::endl; return;
   }

   itShiftX = crntGeometry->GetProjectionOffsetsX().begin();
   itShiftY = crntGeometry->GetProjectionOffsetsY().begin();

   std::ofstream fout;
   fout.open(strPath.toLocal8Bit().constData());

   fout << "MV_Gantry_Angle" << "	" << "PanelShiftX(mm)" << "	" <<
   "PanelShiftY(mm)" << std::endl;

   for (itAng = crntGeometry->GetGantryAngles().begin(); itAng !=
   crntGeometry->GetGantryAngles().end(); itAng++)
   {
   fout << (*itAng) << "	" << (*itShiftX) << "	" << (*itShiftY) <<
   std::endl;

   itShiftX++;
   itShiftY++;
   }

   fout.close();*/
}

void CbctReconWidget::SLTM_FineResolScatterCorrectrionMacro() {
  const auto curResampleF = this->m_cbctrecon->m_fResampleF;
  this->ui.lineEdit_DownResolFactor->setText("1.0");
  SLT_LoadSelectedProjFiles();
  this->m_cbctrecon->m_fResampleF = curResampleF;
  this->ui.lineEdit_DownResolFactor->setText(
      QString("%1").arg(this->m_cbctrecon->m_fResampleF));

  // Scatter correction

  SLT_DoScatterCorrection_APRIORI();
}

void CbctReconWidget::SLTM_FullScatterCorrectionMacroAP() // single. should be
                                                          // called after HIS
                                                          // folder is defined
{
  if (this->m_cbctrecon->m_strPathPatientDir.length() < 2) {
    return;
  }

  const auto enRegImg = REGISTER_DEFORM_FINAL;
  const auto bFullResolForFinalRecon = false;

  const auto bIntensityShift = true;

  FullScatterCorrectionMacroSingle(this->m_cbctrecon->m_strPathPatientDir,
                                   enRegImg, bFullResolForFinalRecon, false,
                                   bIntensityShift);
}

void CbctReconWidget::SLTM_BatchScatterCorrectionMacroAP() {
  // Scatter parameters
  auto batchmodeTime = QTime::currentTime();

  // 1) Get img_ file lists
  const auto dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open IMAGES Directory"), this->m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  auto dirIMAGES = QDir(dirPath);

  if (!dirIMAGES.exists()) {
    return;
  }

  const auto listDir = dirIMAGES.entryInfoList(QDir::Dirs, QDir::Name);

  const auto iCntProjDir = listDir.size();

  std::cout << "Found directory number= " << iCntProjDir - 2 << std::endl;

  if (iCntProjDir <= 2) // only /. and /.. exist
  {
    std::cout << "Error! No projection directory exists." << std::endl;
    return;
  }
  // Several questions to set params
  // 1) Output Dir
  auto strOutDirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Output Directory"), this->m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (strOutDirPath.length() < 1) {
    return;
  }

  // 2) fwd reference: manual or auto-rigid or deformed

  bool ok;
  auto text = QInputDialog::getText(
      this, "Input Dialog",
      "Reference for forward projection([0] Deformed CT, [1] auto-rigid CT, "
      "[2] manual-aligned CT(dcm_plan), [3] DeformedCT_skipAutoRigid",
      QLineEdit::Normal, "0", &ok);

  auto enRegImg = REGISTER_DEFORM_FINAL;

  if (ok && !text.isEmpty()) {
    const auto iRefImgVal = text.toInt();

    if (iRefImgVal == 0) {
      enRegImg = REGISTER_DEFORM_FINAL;
    } else if (iRefImgVal == 1) {
      enRegImg = REGISTER_AUTO_RIGID;
    } else if (iRefImgVal == 2) {
      enRegImg = REGISTER_MANUAL_RIGID;
    } else if (iRefImgVal == 3) {
      enRegImg = REGISTER_DEFORM_SKIP_AUTORIGID;
      //} else {
      //  enRegImg = REGISTER_DEFORM_FINAL; <- initialized value
    }
  }

  // 3) Fine resol option
  const auto bFullResolForFinalRecon = false;
  auto bIntensityShift = false;
  /*QMessageBox msgBox;
  QString strMsg = "Full-resolution reconstruction after scatter generation?";
  msgBox.setText(strMsg);
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  int res = msgBox.exec();
  if (res == QMessageBox::Yes)
  {
  bFullResolForFinalRecon = true;
  }*/

  /* QMessageBox msgBox;
   QString strMsg = "Apply truncation for raw CBCT?";
   msgBox.setText(strMsg);
   msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
   res = msgBox.exec();

   if (res == QMessageBox::Yes)
   {
   bTrancOnlyRaw = true;
   }*/

  QMessageBox msgBox;
  const auto strMsg = QString("Intensity shift for raw CBCT?");
  msgBox.setText(strMsg);
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  auto res = msgBox.exec();

  if (res == QMessageBox::Yes) {
    bIntensityShift = true;
  }

  auto bExportShortImages = false;
  QMessageBox msgBox2;
  const auto strMsg2 = QString("Export short images after correction?");
  msgBox2.setText(strMsg2);
  msgBox2.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  res = msgBox2.exec();
  if (res == QMessageBox::Yes) {
    bExportShortImages = true;
  }

  auto cntHisDir = 0;
  for (auto i = 2; i < iCntProjDir; i++) // start from 2
  {
    auto curProjDirPath = listDir.at(i).absoluteFilePath();

    if (curProjDirPath.contains("img_")) {
      std::cout << "Found projection dir number: " << cntHisDir
                << ", Current Proj Path: "
                << curProjDirPath.toLocal8Bit().constData() << std::endl;
      this->m_cbctrecon->SetProjDir(curProjDirPath);
      FullScatterCorrectionMacroSingle(strOutDirPath, enRegImg,
                                       bFullResolForFinalRecon,
                                       bExportShortImages, bIntensityShift);
      cntHisDir++;
    }
  }

  const auto elapsedSec = batchmodeTime.elapsed() / 1000.0f;

  std::cout << "Batch mode calculation is done! "
            << QString::number(elapsedSec, 'f', 2).toLocal8Bit().constData()
            << " seconds was spent for " << cntHisDir << " cases" << std::endl;
}

// Uses SLT functions heavily:
bool CbctReconWidget::FullScatterCorrectionMacroSingle(
    QString &outputDirPath, const enREGI_IMAGES enFwdRefImg,
    const bool bFullResolRecon, const bool bExportImages,
    const bool bCBCT_IntensityShift) {
  if (this->m_cbctrecon->m_strDCMUID.length() < 1) {
    return false;
  }

  this->m_cbctrecon->m_bMacroContinue = true;

  const auto bFOVCropping = this->ui.checkBox_PostDispObjOn
                                ->isChecked(); // this button is for display,
                                               // but use it for cropping option
                                               // in macro mode
  const auto physPosX = this->ui.lineEdit_PostFOV_X->text().toFloat();
  const auto physPosY = this->ui.lineEdit_PostFOV_Y->text().toFloat();
  const auto physRadius = this->ui.lineEdit_PostFOV_R->text().toFloat();
  const auto physTablePosY = this->ui.lineEdit_PostTablePosY->text().toFloat();

  // Load Pushbutton
  SLT_LoadSelectedProjFiles();

  // float fOldValTruncation =
  // this->ui.lineEdit_Ramp_TruncationCorrection->text().toFloat();;

  SLT_DoReconstruction();
  // this->ui.lineEdit_Ramp_TruncationCorrection->setText(QString("0.0"));

  const auto addingVal = this->ui.lineEdit_AddConstHU->text().toInt();

  if (addingVal != 0) {
    std::cout << "Raw CBCT is being added by HU of: " << addingVal << std::endl;
    SLT_AddConstHUToCurImg();
  }

  if (bFOVCropping) {
    // Crop CBCT with predetermined FOV/ Table
    std::cout << "FOV cropping is under way..." << std::endl;
    this->m_cbctrecon->CropFOV3D(this->m_cbctrecon->m_spRawReconImg, physPosX,
                                 physPosY, physRadius, physTablePosY);
  }

  SLT_ViewRegistration();

  m_dlgRegistration->SLT_PreProcessCT();
  if (!this->m_cbctrecon->m_bMacroContinue) {
    std::cout << "Stopped during MacroSingle due to error in PreProcessCT"
              << std::endl;
    return false;
  }

  QString strSuffix;
  switch (enFwdRefImg) {
  case REGISTER_MANUAL_RIGID:
    m_dlgRegistration->SLT_ManualMoveByDCMPlan();
    m_dlgRegistration->SLT_ConfirmManualRegistration(); // skin cropping for
                                                        // CBCT. only works when
                                                        // CBCT_skin crop is on
    strSuffix = strSuffix + "_man";

    if (bCBCT_IntensityShift) {
      m_dlgRegistration->SLT_IntensityNormCBCT();
    }

    break;
  case REGISTER_AUTO_RIGID:
    m_dlgRegistration->SLT_ManualMoveByDCMPlan();
    m_dlgRegistration->SLT_ConfirmManualRegistration(); // skin cropping

    if (bCBCT_IntensityShift) {
      m_dlgRegistration->SLT_IntensityNormCBCT();
    }

    // OPtional
    if (bFOVCropping) {
      this->m_cbctrecon->CropFOV3D(this->m_cbctrecon->m_spManualRigidCT,
                                   physPosX, physPosY, physRadius,
                                   physTablePosY);
    }

    m_dlgRegistration->SLT_DoRegistrationRigid();
    strSuffix = strSuffix + "_rigid";
    break;
  case REGISTER_DEFORM_FINAL:
    std::cout << "REGISTER_DEFORM_FINAL was chosen." << std::endl;
    m_dlgRegistration->SLT_ManualMoveByDCMPlan();
    m_dlgRegistration->SLT_ConfirmManualRegistration(); // skin cropping

    if (bCBCT_IntensityShift) {
      std::cout << "IntensityShift is underway" << std::endl;
      m_dlgRegistration->SLT_IntensityNormCBCT();
    }

    if (bFOVCropping) {
      this->m_cbctrecon->CropFOV3D(this->m_cbctrecon->m_spManualRigidCT,
                                   physPosX, physPosY, physRadius,
                                   physTablePosY);
    }

    m_dlgRegistration->SLT_DoRegistrationRigid();
    m_dlgRegistration->SLT_DoRegistrationDeform();
    strSuffix = strSuffix + "_defrm";
    break;

  case REGISTER_DEFORM_SKIP_AUTORIGID:
    m_dlgRegistration->SLT_ManualMoveByDCMPlan();
    m_dlgRegistration->SLT_ConfirmManualRegistration(); // skin cropping
    // m_pDlgRegistration->SLT_DoRegistrationRigid();

    if (bCBCT_IntensityShift) {
      m_dlgRegistration->SLT_IntensityNormCBCT();
    }

    if (bFOVCropping) {
      this->m_cbctrecon->CropFOV3D(this->m_cbctrecon->m_spManualRigidCT,
                                   physPosX, physPosY, physRadius,
                                   physTablePosY);
    }

    m_dlgRegistration->SLT_DoRegistrationDeform();
    strSuffix = strSuffix + "_defrm_skipRigid";
    break;

  default:
    break;
  }

  if (bFullResolRecon) // if fullResol recon is on, load original proj files
                       // again
  {
    const auto curResampleF = this->m_cbctrecon->m_fResampleF;
    this->ui.lineEdit_DownResolFactor->setText("1.0");
    SLT_LoadSelectedProjFiles();
    this->m_cbctrecon->m_fResampleF = curResampleF;
    this->ui.lineEdit_DownResolFactor->setText(
        QString("%1").arg(this->m_cbctrecon->m_fResampleF));

    strSuffix = strSuffix + "_HD";
  } else {
    // strSuffix = strSuffix + "_HD"
  }
  SLT_DoScatterCorrection_APRIORI();

  // If there is couch shift information and this cbct is a pre-treatment CBCT,
  // couch shift can be applied  to represent the final treatment position

  if (this->ui.checkBox_CouchShiftAddToMacro->isChecked()) {
    SLT_DoCouchCorrection();
  }

  // 1) Save the corrCBCT image as signed short
  auto outputPath_rawCBCT = outputDirPath + "/" +
                            this->m_cbctrecon->m_strDCMUID + strSuffix +
                            "_rawCBCT.mha";
  auto outputPath_corrCBCT = outputDirPath + "/" +
                             this->m_cbctrecon->m_strDCMUID + strSuffix +
                             "_corrCBCT.mha";
  auto outputPath_manCT = outputDirPath + "/" + this->m_cbctrecon->m_strDCMUID +
                          strSuffix + "_manCT.mha";
  auto outputPath_rigidCT = outputDirPath + "/" +
                            this->m_cbctrecon->m_strDCMUID + strSuffix +
                            "_rigidCT.mha";
  auto outputPath_deformCT = outputDirPath + "/" +
                             this->m_cbctrecon->m_strDCMUID + strSuffix +
                             "_deformCT.mha";

  if (bExportImages) {
    ExportReconSHORT_HU(this->m_cbctrecon->m_spRawReconImg, outputPath_rawCBCT);
    ExportReconSHORT_HU(this->m_cbctrecon->m_spScatCorrReconImg,
                        outputPath_corrCBCT);
    ExportReconSHORT_HU(this->m_cbctrecon->m_spManualRigidCT, outputPath_manCT);
    ExportReconSHORT_HU(this->m_cbctrecon->m_spAutoRigidCT, outputPath_rigidCT);
    ExportReconSHORT_HU(this->m_cbctrecon->m_spDeformedCT_Final,
                        outputPath_deformCT);
  }

  // 2) Calculate batched WEPL points
  if (!this->m_cbctrecon->m_vPOI_DCM.empty()) {
    auto outputTxtPath = outputDirPath + "/" + this->m_cbctrecon->m_strDCMUID +
                         strSuffix + "_WEPL.txt";

    const auto fAngleStart = this->ui.lineEdit_AngStart->text().toDouble();
    const auto fAngleEnd = this->ui.lineEdit_AngEnd->text().toDouble();
    const auto fAngleGap = this->ui.lineEdit_WEPL_AngRes->text().toDouble();

    this->m_cbctrecon->ExportAngularWEPL_byFile(outputTxtPath, fAngleStart,
                                                fAngleEnd, fAngleGap);
  }
  return true;
}

void CbctReconWidget::SLT_OpenPhaseData() {
  if (!this->m_cbctrecon->m_vPhaseFloat.empty()) {
    this->m_cbctrecon->m_vPhaseFloat.clear();
  }

  // Open file
  auto filePath = QFileDialog::getOpenFileName(
      this, "Open phase text", this->m_cbctrecon->m_strPathDirDefault,
      "Phase text file (*.txt)", nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  std::ifstream fin;
  fin.open(filePath.toLocal8Bit().constData(), std::ios::in);
  if (fin.fail()) {
    return;
  }

  this->ui.lineEdit_PhaseTxtPath->setText(filePath);

  char str[MAX_LINE_LENGTH];

  float tmpPhase = 0.0;

  float phaseSum = 0.0;
  auto phaseCnt = 0;
  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString strLine(&str[0]);

    if (strLine.length() < 1) {
      break;
    }

    tmpPhase = strLine.toFloat();
    this->m_cbctrecon->m_vPhaseFloat.push_back(tmpPhase);
    phaseCnt++;
    phaseSum = phaseSum + tmpPhase;
  }
  fin.close();
  std::cout << "NumOfPhaseData[Float]= " << phaseCnt << "  Mean Phase value= "
            << phaseSum / static_cast<double>(phaseCnt) << std::endl;
}

void CbctReconWidget::SLT_Export4DCBCT() const {
  if (this->m_cbctrecon->m_spCustomGeometry == nullptr) {
    std::cout << "Error! no Geometry information loaded yet" << std::endl;
    return;
  }

  const auto NumOfGanAngle =
      this->m_cbctrecon->m_spCustomGeometry->GetGantryAngles().size();
  const auto NumOfPhase = this->m_cbctrecon->m_vPhaseFloat.size();

  if (NumOfGanAngle != NumOfPhase) {
    std::cout << "Size not matched. NumOfProjection= " << NumOfGanAngle
              << " NumOfProjection= " << NumOfPhase << std::endl;
    return;
  }
  // build phase bins
  auto strPhaseTextFull = this->ui.lineEdit_PhaseExportString->text();
  auto strlistPhaseFull = strPhaseTextFull.split(";");

  const auto cntGroup = strlistPhaseFull.count();

  std::vector<int> vPhaseBinsSelected;
  // m_strPatientDirName = tmpDir_PatientFolder.dirName();
  // m_strPathPatientDir = tmpDir_PatientFolder.absolutePath();

  // QDir tmpDir_RootFolder(movingDir.absolutePath()); //root folder

  // if (tmpDir_RootFolder.absolutePath().length() > 1)
  //    m_strPathDirDefault = tmpDir_RootFolder.absolutePath();

  ////option 1: already made rtk xml file
  // QString tmpPathRTKGeometry = tmpDir_RootFolder.absolutePath() + "/" +
  // "ElektaGeom_" + m_strDCMUID + ".xml";  QFileInfo
  // rtkGeomInfo(tmpPathRTKGeometry);

  auto strDirForXML =
      this->m_cbctrecon->m_strPathDirDefault; // where xml file is located
  // QString strUID ;//P00102030P + m_strDCMUID
  auto strDirForProj = this->m_cbctrecon->m_strPathIMAGES;

  for (auto i = 0; i < cntGroup; i++) {
    // for a single group
    auto strListGroup = strlistPhaseFull.at(i).split(",");
    const auto iPhaseCnt = strListGroup.count();
    vPhaseBinsSelected.clear();

    for (auto j = 0; j < iPhaseCnt; j++) {
      vPhaseBinsSelected.push_back(strListGroup.at(j).toInt());
    }
    // m_vSelectedFileNames: full file paths of projections
    // m_spCustomGeometry: full information
    // m_vPhaseFloat: full data of phase

    // Create Dir, xml, etc
    if (!this->m_cbctrecon->ResortCBCTProjection(
            vPhaseBinsSelected, strDirForXML, strDirForProj,
            this->m_cbctrecon->m_strDCMUID, this->m_cbctrecon->m_vPhaseFloat,
            this->m_cbctrecon->m_spCustomGeometry,
            this->m_cbctrecon->m_vSelectedFileNames)) {
      std::cout << "Error in ResortCBCTProjection "
                << strlistPhaseFull.at(i).toLocal8Bit().constData()
                << std::endl;
      return;
    }
  }

  // mkdir
  // QString strCrntDir = m_strPathPatientDir + "/" + "IMAGES"; //current Proj
  // folder  QString strCrntDir = this->ui.lineEdit_HisDirPath->text();

  ////Make a sub directory
  // QDir crntDir(strCrntDir);

  // if (!crntDir.exists())
  //{
  //    std::cout << "File save error: The specified folder does not exist." <<
  //    std::endl; return;
  //}

  ////QString fwdDirName = "fwd_" + m_strDCMUID;

  // QString str4D = "4DCBCT";
  // bool tmpResult = crntDir.mkdir(str4D);

  // QString strSubDir = strCrntDir + "/" + str4D;
  // QDir crntSubDir(strSubDir);
  // if (!crntSubDir.exists())
  //{
  //    std::cout << "File save error" << std::endl;
  //    return;
  //}

  // crntSubDir.mkdir();
}

template <enREGI_IMAGES imagetype> void CbctReconWidget::LoadMHAfileAs() {
  auto fileName = QFileDialog::getOpenFileName(
      this, "Open Image", this->m_cbctrecon->m_strPathDirDefault,
      "Short image file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  this->m_cbctrecon->LoadShort3DImage(fileName, imagetype);
  auto imgDim =
      this->m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize();

  this->ui.lineEdit_Cur3DFileName->setText(fileName);

  this->ui.spinBoxReconImgSliceNo->setMinimum(0);
  this->ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  const auto initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  this->ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called
  this->ui.radioButton_graph_recon->setChecked(true);

  m_dlgRegistration->UpdateListOfComboBox(0); // combo selection
                                              // signalis called
  m_dlgRegistration->UpdateListOfComboBox(1);
  m_dlgRegistration->SelectComboExternal(
      0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  m_dlgRegistration->SelectComboExternal(1, imagetype);
  // std::cout << m_spScatCorrReconImg->GetBufferedRegion().GetSize() <<
  // std::endl;

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_LoadCBCTcorrMHA() {
  LoadMHAfileAs<REGISTER_COR_CBCT>();
}

void CbctReconWidget::SLT_LoadCTrigidMHA() {
  LoadMHAfileAs<REGISTER_AUTO_RIGID>();
}

void CbctReconWidget::SLT_LoadCTdeformMHA() {
  LoadMHAfileAs<REGISTER_DEFORM_FINAL>();
}

void CbctReconWidget::SLT_DoCouchCorrection() {
  auto strTrans = this->ui.lineEdit_CouchTrans->text();
  auto strRot = this->ui.lineEdit_CouchRot->text();

  auto strListTrans = strTrans.split(",");
  auto strListRot = strRot.split(",");

  if (strListTrans.count() != 3 || strListRot.count() != 3) {
    std::cout << "Error! No couch shift data is available!" << std::endl;
    return;
  }

  const auto couchShiftTrans =
      VEC3D{strListTrans.at(0).toDouble(), // mm
            strListTrans.at(1).toDouble(), strListTrans.at(2).toDouble()};

  const auto couchShiftRot =
      VEC3D{strListRot.at(0).toDouble(), strListRot.at(1).toDouble(),
            strListRot.at(2).toDouble()};
  // Images to correct:
  /*m_spRawReconImg;
  m_spScatCorrReconImg;
  m_spDeformedCT_Final;
  m_spAutoRigidCT;*/

  // not manual CT!!!

  ImageTransformUsingCouchCorrection(this->m_cbctrecon->m_spRawReconImg,
                                     this->m_cbctrecon->m_spRawReconImg,
                                     couchShiftTrans, couchShiftRot);
  ImageTransformUsingCouchCorrection(this->m_cbctrecon->m_spScatCorrReconImg,
                                     this->m_cbctrecon->m_spScatCorrReconImg,
                                     couchShiftTrans, couchShiftRot);
  ImageTransformUsingCouchCorrection(this->m_cbctrecon->m_spDeformedCT_Final,
                                     this->m_cbctrecon->m_spDeformedCT_Final,
                                     couchShiftTrans, couchShiftRot);
  ImageTransformUsingCouchCorrection(this->m_cbctrecon->m_spAutoRigidCT,
                                     this->m_cbctrecon->m_spAutoRigidCT,
                                     couchShiftTrans, couchShiftRot);

  m_dlgRegistration->UpdateListOfComboBox(0); // combo selection
                                              // signalis called
  m_dlgRegistration->UpdateListOfComboBox(1);
  m_dlgRegistration->SelectComboExternal(
      0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  m_dlgRegistration->SelectComboExternal(1, REGISTER_COR_CBCT);

  this->m_cbctrecon->m_spCrntReconImg = this->m_cbctrecon->m_spScatCorrReconImg;
  SLT_DrawReconImage();

  std::cout << "Couch shift and rotation was successfully applied."
            << std::endl;
}

// Multiple mha files
void CbctReconWidget::SLTM_WELPCalcMultipleFiles() {
  // Singed short
  auto listFilePath = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open",
      this->m_cbctrecon->m_strPathDirDefault, "signed short 3D images (*.mha)");

  const auto iCntFiles = listFilePath.count();
  if (iCntFiles < 1) {
    return;
  }

  int iCntPOI = this->m_cbctrecon->m_vPOI_DCM.size();

  if (iCntPOI < 1) {
    std::cout << "There is no POI file loaded." << std::endl;
    SLT_LoadPOIData();
  }
  iCntPOI = this->m_cbctrecon->m_vPOI_DCM.size();
  if (iCntPOI < 1) {
    std::cout << "Error! still no POI" << std::endl;
    return;
  }

  auto strPathOutText = QFileDialog::getSaveFileName(
      this, "File path to save", this->m_cbctrecon->m_strPathDirDefault,
      "WEPL_value (*.txt)", nullptr, nullptr); // Filename don't need to exist
  if (strPathOutText.length() <= 1) {
    return;
  }

  std::vector<std::vector<WEPLData>> vArrOutputWEPL;
  vArrOutputWEPL.resize(iCntFiles);

  const auto fAngleStart = this->ui.lineEdit_AngStart->text().toDouble();
  const auto fAngleEnd = this->ui.lineEdit_AngEnd->text().toDouble();
  for (auto i = 0; i < iCntFiles; i++) {
    this->m_cbctrecon->GetWEPLDataFromSingleFile(
        listFilePath.at(i), this->m_cbctrecon->m_vPOI_DCM, vArrOutputWEPL.at(i),
        fAngleStart, fAngleEnd);
  }

  std::ofstream fout;
  fout.open(strPathOutText.toLocal8Bit().constData());

  fout << "Point Index"
       << "\t"
       << "Gantry Angle"
       << "\t"
       << "Sample Number";

  for (auto i = 0; i < iCntFiles; i++) {
    QFileInfo fInfo(listFilePath.at(i));
    auto strFileName = fInfo.fileName();

    fout << "\t" << strFileName.toLocal8Bit().constData();
  }
  fout << std::endl;

  const auto cnt_wepl = vArrOutputWEPL.at(0).size();
  for (auto i = 0; i < iCntFiles; i++) {
    const auto cur_count = vArrOutputWEPL.at(i).size();
    if (cnt_wepl != cur_count) {
      std::cout << "Error! some of the WEPL count doesn't match!" << std::endl;
      return;
    }
  }

  for (auto i = 0U; i < cnt_wepl; i++) {
    fout << vArrOutputWEPL.at(0).at(i).ptIndex << "\t"
         << vArrOutputWEPL.at(0).at(i).fGanAngle << "\t" << i;

    for (auto j = 0; j < iCntFiles; j++) {
      fout << "\t" << vArrOutputWEPL.at(j).at(i).fWEPL;
    }
    fout << std::endl;
  }

  fout.close();

  std::cout << "Saving angular WEPL is completed" << std::endl;

  // delete[] vArrOutputWEPL;
}

void CbctReconWidget::SLTM_ScatterCorPerProjRef() // load text file
{
  const auto scaMedian = this->ui.lineEdit_scaMedian->text().toDouble();
  const auto scaGaussian = this->ui.lineEdit_scaGaussian->text().toDouble();
  const auto postScatMedianSize =
      this->ui.lineEdit_scaPostMedian->text().toInt();

  const auto fdk_options = getFDKoptions();

  this->m_cbctrecon->ScatterCorPerProjRef(
      scaMedian, scaGaussian, postScatMedianSize,
      this->ui.radioButton_UseCUDA->isChecked(),
      this->ui.radioButton_UseOpenCL->isChecked(),
      this->ui.checkBox_ExportVolDICOM->isChecked(),
      fdk_options); // load text file

  SLT_DrawProjImages();
}

void CbctReconWidget::SLTM_LoadPerProjRefList() {
  auto filePath = QFileDialog::getOpenFileName(
      this, "PerProjVol list", this->m_cbctrecon->m_strPathDirDefault,
      "File path list (*.txt)", nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  std::ifstream fin;
  fin.open(filePath.toLocal8Bit().constData(), std::ios::in);
  if (fin.fail()) {
    return;
  }

  this->m_cbctrecon->m_strListPerProjRefVol.clear();

  char str[MAX_LINE_LENGTH];
  // File format:
  // header [Projection]	 [phase] 	 [amplitude] 	 [filename]
  // 0	0.45	1	D:/4DCT_ScatterCor/03_MotionModelA/
  // model_out_N_0000_phase_0.45_amp_1.mha

  memset(&str[0], 0, MAX_LINE_LENGTH);
  fin.getline(&str[0], MAX_LINE_LENGTH);

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString strLine(&str[0]);

    auto strList = strLine.split('\t'); // tab

    if (strList.count() > 0 && strList.count() < 4) {
      std::cout << "Str = " << strLine.toLocal8Bit().constData() << std::endl;
      std::cout << "abnormal file expression." << std::endl;
      break;
    }

    this->m_cbctrecon->m_strListPerProjRefVol.push_back(strList.at(3));
  }

  std::cout << this->m_cbctrecon->m_strListPerProjRefVol.count()
            << " image paths were found" << std::endl;

  fin.close();
}

void CbctReconWidget::SLTM_CropMaskBatch() {
  // Specify mask file (USHORT)
  auto maskFilePath = QFileDialog::getOpenFileName(
      this, "Mask image (Ushort)", this->m_cbctrecon->m_strPathDirDefault,
      "3D mask file (*.mha)", nullptr, nullptr);

  if (maskFilePath.length() < 1) {
    return;
  }

  // QString strPath_mskSkinCT_final;
  // QString strPath_mskSkinCT_autoRegi_exp = m_strPathPlastimatch +
  // "/msk_skin_CT_autoRegi_exp.mha";  QFileInfo
  // maskInfoAuto(strPath_mskSkinCT_autoRegi_exp);

  // QString strPath_mskSkinCT_manualRegi_exp = m_strPathPlastimatch +
  // "/msk_skin_CT_manRegi_exp.mha";  QFileInfo
  // maskInfoManual(strPath_mskSkinCT_manualRegi_exp);

  // if (maskInfoAuto.exists()) //if the mask file is not prepared, give up the
  // skin removal
  //{
  //    strPath_mskSkinCT_final = strPath_mskSkinCT_autoRegi_exp;
  //}
  // else
  //{
  //    std::cout << "Mask file of auto-registration is not prepared. Use manual
  //    regi-mask instead" << std::endl;

  //    if (maskInfoManual.exists())
  //    {
  //        strPath_mskSkinCT_final = strPath_mskSkinCT_manualRegi_exp;
  //    }
  //    else
  //    {
  //        std::cout << "Mask file of manual registration is not prepared. Skip
  //        skin removal!" << std::endl; return;
  //    }
  //}

  // Get File names (SHORT) where to apply Mask cropping
  auto targetFilePaths = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open",
      this->m_cbctrecon->m_strPathDirDefault, "target files (*.mha)");

  const auto iCnt = targetFilePaths.size();

  if (iCnt < 1) {
    return;
  }

  for (auto i = 0; i < iCnt; i++) {
    const auto &curPath = targetFilePaths.at(i);

    // Overritting
    const auto mask_option = MASK_OPERATION_MASK;
    QString input_fn = curPath.toLocal8Bit().constData();
    // QString mask_fn = strPath_mskSkinCT_final.toLocal8Bit().constData();
    QString mask_fn = maskFilePath.toLocal8Bit().constData();
    QString output_fn = curPath.toLocal8Bit().constData();
    const auto mask_value = -1024.0; // unsigned short
    m_cbctregistration->plm_mask_main(mask_option, input_fn, mask_fn, output_fn,
                                      mask_value);
    std::cout << i + 1 << "/" << iCnt << std::endl;
  }

  //  QString strPath_mskSkinCT_final;
  // QString strPath_mskSkinCT_autoRegi_exp = m_strPathPlastimatch +
  // "/msk_skin_CT_autoRegi_exp.mha";  QFileInfo
  // maskInfoAuto(strPath_mskSkinCT_autoRegi_exp);

  // QString strPath_mskSkinCT_manualRegi_exp = m_strPathPlastimatch +
  // "/msk_skin_CT_manRegi_exp.mha";  QFileInfo
  // maskInfoManual(strPath_mskSkinCT_manualRegi_exp);

  // if (maskInfoAuto.exists()) //if the mask file is not prepared, give up the
  // skin removal
  //{
  //    strPath_mskSkinCT_final = strPath_mskSkinCT_autoRegi_exp;
  //}
  // else
  //{
  //    std::cout << "Mask file of auto-registration is not prepared. Use manual
  //    regi-mask instead" << std::endl;

  //    if (maskInfoManual.exists())
  //    {
  //        strPath_mskSkinCT_final = strPath_mskSkinCT_manualRegi_exp;
  //    }
  //    else
  //    {
  //        std::cout << "Mask file of manual registration is not prepared. Skip
  //        skin removal!" << std::endl; return;
  //    }
  //}
}

void CbctReconWidget::SLT_CropSupInf() {
  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  auto bRaw = 0;

  if (this->m_cbctrecon->m_spCrntReconImg ==
      this->m_cbctrecon->m_spRawReconImg) {
    bRaw = 1;
  }

  const auto dcmPosCutSup = this->ui.lineEdit_SupCutPos->text().toFloat(); // mm
  const auto dcmPosCutInf = this->ui.lineEdit_InfCutPos->text().toFloat(); // mm

  // CropFOV3D(m_spCrntReconImg, physPosX, physPosY, physRadius, physTablePosY);
  this->m_cbctrecon->CropSupInf(this->m_cbctrecon->m_spCrntReconImg,
                                dcmPosCutInf, dcmPosCutSup);
  // QString strPath = m_strPathDirDefault + "/" + "TempSI_Cropped.mha";
  // QString strTmpFile = "C:/TmpSI_Cropped.mha";

  auto strPath =
      m_cbctregistration->m_strPathPlastimatch + "/" + "tmp_SI_cropped.mha";
  ExportReconSHORT_HU(this->m_cbctrecon->m_spCrntReconImg, strPath);

  QString strName = "SI_Cropped";
  if (bRaw != 0) {
    if (!LoadShortImageToUshort(strPath, this->m_cbctrecon->m_spRawReconImg)) {
      std::cout << "error! in LoadShortImageToUshort" << std::endl;
    }
    UpdateReconImage(this->m_cbctrecon->m_spRawReconImg, strName);
  } else {
    if (LoadShortImageToUshort(strPath, this->m_cbctrecon->m_spRefCTImg)) {
      std::cout << "error! in LoadShortImageToUshort" << std::endl;
    }
    UpdateReconImage(this->m_cbctrecon->m_spRefCTImg, strName);
  }

  ///*So buggy*/

  /*QString strName = "SI_Cropped";
  UpdateReconImage(m_spCrntReconImg, strName);*/

  // SLT_DrawReconImage();
}

// IO_START
void CbctReconWidget::SLT_LoadImageFloat3D() // Dose image for JPhillips
{
  using ReaderType = itk::ImageFileReader<FloatImageType>;
  auto reader = ReaderType::New();

  auto fileName = QFileDialog::getOpenFileName(
      this, "Open Image", this->m_cbctrecon->m_strPathDirDefault,
      "3D dose float file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  reader->SetFileName(fileName.toLocal8Bit().constData());
  reader->Update();

  // Multiply: Gy to mGy
  using MultiplyImageFilterType =
      itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  auto multiplyImageFilter = MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(reader->GetOutput());
  multiplyImageFilter->SetConstant(100.0); // calculated already //Gy to cGy

  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
  auto castFilter = CastFilterType::New();
  castFilter->SetInput(multiplyImageFilter->GetOutput());

  castFilter->Update();

  this->m_cbctrecon->m_spRawReconImg = castFilter->GetOutput();
  this->m_cbctrecon->m_spCrntReconImg = this->m_cbctrecon->m_spRawReconImg;

  // Update UI
  auto imgDim =
      this->m_cbctrecon->m_spRawReconImg->GetBufferedRegion().GetSize();
  auto spacing = this->m_cbctrecon->m_spRawReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  this->ui.lineEdit_Cur3DFileName->setText(fileName);

  this->m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  this->ui.spinBoxReconImgSliceNo->setMinimum(0);
  this->ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  const auto initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  this->ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  this->ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_Load3DImage() // mha reconstructed file, from external
                                        // source
{
  using ReaderType = itk::ImageFileReader<UShortImageType>;
  auto reader = ReaderType::New();

  auto fileName = QFileDialog::getOpenFileName(
      this, "Open Image", this->m_cbctrecon->m_strPathDirDefault,
      "Projection file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  reader->SetFileName(fileName.toLocal8Bit().constData());
  reader->Update();

  this->m_cbctrecon->m_spRawReconImg = reader->GetOutput();
  this->m_cbctrecon->m_spCrntReconImg = this->m_cbctrecon->m_spRawReconImg;

  auto imgDim =
      this->m_cbctrecon->m_spRawReconImg->GetBufferedRegion().GetSize();
  auto spacing = this->m_cbctrecon->m_spRawReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  this->ui.lineEdit_Cur3DFileName->setText(fileName);

  this->m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  this->ui.spinBoxReconImgSliceNo->setMinimum(0);
  this->ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  const auto initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  this->ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  this->ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_Load3DImageShort() {
  auto fileName = QFileDialog::getOpenFileName(
      this, "Open Image", this->m_cbctrecon->m_strPathDirDefault,
      "short mha file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  if (!LoadShortImageToUshort(fileName, this->m_cbctrecon->m_spRawReconImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  auto imageCalculatorFilter2 = ImageCalculatorFilterType2::New();
  // imageCalculatorFilter2->SetImage(m_spReconImg);
  imageCalculatorFilter2->SetImage(this->m_cbctrecon->m_spRawReconImg);
  imageCalculatorFilter2->Compute();

  const auto minVal2 =
      static_cast<double>(imageCalculatorFilter2->GetMinimum());
  const auto maxVal2 =
      static_cast<double>(imageCalculatorFilter2->GetMaximum());
  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  this->m_cbctrecon->m_spCrntReconImg = this->m_cbctrecon->m_spRawReconImg;

  auto imgDim =
      this->m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize();
  auto spacing = this->m_cbctrecon->m_spCrntReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  this->ui.lineEdit_Cur3DFileName->setText(fileName);

  this->m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  this->ui.spinBoxReconImgSliceNo->setMinimum(0);
  this->ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  const auto initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  this->ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  this->ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_LoadNKIImage() {
  auto filePath = QFileDialog::getOpenFileName(
      this, "Open Image", this->m_cbctrecon->m_strPathDirDefault,
      "NKI file (*.SCAN)", nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  auto v =
      nki_load(filePath.toLocal8Bit().constData()); // NKI is unsigned short!!!
  if (v == nullptr) {
    std::cerr << "file reading error" << std::endl;
    return;
  }

  const auto endFix = QString("_conv");
  auto srcFileInfo = QFileInfo(filePath);
  auto dir = srcFileInfo.absoluteDir();
  auto baseName = srcFileInfo.completeBaseName();
  const auto extName = QString("mha");

  const auto newFileName = baseName.append(endFix).append(".").append(extName);
  auto newPath = dir.absolutePath() + "/" + newFileName;

  write_mha(newPath.toLocal8Bit().constData(), v);
  std::cout << "File conversion is done. Trying to read mha file.."
            << std::endl;
  // corrImg.ReleaseBuffer();
  // NKI to mha

  if (!LoadShortImageToUshort(newPath, this->m_cbctrecon->m_spRawReconImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  auto imageCalculatorFilter2 = ImageCalculatorFilterType2::New();
  // imageCalculatorFilter2->SetImage(m_spReconImg);
  imageCalculatorFilter2->SetImage(this->m_cbctrecon->m_spRawReconImg);
  imageCalculatorFilter2->Compute();

  const auto minVal2 =
      static_cast<double>(imageCalculatorFilter2->GetMinimum());
  const auto maxVal2 =
      static_cast<double>(imageCalculatorFilter2->GetMaximum());
  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  this->m_cbctrecon->m_spCrntReconImg = this->m_cbctrecon->m_spRawReconImg;

  auto imgDim =
      this->m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize();
  auto spacing = this->m_cbctrecon->m_spCrntReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  this->ui.lineEdit_Cur3DFileName->setText(newFileName);

  this->m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  this->ui.spinBoxReconImgSliceNo->setMinimum(0);
  this->ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  const auto initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  this->ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  this->ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_ExportHis() {
  if (this->m_cbctrecon->m_iImgCnt < 1) {
    std::cout << "Error: Load raw his images first" << std::endl;
    return;
  }

  // Get Folder Name!

  // For displaying Dir only..
  const auto dir = QFileDialog::getExistingDirectory(
      this, "Open Directory", "/home",
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  // FileName should be same, only selected folder

  for (auto i = 0; i < this->m_cbctrecon->m_iImgCnt; i++) {
    auto tmpInfo = QFileInfo(this->m_cbctrecon->m_arrYKImage[i].m_strFilePath);
    auto newPath = dir + "/" + tmpInfo.fileName();
    this->m_cbctrecon->m_arrYKImage[i].SaveDataAsHis(
        newPath.toLocal8Bit().constData(), false);
  }

  std::cout << "File export was done successfully" << std::endl;
}

void CbctReconWidget::SLT_ReloadProjections() {
  QFile projFile("Projections.mha");
  if (!projFile.exists()) {
    std::cerr
        << "Projections were never saved! i.e. Projections.mha doesn't exist."
        << std::endl;
    return;
  }
  std::cout << "Reading: " << projFile.fileName().toStdString() << std::endl;
  using ImageReaderType = itk::ImageFileReader<FloatImageType>;
  auto ImageReader = ImageReaderType::New();
  ImageReader->SetFileName(projFile.fileName().toStdString());
  ImageReader->Update();
  this->m_cbctrecon->m_spProjImg3DFloat = ImageReader->GetOutput();

  // Copied from SLT_LoadSelectedFiles:

  if (this->m_cbctrecon->m_fResampleF != 1.0) {
    this->m_cbctrecon->ResampleItkImage(
        this->m_cbctrecon->m_spProjImg3DFloat,
        this->m_cbctrecon->m_spProjImg3DFloat,
        this->m_cbctrecon
            ->m_fResampleF); // was! BROKEN AF for .his where input size
                             // != 1024 (tested with 1016) -> outputs
                             // offset -inputoffset/refactor^2 and 4
                             // pixels too few in x and y
  }

  if (this->m_cbctrecon->m_projFormat == HND_FORMAT) { // -> hnd
    std::cout << "Fitted bowtie-filter correction ongoing..." << std::endl;
    SLT_DoBowtieCorrection();
  }

  this->m_cbctrecon->ConvertLineInt2Intensity(
      this->m_cbctrecon->m_spProjImg3DFloat,
      this->m_cbctrecon->m_spProjImgRaw3D,
      65535); // if X not 1024 == input size: out_offset =
              // in_offset + (1024*res_f -
              // X*res_f)*out_spacing     <- will still
              // break down at fw_projection

  auto originPt = this->m_cbctrecon->m_spProjImg3DFloat->GetOrigin();
  auto FloatImgSize =
      this->m_cbctrecon->m_spProjImg3DFloat->GetBufferedRegion().GetSize();
  auto FloatImgSpacing = this->m_cbctrecon->m_spProjImg3DFloat->GetSpacing();

  std::cout << "YKDEBUG: Origin" << originPt[0] << ", " << originPt[1] << ", "
            << originPt[2] << std::endl;
  std::cout << "YKDEBUG: Size" << FloatImgSize[0] << ", " << FloatImgSize[1]
            << ", " << FloatImgSize[2] << std::endl;
  std::cout << "YKDEBUG: Spacing" << FloatImgSpacing[0] << ", "
            << FloatImgSpacing[1] << ", " << FloatImgSpacing[2] << std::endl;

  std::cout
      << "Raw3DProj dimension "
      << this->m_cbctrecon->m_spProjImgRaw3D->GetRequestedRegion().GetSize()
      << std::endl;

  // m_spProjImgRaw3D is Ushort

  std::cout << "Projection reading succeeded."
            << this->m_cbctrecon->m_vSelectedFileNames.size()
            << " files were read" << std::endl;

  // Because you can load projections from previous run:
  this->ui.pushButton_DoRecon->setEnabled(true);

  this->ui.spinBoxImgIdx->setMinimum(0);
  this->ui.spinBoxImgIdx->setMaximum(
      this->m_cbctrecon->m_vSelectedFileNames.size() - 1);
  this->ui.spinBoxImgIdx->setValue(0);

  this->m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();

  this->SLT_DrawProjImages(); // Update Table is called
}

void CbctReconWidget::SLT_LoadPlanCT_mha() // m_spRecon -->m_spRefCT
{
  // typedef itk::ImageFileReader<ShortImageType> ReaderType;
  // ReaderType::Pointer reader = ReaderType::New();

  auto fileName = QFileDialog::getOpenFileName(
      this, "Open Image", this->m_cbctrecon->m_strPathDirDefault,
      "Plan CT file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  if (!LoadShortImageToUshort(fileName, this->m_cbctrecon->m_spRefCTImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  auto imageCalculatorFilter2 = ImageCalculatorFilterType2::New();
  imageCalculatorFilter2->SetImage(this->m_cbctrecon->m_spRefCTImg);
  imageCalculatorFilter2->Compute();

  const auto minVal2 =
      static_cast<double>(imageCalculatorFilter2->GetMinimum());
  const auto maxVal2 =
      static_cast<double>(imageCalculatorFilter2->GetMaximum());

  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  auto imgDim = this->m_cbctrecon->m_spRefCTImg->GetBufferedRegion().GetSize();
  auto spacing = this->m_cbctrecon->m_spRefCTImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  this->m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT,
                                            REGISTER_MANUAL_RIGID);
  this->m_cbctrecon->m_spCrntReconImg = this->m_cbctrecon->m_spRefCTImg;

  this->ui.lineEdit_Cur3DFileName->setText(fileName);
  this->m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  this->ui.spinBoxReconImgSliceNo->setMinimum(0);
  this->ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  const auto initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  this->ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called
  this->ui.radioButton_graph_recon->setChecked(true);
}

void CbctReconWidget::SLT_ExportReconUSHORT() {
  if (this->m_cbctrecon->m_spCrntReconImg == nullptr) {
    std::cout << " no image to export" << std::endl;
    return;
  }

  auto strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "unsigned short meta image (*.mha)", nullptr,
      nullptr);
  if (strPath.length() <= 1) {
    return;
  }

  using WriterType = itk::ImageFileWriter<UShortImageType>;
  auto writer = WriterType::New();
  writer->SetFileName(strPath.toLocal8Bit().constData());
  writer->SetUseCompression(true); // not exist in original code (rtkfdk)
  writer->SetInput(this->m_cbctrecon->m_spCrntReconImg);

  std::cout << "Writing is under progress...: "
            << strPath.toLocal8Bit().constData() << std::endl;
  writer->Update();
  std::cout << "Writing was successfully done" << std::endl;

  const auto msgStr = QString("USHORT File Writing was successfully done");
  QMessageBox::information(this, "Procedure Done", msgStr);
}

// Function for independent projection his images
void CbctReconWidget::LoadRawHisImages() {

  auto files =
      QFileDialog::getOpenFileNames(this, "Select one or more files to open",
                                    this->m_cbctrecon->m_strPathDirDefault,
                                    "projection images (*.his,*.hnd,*.xim)");

  this->m_cbctrecon->m_iImgCnt = files.size();
  std::vector<std::string> fileVector;

  for (auto &cur_file : files) {
    fileVector.push_back(cur_file.toStdString());
  }

  if (this->m_cbctrecon->m_iImgCnt < 1) {
    return;
  }
  m_pTableModel.reset();
  this->m_cbctrecon->ReleaseMemory();

  this->m_cbctrecon->m_arrYKImage.resize(this->m_cbctrecon->m_iImgCnt);
  using ReaderType = rtk::ProjectionsReader<FloatImageType>;
  auto reader = ReaderType::New();
  reader->SetFileNames(fileVector);
  reader->UpdateOutputInformation();
  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
  auto castFilter = CastFilterType::New();
  castFilter->SetInput(reader->GetOutput());
  castFilter->Update();

  const auto width =
      castFilter->GetOutput()->GetLargestPossibleRegion().GetSize()[0]; // width
  const auto height = castFilter->GetOutput()
                          ->GetLargestPossibleRegion()
                          .GetSize()[1]; // height
  const auto sizePix = width * height;
  const auto sizeBuf = sizePix * sizeof(FloatImageType::PixelType);
  const auto bytesPerPix = qRound(sizeBuf / static_cast<double>(sizePix));

  size_t index = 0;
  for (auto &it : this->m_cbctrecon->m_arrYKImage) {
    const auto &strFile = files.at(index);

    if (bytesPerPix != 2) {
      break;
    }

    it.CreateImage(width, height, 0);

    // reader->Read(it.m_pData);
    it.m_pData = &castFilter->GetOutput()->GetBufferPointer()[sizePix * index];

    it.m_strFilePath = strFile;
    // Copy his header - 100 bytes
    it.CopyHisHeader(strFile.toLocal8Bit().constData());

    index++;
  }

  this->m_cbctrecon->m_multiplyFactor = 1.0;

  this->ui.spinBoxImgIdx->setMinimum(0);
  this->ui.spinBoxImgIdx->setMaximum(this->m_cbctrecon->m_iImgCnt - 1);
  this->ui.spinBoxImgIdx->setValue(0);

  this->m_cbctrecon->SetMaxAndMinValueOfProjectionImage();
  SLT_InitializeGraphLim();

  SLT_DrawRawImages(); // Change FileName as well.. read spinbox value and draw
                       // image
}
// IO_END

void CbctReconWidget::SLT_SaveCurrentSetting() const {
  if (!SaveCurrentSetting(this->m_cbctrecon->m_strPathDefaultConfigFile)) {
    std::cout << "Error! in SaveCurrentSetting" << std::endl;
  }
}

bool CbctReconWidget::SaveCurrentSetting(QString &strPathConfigFile) const {
  QFileInfo fInfo(strPathConfigFile);
  if (!fInfo.exists()) {
    std::cout << "Config file not exist. will be created now" << std::endl;
  } else {
    std::cout << "Config file is found. it will be overwritten now"
              << std::endl;
  }
  auto &p_dlgreg_ui = m_dlgRegistration->ui;

  std::ofstream fout;
  fout.open(strPathConfigFile.toLocal8Bit().constData());

  auto strRefmAs = this->ui.lineEdit_RefmAs->text();
  auto PostFOV_R = this->ui.lineEdit_PostFOV_R->text();
  auto PostTablePosY = this->ui.lineEdit_PostTablePosY->text();

  auto strBkFillCT = p_dlgreg_ui.lineEditBkFillCT->text();
  auto strBkDetectCT = p_dlgreg_ui.lineEditBkDetectCT->text();
  auto strBubFillCT = p_dlgreg_ui.lineEditBubFillCT->text();

  auto strBkFillCBCT = p_dlgreg_ui.lineEditBkFillCBCT->text();
  auto strBkDetectCBCT = p_dlgreg_ui.lineEditBubDetectCBCT->text();
  auto strBubFillCBCT = p_dlgreg_ui.lineEditBubFillCBCT->text();

  auto strCropContourName = p_dlgreg_ui.lineEditCropContourName->text();

  auto strFOVPos = p_dlgreg_ui.lineEditFOVPos->text();

  auto strArgument1 = p_dlgreg_ui.lineEditArgument1->text();
  auto strArgument2 = p_dlgreg_ui.lineEditArgument2->text();
  auto strArgument3 = p_dlgreg_ui.lineEditArgument3->text();

  const auto bExportFwd = this->ui.checkBox_ExportFwd->isChecked();
  const auto bExportScat = this->ui.checkBox_ExportScat->isChecked();
  const auto bExportCor = this->ui.checkBox_ExportCor->isChecked();
  const auto bExportVolDICOM = this->ui.checkBox_ExportVolDICOM->isChecked();
  const auto bCouchShiftAddToMacro =
      this->ui.checkBox_CouchShiftAddToMacro->isChecked();

  // From Registration GUI
  const auto bCropBkgroundCT = p_dlgreg_ui.checkBoxCropBkgroundCT->isChecked();
  const auto bCropBkgroundCBCT =
      p_dlgreg_ui.checkBoxCropBkgroundCBCT->isChecked();

  const auto bFillBubbleCT = p_dlgreg_ui.checkBoxFillBubbleCT->isChecked();
  const auto bFillBubbleCBCT = p_dlgreg_ui.checkBoxFillBubbleCBCT->isChecked();

  const auto bUseROIForRigid = p_dlgreg_ui.checkBoxUseROIForRigid->isChecked();
  const auto bUseROIForDIR = p_dlgreg_ui.checkBoxUseROIForDIR->isChecked();

  const auto bRadioButton_mse = p_dlgreg_ui.radioButton_mse->isChecked();
  const auto bRadioButton_mi = p_dlgreg_ui.radioButton_mi->isChecked();

  fout << "strRefmAs"
       << "\t" << strRefmAs.toLocal8Bit().constData() << "\n";
  fout << "PostFOV_R"
       << "\t" << PostFOV_R.toLocal8Bit().constData() << "\n";
  fout << "PostTablePosY"
       << "\t" << PostTablePosY.toLocal8Bit().constData() << "\n";

  fout << "strBkFillCT"
       << "\t" << strBkFillCT.toLocal8Bit().constData() << "\n";
  fout << "strBkDetectCBCT"
       << "\t" << strBkDetectCBCT.toLocal8Bit().constData() << "\n";
  fout << "strBubFillCBCT"
       << "\t" << strBubFillCBCT.toLocal8Bit().constData() << "\n";

  fout << "strBkFillCBCT"
       << "\t" << strBkFillCBCT.toLocal8Bit().constData() << "\n";
  fout << "strBkDetectCBCT"
       << "\t" << strBkDetectCBCT.toLocal8Bit().constData() << "\n";
  fout << "strBubFillCBCT"
       << "\t" << strBubFillCBCT.toLocal8Bit().constData() << "\n";

  fout << "strCropContourName"
       << "\t" << strCropContourName.toLocal8Bit().constData() << "\n";

  fout << "strFOVPos"
       << "\t" << strFOVPos.toLocal8Bit().constData() << "\n";
  fout << "strArgument1"
       << "\t" << strArgument1.toLocal8Bit().constData() << "\n";
  fout << "strArgument2"
       << "\t" << strArgument2.toLocal8Bit().constData() << "\n";
  fout << "strArgument3"
       << "\t" << strArgument3.toLocal8Bit().constData() << "\n";

  fout << "bExportFwd"
       << "\t" << bExportFwd << "\n";
  fout << "bExportScat"
       << "\t" << bExportScat << "\n";
  fout << "bExportCor"
       << "\t" << bExportCor << "\n";
  fout << "bExportVolDICOM"
       << "\t" << bExportVolDICOM << "\n";
  fout << "bCouchShiftAddToMacro"
       << "\t" << bCouchShiftAddToMacro << "\n";

  fout << "bCropBkgroundCT"
       << "\t" << bCropBkgroundCT << "\n";
  fout << "bCropBkgroundCBCT"
       << "\t" << bCropBkgroundCBCT << "\n";

  fout << "bFillBubbleCT"
       << "\t" << bFillBubbleCT << "\n";
  fout << "bFillBubbleCBCT"
       << "\t" << bFillBubbleCBCT << "\n";

  fout << "bUseROIForRigid"
       << "\t" << bUseROIForRigid << "\n";
  fout << "bUseROIForDIR"
       << "\t" << bUseROIForDIR << "\n";

  fout << "radioButton_mse"
       << "\t" << bRadioButton_mse << "\n";
  fout << "radioButton_mi"
       << "\t" << bRadioButton_mi << "\n";
  fout << std::flush;
  fout.close();
  return true;
}

bool CbctReconWidget::LoadCurrentSetting(QString &strPathConfigFile) const {
  auto fInfo = QFileInfo(strPathConfigFile);

  if (!fInfo.exists()) {
    //  std::cout << "Config file doesn't exist" << std::endl;
    return false;
  }

  std::ifstream fin;
  char str[MAX_LINE_LENGTH];

  fin.open(strPathConfigFile.toLocal8Bit().constData());

  if (fin.fail()) {
    return false;
  }
  auto &p_dlgreg_ui = m_dlgRegistration->ui;

  // int cnt = 0;
  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH); // Read out header
    auto tmpStr = QString(&str[0]);
    auto strList = tmpStr.split("\t"); // tab

    QString strHeader, strContent;
    if (strList.count() == 2) {
      strHeader = strList.at(0);
      strContent = strList.at(1);

      const auto bFlagContent = strContent.toInt() == 1;

      if (strHeader == "strRefmAs") {
        this->ui.lineEdit_RefmAs->setText(strContent);
      } else if (strHeader == "PostFOV_R") {
        this->ui.lineEdit_PostFOV_R->setText(strContent);
      } else if (strHeader == "PostTablePosY") {
        this->ui.lineEdit_PostTablePosY->setText(strContent);

      } else if (strHeader == "strBkFillCT") {
        p_dlgreg_ui.lineEditBkFillCT->setText(strContent);
      } else if (strHeader == "strBkDetectCT") {
        p_dlgreg_ui.lineEditBkDetectCT->setText(strContent);
      } else if (strHeader == "strBubFillCT") {
        p_dlgreg_ui.lineEditBubFillCT->setText(strContent);

      } else if (strHeader == "strBkFillCBCT") {
        p_dlgreg_ui.lineEditBkFillCBCT->setText(strContent);
      } else if (strHeader == "strBkDetectCBCT") {
        p_dlgreg_ui.lineEditBubDetectCBCT->setText(strContent);
      } else if (strHeader == "strBubFillCBCT") {
        p_dlgreg_ui.lineEditBubFillCBCT->setText(strContent);
      }

      if (strHeader == "strCropContourName") {
        p_dlgreg_ui.lineEditCropContourName->setText(strContent);
      } else if (strHeader == "strFOVPos") {
        p_dlgreg_ui.lineEditFOVPos->setText(strContent);

      } else if (strHeader == "strArgument1") {
        p_dlgreg_ui.lineEditArgument1->setText(strContent);
      } else if (strHeader == "strArgument2") {
        p_dlgreg_ui.lineEditArgument2->setText(strContent);
      } else if (strHeader == "strArgument3") {
        p_dlgreg_ui.lineEditArgument3->setText(strContent);

      } else if (strHeader == "bExportFwd") {
        this->ui.checkBox_ExportFwd->setChecked(bFlagContent);
      } else if (strHeader == "bExportScat") {
        this->ui.checkBox_ExportScat->setChecked(bFlagContent);
      } else if (strHeader == "bExportCor") {
        this->ui.checkBox_ExportCor->setChecked(bFlagContent);

      } else if (strHeader == "bExportVolDICOM") {
        this->ui.checkBox_ExportVolDICOM->setChecked(bFlagContent);
      } else if (strHeader == "bCouchShiftAddToMacro") {
        this->ui.checkBox_CouchShiftAddToMacro->setChecked(bFlagContent);

      } else if (strHeader == "bCropBkgroundCT") {
        p_dlgreg_ui.checkBoxCropBkgroundCT->setChecked(bFlagContent);
      } else if (strHeader == "bCropBkgroundCBCT") {
        p_dlgreg_ui.checkBoxCropBkgroundCBCT->setChecked(bFlagContent);

      } else if (strHeader == "bFillBubbleCT") {
        p_dlgreg_ui.checkBoxFillBubbleCT->setChecked(bFlagContent);
      } else if (strHeader == "bFillBubbleCBCT") {
        p_dlgreg_ui.checkBoxFillBubbleCBCT->setChecked(bFlagContent);

      } else if (strHeader == "bUseROIForRigid") {
        p_dlgreg_ui.checkBoxUseROIForRigid->setChecked(bFlagContent);
      } else if (strHeader == "bUseROIForDIR") {
        p_dlgreg_ui.checkBoxUseROIForDIR->setChecked(bFlagContent);

      } else if (strHeader == "radioButton_mse") {
        p_dlgreg_ui.radioButton_mse->setChecked(bFlagContent);
      } else if (strHeader == "radioButton_mi") {
        p_dlgreg_ui.radioButton_mi->setChecked(bFlagContent);
      }
    }
  }
  fin.close();

  return true;
}
