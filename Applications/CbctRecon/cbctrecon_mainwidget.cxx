#include "cbctrecon_mainwidget.h"

// Std
#include <iostream>
#include <thread>
#include <vector>

// Qt
#include <qclipboard.h>
#include <qfiledialog.h>
#include <qinputdialog.h>
#include <qmessagebox.h>
#include <qstandarditemmodel.h>
#include <qstring.h>

// ITK
#include "itkCastImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageDuplicator.h"
#include "itkImageSliceIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkMultiplyImageFilter.h"

// PLM
#include "mha_io.h"
#include "nki_io.h"

// Local
#include "cbctrecon.h"
#include "cbctregistration.h"
#include "DlgExternalCommand.h"
#include "DlgRegistration.h"
#include "qcustomplot.h"

#pragma GCC poison new

CbctReconWidget::CbctReconWidget(QWidget *parent, Qt::WindowFlags flags)
    : QMainWindow(parent, flags) {

  ui.setupUi(this);

  // Disable cuda & opencl as defaults
  ui.radioButton_UseCPU->setChecked(true);
  ui.radioButton_UseCUDA->setDisabled(true);
  ui.radioButton_UseOpenCL->setDisabled(true);

  ui.pushButton_DoRecon->setDisabled(true);

  connect(ui.labelImageRaw, SIGNAL(Mouse_Move()), this,
          SLOT(SLT_DataProbeProj()));
  connect(ui.labelImageRaw, SIGNAL(Mouse_Pressed_Left()), this,
          SLOT(SLT_CalculateROI_Proj())); // added

  connect(ui.labelReconImage, SIGNAL(Mouse_Move()), this,
          SLOT(SLT_DataProbeRecon()));
  connect(ui.labelReconImage, SIGNAL(Mouse_Pressed_Left()), this,
          SLOT(SLT_CalculateROI_Recon()));

  // Mouse_Left

  // connect(ui.MyLabel, SIGNAL(Mouse_Pos()), this, SLOT(MouseCurrPos()));

#if USE_CUDA
  ui.radioButton_UseCUDA->setDisabled(false);
  ui.radioButton_UseCUDA->setChecked(true);
#endif

#if USE_OPENCL_PLM
  ui.radioButton_UseOpenCL->setDisabled(false);
#endif

  m_cbctrecon = std::make_unique<CbctRecon>();
  m_dlgRegistration = std::make_unique<DlgRegistration>(this);
  m_cbctregistration = m_dlgRegistration->m_cbctregistration.get();

  QString tmp_folder("tmp");
  init_DlgRegistration(tmp_folder); // to Setup plastimatch folder. this is
                                    // useful if registration will be only done
  // m_pDlgHistogram = new DlgHistogram(this);
  m_dlgExternalCommand = std::make_unique<DlgExternalCommand>(this);

  // 20141017 QTIMER for sync
  m_Timer = std::make_unique<QTimer>(this);
  connect(std::move(m_Timer.get()), SIGNAL(timeout()), this,
          SLOT(SLT_TimerEvent()));
  m_busyTimer = false;

  m_cbctrecon->m_strPathDirDefault =
      R"(D:\Program_data\01_20140827_CBCT_All\04_patient_phan_pelvis_M\IMAGES\img_1.3.46.423632.135786.1409186054.9_M20mAs6440)";

  QString strPathCurAppDir = QDir::currentPath(); // should be same as .exe file
  std::cout << "Current app path= "
            << strPathCurAppDir.toLocal8Bit().constData() << std::endl;
  m_cbctrecon->m_strPathDefaultConfigFile =
      strPathCurAppDir + "/" + "DefaultConfig.cfg";

  if (!LoadCurrentSetting(
          m_cbctrecon->m_strPathDefaultConfigFile)) // Update GUI
  {
    std::cout << "DefaultConfig.cfg is not found in the application folder. A "
                 "new one will be created"
              << std::endl;
    if (!SaveCurrentSetting(m_cbctrecon->m_strPathDefaultConfigFile)) {
      std::cout << "Error in SaveCurrentSetting" << std::endl;
    }
  }
}

void CbctReconWidget::init_DlgRegistration(
  QString &strDCM_UID) // init dlgRegistrations
{
  m_dlgRegistration->initDlgRegistration(
    strDCM_UID); // NULLing all temporary spImage
}

void CbctReconWidget::SLT_LoadRawImages() { LoadRawHisImages(); }

void CbctReconWidget::SLT_DrawRawImages() {
  int crntIdx = ui.spinBoxImgIdx->value();

  if (crntIdx >= m_cbctrecon->m_iImgCnt) {
    return;
  }

  int windowMin = ui.sliderRawMin->value();
  int windowMax = ui.sliderRawMax->value();

  QFileInfo tmpInfo =
      QFileInfo(m_cbctrecon->m_arrYKImage[crntIdx].m_strFilePath);
  ui.lineEditFileName->setText(tmpInfo.fileName());

  int width = m_cbctrecon->m_arrYKImage[crntIdx].m_iWidth;
  int height = m_cbctrecon->m_arrYKImage[crntIdx].m_iHeight;
  m_cbctrecon->m_dspYKImgProj->CreateImage(width, height, 0);
  m_cbctrecon->m_dspYKImgProj->CopyFromBuffer(
      m_cbctrecon->m_arrYKImage[crntIdx].m_pData, width, height);

  m_cbctrecon->m_dspYKImgProj->FillPixMapMinMax(windowMin, windowMax);
  ui.labelImageRaw->SetBaseImage(m_cbctrecon->m_dspYKImgProj.get());
  ui.labelImageRaw->update();
}

void CbctReconWidget::SLT_DrawProjImages() {
  if (m_cbctrecon->m_dspYKImgProj == nullptr) {
    return;
  }

  if (m_cbctrecon->m_iImgCnt > 0) {
    SLT_DrawRawImages();
    //		SLT_DrawGraph();
    SLT_UpdateTable();
    return;
  }

  int iReqSlice = ui.spinBoxImgIdx->value();

  if (!m_cbctrecon->FillProjForDisplay(ui.spinBoxImgIdx->value())) {
    return;
  }

  m_cbctrecon->m_dspYKImgProj->FillPixMapMinMax(ui.sliderRawMin->value(),
                                                ui.sliderRawMax->value());

  ui.labelImageRaw->SetBaseImage(m_cbctrecon->m_dspYKImgProj.get());
  ui.labelImageRaw->update();

  SLT_UpdateTable();
}

void CbctReconWidget::SLT_FileNameHex2Dec() {
  QStringList files = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open",
      m_cbctrecon->m_strPathDirDefault, "projection images (*.his)");

  int cnt = files.size();
  if (cnt <= 0) {
    return;
  }

  QString strMsg = "Original file names will be gone. Backup is strongly "
                   "recommended. Continue?";

  QMessageBox msgBox;
  msgBox.setText(strMsg);
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);

  int res = msgBox.exec();

  if (res == QMessageBox::Yes) {
    m_cbctrecon->RenameFromHexToDecimal(files);
  }
}

void CbctReconWidget::SLT_MakeElektaXML() {
  // Define IMAGE.DBF path
  QString filePath_ImageDBF = QFileDialog::getOpenFileName(
      this, "SelectIMAGE.DBF file", m_cbctrecon->m_strPathDirDefault,
      "Elekta DB file (*.dbf)", nullptr, nullptr);

  if (filePath_ImageDBF.length() < 2) {
    return;
  }

  QString filePath_FrameDBF = QFileDialog::getOpenFileName(
      this, "Select FRAME.DBF file", m_cbctrecon->m_strPathDirDefault,
      "Elekta DB file (*.dbf)", nullptr, nullptr);

  if (filePath_FrameDBF.length() < 2) {
    return;
  }

  QString DICOM_UID;
  QInputDialog inputDlg;

  bool ok;
  QString text = QInputDialog::getText(
      this, "Input Dialog", "DICOM_UID:", QLineEdit::Normal, "DICOM_UID", &ok);

  if (ok && !text.isEmpty()) {
    DICOM_UID = text;
  }

  if (DICOM_UID.length() < 2) {
    return;
  }

  QString genFilePath = m_cbctrecon->MakeElektaXML(
      filePath_ImageDBF, filePath_FrameDBF, DICOM_UID);
  std::cout << "Generated ElektaXML path: "
            << genFilePath.toLocal8Bit().constData() << std::endl;
}

void CbctReconWidget::SLT_OpenOffsetFile() {
  // QString strPath = QFileDialog::getOpenFileNames(this,"Select one or more
  // files to open","/home","Images (*.raw)");
  QString strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open", m_cbctrecon->m_strPathDirDefault,
      "raw image (*.raw)");

  if (strPath.length() <= 1) {
    return;
  }

  m_cbctrecon->LoadCalibData(strPath.toLocal8Bit().constData(), OFFSET_CALIB);

  ui.lineEdit_offsetPath->setText(strPath);
}

void CbctReconWidget::SLT_OpenGainFile() {
  QString strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open", m_cbctrecon->m_strPathDirDefault,
      "raw image (*.raw)");

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_gainPath->setText(strPath);

  m_cbctrecon->LoadCalibData(strPath.toLocal8Bit().constData(), GAIN_CALIB);
}

void CbctReconWidget::SLT_OpenBadpixelFile() {
  QString strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open", m_cbctrecon->m_strPathDirDefault,
      "bad pixel map file (*.pmf)");

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_badpixelPath->setText(strPath);
  m_cbctrecon->LoadCalibData(strPath.toLocal8Bit().constData(), BADPIXEL_CALIB);
  // m_pImgGain->LoadRawImage(strPath.toLocal8Bit(),IMG_WIDTH, IMG_HEIGHT);
}

void CbctReconWidget::SLT_ApplyCalibration() {
  if (m_cbctrecon->m_iImgCnt < 1) {
    return;
  }

  bool bDarkCorrApply = ui.checkBox_offsetOn->isChecked();
  bool bGainCorrApply = ui.checkBox_gainOn->isChecked();
  bool bDefectMapApply = ui.checkBox_badpixelOn->isChecked();
  for (int i = 0; i < m_cbctrecon->m_iImgCnt; i++) {
    m_cbctrecon->CorrectSingleFile(
        &m_cbctrecon->m_arrYKImage[i], bDarkCorrApply, bGainCorrApply,
        bDefectMapApply); // pixel value will be changed
  }
  SLT_DrawRawImages();
}

void CbctReconWidget::SLT_DrawReconImage() {
  if (m_cbctrecon->m_dspYKReconImage == nullptr) {
    return;
  }

  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    std::cout << "no recon image to be displayed" << std::endl;
    return;
  }

  //  The ExtractImageFilter type is instantiated using the input and
  //  output image types. A filter object is created with the New()
  //  method and assigned to a SmartPointer.

  using ExtractFilterType =
      itk::ExtractImageFilter<UShortImageType, UShortImage2DType>;
  ExtractFilterType::Pointer extractFilter = ExtractFilterType::New();

  using DuplicatorType = itk::ImageDuplicator<UShortImageType>;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(m_cbctrecon->m_spCrntReconImg);
  duplicator->Update();
  UShortImageType::Pointer clonedImage = duplicator->GetOutput();

  extractFilter->SetDirectionCollapseToSubmatrix();

  UShortImageType::RegionType crntRegion3D = clonedImage->GetBufferedRegion();

  crntRegion3D = clonedImage->GetBufferedRegion();

  //  We take the size from the region and collapse the size in the $Z$
  //  component by setting its value to $1$.

  // Get Image Size and Extraction Index info.
  UShortImageType::SizeType size = crntRegion3D.GetSize();
  size[2] = 0; // z size number = 0 --> should not be 1

  UShortImageType::IndexType start = crntRegion3D.GetIndex();
  const int iSliceNumber = ui.spinBoxReconImgSliceNo->value();
  start[2] = iSliceNumber; // 60

  double originZ = m_cbctrecon->m_spCrntReconImg->GetOrigin()[2];
  double spacingZ = m_cbctrecon->m_spCrntReconImg->GetSpacing()[2];
  double posZ = originZ + iSliceNumber * spacingZ;

  QString strPosZ = QString("%1").arg(posZ, 0, 'f', 2);
  // strPosZ.sprintf("%4.2f", posZ);
  ui.lineEdit_CurrentPosZ->setText(strPosZ);

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
  m_cbctrecon->m_dspYKReconImage = YK16GrayImage::CopyItkImage2YKImage(
      pCrnt2D,
      std::move(m_cbctrecon->m_dspYKReconImage)); // dimension should be same
                                                  // automatically.

  // m_dspYKReconImage->SaveDataAsRaw("D:\\RawFile.raw"); //410 410 OK

  float physPosX = ui.lineEdit_PostFOV_X->text().toFloat();
  float physPosY = ui.lineEdit_PostFOV_Y->text().toFloat();
  float physRadius = ui.lineEdit_PostFOV_R->text().toFloat();
  float physTablePosY = ui.lineEdit_PostTablePosY->text().toFloat();
  m_cbctrecon->PostApplyFOVDispParam(physPosX, physPosY, physRadius,
                                     physTablePosY);
  // SLT_UpdatePostProcDispObj();

  if (ui.checkBox_PostDispObjOn->isChecked()) {
    m_cbctrecon->m_dspYKReconImage->m_bDrawFOVCircle = true;
    m_cbctrecon->m_dspYKReconImage->m_bDrawTableLine = true;
  }

  else {
    m_cbctrecon->m_dspYKReconImage->m_bDrawFOVCircle = false;
    m_cbctrecon->m_dspYKReconImage->m_bDrawTableLine = false;
  }

  m_cbctrecon->m_dspYKReconImage->FillPixMapMinMax(
      ui.sliderReconImgMin->value(), ui.sliderReconImgMax->value());
  ui.labelReconImage->SetBaseImage(m_cbctrecon->m_dspYKReconImage.get());
  ui.labelReconImage->update();

  // SLT_DrawGraph();
  SLT_UpdateTable();
}

void CbctReconWidget::SLT_OpenElektaGeomFile() {
  QString strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open", m_cbctrecon->m_strPathDirDefault,
      "Geometry file (*.xml)");

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_ElektaGeomPath->setText(strPath);
}

void CbctReconWidget::SLT_SetOutputPath() {
  QString strPath = QFileDialog::getSaveFileName(
      this, "File path to save", "D:\\", "meta 3D image data (*.mha)", nullptr,
      nullptr); // Filename don't need to exist

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_OutputFilePath->setText(strPath);
}

void CbctReconWidget::SLT_DoReconstruction() {
  if (ui.radioButton_UseCUDA->isChecked()) {
    m_cbctrecon->DoReconstructionFDK<CUDA_DEVT>(REGISTER_RAW_CBCT);
  } else if (ui.radioButton_UseOpenCL->isChecked()) {
    m_cbctrecon->OpenCLDoReconstructionFDK(REGISTER_RAW_CBCT);
  } else {
    m_cbctrecon->DoReconstructionFDK<CPU_DEVT>(REGISTER_RAW_CBCT);
  }

  QString update_text("RAW_CBCT");
  UpdateReconImage(m_cbctrecon->m_spCrntReconImg, update_text);

  m_dlgRegistration->UpdateListOfComboBox(0); // combo selection
                                                            // signalis called
  m_dlgRegistration->UpdateListOfComboBox(1);
  // m_pDlgRegistration->SelectComboExternal(0, REGISTER_RAW_CBCT); // will call
  // fixedImageSelected  m_pDlgRegistration->SelectComboExternal(1,
  // REGISTER_RAW_CBCT );

  // After first reconstruction, set Median size to 0 0 1 for scatter corrected
  // solution
  /* ui.lineEdit_PostMedSizeX->setText(QString("%1").arg(0.0));
  ui.lineEdit_PostMedSizeY->setText(QString("%1").arg(0.0));
  ui.lineEdit_PostMedSizeZ->setText(QString("%1").arg(1.0));*/
}

void CbctReconWidget::SLT_InitializeGraphLim() {
  // Set Max Min at graph
  if (ui.radioButton_graph_proj->isChecked()) {
    if (m_cbctrecon->m_iImgCnt > 0) // if indep raw his images are loaded
    {
      int horLen = m_cbctrecon->m_dspYKImgProj->m_iWidth;
      // int verLen = m_dspYKImgProj->m_iHeight;

      // set edit maxium min
      QString strXMin = QString("%1").arg(horLen);
      ui.lineEditXMin->setText("0");
      ui.lineEditXMax->setText(strXMin);

      QString strYMin, strYMax;
      strYMin = QString("%1").arg(m_cbctrecon->m_fProjImgValueMin, 0, 'f', 1);
      strYMax = QString("%1").arg(m_cbctrecon->m_fProjImgValueMax, 0, 'f', 1);

      ui.lineEditYMin->setText(strYMin);
      ui.lineEditYMax->setText(strYMax);
    }

    if (m_cbctrecon->m_spProjImg3DFloat == nullptr) {
      return;
    }

    int horLen =
        m_cbctrecon->m_spProjImg3DFloat->GetBufferedRegion().GetSize()[0];
    // int verLen = m_spProjImg3DFloat->GetBufferedRegion().GetSize()[1];

    // set edit maxium min
    QString strXMin = QString("%1").arg(horLen);
    ui.lineEditXMin->setText("0");
    ui.lineEditXMax->setText(strXMin);

    QString strYMin, strYMax;
    strYMin = QString("%1").arg(m_cbctrecon->m_fProjImgValueMin, 0, 'f', 1);
    strYMax = QString("%1").arg(m_cbctrecon->m_fProjImgValueMax, 0, 'f', 1);

    ui.lineEditYMin->setText(strYMin);
    ui.lineEditYMax->setText(strYMax);
  } else if (ui.radioButton_graph_recon->isChecked()) {
    if (m_cbctrecon->m_spCrntReconImg == nullptr) {
      return;
    }

    int horLen =
        m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize()[0];
    // int verLen = m_spCrntReconImg->GetBufferedRegion().GetSize()[1];

    // set edit maxium min

    QString strXMax = QString("%1").arg(horLen);
    ui.lineEditXMin->setText("0");
    ui.lineEditXMax->setText(strXMax);

    QString strYMin, strYMax;
    strYMin = QString("%1").arg(0.0, 0, 'f', 1);
    strYMax = QString("%1").arg(2000.0, 0, 'f', 1);

    ui.lineEditYMin->setText(strYMin);
    ui.lineEditYMax->setText(strYMax);
  }
}

void CbctReconWidget::SLT_CopyTableToClipBoard() {
  qApp->clipboard()->clear();

  QStringList list;

  int rowCnt = m_cbctrecon->m_pTableModel->rowCount();
  int columnCnt = m_cbctrecon->m_pTableModel->columnCount();

  list << "\n";
  // for (int i = 0 ; i < columnCnt ; i++)
  //{
  QFileInfo tmpInfo = QFileInfo(ui.lineEdit_Cur3DFileName->text());
  // list << "Index";
  list << tmpInfo.baseName();
  list << "\n";

  list << "Pos(mm)";
  list << "Intensity";
  list << "\n";

  for (int j = 0; j < rowCnt; j++) {
    for (int i = 0; i < columnCnt; i++) {
      QStandardItem *item = m_cbctrecon->m_pTableModel->item(j, i);
      list << item->text();
    }
    list << "\n";
  }

  qApp->clipboard()->setText(list.join("\t"));
}

void CbctReconWidget::SLT_SetHisDir() // Initialize all image buffer
{
  // Initializing..

  // Set folder --> then use RTK HIS Reader
  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.length() <= 1) {
    return;
  }

  ui.lineEdit_HisDirPath->setText(dirPath);

  m_cbctrecon->SetProjDir(dirPath);
  init_DlgRegistration(m_cbctrecon->m_strDCMUID);

  ui.lineEdit_ElektaGeomPath->setText(m_cbctrecon->m_strPathGeomXML);
  ui.lineEdit_PathCBCTSkinPath->setText(m_cbctrecon->m_strPathRS_CBCT);

  float kVp = 0.0;
  float mA = 0.0;
  float ms = 0.0;
  m_cbctrecon->GetXrayParamFromINI(m_cbctrecon->m_strPathElektaINI, kVp, mA,
                                   ms);

  if (kVp * mA * ms != 0) {
    // update GUI
    std::cout << "Updating current mAs setting from INI file: "
              << "kVp= " << kVp << ", mA= " << mA << ", ms= " << ms
              << std::endl;
  }
  ui.lineEdit_CurmAs->setText(QString("%1, %2").arg(mA).arg(ms));

  VEC3D couch_trans = {-999, -999,
                       -999}; // mm. In the text file, these values are in cm.
  VEC3D couch_rot = {-999, -999,
                     -999}; // mm. In the text file, these values are in cm.

  bool res = m_cbctrecon->GetCouchShiftFromINIXVI(
      m_cbctrecon->m_strPathElektaINIXVI2, &couch_trans, &couch_rot);

  if (res) {
    QString strTransX = QString::number(couch_trans.x, 'f', 1);
    QString strTransY = QString::number(couch_trans.y, 'f', 1);
    QString strTransZ = QString::number(couch_trans.z, 'f', 1);
    QString strTransAll = strTransX + "," + strTransY + "," + strTransZ;

    QString strRotX = QString::number(couch_rot.x, 'f', 1);
    QString strRotY = QString::number(couch_rot.y, 'f', 1);
    QString strRotZ = QString::number(couch_rot.z, 'f', 1);

    QString strRotAll = strRotX + "," + strRotY + "," + strRotZ;

    ui.lineEdit_CouchTrans->setText(strTransAll);
    ui.lineEdit_CouchRot->setText(strRotAll);
  } else {
    ui.lineEdit_CouchTrans->setText("Not available");
    ui.lineEdit_CouchRot->setText("Not available");
  }

  m_cbctrecon->m_vSelectedFileNames.clear();

  std::cout << "Push Load button to load projection images" << std::endl;
}

QString getBowtiePath(QWidget *parent, const QDir &calDir) {
  return QFileDialog::getOpenFileName(
      parent, "Find air(+bowtie) filter image for subtraction",
      calDir.absolutePath(), "Projection (*.xim)", nullptr, nullptr);
}

std::tuple<bool, bool> CbctReconWidget::probeUser(const QString &guessDir) {

  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open CT DICOM Directory"), guessDir,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  bool dcm_success = false;
  if (!(dirPath.length() <= 1)) {

    if (m_cbctrecon->ReadDicomDir(dirPath)) {

      m_dlgRegistration->UpdateVOICombobox(PLAN_CT);
      // UpdateReconImage(m_spRefCTImg, QString("DICOM reference image"));

      m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT,
                                          REGISTER_MANUAL_RIGID);
      dcm_success = true;
    }
  }
  bool instaRecon = false;
  QMessageBox::StandardButton reply;
  reply =
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

  FilterReaderType::Pointer bowtiereader =
      FilterReaderType::New(); // we use is because we need the projections to
                               // be in the same unit (order of magnitude)

  QDir guessDir(proj_path + QString("/../"));

  QDir calDir(proj_path + QString("/Calibrations/"));

  QString bowtiePath;

  switch (m_cbctrecon->m_projFormat) {
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
  ui.pushButton_DoRecon->setDisabled(true);
  // 1) Get all projection file names
  QString dirPath = ui.lineEdit_HisDirPath->text();
  //.toLocal8Bit().constData();

  if (!QFile::exists(dirPath)) {
    std::cout << "Projection file directory was not found. Retry." << std::endl;
    return;
  }

  auto names = m_cbctrecon->GetProjFileNames(dirPath);

  if (!m_cbctrecon->IsFileNameOrderCorrect(names) &&
      (m_cbctrecon->m_projFormat != XIM_FORMAT)) {
    std::cout << "Check the file name order" << std::endl;
    QMessageBox::warning(this, "warning", "Error on File Name Sorting!");
    return;
  }

  std::cout << "File name order was cross-checked and found to be OK!"
            << std::endl;

  int fullCnt = names.size();
  if (fullCnt <= 0) {
    std::cout << "No projection file was found. Retry." << std::endl;
    return;
  }

  std::cout << fullCnt << "  projection files were found." << std::endl;

  // 2) Elekta Geometry file
  QString geomPath = ui.lineEdit_ElektaGeomPath->text();
  QFileInfo geomFileInfo(geomPath);
  //! QFile::exists(geomPath)

  if (!m_cbctrecon->LoadGeometry(geomFileInfo, names)) {
    if (!m_cbctrecon->m_strError.isEmpty()) {
      QMessageBox::critical(this, "LoadXVIGeometryFile",
                            m_cbctrecon->m_strError, QMessageBox::Ok);
    }
  }

  int iFullGeoDataSize =
      m_cbctrecon->m_spFullGeometry->GetGantryAngles().size();
  if (iFullGeoDataSize < 1) {
    std::cout << "Not enough projection image (should be > 0)" << std::endl;
    return;
  }

  if (iFullGeoDataSize != fullCnt) {
    if (m_cbctrecon->m_projFormat != XIM_FORMAT) {
      std::cout << "Size of geometry data and file numbers are not same! Check "
                   "and retry"
                << std::endl;
      return;
    }

    QMessageBox::StandardButton reply;
    reply = QMessageBox::question(
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

  auto angle_gaps = m_cbctrecon->m_spFullGeometry->GetAngularGaps(
      m_cbctrecon->m_spFullGeometry->GetSourceAngles());

  auto sum_gap =
      std::accumulate(std::begin(angle_gaps), std::end(angle_gaps), 0.0);
  sum_gap /= itk::Math::pi * 180.0;
  auto mean_gap = sum_gap / angle_gaps.size();

  std::cout << "AngularGaps Sum (deg):" << sum_gap
            << ", Mean (deg): " << mean_gap << std::endl;

  double gantryAngleInterval =
      ui.lineEdit_ManualProjAngleGap->text().toDouble();

  // if (ui.Radio_KeepOriginalAngles->isChecked())
  if (ui.Radio_ManualProjAngleGap->isChecked()) {
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

  auto exclude_ids = m_cbctrecon->GetExcludeProjFiles(
      ui.Radio_ManualProjAngleGap->isChecked(), gantryAngleInterval);

  m_cbctrecon->LoadSelectedProj(exclude_ids, names);

  // Reads the cone beam projections
  using ReaderType = rtk::ProjectionsReader<FloatImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames(m_cbctrecon->m_vSelectedFileNames);
  // TRY_AND_EXIT_ON_ITK_EXCEPTION(
  // std::thread calc_thread(read_projections, reader);
  std::thread calc_thread([&reader]() { reader->Update(); });
  // calc_thread.detach();

  std::cout << "Reader detached from main thread" << std::endl;

  // After reading the whole file,
  // HIS header should be saved
  m_cbctrecon->saveHisHeader();

  //  Insta Recon, Dcm read
  std::tuple<bool, bool> answers;
  auto &bowtie_reader =
      ReadBowtieFileWhileProbing(geomFileInfo.absolutePath(), answers);

  calc_thread.join();
  std::cout << "Reader re-attached to main thread" << std::endl;

  if (bowtie_reader != nullptr) {
    m_cbctrecon->ApplyBowtie(reader, bowtie_reader);
  }
  if (m_cbctrecon->m_projFormat == HND_FORMAT) {
    std::cout << "Fitted bowtie-filter correction ongoing..." << std::endl;
    SLT_DoBowtieCorrection();
  }

  std::thread save_thread(saveImageAsMHA<FloatImageType>,
                          m_cbctrecon->m_spProjImg3DFloat);
  // Make sure the projections are saved before going out of scope.
  save_thread.join();

  if (!m_cbctrecon->ResampleProjections(
          ui.lineEdit_DownResolFactor->text().toDouble())) { // 0.5
    // reset factor if image was not resampled
    ui.lineEdit_DownResolFactor->setText("1.0");
  }

  m_cbctrecon->ConvertLineInt2Intensity(
      m_cbctrecon->m_spProjImg3DFloat, m_cbctrecon->m_spProjImgRaw3D,
      65535); // if X not 1024 == input size: out_offset =
              // in_offset + (1024*res_f -
              // X*res_f)*out_spacing     <- will still
              // break down at fw_projection

  ui.pushButton_DoRecon->setEnabled(true);

  ui.spinBoxImgIdx->setMinimum(0);
  ui.spinBoxImgIdx->setMaximum(m_cbctrecon->m_vSelectedFileNames.size() - 1);
  ui.spinBoxImgIdx->setValue(0); // it doesn't call Draw Event .. don't know
                                 // why.

  m_cbctrecon
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

void CbctReconWidget::SLT_DataProbeProj() {
  double dspWidth = ui.labelImageRaw->width();
  double dspHeight = ui.labelImageRaw->height();
  int dataWidth = 0;
  int dataHeight = 0;
  int dataX = 0;
  int dataY = 0;
  int dataZ = 0;
  double fProbeValue = 0.0;

  if (m_cbctrecon->m_iImgCnt > 0) // there is indep loaded projection files
  {
    dataWidth = m_cbctrecon->m_dspYKImgProj->m_iWidth;
    dataHeight = m_cbctrecon->m_dspYKImgProj->m_iHeight;

    dataX = qRound(ui.labelImageRaw->x / dspWidth * dataWidth);
    dataY = qRound(ui.labelImageRaw->y / dspHeight * dataHeight);
    dataZ = ui.spinBoxImgIdx->value();
    fProbeValue =
        m_cbctrecon->m_dspYKImgProj->m_pData[dataWidth * dataY + dataX];
  } else {
    if (m_cbctrecon->m_spProjImg3DFloat == nullptr) {
      return;
    }

    dataWidth = static_cast<int>(
        m_cbctrecon->m_spProjImg3DFloat->GetBufferedRegion().GetSize()[0]);
    dataHeight = static_cast<int>(
        m_cbctrecon->m_spProjImg3DFloat->GetBufferedRegion().GetSize()[1]);

    // int crntIdx = ui.spinBoxImgIdx->value();
    // These are displayed data (just index data)
    dataX = qRound(ui.labelImageRaw->x / dspWidth * dataWidth);
    dataY = qRound(ui.labelImageRaw->y / dspHeight * dataHeight);
    dataZ = ui.spinBoxImgIdx->value();

    // fProbeValue = m_dspYKImgProj->m_pData[dataWidth*dataY +
    // dataX]/m_multiplyFactor;
    fProbeValue =
        m_cbctrecon->m_dspYKImgProj->m_pData[dataWidth * dataY + dataX] /
            m_cbctrecon->m_multiplyFactor +
        m_cbctrecon->m_fProjImgValueMin;
  }

  QString dspText;
  dspText = QString("(%1, %2, %3): %4")
                .arg(dataX)
                .arg(dataY)
                .arg(dataZ)
                .arg(fProbeValue, 0, 'f', 2);
  ui.lineEdit_DataProbe_Proj->setText(dspText);
}

void CbctReconWidget::SLT_DataProbeRecon() {
  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  double dspWidth = ui.labelReconImage->width();
  double dspHeight = ui.labelReconImage->height();

  auto dataWidth =
      m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize()[0];
  auto dataHeight =
      m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize()[1];

  // int crntIdx = ui.spinBoxImgIdx->value();
  // These are displayed data (just index data)

  auto dataX = std::lround(ui.labelReconImage->x / dspWidth * dataWidth);
  auto dataY = std::lround(ui.labelReconImage->y / dspHeight * dataHeight);
  auto dataZ = ui.spinBoxReconImgSliceNo->value();

  auto iProbeValue =
      m_cbctrecon->m_dspYKReconImage->m_pData[dataWidth * dataY + dataX];
  // unsigned short iProbeValue = GetValueFrom3DImageUshort(dataX, dataY, dataZ,
  // m_spReconImg);

  QChar zero('0');
  QString dspText;
  dspText = QString("(%1, %2, %3): %4")
                .arg(dataX, 3, 10, zero)
                .arg(dataY, 3, 10, zero)
                .arg(dataZ, 3, 10, zero)
                .arg(iProbeValue);
  // dspText.sprintf("(%03d, %03d, %03d): %d", dataX, dataY, dataZ,
  // iProbeValue);
  ui.lineEdit_DataProbe_Recon->setText(dspText);
}

void CbctReconWidget::SLT_DrawGraph() // based on profile
{
  if (m_cbctrecon->m_pTableModel == nullptr) {
    return;
  }

  // Draw only horizontal, center

  QVector<double> vAxisX; // can be rows or columns
  QVector<double> vAxisY;

  // QStandardItemModel 	m_pTableModel.item()
  int dataLen = m_cbctrecon->m_pTableModel->rowCount();

  if (dataLen < 1) {
    return;
  }

  // std::cout << "check graph 1" << std::endl;
  ui.customPlot->clearGraphs();

  double minX = 9999.0;
  double maxX = -1.0;

  for (int i = 0; i < dataLen; i++) {
    QStandardItem *tableItem1 = m_cbctrecon->m_pTableModel->item(i, 0);
    QStandardItem *tableItem2 = m_cbctrecon->m_pTableModel->item(i, 1);
    double tableVal1 = tableItem1->text().toDouble();
    double tableVal2 = tableItem2->text().toDouble();

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

  ui.customPlot->addGraph();
  ui.customPlot->graph(0)->setData(vAxisX, vAxisY);
  ui.customPlot->graph(0)->setPen(QPen(Qt::blue));
  ui.customPlot->graph(0)->setName("Image profile");

  ui.lineEditXMin->setText(QString("%1").arg(minX));
  ui.lineEditXMax->setText(QString("%1").arg(maxX));

  double tmpXMin = ui.lineEditXMin->text().toDouble();
  double tmpXMax = ui.lineEditXMax->text().toDouble();
  double tmpYMin = ui.lineEditYMin->text().toDouble();
  double tmpYMax = ui.lineEditYMax->text().toDouble();

  // std::cout << "check graph 3" << std::endl;

  ui.customPlot->xAxis->setRange(tmpXMin, tmpXMax);
  ui.customPlot->yAxis->setRange(tmpYMin, tmpYMax);

  ui.customPlot->xAxis->setLabel("mm");
  ui.customPlot->yAxis->setLabel("Intensity");
  ui.customPlot->setTitle("Image Profile");
  QFont titleFont = font();
  titleFont.setPointSize(10);
  ui.customPlot->setTitleFont(titleFont);

  // std::cout << "check graph 4" << std::endl;

  ui.customPlot->legend->setVisible(false);
  QFont legendFont = font();  // start out with MainWindow's font..
  legendFont.setPointSize(9); // and make a bit smaller for legend
  ui.customPlot->legend->setFont(legendFont);
  ui.customPlot->legend->setPositionStyle(QCPLegend::psTopRight);
  ui.customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 200)));

  // std::cout << "check graph 5" << std::endl;
  ui.customPlot->replot();

  // SLT_UpdateTable();
}

void CbctReconWidget::SLT_UpdateTable() {
  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    ui.radioButton_graph_proj->setChecked(true);
  }

  // std::cout << "check 1" << std::endl;
  YK16GrayImage *pYKImg = nullptr;
  double fMultiPlyFactor = 1.0;
  double fMinValue = 0.0;

  if (ui.radioButton_graph_proj->isChecked()) {
    pYKImg =
        m_cbctrecon->m_dspYKImgProj.get(); // you may look, but no touching!

    if (m_cbctrecon->m_iImgCnt > 0) { // if indep image
      fMultiPlyFactor = 1.0;
    } else {
      fMultiPlyFactor = m_cbctrecon->m_multiplyFactor;
      fMinValue = m_cbctrecon->m_fProjImgValueMin;
    }
  } else {
    pYKImg = m_cbctrecon->m_dspYKReconImage.get();
    fMultiPlyFactor = 1.0;
    fMinValue = 0.0;
  }
  if (pYKImg == nullptr) {
    return;
  }

  // std::cout << "check 2" << std::endl;

  // std::cout << "check 3" << std::endl;
  int columnSize = 1;
  int rowSize = 0;

  /// int rowSize = pYKImg->m_iWidth;

  if (ui.radioButton_Profile_Hor->isChecked()) {
    columnSize = 2;
    rowSize = pYKImg->m_iWidth;
  } else {
    columnSize = 2;
    rowSize = pYKImg->m_iHeight;
  }

  // std::cout << "check 4" << std::endl;
  m_cbctrecon->m_pTableModel.reset();
  m_cbctrecon->m_pTableModel = std::make_unique<QStandardItemModel>(
      rowSize, columnSize, this); // 2 Rows and 3 Columns

  // for (int i = 0 ; i<columnSize ; i++)
  //{
  // QFileInfo tmpInfo = QFileInfo(m_arrYKImage[i].m_strFilePath);
  // m_pTableModel->setHorizontalHeaderItem(0, new
  // QStandardItem(QString("Index"))); m_pTableModel->setHorizontalHeaderItem(0,
  // new QStandardItem(QString("Profile")));
  auto pos_item = std::make_unique<QStandardItem>(QString("Position(mm)"));
  auto val_item = std::make_unique<QStandardItem>(QString("Value"));

  m_cbctrecon->m_pTableModel->setHorizontalHeaderItem(
      0, std::move(pos_item.get()));
  m_cbctrecon->m_pTableModel->setHorizontalHeaderItem(
      1, std::move(val_item.get()));
  //}

  // std::cout << "check 5" << std::endl;
  // int width = pYKImg->m_iWidth;
  // int height = pYKImg->m_iHeight;
  // int fixedY = qRound(height / 2.0);

  double originX, originY;
  double spacingX, spacingY;
  originX = 0.0;
  originY = 0.0;
  spacingX = 1.0;
  spacingY = 1.0;

  if (!ui.radioButton_graph_proj->isChecked()) {
    if (m_cbctrecon->m_spCrntReconImg != nullptr) {
      UShortImageType::PointType tmpOrigin =
          m_cbctrecon->m_spCrntReconImg->GetOrigin();
      UShortImageType::SpacingType tmpSpacing =
          m_cbctrecon->m_spCrntReconImg->GetSpacing();
      originX = tmpOrigin[0];
      originY = tmpOrigin[1];
      spacingX = tmpSpacing[0];
      spacingY = tmpSpacing[1];
    }
  }

  // std::cout << "check 6" << std::endl;

  QVector<qreal> vPos;
  if (ui.radioButton_Profile_Hor->isChecked()) {
    for (int i = 0; i < rowSize; i++) {
      vPos.push_back(originX + i * spacingX);
    }
  } else {
    for (int i = 0; i < rowSize; i++) {
      vPos.push_back(originY + i * spacingY);
    }
  }

  QVector<qreal> vProfile;
  if (ui.radioButton_Profile_Hor->isChecked()) {
    pYKImg->GetProfileData(vProfile, DIRECTION_HOR);
  } else {
    pYKImg->GetProfileData(vProfile, DIRECTION_VER);
  }

  // int i = fixedY;
  for (int i = 0; i < rowSize; i++) {
    auto tmpVal1 = vPos[i];
    auto xpos_item =
        std::make_unique<QStandardItem>(QString("%1").arg(tmpVal1));
    m_cbctrecon->m_pTableModel->setItem(i, 0, std::move(xpos_item.get()));

    auto tmpVal2 = vProfile[i] / fMultiPlyFactor + fMinValue;
    auto profval_item =
        std::make_unique<QStandardItem>(QString("%1").arg(tmpVal2));
    m_cbctrecon->m_pTableModel->setItem(i, 1, std::move(profval_item.get()));
  }

  ui.tableViewReconImgProfile->setModel(
      m_cbctrecon->m_pTableModel.get()); // also for proj

  // std::cout << "check 7" << std::endl;
  SLT_DrawGraph();
}

// Mouse Left Click
void CbctReconWidget::SLT_CalculateROI_Recon() {
  if (m_cbctrecon->m_dspYKReconImage == nullptr) {
    return;
  }

  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }
  auto &CrntReconImg = m_cbctrecon->m_spCrntReconImg;
  auto &dspReconImg = m_cbctrecon->m_dspYKReconImage;

  double dspWidth = ui.labelReconImage->width();
  double dspHeight = ui.labelReconImage->height();

  int dataWidth = dspReconImg->m_iWidth;
  int dataHeight = dspReconImg->m_iHeight;
  if (dataWidth * dataHeight == 0) {
    return;
  }

  // int crntIdx = ui.spinBoxImgIdx->value();
  // These are displayed data (just index data)

  int dataX = qRound(ui.labelReconImage->x / dspWidth * dataWidth);
  int dataY = qRound(ui.labelReconImage->y / dspHeight * dataHeight);
  int dataZ = ui.spinBoxReconImgSliceNo->value();

  auto originX = static_cast<double>(CrntReconImg->GetOrigin()[0]);
  auto originY = static_cast<double>(CrntReconImg->GetOrigin()[1]);
  auto originZ = static_cast<double>(CrntReconImg->GetOrigin()[2]);

  auto spacingX = static_cast<double>(CrntReconImg->GetSpacing()[0]);
  auto spacingY = static_cast<double>(CrntReconImg->GetSpacing()[1]);
  auto spacingZ = static_cast<double>(CrntReconImg->GetSpacing()[2]);

  double posX = originX + dataX * spacingX;
  double posY = originY + dataY * spacingY;
  double posZ = originZ + dataZ * spacingZ;

  QString tmpStr1, tmpStr2, tmpStr3;
  tmpStr1 = QString("%1").arg(posX, 0, 'f', 2);
  tmpStr2 = QString("%1").arg(posY, 0, 'f', 2);
  tmpStr3 = QString("%1").arg(posZ, 0, 'f', 2);
  ui.lineEdit_ForcedProbePosX->setText(tmpStr1);
  ui.lineEdit_ForcedProbePosY->setText(tmpStr2);
  ui.lineEdit_ForcedProbePosZ->setText(tmpStr3);

  dspReconImg->SetProfileProbePos(dataX, dataY);
  if (ui.radioButton_Profile_Hor->isChecked()) {
    dspReconImg->m_bDrawProfileX = true;
    dspReconImg->m_bDrawProfileY = false;
  } else {
    dspReconImg->m_bDrawProfileX = false;
    dspReconImg->m_bDrawProfileY = true;
  }

  // dspReconImg value itself
  int ROI_size = ui.lineEdit_ROI_size->text().toInt();
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

    QString strMean;
    strMean = QString("%1").arg(dspReconImg->m_fPixelMean_ROI, 0, 'f', 2);
    QString strSD;
    strSD = QString("%1").arg(dspReconImg->m_fPixelSD_ROI, 0, 'f', 2);
    ui.lineEdit_ROI_mean->setText(strMean);
    ui.lineEdit_ROI_SD->setText(strSD);
  } else {
    dspReconImg->DrawROIOn(false);
  }

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_CalculateROI_Proj() {
  if (m_cbctrecon->m_dspYKImgProj == nullptr) {
    return;
  }
  auto &dspProjImg = m_cbctrecon->m_dspYKImgProj;

  double dspWidth = ui.labelImageRaw->width();
  double dspHeight = ui.labelImageRaw->height();

  int dataWidth = dspProjImg->m_iWidth;
  int dataHeight = dspProjImg->m_iHeight;
  if (dataWidth * dataHeight == 0) {
    return;
  }

  // int crntIdx = ui.spinBoxImgIdx->value();
  // These are displayed data (just index data)
  int dataX = qRound(ui.labelImageRaw->x / dspWidth * dataWidth);
  int dataY = qRound(ui.labelImageRaw->y / dspHeight * dataHeight);
  int dataZ = ui.spinBoxImgIdx->value();

  double originX = 0.0; // 0
  double originY = 0.0; // 0
  double originZ = 0.0; // 0

  double spacingX = 1.0; // 1
  double spacingY = 1.0; // 1
  double spacingZ = 1.0; // 1

  double posX = originX + dataX * spacingX;
  double posY = originY + dataY * spacingY;
  double posZ = originZ + dataZ * spacingZ;

  ui.lineEdit_ForcedProbePosX->setText(QString("%1").arg(posX));
  ui.lineEdit_ForcedProbePosY->setText(QString("%1").arg(posY));
  ui.lineEdit_ForcedProbePosZ->setText(QString("%1").arg(posZ));

  dspProjImg->SetProfileProbePos(dataX, dataY);

  if (ui.radioButton_Profile_Hor->isChecked()) {
    dspProjImg->m_bDrawProfileX = true;
    dspProjImg->m_bDrawProfileY = false;
  } else {
    dspProjImg->m_bDrawProfileX = false;
    dspProjImg->m_bDrawProfileY = true;
  }

  // m_dspYKReconImage value itself
  int ROI_size = ui.lineEdit_ROI_size->text().toInt();
  if (ROI_size < 0) {
    return;
  }

  if (ROI_size > 0) {
    dspProjImg->setROI(
        qRound(dataX - ROI_size / 2.0), qRound(dataY - ROI_size / 2.0),
        qRound(dataX + ROI_size / 2.0), qRound(dataY + ROI_size / 2.0));
    dspProjImg->CalcImageInfo_ROI();
    dspProjImg->DrawROIOn(true);
    QString strMean;
    // strMean.sprintf("%5.1f", dspProjImg->m_fPixelMean_ROI);
    strMean = QString("%1").arg(
        (dspProjImg->m_fPixelMean_ROI / m_cbctrecon->m_multiplyFactor) +
            m_cbctrecon->m_fProjImgValueMin,
        0, 'f', 2);

    QString strSD;
    strSD = QString("%1").arg(
        dspProjImg->m_fPixelSD_ROI / m_cbctrecon->m_multiplyFactor, 0, 'f', 2);
    ui.lineEdit_ROI_mean->setText(strMean);
    ui.lineEdit_ROI_SD->setText(strSD);
  } else {
    dspProjImg->DrawROIOn(false);
  }

  SLT_DrawProjImages();
}

void CbctReconWidget::SLT_GoForcedProbePos() // when forced probe button was
                                             // clicked
{
  double fForcedProbePosX =
      ui.lineEdit_ForcedProbePosX->text().toDouble(); // data is the reference
  double fForcedProbePosY = ui.lineEdit_ForcedProbePosY->text().toDouble();
  double fForcedProbePosZ = ui.lineEdit_ForcedProbePosZ->text().toDouble();

  double dspWidth = 0.0;
  double dspHeight = 0.0;
  int dataWidth = 0;
  int dataHeight = 0;

  // First change the scene acc to Z value
  double originX, originY, originZ;
  double spacingX, spacingY, spacingZ;
  int sliceIdx = 0;

  int dataX, dataY;

  if (ui.radioButton_graph_proj->isChecked()) {
    if (m_cbctrecon->m_spProjImg3DFloat == nullptr) {
      return;
    }
    auto &ProjImg3D = m_cbctrecon->m_spProjImg3DFloat;

    originX = ProjImg3D->GetOrigin()[0];
    originY = ProjImg3D->GetOrigin()[1];
    originZ = ProjImg3D->GetOrigin()[2];

    spacingX = ProjImg3D->GetSpacing()[0];
    spacingY = ProjImg3D->GetSpacing()[1];
    spacingZ = ProjImg3D->GetSpacing()[2];

    sliceIdx = qRound((fForcedProbePosZ - originZ) / spacingZ);

    if (sliceIdx < 0 || sliceIdx >= m_cbctrecon->m_iImgCnt) {
      return;
    }

    ui.spinBoxImgIdx->setValue(sliceIdx); // Draw function is called

    dspWidth = ui.labelImageRaw->width();
    dspHeight = ui.labelImageRaw->height();

    dataWidth = m_cbctrecon->m_dspYKImgProj->m_iWidth;
    dataHeight = m_cbctrecon->m_dspYKImgProj->m_iHeight;

    dataX = qRound((fForcedProbePosX - originX) / spacingX);
    dataY = qRound((fForcedProbePosY - originY) / spacingY);

    if (dataX < 0 || dataX >= dataWidth || dataY < 0 || dataY >= dataHeight) {
      return;
    }

    ui.labelImageRaw->x =
        qRound(dataX / static_cast<double>(dataWidth) * dspWidth);
    ui.labelImageRaw->y =
        qRound(dataY / static_cast<double>(dataHeight) * dspHeight);

    SLT_CalculateROI_Proj();
  } else if (ui.radioButton_graph_recon->isChecked()) {
    if (m_cbctrecon->m_spCrntReconImg == nullptr) {
      return;
    }
    auto &CrntReconImg = m_cbctrecon->m_spCrntReconImg;
    originX = CrntReconImg->GetOrigin()[0];
    originY = CrntReconImg->GetOrigin()[1];
    originZ = CrntReconImg->GetOrigin()[2];

    spacingX = CrntReconImg->GetSpacing()[0];
    spacingY = CrntReconImg->GetSpacing()[1];
    spacingZ = CrntReconImg->GetSpacing()[2];

    sliceIdx = qRound((fForcedProbePosZ - originZ) / spacingZ);

    if (sliceIdx < 0 ||
        sliceIdx >=
            static_cast<int>(CrntReconImg->GetBufferedRegion().GetSize()[2])) {
      return;
    }

    ui.spinBoxReconImgSliceNo->setValue(sliceIdx); // Draw function is called

    dspWidth = ui.labelReconImage->width();
    dspHeight = ui.labelReconImage->height();

    dataWidth = m_cbctrecon->m_dspYKReconImage->m_iWidth;
    dataHeight = m_cbctrecon->m_dspYKReconImage->m_iHeight;

    dataX = qRound((fForcedProbePosX - originX) / spacingX);
    dataY = qRound((fForcedProbePosY - originY) / spacingY);

    if (dataX < 0 || dataX >= dataWidth || dataY < 0 || dataY >= dataHeight) {
      return;
    }

    ui.labelReconImage->x =
        qRound(dataX / static_cast<double>(dataWidth) * dspWidth);
    ui.labelReconImage->y =
        qRound(dataY / static_cast<double>(dataHeight) * dspHeight);

    SLT_CalculateROI_Recon();
  }
}

void CbctReconWidget::SLT_PostApplyFOVDispParam() {
  // m_cbctrecon->PostApplyFOVDispParam();
  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_DoPostProcessing() {
  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position

  float physPosX = ui.lineEdit_PostFOV_X->text().toFloat();
  float physPosY = ui.lineEdit_PostFOV_Y->text().toFloat();

  float physRadius = ui.lineEdit_PostFOV_R->text().toFloat();
  float physTablePosY = ui.lineEdit_PostTablePosY->text().toFloat();

  std::cout << "YKDEBUG " << physPosX << "," << physPosY << "," << physRadius
            << "," << physTablePosY << std::endl;

  m_cbctrecon->CropFOV3D(m_cbctrecon->m_spCrntReconImg, physPosX, physPosY,
                         physRadius, physTablePosY);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_PostProcCropInv() {
  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position

  double physPosX = ui.lineEdit_PostFOV_X->text().toDouble();
  double physPosY = ui.lineEdit_PostFOV_Y->text().toDouble();

  double physRadius = ui.lineEdit_PostFOV_R->text().toDouble();
  // double physTablePosY = ui.lineEdit_PostTablePosY->text().toDouble();

  UShortImageType::PointType origin =
      m_cbctrecon->m_spCrntReconImg->GetOrigin();
  UShortImageType::SpacingType spacing =
      m_cbctrecon->m_spCrntReconImg->GetSpacing();
  // UShortImageType::SizeType size =
  // m_spCrntReconImg->GetBufferedRegion().GetSize();

  // itk::ImageSliceConstIteratorWithIndex<FloatImageType> it (m_spReconImg,
  // m_spReconImg->GetRequestedRegion());
  itk::ImageSliceIteratorWithIndex<UShortImageType> it(
      m_cbctrecon->m_spCrntReconImg,
      m_cbctrecon->m_spCrntReconImg->GetRequestedRegion());

  // ImageSliceConstIteratorWithIndex<ImageType> it( image,
  // image->GetRequestedRegion() );
  // UShortImageType::SizeType imgSize =
  // m_spCrntReconImg->GetRequestedRegion().GetSize(); //1016x1016 x z

  // int width = imgSize[0];
  // int height = imgSize[1];

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  int iNumSlice = 0;
  int iPosX = 0;
  int iPosY = 0;

  // int i = 0;//height
  // int j = 0; // width

  double crntPhysX = 0.0;
  double crntPhysY = 0.0;

  while (!it.IsAtEnd()) {
    iPosY = 0;
    while (!it.IsAtEndOfSlice()) {
      iPosX = 0;
      while (!it.IsAtEndOfLine()) {
        // Calculate physical position

        crntPhysX = iPosX * static_cast<double>(spacing[0]) +
                    static_cast<double>(origin[0]);
        crntPhysY = iPosY * static_cast<double>(spacing[1]) +
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

  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.isEmpty()) {
    return;
  }

  // Get current folder
  QDir crntDir(dirPath);

  QInputDialog inputDlg;

  bool ok;
  QString textInput = QInputDialog::getText(
      this, "Input Dialog", "Set Patient ID and Name", QLineEdit::Normal,
      "PatientID_LastName_FirstName", &ok);

  // QString strEndFix = "YKP";
  QString strPatientID;
  QString strLastName;
  QString strFirstName;

  if (ok && !textInput.isEmpty()) {
    QStringList strListPtInfo = textInput.split("_");

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
      strPatientID = m_cbctrecon->m_strDCMUID;
    }
    // strPatientID = m_strDCMUID + "_" + strEndFix;
  } else {
    strPatientID = m_cbctrecon->m_strDCMUID;
  }

  if (strPatientID.isEmpty()) {
    return;
  }

  for (int i = 0;
       i < m_dlgRegistration->ui.comboBoxImgFixed->count(); i++) {
    m_dlgRegistration->ui.comboBoxImgFixed->setCurrentIndex(i);
    QString strDirName =
        m_dlgRegistration->ui.comboBoxImgFixed->currentText();
    bool tmpResult = crntDir.mkdir(strDirName); // what if the directory exists?
    if (!tmpResult) {
      std::cout
          << "DICOM dir seems to exist already. Files will be overwritten."
          << std::endl;
    }

    QString strSavingFolder = dirPath + "/" + strDirName;
    QString strFullName = strLastName + ", " + strFirstName;
    m_dlgRegistration->LoadImgFromComboBox(0, strDirName);

    m_cbctrecon->SaveUSHORTAsSHORT_DICOM(
        m_cbctregistration->m_spFixed, strPatientID, strFullName,
        strSavingFolder);
    QString mhaFileName = strSavingFolder + "/" + strDirName + ".mha";
    m_cbctrecon->ExportReconSHORT_HU(m_cbctregistration->m_spFixed,
                                     mhaFileName);
  }
  SLT_GeneratePOIData();
  QString angle_end_one("1");
  ui.lineEdit_AngEnd->setText(angle_end_one);
  SLT_ExportAngularWEPL_byFile();
}

void CbctReconWidget::SLT_ExportReconSHORT_HU() {
  QString strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "signed short meta image (*.mha)", nullptr,
      nullptr);
  if (strPath.length() <= 1) {
    return;
  }
  m_cbctrecon->ExportReconSHORT_HU(m_cbctrecon->m_spCrntReconImg, strPath);
}

void CbctReconWidget::SLT_DoBHC() {

  std::cout << "Beam hardening correction is under progress.." << std::endl;
  m_cbctrecon->DoBeamHardeningCorrection(); // only for m_spProjImg3D
  m_cbctrecon->SetMaxAndMinValueOfProjectionImage();

  SLT_DrawProjImages();
}

void CbctReconWidget::SLT_DoBowtieCorrection() {
  if (m_cbctrecon->m_spProjImg3DFloat == nullptr) {
    return;
  }

  if (m_cbctrecon->m_projFormat != HND_FORMAT) {
    std::cout
        << "Bow tie filtering should not be used for His data or Xim data!!"
        << std::endl;
    return;
  }

  QStringList strList = ui.comboBox_fBTcor->currentText().split(';');

  m_cbctrecon->BowtieByFit(ui.checkBox_Fullfan->isChecked(), strList);

  m_cbctrecon->SetMaxAndMinValueOfProjectionImage();
  SLT_DrawProjImages();
  std::cout << "Bow-tie correction done." << std::endl;
}

void CbctReconWidget::SLT_ViewRegistration() // default showing function
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
                                        bool bSave, bool use_cuda) {
  if (spVolImg3D == nullptr) {
    std::cout << "ERROR! No 3D-CT file. Load 3D CT file first" << std::endl;
    return;
  }

  if (m_cbctrecon->m_iCntSelectedProj < 1 && bSave) {
    std::cout << "Error! No projection image is loaded" << std::endl;
    return;
  }

  if (spGeometry->GetGantryAngles().empty()) {
    std::cout << "No geometry!" << std::endl;
    return;
  }

#ifndef USE_CUDA
  m_cbctrecon->CPU_ForwardProjection(spVolImg3D, spGeometry,
                                     spProjCT3D); // final moving image
#else
  if (use_cuda) {
    m_cbctrecon->CUDA_ForwardProjection(spVolImg3D, spGeometry,
                                        spProjCT3D); // final moving image
  } else {
    m_cbctrecon->CPU_ForwardProjection(spVolImg3D, spGeometry,
                                       spProjCT3D); // final moving image
  }
#endif // !USE_CUDA
  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout << " Patient DIR Path: "
              << m_cbctrecon->m_strPathPatientDir.toLocal8Bit().constData()
              << std::endl;

    bool manuallySelectedDir = false; // <- just to make sure I don't break
                                      // usecases of the older version.
    if (m_cbctrecon->m_strPathPatientDir.isEmpty()) {
      std::cout << "File save error!: No patient DIR name" << std::endl;

      m_cbctrecon->m_strPathPatientDir = QFileDialog::getExistingDirectory(
          this, tr("Open Directory"), ".",
          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
      if (m_cbctrecon->m_strPathPatientDir.length() <= 1) {
        return;
      }
      manuallySelectedDir = true;
    }

    // Get current folder
    QString subdir_images("IMAGES");
    QString strCrntDir = m_cbctrecon->m_strPathPatientDir + "/" +
                         subdir_images; // current Proj folder

    // Make a sub directory
    QDir crntDir(strCrntDir);

    if (!crntDir.exists()) {
      if (manuallySelectedDir) {
        QDir current_dir(m_cbctrecon->m_strPathPatientDir);
        bool success = current_dir.mkdir(subdir_images);
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

    QString fwdDirName = "fwd_" + m_cbctrecon->m_strDCMUID;

    bool tmpResult = crntDir.mkdir(fwdDirName); // what if the directory exists?

    if (!tmpResult) {
      std::cout << "FwdProj directory seems to exist already. Files will be "
                   "overwritten."
                << std::endl;
    }

    QString strSavingFolder = strCrntDir + "/" + fwdDirName;
    m_cbctrecon->SaveProjImageAsHIS(spProjCT3D, m_cbctrecon->m_arrYKBufProj,
                                    strSavingFolder, m_cbctrecon->m_fResampleF);
  }
}

void CbctReconWidget::SLT_DoScatterCorrection_APRIORI() {

  bool bExportProj_Fwd = ui.checkBox_ExportFwd->isChecked();
  bool bExportProj_Scat = ui.checkBox_ExportScat->isChecked();
  bool bExportProj_Cor = ui.checkBox_ExportCor->isChecked();

  // ForwardProjection(m_spRefCTImg, m_spCustomGeometry, m_spProjImgCT3D,
  // false); //final moving image
  if (m_cbctregistration->m_spMoving != nullptr) {
    ForwardProjection(
        m_cbctregistration->m_spMoving,
        m_cbctrecon->m_spCustomGeometry, m_cbctrecon->m_spProjImgCT3D,
        bExportProj_Fwd,
        ui.radioButton_UseCUDA->isChecked()); // final moving image
  } else if (m_cbctrecon->m_spRefCTImg != nullptr) {
    std::cout << "No Moving image in Registration is found. Ref CT image will "
                 "be used instead"
              << std::endl;
    ForwardProjection(
        m_cbctrecon->m_spRefCTImg, m_cbctrecon->m_spCustomGeometry,
        m_cbctrecon->m_spProjImgCT3D, bExportProj_Fwd,
        ui.radioButton_UseCUDA->isChecked()); // final moving image
  } else {
    std::cout << "Error!: No ref image for forward projection is found."
              << std::endl;
    return;
  }

  // YKTEMP
  std::cout << "ProjImgCT Size = "
            << m_cbctrecon->m_spProjImgCT3D->GetBufferedRegion().GetSize()[0]
            << ", "
            << m_cbctrecon->m_spProjImgCT3D->GetBufferedRegion().GetSize()[1]
            << ", "
            << m_cbctrecon->m_spProjImgCT3D->GetBufferedRegion().GetSize()[2]
            << std::endl;
  std::cout << "ProjImgCT origin = "
            << m_cbctrecon->m_spProjImgCT3D->GetOrigin()[0] << ", "
            << m_cbctrecon->m_spProjImgCT3D->GetOrigin()[1] << ", "
            << m_cbctrecon->m_spProjImgCT3D->GetOrigin()[2] << std::endl;
  std::cout << "ProjImgCT spacing = "
            << m_cbctrecon->m_spProjImgCT3D->GetSpacing()[0] << ", "
            << m_cbctrecon->m_spProjImgCT3D->GetSpacing()[1] << ", "
            << m_cbctrecon->m_spProjImgCT3D->GetSpacing()[2] << std::endl;

  // double scaResam = ui.lineEdit_scaResam->text().toDouble();
  double scaMedian = ui.lineEdit_scaMedian->text().toDouble();
  double scaGaussian = ui.lineEdit_scaGaussian->text().toDouble();

  std::cout << "Generating scatter map is ongoing..." << std::endl;

  m_cbctrecon->GenScatterMap_PriorCT(
      m_cbctrecon->m_spProjImgRaw3D, m_cbctrecon->m_spProjImgCT3D,
      m_cbctrecon->m_spProjImgScat3D, scaMedian, scaGaussian,
      m_cbctrecon->m_iFixedOffset_ScatterMap,
      bExportProj_Scat); // void GenScatterMap2D_PriorCT()

  std::cout << "To account for the mAs values, the intensity scale factor of "
            << m_cbctrecon->GetRawIntensityScaleFactor(
                   m_cbctrecon->m_strRef_mAs, m_cbctrecon->m_strCur_mAs)
            << "was multiplied during scatter correction to avoid negative "
               "scatter"
            << std::endl;

  ui.lineEdit_CurmAs->setText(m_cbctrecon->m_strCur_mAs);
  ui.lineEdit_RefmAs->setText(m_cbctrecon->m_strRef_mAs);

  m_cbctrecon->m_spProjImgCT3D->Initialize(); // memory saving

  std::cout << "Scatter correction is in progress..." << std::endl;

  int postScatMedianSize = ui.lineEdit_scaPostMedian->text().toInt();
  m_cbctrecon->ScatterCorr_PrioriCT(
      m_cbctrecon->m_spProjImgRaw3D, m_cbctrecon->m_spProjImgScat3D,
      m_cbctrecon->m_spProjImgCorr3D, m_cbctrecon->m_iFixedOffset_ScatterMap,
      postScatMedianSize, bExportProj_Cor);
  m_cbctrecon->m_spProjImgScat3D->Initialize(); // memory saving

  std::cout << "AfterCorrectionMacro is ongoing..." << std::endl;

  // Update UI
  ui.pushButton_DoRecon->setEnabled(true);
  ui.spinBoxImgIdx->setMinimum(0);
  int iSizeZ =
      m_cbctrecon->m_spProjImg3DFloat->GetRequestedRegion().GetSize()[2];
  ui.spinBoxImgIdx->setMaximum(iSizeZ - 1);
  ui.spinBoxImgIdx->setValue(0);
  m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();
  SLT_DrawProjImages(); // Update Table is called

  m_cbctrecon->AfterScatCorrectionMacro(
      ui.radioButton_UseCUDA->isChecked(),
      ui.checkBox_ExportVolDICOM->isChecked());

  // Skin removal (using CT contour w/ big margin)
  std::cout
    << "Post  FDK reconstruction is done. Moving on to post skin removal"
    << std::endl;

  m_cbctregistration->PostSkinRemovingCBCT(m_cbctrecon->m_spRawReconImg);
  m_cbctregistration->PostSkinRemovingCBCT(m_cbctrecon->m_spScatCorrReconImg);

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

  std::cout << "Updating ReconImage..";
  QString updated_text = QString("Scatter corrected CBCT");
  UpdateReconImage(m_cbctrecon->m_spScatCorrReconImg,
                   updated_text); // main GUI update

  std::cout << "FINISHED!Scatter correction: CBCT DICOM files are saved"
            << std::endl;
}

// called whenver recon 3D image for display changes.
void CbctReconWidget::UpdateReconImage(UShortImageType::Pointer &spNewImg,
                                       QString &fileName) {
  m_cbctrecon->m_spCrntReconImg = spNewImg;

  UShortImageType::PointType origin_new =
      m_cbctrecon->m_spCrntReconImg->GetOrigin();
  UShortImageType::SpacingType spacing_new =
      m_cbctrecon->m_spCrntReconImg->GetSpacing();
  UShortImageType::SizeType size_new =
      m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize();

  std::cout << "New Origin" << origin_new << std::endl;
  std::cout << "New spacing" << spacing_new << std::endl;
  std::cout << "New size" << size_new << std::endl;

  ui.lineEdit_Cur3DFileName->setText(fileName);

  UShortImageType::SizeType size =
      m_cbctrecon->m_spCrntReconImg->GetRequestedRegion().GetSize();

  m_cbctrecon->m_dspYKReconImage->CreateImage(size[0], size[1], 0);

  disconnect(ui.spinBoxReconImgSliceNo, SIGNAL(valueChanged(int)), this,
             SLOT(SLT_DrawReconImage()));

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(size[2] - 1);

  int initVal = qRound((size[2] - 1) / 2.0);
  // SLT_DrawReconImage(); //Update Table, Update Graph

  // m_dspYKReconImage->CreateImage(size_trans[0], size_trans[1],0);
  SLT_InitializeGraphLim();

  ui.spinBoxReconImgSliceNo->setValue(initVal);
  ui.radioButton_graph_recon->setChecked(true);

  connect(ui.spinBoxReconImgSliceNo, SIGNAL(valueChanged(int)), this,
          SLOT(SLT_DrawReconImage()));

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_TempAudit() {
  if (m_cbctrecon->m_spRawReconImg != nullptr) {
    std::cout << "m_spRawReconImg " << m_cbctrecon->m_spRawReconImg
              << std::endl;
  }

  if (m_cbctrecon->m_spRefCTImg != nullptr) {
    std::cout << "m_spRefCTImg " << m_cbctrecon->m_spRefCTImg << std::endl;
  }

  if (m_cbctrecon->m_spCrntReconImg != nullptr) {
    std::cout << "m_spCrntReconImg " << m_cbctrecon->m_spCrntReconImg
              << std::endl;
  }
}

void CbctReconWidget::SLT_LoadPlanCT_USHORT() {
  // typedef itk::ImageFileWriter<FloatImageType> WriterType;
  using ReaderType = itk::ImageFileReader<UShortImageType>;
  ReaderType::Pointer reader = ReaderType::New();

  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_cbctrecon->m_strPathDirDefault,
      "Projection file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  reader->SetFileName(fileName.toLocal8Bit().constData());
  reader->Update();

  m_cbctrecon->m_spRefCTImg = reader->GetOutput();
  QString ref_ct = QString("RefCT");
  UpdateReconImage(m_cbctrecon->m_spRefCTImg, ref_ct);

  m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT, REGISTER_MANUAL_RIGID);
}

void CbctReconWidget::SLT_CalcAndSaveAngularWEPL() // single point
{
  std::vector<WEPLData> vOutputWEPL;

  double fAngleGap = ui.lineEdit_WEPL_AngRes->text().toDouble();
  double fAngleStart = ui.lineEdit_AngStart->text().toDouble();
  double fAngleEnd = ui.lineEdit_AngEnd->text().toDouble();

  VEC3D curPOI{};
  curPOI.x = ui.lineEdit_ForcedProbePosX->text().toDouble(); // in mm
  curPOI.y = ui.lineEdit_ForcedProbePosY->text().toDouble();
  curPOI.z = ui.lineEdit_ForcedProbePosZ->text().toDouble();

  m_cbctrecon->GetAngularWEPL_SinglePoint(m_cbctrecon->m_spCrntReconImg,
                                          fAngleGap, fAngleStart, fAngleEnd,
                                          curPOI, 0, vOutputWEPL, true);
  std::cout << "Computed WEPL points: " << vOutputWEPL.size() << std::endl;

  // export arrWEPL
  QString filePath = QFileDialog::getSaveFileName(
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
  if (m_cbctrecon->m_spProjImg3DFloat == nullptr) {
    return;
  }

  UShortImageType::Pointer spIntensityRaw;
  m_cbctrecon->ConvertLineInt2Intensity(m_cbctrecon->m_spProjImg3DFloat,
                                        spIntensityRaw, 65535);

  using ScatterFilterType =
      rtk::BoellaardScatterCorrectionImageFilter<UShortImageType,
                                                 UShortImageType>;

  ScatterFilterType::Pointer spScatFilter = ScatterFilterType::New();

  double airThre = ui.lineEdit_uniAirThre->text().toDouble();
  double scat2PrimRatio = ui.lineEdit_uniSPR->text().toDouble();
  double nonNagativity = ui.lineEdit_uniNegativity->text().toDouble();

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

  m_cbctrecon->ConvertIntensity2LineInt(spIntensityUniformCorr,
                                        m_cbctrecon->m_spProjImg3DFloat, 65535);

  // ConvertLineInt2Intensity(m_spProjImg3DFloat, m_spProjImgRaw3D, 65535);

  m_cbctrecon->m_spProjImgRaw3D = spIntensityUniformCorr;

  ui.spinBoxImgIdx->setValue(0);
  m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();

  this->SLT_DrawProjImages(); // Update Table is called

  std::cout << "FINISHED!: Uniform Scatter correction for raw projection "
               "images (Boallaard method) is completed. Proceed to "
               "reconstruction"
            << std::endl;
}

void CbctReconWidget::SLT_FileExportShortDICOM_CurrentImg() {

  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.isEmpty()) {
    return;
  }

  // Get current folder
  QDir crntDir(dirPath);

  QInputDialog inputDlg;

  bool ok;
  QString textInput = QInputDialog::getText(
      this, "Input Dialog", "Set Patient ID and Name", QLineEdit::Normal,
      "PatientID_LastName_FirstName", &ok);

  // QString strEndFix = "YKP";
  QString strPatientID;
  QString strLastName;
  QString strFirstName;

  if (ok && !textInput.isEmpty()) {
    QStringList strListPtInfo = textInput.split("_");

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
      strPatientID = m_cbctrecon->m_strDCMUID;
    }
    // strPatientID = m_strDCMUID + "_" + strEndFix;
  } else {
    strPatientID = m_cbctrecon->m_strDCMUID;
  }

  if (strPatientID.isEmpty()) {
    return;
  }

  QString strDirName = strPatientID + "_DCM";
  bool tmpResult = crntDir.mkdir(strDirName); // what if the directory exists?
  if (!tmpResult) {
    std::cout << "DICOM dir seems to exist already. Files will be overwritten."
              << std::endl;
  }

  QString strSavingFolder = dirPath + "/" + strDirName;
  QString strFullName = strLastName + ", " + strFirstName;
  m_cbctrecon->SaveUSHORTAsSHORT_DICOM(m_cbctrecon->m_spCrntReconImg,
                                       strPatientID, strFullName,
                                       strSavingFolder);
}

void CbctReconWidget::SLT_AddConstHUToCurImg() {
  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }
  int addingVal = ui.lineEdit_AddConstHU->text().toInt();
  m_cbctrecon->AddConstHU(m_cbctrecon->m_spCrntReconImg, addingVal);
  QString updated_text = QString("Added%1").arg(addingVal);
  UpdateReconImage(m_cbctrecon->m_spCrntReconImg, updated_text);
}

void CbctReconWidget::SLT_SetCBCTSkinRSPath() {
  QString strPath = QFileDialog::getOpenFileName(
      this, "Open RS file", m_cbctrecon->m_strPathDirDefault,
      "DICOM RS (*.dcm)", nullptr, nullptr);

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_PathCBCTSkinPath->setText(strPath);
}

void CbctReconWidget::SLT_CropSkinUsingThreshold() {

  QString update_text = QString("Thresh-based skin cropped image");
  auto thresh = ui.lineEdit_Threshold->text().toInt();
  auto erode = ui.lineEdit_ErodeRadius->text().toInt();
  auto dilate = ui.lineEdit_DilateRadius->text().toInt();

  auto imgType = m_cbctrecon->CropSkinUsingThreshold(thresh, erode, dilate);

  switch (imgType) {
  case 1:
    UpdateReconImage(m_cbctrecon->m_spRawReconImg, update_text);
    break;
  case 2:
    UpdateReconImage(m_cbctrecon->m_spRefCTImg, update_text);
    break;
  case 3:
    UpdateReconImage(m_cbctrecon->m_spScatCorrReconImg, update_text);
    break;
  default:
    std::cerr << "WTF!?" << std::endl;
    break;
  }
}

void CbctReconWidget::SLT_CropSkinUsingRS() {
  QString strPathRS = ui.lineEdit_PathCBCTSkinPath->text();
  if (strPathRS.length() < 1) {
    return;
  }

  double croppingMargin = ui.lineEdit_SkinMargin->text().toDouble();
  QString update_text = QString("RS-based skin cropped image");
  if (m_cbctrecon->m_spCrntReconImg == m_cbctrecon->m_spRawReconImg) {
    m_cbctregistration->CropSkinUsingRS(
        m_cbctrecon->m_spRawReconImg, strPathRS, croppingMargin);
    UpdateReconImage(m_cbctrecon->m_spRawReconImg, update_text);
  } else if (m_cbctrecon->m_spCrntReconImg == m_cbctrecon->m_spRefCTImg) {
    m_cbctregistration->CropSkinUsingRS(m_cbctrecon->m_spRefCTImg,
                                                     strPathRS, croppingMargin);
    UpdateReconImage(m_cbctrecon->m_spRefCTImg, update_text);
  } else if (m_cbctrecon->m_spCrntReconImg ==
             m_cbctrecon->m_spScatCorrReconImg) {
    m_cbctregistration->CropSkinUsingRS(
        m_cbctrecon->m_spScatCorrReconImg, strPathRS, croppingMargin);
    UpdateReconImage(m_cbctrecon->m_spScatCorrReconImg, update_text);
  }
}

void CbctReconWidget::SLT_ExportAngularWEPL_byFile() {
  // export arrWEPL
  QString filePath = QFileDialog::getSaveFileName(
      this, "Save data", m_cbctrecon->m_strPathDirDefault,
      "txt image file (*.txt)", nullptr,
      nullptr); // Filename don't need to exist

  if (filePath.length() < 1) {
    return;
  }

  double fAngleStart = ui.lineEdit_AngStart->text().toDouble();
  double fAngleEnd = ui.lineEdit_AngEnd->text().toDouble();
  double fAngleGap = ui.lineEdit_WEPL_AngRes->text().toDouble();

  m_cbctrecon->ExportAngularWEPL_byFile(filePath, fAngleStart, fAngleEnd,
                                        fAngleGap);
}

void CbctReconWidget::SLT_GeneratePOIData() // it fills m_vPOI_DCM
{
  m_cbctrecon->GeneratePOIData(ui.checkBox_AP->isChecked(),
                               ui.lineEdit_PostTablePosY->text().toDouble());
}

void CbctReconWidget::SLT_LoadPOIData() // it fills m_vPOI_DCM
{
  if (!m_cbctrecon->m_vPOI_DCM.empty()) {
    m_cbctrecon->m_vPOI_DCM.clear();
  }

  QString filePath = QFileDialog::getOpenFileName(
      this, "POI data file", m_cbctrecon->m_strPathDirDefault,
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
    QStringList strList = strLine.split('\t'); // tab
    // third one is the organ name
    if (strList.length() < 3) {
      std::cout << "abnormal file expression." << std::endl;
      break;
    }
    fPOI.x = strList.at(0).toDouble();
    fPOI.y = strList.at(1).toDouble();
    fPOI.z = strList.at(2).toDouble();

    m_cbctrecon->m_vPOI_DCM.push_back(fPOI);
  }
  for (int i = 0; i < static_cast<int>(m_cbctrecon->m_vPOI_DCM.size()); i++) {
    std::cout << "Data " << i << "	" << m_cbctrecon->m_vPOI_DCM.at(i).x
              << ", " << m_cbctrecon->m_vPOI_DCM.at(i).y << ", "
              << m_cbctrecon->m_vPOI_DCM.at(i).z << std::endl;
  }
  std::cout << "POI data has been loaded. " << m_cbctrecon->m_vPOI_DCM.size()
            << " data points are read" << std::endl;
  fin.close();
}

void CbctReconWidget::SLT_StartSyncFromSharedMem() {
  // int msInterval = ui.lineEditTimerInterval->text().toInt();

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
  HANDLE hSemaphore = OpenSemaphore(SEMAPHORE_ALL_ACCESS, FALSE, "YKSemaphore");
  // increase counter
  //  LONG prev_counter;
  // ReleaseSemaphore(hSemaphore, 1, &prev_counter);
  // decrease counter
  // ReleaseSemaphore(hSemaphore, 1, &prev_counter);

  std::ofstream fout;
  fout.open("E:\\SemphoreLogC++.txt");

  int max_cnt = 200;

  int cnt = 0;
  while (cnt < max_cnt) {
    cnt++;

    DWORD dwWaitResult = WaitForSingleObject(hSemaphore, 5);

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

  if (m_cbctrecon->m_arrYKImage.empty() || m_cbctrecon->m_iImgCnt != 1) {
    return;
  }

  m_busyTimer = true;

  // Look into the shared mem
  TCHAR szName[] = TEXT("YKSharedMemory");
  HANDLE handle = OpenFileMapping(FILE_MAP_READ, FALSE, szName);

  if (handle == nullptr) {
    std::cout << "Cannot open Mapped file" << std::endl;
    SLT_StopSyncFromSharedMem();
    return;
  }

  int size = 1024 * 1024 * 2;
  auto pix_size = static_cast<int>(size / 2.0);
  unsigned char *charBuf = nullptr;
  charBuf = static_cast<unsigned char *>(
      MapViewOfFile(handle, FILE_MAP_READ, 0, 0, size));

  if (charBuf == nullptr) {
    std::cout << "Shared memory was not read. Timer will be stopped"
              << std::endl;
    SLT_StopSyncFromSharedMem();
    CloseHandle(handle);
    delete[] charBuf;
    return;
  }

  // byte array to unsigned short
  // assuming little endian

  // unsigned short* imgBuf = new unsigned short [pix_size];
  int idxA = 0;
  int idxB = 0;

  for (int i = 0; i < pix_size; i++) {
    idxA = i * 2 + 1;
    idxB = i * 2;
    // 0: 1,0  1: 3,2 ...
    m_cbctrecon->m_arrYKImage.at(0).m_pData[i] =
        ((charBuf[idxA] << 8) | charBuf[idxB]); // little endian
  }

  ui.spinBoxImgIdx->setValue(0);
  SLT_DrawRawImages();

  CloseHandle(handle);

  m_busyTimer = false;
#endif
}

void CbctReconWidget::SLTM_ViewExternalCommand() {
  m_dlgExternalCommand->show();
}

void CbctReconWidget::SLTM_LoadDICOMdir() {
  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.length() <= 1) {
    return;
  }

  if (m_cbctrecon->ReadDicomDir(dirPath)) {

    m_dlgRegistration->UpdateVOICombobox(PLAN_CT);
    QString update_text = QString("DICOM reference image");
    UpdateReconImage(m_cbctrecon->m_spRefCTImg, update_text);

    m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT, REGISTER_MANUAL_RIGID);
  }
}

void CbctReconWidget::SLTM_LoadRTKoutput() {
  QString filePath = QFileDialog::getOpenFileName(
      this, "Open Image", m_cbctrecon->m_strPathDirDefault,
      "rtk output float image (*.mha)", nullptr, nullptr);
  m_cbctrecon->LoadExternalFloatImage(filePath, true);
  QString strCrntFileName;
  QFileInfo outFileInfo(filePath);
  strCrntFileName = outFileInfo.fileName();
  UpdateReconImage(m_cbctrecon->m_spRawReconImg, strCrntFileName);
}

// Only can be used for m_spRawRecon // NOT USED AT ALL?
void CbctReconWidget::FileExportByGUI() // USHORT
{
  QString outputFilePath = ui.lineEdit_OutputFilePath->text();
  QFileInfo outFileInfo(outputFilePath);
  QDir outFileDir = outFileInfo.absoluteDir();

  // bool b = outFileDir.exists();
  // QString tmpPath = outFileDir.absolutePath();

  if (outputFilePath.length() < 2 || !outFileDir.exists()) {
    std::cout << "No available output path. Should be exported later"
              << std::endl;
  } else {
    using WriterType = itk::ImageFileWriter<UShortImageType>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputFilePath.toLocal8Bit().constData());
    writer->SetUseCompression(true); // not exist in original code (rtkfdk)
    writer->SetInput(m_cbctrecon->m_spRawReconImg);

    std::cout << "Writing the image to: "
              << outputFilePath.toLocal8Bit().constData() << std::endl;

    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());

    std::cout << std::endl;
    std::cout << "The output image was successfully saved" << std::endl;
  }
}

void CbctReconWidget::SLT_OutPathEdited() {
  if (!ui.lineEdit_OutputFilePath->text().isEmpty()) {
    ui.lineEdit_outImgDim_LR->setEnabled(true);
    ui.lineEdit_outImgDim_AP->setEnabled(true);
    ui.lineEdit_outImgDim_SI->setEnabled(true);
    ui.lineEdit_outImgSp_LR->setEnabled(true);
    ui.lineEdit_outImgSp_AP->setEnabled(true);
    ui.lineEdit_outImgSp_SI->setEnabled(true);
  } else {
    ui.lineEdit_outImgDim_LR->setEnabled(false);
    ui.lineEdit_outImgDim_AP->setEnabled(false);
    ui.lineEdit_outImgDim_SI->setEnabled(false);
    ui.lineEdit_outImgSp_LR->setEnabled(false);
    ui.lineEdit_outImgSp_AP->setEnabled(false);
    ui.lineEdit_outImgSp_SI->setEnabled(false);
  }
}

void CbctReconWidget::SLT_MedianFilterDoNow() {
  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  /*QString strCrntFileName;
  QFileInfo outFileInfo(strPath);
  strCrntFileName = outFileInfo.fileName();*/

  UShortImageType::SizeType indexRadius{};
  indexRadius[0] = ui.lineEdit_PostMedSizeX->text().toInt(); // radius along x
  indexRadius[1] = ui.lineEdit_PostMedSizeY->text().toInt(); // radius along y
  indexRadius[2] = ui.lineEdit_PostMedSizeZ->text().toInt(); // radius along y

  if (indexRadius[0] != 0 || indexRadius[1] != 0 || indexRadius[2] != 0) {
    m_cbctrecon->MedianFilterByGUI(indexRadius);

    QString prevFileName = ui.lineEdit_Cur3DFileName->text();
    UpdateReconImage(m_cbctrecon->m_spCrntReconImg,
                     prevFileName.append("_med"));
  } else {
    std::cout << "Not valid median window" << std::endl;
  }
}

void CbctReconWidget::SLT_Export2DDose_TIF() // 2D dose from current displayed
                                             // image of reconstruction
{
  if (m_cbctrecon->m_dspYKReconImage == nullptr) {
    return;
  }

  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  QString strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "signed short meta image (*.tif)", nullptr,
      nullptr);
  if (strPath.length() <= 1) {
    return;
  }

  auto originLeft =
      static_cast<double>(m_cbctrecon->m_spCrntReconImg->GetOrigin()[0]);
  auto originTop = static_cast<double>(
      m_cbctrecon->m_spCrntReconImg->GetOrigin()[1]); // not sure...

  auto spacingX =
      static_cast<double>(m_cbctrecon->m_spCrntReconImg->GetSpacing()[0]);
  auto spacingY = static_cast<double>(
      m_cbctrecon->m_spCrntReconImg->GetSpacing()[1]); // not sure...

  if (!SaveDoseGrayImage(strPath.toLocal8Bit().constData(),
                         m_cbctrecon->m_dspYKReconImage->m_iWidth,
                         m_cbctrecon->m_dspYKReconImage->m_iHeight, spacingX,
                         spacingY, originLeft, originTop,
                         m_cbctrecon->m_dspYKReconImage->m_pData)) {
    std::cout << "Failed in save gray dose file" << std::endl;
  } else {
    std::cout << "image exported successfully." << std::endl;
  }
}
void CbctReconWidget::SLTM_Export2DDoseMapAsMHA() {

  QString strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "itk compatible meta image (*.mha)", nullptr,
      nullptr);

  m_cbctrecon->Export2DDoseMapAsMHA(strPath);
}

void CbctReconWidget::SLTM_ExportProjGeometryTXT() {

  QString strPath = QFileDialog::getSaveFileName(
      this, "Save text file", "", "text (*.txt)", nullptr, nullptr);

  m_cbctrecon->ExportProjGeometryTXT(strPath);
}

void CbctReconWidget::SLTM_ForwardProjection() {
  if (m_cbctrecon->m_spRawReconImg == nullptr) {
    return;
  }

  GeometryType::Pointer crntGeometry = GeometryType::New();

  if (m_cbctrecon->m_spCustomGeometry == nullptr) {
    std::cout << "No geometry is ready. moving on to 360 projection"
              << std::endl;

    double curSID = 1000.0;
    double curSDD = 1536.0;
    double curProjOffsetX = 0.0;
    double curProjOffsetY = 0.0;
    double curOutOfPlaneAngles = 0.0;
    double curInPlaneAngles = 0.0;
    double curSrcOffsetX = 0.0;
    double curSrcOffsetY = 0.0;

    double curMVGantryAngle = 0.0;

    // double startAngle = 180.0; //kV = 270.0, CW
    double startAngle = 270; // kV = 360.0, CW
    // int NumOfProj = 360;
    int NumOfProj = 1;

    for (int i = 0; i < NumOfProj; i++) {
      curMVGantryAngle = startAngle + i;
      if (curMVGantryAngle > 360.0) {
        curMVGantryAngle = curMVGantryAngle - 360.0;
      }
      // AddProjection: current CBCT software version only requires MV gantry
      // angle!!!
      crntGeometry->AddProjection(
          curSID, curSDD, curMVGantryAngle, curProjOffsetX,
          curProjOffsetY,                        // Flexmap
          curOutOfPlaneAngles, curInPlaneAngles, // In elekta, these are 0
          curSrcOffsetX, curSrcOffsetY);         // In elekta, these are 0
    }

    ForwardProjection(m_cbctrecon->m_spRawReconImg, crntGeometry,
                      m_cbctrecon->m_spProjImgRaw3D, false,
                      ui.radioButton_UseCUDA->isChecked());
    // Save proj3D;

    // QString outputPath = "D:/ProjTemplate.mha";
    // QString outputPath = "D:/2D3DRegi/FwdProj_0.mha";
    QString outputPath = QFileDialog::getSaveFileName(
        this, "File path to save", m_cbctrecon->m_strPathDirDefault,
        "Projection stack (*.mha)", nullptr,
        nullptr); // Filename don't need to exist
    if (outputPath.length() <= 1) {
      return;
    }

    using WriterType = itk::ImageFileWriter<UShortImageType>;
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName(outputPath.toLocal8Bit().constData());
    // writer->SetUseCompression(true);
    writer->SetUseCompression(true); // for plastimatch
    writer->SetInput(m_cbctrecon->m_spProjImgRaw3D);
    writer->Update();

    return;
  }
  // if there is a geometry

  int cntProj = m_cbctrecon->m_spCustomGeometry->GetGantryAngles().size();

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

  for (int i = 0; i < cntProj; i++) {
    double curSID =
        m_cbctrecon->m_spCustomGeometry->GetSourceToIsocenterDistances().at(i);
    double curSDD =
        m_cbctrecon->m_spCustomGeometry->GetSourceToDetectorDistances().at(i);
    double curGantryAngle =
        m_cbctrecon->m_spCustomGeometry->GetGantryAngles().at(i);
    double kVAng = curGantryAngle * 360. / (2. * itk::Math::pi);
    double MVAng =
        kVAng - (m_cbctrecon->m_projFormat == HIS_FORMAT ? 0.0 : 90.0);
    if (MVAng < 0.0) {
      MVAng = MVAng + 360.0;
    }
    curGantryAngle = MVAng;

    double curProjOffsetX =
        m_cbctrecon->m_spCustomGeometry->GetProjectionOffsetsX().at(i);
    double curProjOffsetY =
        m_cbctrecon->m_spCustomGeometry->GetProjectionOffsetsY().at(i);

    double curOutOfPlaneAngles =
        m_cbctrecon->m_spCustomGeometry->GetOutOfPlaneAngles().at(i);
    double curInPlaneAngles =
        m_cbctrecon->m_spCustomGeometry->GetInPlaneAngles().at(i);

    double curSrcOffsetX =
        m_cbctrecon->m_spCustomGeometry->GetSourceOffsetsX().at(i);
    double curSrcOffsetY =
        m_cbctrecon->m_spCustomGeometry->GetSourceOffsetsY().at(i);

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

  ForwardProjection(m_cbctrecon->m_spRawReconImg, crntGeometry,
                    m_cbctrecon->m_spProjImgRaw3D, true,
                    ui.radioButton_UseCUDA->isChecked());

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
  float curResampleF = m_cbctrecon->m_fResampleF;
  ui.lineEdit_DownResolFactor->setText("1.0");
  SLT_LoadSelectedProjFiles();
  m_cbctrecon->m_fResampleF = curResampleF;
  ui.lineEdit_DownResolFactor->setText(
      QString("%1").arg(m_cbctrecon->m_fResampleF));

  // Scatter correction

  SLT_DoScatterCorrection_APRIORI();
}

void CbctReconWidget::SLTM_FullScatterCorrectionMacroAP() // single. should be
                                                          // called after HIS
                                                          // folder is defined
{
  if (m_cbctrecon->m_strPathPatientDir.length() < 2) {
    return;
  }

  enREGI_IMAGES enRegImg = REGISTER_DEFORM_FINAL;
  bool bFullResolForFinalRecon = false;

  bool bIntensityShift = true;

  FullScatterCorrectionMacroSingle(m_cbctrecon->m_strPathPatientDir, enRegImg,
                                   bFullResolForFinalRecon, false,
                                   bIntensityShift);
}

void CbctReconWidget::SLTM_BatchScatterCorrectionMacroAP() {
  // Scatter parameters
  QTime batchmodeTime = QTime::currentTime();

  // 1) Get img_ file lists
  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open IMAGES Directory"), m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  QDir dirIMAGES = QDir(dirPath);

  if (!dirIMAGES.exists()) {
    return;
  }

  QFileInfoList listDir;
  listDir = dirIMAGES.entryInfoList(QDir::Dirs, QDir::Name);

  int iCntProjDir = listDir.size();

  std::cout << "Found directory number= " << iCntProjDir - 2 << std::endl;

  QString curProjDirPath;
  if (iCntProjDir <= 2) // only /. and /.. exist
  {
    std::cout << "Error! No projection directory exists." << std::endl;
    return;
  }
  // Several questions to set params
  // 1) Output Dir
  QString strOutDirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Output Directory"), m_cbctrecon->m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (strOutDirPath.length() < 1) {
    return;
  }

  // 2) fwd reference: manual or auto-rigid or deformed

  bool ok;
  QString text = QInputDialog::getText(
      this, "Input Dialog",
      "Reference for forward projection([0] Deformed CT, [1] auto-rigid CT, "
      "[2] manual-aligned CT(dcm_plan), [3] DeformedCT_skipAutoRigid",
      QLineEdit::Normal, "0", &ok);

  enREGI_IMAGES enRegImg = REGISTER_DEFORM_FINAL;

  if (ok && !text.isEmpty()) {
    int iRefImgVal = text.toInt();

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

  int res;
  // 3) Fine resol option
  bool bFullResolForFinalRecon = false;
  bool bIntensityShift = false;
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
  QString strMsg = "Intensity shift for raw CBCT?";
  msgBox.setText(strMsg);
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  res = msgBox.exec();

  if (res == QMessageBox::Yes) {
    bIntensityShift = true;
  }

  bool bExportShortImages = false;
  QMessageBox msgBox2;
  QString strMsg2 = "Export short images after correction?";
  msgBox2.setText(strMsg2);
  msgBox2.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  res = msgBox2.exec();
  if (res == QMessageBox::Yes) {
    bExportShortImages = true;
  }

  int cntHisDir = 0;
  for (int i = 2; i < iCntProjDir; i++) // start from 2
  {
    curProjDirPath = listDir.at(i).absoluteFilePath();

    if (curProjDirPath.contains("img_")) {
      std::cout << "Found projection dir number: " << cntHisDir
                << ", Current Proj Path: "
                << curProjDirPath.toLocal8Bit().constData() << std::endl;
      m_cbctrecon->SetProjDir(curProjDirPath);
      FullScatterCorrectionMacroSingle(strOutDirPath, enRegImg,
                                       bFullResolForFinalRecon,
                                       bExportShortImages, bIntensityShift);
      cntHisDir++;
    }
  }

  float elapsedSec = batchmodeTime.elapsed() / 1000.0;

  std::cout << "Batch mode calculation is done! "
            << QString::number(elapsedSec, 'f', 2).toLocal8Bit().constData()
            << " seconds was spent for " << cntHisDir << " cases" << std::endl;
}

// Uses SLT functions heavily:
bool CbctReconWidget::FullScatterCorrectionMacroSingle(
    QString &outputDirPath, enREGI_IMAGES enFwdRefImg, bool bFullResolRecon,
    bool bExportImages, bool bCBCT_IntensityShift) {
  if (m_cbctrecon->m_strDCMUID.length() < 1) {
    return false;
  }

  m_cbctrecon->m_bMacroContinue = true;

  bool bFOVCropping =
      ui.checkBox_PostDispObjOn->isChecked(); // this button is for display, but
                                              // use it for cropping option in
                                              // macro mode
  float physPosX = ui.lineEdit_PostFOV_X->text().toFloat();
  float physPosY = ui.lineEdit_PostFOV_Y->text().toFloat();
  float physRadius = ui.lineEdit_PostFOV_R->text().toFloat();
  float physTablePosY = ui.lineEdit_PostTablePosY->text().toFloat();

  // Load Pushbutton
  SLT_LoadSelectedProjFiles();

  // float fOldValTruncation =
  // ui.lineEdit_Ramp_TruncationCorrection->text().toFloat();;

  SLT_DoReconstruction();
  // ui.lineEdit_Ramp_TruncationCorrection->setText(QString("0.0"));

  int addingVal = ui.lineEdit_AddConstHU->text().toInt();

  if (addingVal != 0) {
    std::cout << "Raw CBCT is being added by HU of: " << addingVal << std::endl;
    SLT_AddConstHUToCurImg();
  }

  if (bFOVCropping) {
    // Crop CBCT with predetermined FOV/ Table
    std::cout << "FOV cropping is under way..." << std::endl;
    m_cbctrecon->CropFOV3D(m_cbctrecon->m_spRawReconImg, physPosX, physPosY,
                           physRadius, physTablePosY);
  }

  SLT_ViewRegistration();

  m_dlgRegistration->SLT_PreProcessCT();
  if (!m_cbctrecon->m_bMacroContinue) {
    std::cout << "Stopped during MacroSingle due to error in PreProcessCT"
              << std::endl;
    return false;
  }

  QString strSuffix;
  switch (enFwdRefImg) {
  case REGISTER_MANUAL_RIGID:
    m_dlgRegistration->SLT_ManualMoveByDCMPlan();
    m_dlgRegistration
        ->SLT_ConfirmManualRegistration(); // skin cropping for
                                           // CBCT. only works when
                                           // CBCT_skin crop is on
    strSuffix = strSuffix + "_man";

    if (bCBCT_IntensityShift) {
      m_dlgRegistration->SLT_IntensityNormCBCT();
    }

    break;
  case REGISTER_AUTO_RIGID:
    m_dlgRegistration->SLT_ManualMoveByDCMPlan();
    m_dlgRegistration
        ->SLT_ConfirmManualRegistration(); // skin cropping

    if (bCBCT_IntensityShift) {
      m_dlgRegistration->SLT_IntensityNormCBCT();
    }

    // OPtional
    if (bFOVCropping) {
      m_cbctrecon->CropFOV3D(m_cbctrecon->m_spManualRigidCT, physPosX, physPosY,
                             physRadius, physTablePosY);
    }

    m_dlgRegistration->SLT_DoRegistrationRigid();
    strSuffix = strSuffix + "_rigid";
    break;
  case REGISTER_DEFORM_FINAL:
    std::cout << "REGISTER_DEFORM_FINAL was chosen." << std::endl;
    m_dlgRegistration->SLT_ManualMoveByDCMPlan();
    m_dlgRegistration
        ->SLT_ConfirmManualRegistration(); // skin cropping

    if (bCBCT_IntensityShift) {
      std::cout << "IntensityShift is underway" << std::endl;
      m_dlgRegistration->SLT_IntensityNormCBCT();
    }

    if (bFOVCropping) {
      m_cbctrecon->CropFOV3D(m_cbctrecon->m_spManualRigidCT, physPosX, physPosY,
                             physRadius, physTablePosY);
    }

    m_dlgRegistration->SLT_DoRegistrationRigid();
    m_dlgRegistration->SLT_DoRegistrationDeform();
    strSuffix = strSuffix + "_defrm";
    break;

  case REGISTER_DEFORM_SKIP_AUTORIGID:
    m_dlgRegistration->SLT_ManualMoveByDCMPlan();
    m_dlgRegistration
        ->SLT_ConfirmManualRegistration(); // skin cropping
    // m_pDlgRegistration->SLT_DoRegistrationRigid();

    if (bCBCT_IntensityShift) {
      m_dlgRegistration->SLT_IntensityNormCBCT();
    }

    if (bFOVCropping) {
      m_cbctrecon->CropFOV3D(m_cbctrecon->m_spManualRigidCT, physPosX, physPosY,
                             physRadius, physTablePosY);
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
    float curResampleF = m_cbctrecon->m_fResampleF;
    ui.lineEdit_DownResolFactor->setText("1.0");
    SLT_LoadSelectedProjFiles();
    m_cbctrecon->m_fResampleF = curResampleF;
    ui.lineEdit_DownResolFactor->setText(
        QString("%1").arg(m_cbctrecon->m_fResampleF));

    strSuffix = strSuffix + "_HD";
  } else {
    // strSuffix = strSuffix + "_HD"
  }
  SLT_DoScatterCorrection_APRIORI();

  // If there is couch shift information and this cbct is a pre-treatment CBCT,
  // couch shift can be applied  to represent the final treatment position

  if (ui.checkBox_CouchShiftAddToMacro->isChecked()) {
    SLT_DoCouchCorrection();
  }

  // 1) Save the corrCBCT image as signed short
  QString outputPath_rawCBCT = outputDirPath + "/" + m_cbctrecon->m_strDCMUID +
                               strSuffix + "_rawCBCT.mha";
  QString outputPath_corrCBCT = outputDirPath + "/" + m_cbctrecon->m_strDCMUID +
                                strSuffix + "_corrCBCT.mha";
  QString outputPath_manCT =
      outputDirPath + "/" + m_cbctrecon->m_strDCMUID + strSuffix + "_manCT.mha";
  QString outputPath_rigidCT = outputDirPath + "/" + m_cbctrecon->m_strDCMUID +
                               strSuffix + "_rigidCT.mha";
  QString outputPath_deformCT = outputDirPath + "/" + m_cbctrecon->m_strDCMUID +
                                strSuffix + "_deformCT.mha";

  if (bExportImages) {
    m_cbctrecon->ExportReconSHORT_HU(m_cbctrecon->m_spRawReconImg,
                                     outputPath_rawCBCT);
    m_cbctrecon->ExportReconSHORT_HU(m_cbctrecon->m_spScatCorrReconImg,
                                     outputPath_corrCBCT);
    m_cbctrecon->ExportReconSHORT_HU(m_cbctrecon->m_spManualRigidCT,
                                     outputPath_manCT);
    m_cbctrecon->ExportReconSHORT_HU(m_cbctrecon->m_spAutoRigidCT,
                                     outputPath_rigidCT);
    m_cbctrecon->ExportReconSHORT_HU(m_cbctrecon->m_spDeformedCT_Final,
                                     outputPath_deformCT);
  }

  // 2) Calculate batched WEPL points
  if (!m_cbctrecon->m_vPOI_DCM.empty()) {
    QString outputTxtPath = outputDirPath + "/" + m_cbctrecon->m_strDCMUID +
                            strSuffix + "_WEPL.txt";

    double fAngleStart = ui.lineEdit_AngStart->text().toDouble();
    double fAngleEnd = ui.lineEdit_AngEnd->text().toDouble();
    double fAngleGap = ui.lineEdit_WEPL_AngRes->text().toDouble();

    m_cbctrecon->ExportAngularWEPL_byFile(outputTxtPath, fAngleStart, fAngleEnd,
                                          fAngleGap);
  }
  return true;
}

void CbctReconWidget::SLT_OpenPhaseData() {
  if (!m_cbctrecon->m_vPhaseFloat.empty()) {
    m_cbctrecon->m_vPhaseFloat.clear();
  }

  // Open file
  QString filePath = QFileDialog::getOpenFileName(
      this, "Open phase text", m_cbctrecon->m_strPathDirDefault,
      "Phase text file (*.txt)", nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  std::ifstream fin;
  fin.open(filePath.toLocal8Bit().constData(), std::ios::in);
  if (fin.fail()) {
    return;
  }

  ui.lineEdit_PhaseTxtPath->setText(filePath);

  char str[MAX_LINE_LENGTH];

  float tmpPhase = 0.0;

  float phaseSum = 0.0;
  int phaseCnt = 0;
  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString strLine(&str[0]);

    if (strLine.length() < 1) {
      break;
    }

    tmpPhase = strLine.toFloat();
    m_cbctrecon->m_vPhaseFloat.push_back(tmpPhase);
    phaseCnt++;
    phaseSum = phaseSum + tmpPhase;
  }
  fin.close();
  std::cout << "NumOfPhaseData[Float]= " << phaseCnt << "  Mean Phase value= "
            << phaseSum / static_cast<double>(phaseCnt) << std::endl;
}

void CbctReconWidget::SLT_Export4DCBCT() {
  if (m_cbctrecon->m_spCustomGeometry == nullptr) {
    std::cout << "Error! no Geometry information loaded yet" << std::endl;
    return;
  }

  int NumOfGanAngle = m_cbctrecon->m_spCustomGeometry->GetGantryAngles().size();
  int NumOfPhase = m_cbctrecon->m_vPhaseFloat.size();

  if (NumOfGanAngle != NumOfPhase) {
    std::cout << "Size not matched. NumOfProjection= " << NumOfGanAngle
              << " NumOfProjection= " << NumOfPhase << std::endl;
    return;
  }
  // build phase bins
  QString strPhaseTextFull = ui.lineEdit_PhaseExportString->text();
  QStringList strlistPhaseFull = strPhaseTextFull.split(";");

  int cntGroup = strlistPhaseFull.count();

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

  QString strDirForXML =
      m_cbctrecon->m_strPathDirDefault; // where xml file is located
  // QString strUID ;//P00102030P + m_strDCMUID
  QString strDirForProj = m_cbctrecon->m_strPathIMAGES;

  for (int i = 0; i < cntGroup; i++) {
    // for a single group
    QStringList strListGroup = strlistPhaseFull.at(i).split(",");
    int iPhaseCnt = strListGroup.count();
    vPhaseBinsSelected.clear();

    for (int j = 0; j < iPhaseCnt; j++) {
      vPhaseBinsSelected.push_back(strListGroup.at(j).toInt());
    }
    // m_vSelectedFileNames: full file paths of projections
    // m_spCustomGeometry: full information
    // m_vPhaseFloat: full data of phase

    // Create Dir, xml, etc
    if (!m_cbctrecon->ResortCBCTProjection(
            vPhaseBinsSelected, strDirForXML, strDirForProj,
            m_cbctrecon->m_strDCMUID, m_cbctrecon->m_vPhaseFloat,
            m_cbctrecon->m_spCustomGeometry,
            m_cbctrecon->m_vSelectedFileNames)) {
      std::cout << "Error in ResortCBCTProjection "
                << strlistPhaseFull.at(i).toLocal8Bit().constData()
                << std::endl;
      return;
    }
  }

  // mkdir
  // QString strCrntDir = m_strPathPatientDir + "/" + "IMAGES"; //current Proj
  // folder  QString strCrntDir = ui.lineEdit_HisDirPath->text();

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
  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_cbctrecon->m_strPathDirDefault,
      "Short image file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  m_cbctrecon->LoadShort3DImage(fileName, imagetype);
  UShortImageType::SizeType imgDim =
      m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize();

  ui.lineEdit_Cur3DFileName->setText(fileName);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(initVal); // DrawRecon Imge is called
  ui.radioButton_graph_recon->setChecked(true);

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
  QString strTrans = ui.lineEdit_CouchTrans->text();
  QString strRot = ui.lineEdit_CouchRot->text();

  QStringList strListTrans = strTrans.split(",");
  QStringList strListRot = strRot.split(",");

  if (strListTrans.count() != 3 || strListRot.count() != 3) {
    std::cout << "Error! No couch shift data is available!" << std::endl;
    return;
  }

  VEC3D couchShiftTrans{}, couchShiftRot{};

  couchShiftTrans.x = strListTrans.at(0).toDouble(); // mm
  couchShiftTrans.y = strListTrans.at(1).toDouble();
  couchShiftTrans.z = strListTrans.at(2).toDouble();

  couchShiftRot.x = strListRot.at(0).toDouble();
  couchShiftRot.y = strListRot.at(1).toDouble();
  couchShiftRot.z = strListRot.at(2).toDouble();
  // Images to correct:
  /*m_spRawReconImg;
  m_spScatCorrReconImg;
  m_spDeformedCT_Final;
  m_spAutoRigidCT;*/

  // not manual CT!!!

  m_cbctrecon->ImageTransformUsingCouchCorrection(
      m_cbctrecon->m_spRawReconImg, m_cbctrecon->m_spRawReconImg,
      couchShiftTrans, couchShiftRot);
  m_cbctrecon->ImageTransformUsingCouchCorrection(
      m_cbctrecon->m_spScatCorrReconImg, m_cbctrecon->m_spScatCorrReconImg,
      couchShiftTrans, couchShiftRot);
  m_cbctrecon->ImageTransformUsingCouchCorrection(
      m_cbctrecon->m_spDeformedCT_Final, m_cbctrecon->m_spDeformedCT_Final,
      couchShiftTrans, couchShiftRot);
  m_cbctrecon->ImageTransformUsingCouchCorrection(
      m_cbctrecon->m_spAutoRigidCT, m_cbctrecon->m_spAutoRigidCT,
      couchShiftTrans, couchShiftRot);

  m_dlgRegistration->UpdateListOfComboBox(0); // combo selection
                                                            // signalis called
  m_dlgRegistration->UpdateListOfComboBox(1);
  m_dlgRegistration->SelectComboExternal(
      0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  m_dlgRegistration->SelectComboExternal(1, REGISTER_COR_CBCT);

  m_cbctrecon->m_spCrntReconImg = m_cbctrecon->m_spScatCorrReconImg;
  SLT_DrawReconImage();

  std::cout << "Couch shift and rotation was successfully applied."
            << std::endl;
}

// Multiple mha files
void CbctReconWidget::SLTM_WELPCalcMultipleFiles() {
  // Singed short
  QStringList listFilePath = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open",
      m_cbctrecon->m_strPathDirDefault, "signed short 3D images (*.mha)");

  int iCntFiles = listFilePath.count();
  if (iCntFiles < 1) {
    return;
  }

  int iCntPOI = m_cbctrecon->m_vPOI_DCM.size();

  if (iCntPOI < 1) {
    std::cout << "There is no POI file loaded." << std::endl;
    SLT_LoadPOIData();
  }
  iCntPOI = m_cbctrecon->m_vPOI_DCM.size();
  if (iCntPOI < 1) {
    std::cout << "Error! still no POI" << std::endl;
    return;
  }

  QString strPathOutText = QFileDialog::getSaveFileName(
      this, "File path to save", m_cbctrecon->m_strPathDirDefault,
      "WEPL_value (*.txt)", nullptr, nullptr); // Filename don't need to exist
  if (strPathOutText.length() <= 1) {
    return;
  }

  std::vector<std::vector<WEPLData>> vArrOutputWEPL;
  vArrOutputWEPL.resize(iCntFiles);

  double fAngleStart = ui.lineEdit_AngStart->text().toDouble();
  double fAngleEnd = ui.lineEdit_AngEnd->text().toDouble();
  for (int i = 0; i < iCntFiles; i++) {
    m_cbctrecon->GetWEPLDataFromSingleFile(
        listFilePath.at(i), m_cbctrecon->m_vPOI_DCM, vArrOutputWEPL.at(i),
        fAngleStart, fAngleEnd);
  }

  std::ofstream fout;
  fout.open(strPathOutText.toLocal8Bit().constData());

  fout << "Point Index"
       << "\t"
       << "Gantry Angle"
       << "\t"
       << "Sample Number";

  for (int i = 0; i < iCntFiles; i++) {
    QFileInfo fInfo(listFilePath.at(i));
    QString strFileName = fInfo.fileName();

    fout << "\t" << strFileName.toLocal8Bit().constData();
  }
  fout << std::endl;

  int cntWEPL = vArrOutputWEPL.at(0).size();
  int curCount = 0;
  for (int i = 0; i < iCntFiles; i++) {
    curCount = vArrOutputWEPL.at(i).size();
    if (cntWEPL != curCount) {
      std::cout << "Error! some of the WEPL count doesn't match!" << std::endl;
      return;
    }
  }

  for (int i = 0; i < cntWEPL; i++) {
    fout << vArrOutputWEPL.at(0).at(i).ptIndex << "\t"
         << vArrOutputWEPL.at(0).at(i).fGanAngle << "\t" << i;

    for (int j = 0; j < iCntFiles; j++) {
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
  double scaMedian = ui.lineEdit_scaMedian->text().toDouble();
  double scaGaussian = ui.lineEdit_scaGaussian->text().toDouble();
  int postScatMedianSize = ui.lineEdit_scaPostMedian->text().toInt();

  m_cbctrecon->ScatterCorPerProjRef(
      scaMedian, scaGaussian, postScatMedianSize,
      ui.radioButton_UseCUDA->isChecked(),
      ui.checkBox_ExportVolDICOM->isChecked()); // load text file
}

void CbctReconWidget::SLTM_LoadPerProjRefList() {
  QString filePath = QFileDialog::getOpenFileName(
      this, "PerProjVol list", m_cbctrecon->m_strPathDirDefault,
      "File path list (*.txt)", nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  std::ifstream fin;
  fin.open(filePath.toLocal8Bit().constData(), std::ios::in);
  if (fin.fail()) {
    return;
  }

  m_cbctrecon->m_strListPerProjRefVol.clear();

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

    QStringList strList = strLine.split('\t'); // tab

    if (strList.count() > 0 && strList.count() < 4) {
      std::cout << "Str = " << strLine.toLocal8Bit().constData() << std::endl;
      std::cout << "abnormal file expression." << std::endl;
      break;
    }

    m_cbctrecon->m_strListPerProjRefVol.push_back(strList.at(3));
  }

  std::cout << m_cbctrecon->m_strListPerProjRefVol.count()
            << " image paths were found" << std::endl;

  fin.close();
}

void CbctReconWidget::SLTM_CropMaskBatch() {
  // Specify mask file (USHORT)
  QString maskFilePath = QFileDialog::getOpenFileName(
      this, "Mask image (Ushort)", m_cbctrecon->m_strPathDirDefault,
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
  QStringList targetFilePaths = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open",
      m_cbctrecon->m_strPathDirDefault, "target files (*.mha)");

  int iCnt = targetFilePaths.size();

  if (iCnt < 1) {
    return;
  }

  for (int i = 0; i < iCnt; i++) {
    const QString &curPath = targetFilePaths.at(i);

    // Overritting
    Mask_operation mask_option = MASK_OPERATION_MASK;
    QString input_fn = curPath.toLocal8Bit().constData();
    // QString mask_fn = strPath_mskSkinCT_final.toLocal8Bit().constData();
    QString mask_fn = maskFilePath.toLocal8Bit().constData();
    QString output_fn = curPath.toLocal8Bit().constData();
    float mask_value = -1024.0; // unsigned short
    m_cbctregistration->plm_mask_main(
        mask_option, input_fn, mask_fn, output_fn, mask_value);
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
  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    return;
  }

  int bRaw = 0;

  if (m_cbctrecon->m_spCrntReconImg == m_cbctrecon->m_spRawReconImg) {
    bRaw = 1;
  }

  float dcmPosCutSup = ui.lineEdit_SupCutPos->text().toFloat(); // mm
  float dcmPosCutInf = ui.lineEdit_InfCutPos->text().toFloat(); // mm

  // CropFOV3D(m_spCrntReconImg, physPosX, physPosY, physRadius, physTablePosY);
  m_cbctrecon->CropSupInf(m_cbctrecon->m_spCrntReconImg, dcmPosCutInf,
                          dcmPosCutSup);
  // QString strPath = m_strPathDirDefault + "/" + "TempSI_Cropped.mha";
  // QString strTmpFile = "C:/TmpSI_Cropped.mha";

  QString strPath = m_cbctregistration->m_strPathPlastimatch +
                    "/" + "tmp_SI_cropped.mha";
  m_cbctrecon->ExportReconSHORT_HU(m_cbctrecon->m_spCrntReconImg, strPath);

  QString strName = "SI_Cropped";
  if (bRaw != 0) {
    if (!m_cbctrecon->LoadShortImageToUshort(strPath,
                                             m_cbctrecon->m_spRawReconImg)) {
      std::cout << "error! in LoadShortImageToUshort" << std::endl;
    }
    UpdateReconImage(m_cbctrecon->m_spRawReconImg, strName);
  } else {
    if (!m_cbctrecon->LoadShortImageToUshort(strPath,
                                             m_cbctrecon->m_spRefCTImg)) {
      std::cout << "error! in LoadShortImageToUshort" << std::endl;
    }
    UpdateReconImage(m_cbctrecon->m_spRefCTImg, strName);
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
  ReaderType::Pointer reader = ReaderType::New();

  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_cbctrecon->m_strPathDirDefault,
      "3D dose float file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  reader->SetFileName(fileName.toLocal8Bit().constData());
  reader->Update();

  // Multiply: Gy to mGy
  using MultiplyImageFilterType =
      itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  MultiplyImageFilterType::Pointer multiplyImageFilter =
      MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(reader->GetOutput());
  multiplyImageFilter->SetConstant(100.0); // calculated already //Gy to cGy

  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(multiplyImageFilter->GetOutput());

  castFilter->Update();

  m_cbctrecon->m_spRawReconImg = castFilter->GetOutput();
  m_cbctrecon->m_spCrntReconImg = m_cbctrecon->m_spRawReconImg;

  // Update UI
  UShortImageType::SizeType imgDim =
      m_cbctrecon->m_spRawReconImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing =
      m_cbctrecon->m_spRawReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  ui.lineEdit_Cur3DFileName->setText(fileName);

  m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_Load3DImage() // mha reconstructed file, from external
                                        // source
{
  using ReaderType = itk::ImageFileReader<UShortImageType>;
  ReaderType::Pointer reader = ReaderType::New();

  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_cbctrecon->m_strPathDirDefault,
      "Projection file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  reader->SetFileName(fileName.toLocal8Bit().constData());
  reader->Update();

  m_cbctrecon->m_spRawReconImg = reader->GetOutput();
  m_cbctrecon->m_spCrntReconImg = m_cbctrecon->m_spRawReconImg;

  UShortImageType::SizeType imgDim =
      m_cbctrecon->m_spRawReconImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing =
      m_cbctrecon->m_spRawReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  ui.lineEdit_Cur3DFileName->setText(fileName);

  m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_Load3DImageShort() {
  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_cbctrecon->m_strPathDirDefault,
      "short mha file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  if (!m_cbctrecon->LoadShortImageToUshort(fileName,
                                           m_cbctrecon->m_spRawReconImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  ImageCalculatorFilterType2::Pointer imageCalculatorFilter2 =
      ImageCalculatorFilterType2::New();
  // imageCalculatorFilter2->SetImage(m_spReconImg);
  imageCalculatorFilter2->SetImage(m_cbctrecon->m_spRawReconImg);
  imageCalculatorFilter2->Compute();

  auto minVal2 = static_cast<double>(imageCalculatorFilter2->GetMinimum());
  auto maxVal2 = static_cast<double>(imageCalculatorFilter2->GetMaximum());
  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  m_cbctrecon->m_spCrntReconImg = m_cbctrecon->m_spRawReconImg;

  UShortImageType::SizeType imgDim =
      m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing =
      m_cbctrecon->m_spCrntReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  ui.lineEdit_Cur3DFileName->setText(fileName);

  m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_LoadNKIImage() {
  QString filePath = QFileDialog::getOpenFileName(
      this, "Open Image", m_cbctrecon->m_strPathDirDefault, "NKI file (*.SCAN)",
      nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  Volume *v =
      nki_load(filePath.toLocal8Bit().constData()); // NKI is unsigned short!!!
  if (v == nullptr) {
    std::cerr << "file reading error" << std::endl;
    return;
  }

  QString endFix = "_conv";
  QFileInfo srcFileInfo = QFileInfo(filePath);
  QDir dir = srcFileInfo.absoluteDir();
  QString baseName = srcFileInfo.completeBaseName();
  QString extName = "mha";

  QString newFileName = baseName.append(endFix).append(".").append(extName);
  QString newPath = dir.absolutePath() + "/" + newFileName;

  write_mha(newPath.toLocal8Bit().constData(), v);
  std::cout << "File conversion is done. Trying to read mha file.."
            << std::endl;
  // corrImg.ReleaseBuffer();
  // NKI to mha

  if (!m_cbctrecon->LoadShortImageToUshort(newPath,
                                           m_cbctrecon->m_spRawReconImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  ImageCalculatorFilterType2::Pointer imageCalculatorFilter2 =
      ImageCalculatorFilterType2::New();
  // imageCalculatorFilter2->SetImage(m_spReconImg);
  imageCalculatorFilter2->SetImage(m_cbctrecon->m_spRawReconImg);
  imageCalculatorFilter2->Compute();

  auto minVal2 = static_cast<double>(imageCalculatorFilter2->GetMinimum());
  auto maxVal2 = static_cast<double>(imageCalculatorFilter2->GetMaximum());
  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  m_cbctrecon->m_spCrntReconImg = m_cbctrecon->m_spRawReconImg;

  UShortImageType::SizeType imgDim =
      m_cbctrecon->m_spCrntReconImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing =
      m_cbctrecon->m_spCrntReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  ui.lineEdit_Cur3DFileName->setText(newFileName);

  m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

void CbctReconWidget::SLT_ExportHis() {
  if (m_cbctrecon->m_iImgCnt < 1) {
    std::cout << "Error: Load raw his images first" << std::endl;
    return;
  }

  // Get Folder Name!

  // For displaying Dir only..
  QString dir = QFileDialog::getExistingDirectory(
      this, "Open Directory", "/home",
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  // FileName should be same, only selected folder

  for (int i = 0; i < m_cbctrecon->m_iImgCnt; i++) {
    QFileInfo tmpInfo = QFileInfo(m_cbctrecon->m_arrYKImage[i].m_strFilePath);
    QString newPath = dir + "/" + tmpInfo.fileName();
    m_cbctrecon->m_arrYKImage[i].SaveDataAsHis(
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
  ImageReaderType::Pointer ImageReader = ImageReaderType::New();
  ImageReader->SetFileName(projFile.fileName().toStdString());
  ImageReader->Update();
  m_cbctrecon->m_spProjImg3DFloat = ImageReader->GetOutput();

  // Copied from SLT_LoadSelectedFiles:

  if (m_cbctrecon->m_fResampleF != 1.0) {
    m_cbctrecon->ResampleItkImage(
        m_cbctrecon->m_spProjImg3DFloat, m_cbctrecon->m_spProjImg3DFloat,
        m_cbctrecon->m_fResampleF); // was! BROKEN AF for .his where input size
                                    // != 1024 (tested with 1016) -> outputs
                                    // offset -inputoffset/refactor^2 and 4
                                    // pixels too few in x and y
  }

  if (m_cbctrecon->m_projFormat == HND_FORMAT) { // -> hnd
    std::cout << "Fitted bowtie-filter correction ongoing..." << std::endl;
    SLT_DoBowtieCorrection();
  }

  m_cbctrecon->ConvertLineInt2Intensity(
      m_cbctrecon->m_spProjImg3DFloat, m_cbctrecon->m_spProjImgRaw3D,
      65535); // if X not 1024 == input size: out_offset =
              // in_offset + (1024*res_f -
              // X*res_f)*out_spacing     <- will still
              // break down at fw_projection

  FloatImageType::PointType originPt =
      m_cbctrecon->m_spProjImg3DFloat->GetOrigin();
  FloatImageType::SizeType FloatImgSize =
      m_cbctrecon->m_spProjImg3DFloat->GetBufferedRegion().GetSize();
  FloatImageType::SpacingType FloatImgSpacing =
      m_cbctrecon->m_spProjImg3DFloat->GetSpacing();

  std::cout << "YKDEBUG: Origin" << originPt[0] << ", " << originPt[1] << ", "
            << originPt[2] << std::endl;
  std::cout << "YKDEBUG: Size" << FloatImgSize[0] << ", " << FloatImgSize[1]
            << ", " << FloatImgSize[2] << std::endl;
  std::cout << "YKDEBUG: Spacing" << FloatImgSpacing[0] << ", "
            << FloatImgSpacing[1] << ", " << FloatImgSpacing[2] << std::endl;

  std::cout << "Raw3DProj dimension "
            << m_cbctrecon->m_spProjImgRaw3D->GetRequestedRegion().GetSize()
            << std::endl;

  // m_spProjImgRaw3D is Ushort

  std::cout << "Projection reading succeeded."
            << m_cbctrecon->m_vSelectedFileNames.size() << " files were read"
            << std::endl;

  // Because you can load projections from previous run:
  ui.pushButton_DoRecon->setEnabled(true);

  ui.spinBoxImgIdx->setMinimum(0);
  ui.spinBoxImgIdx->setMaximum(m_cbctrecon->m_vSelectedFileNames.size() - 1);
  ui.spinBoxImgIdx->setValue(0);

  m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();

  this->SLT_DrawProjImages(); // Update Table is called
}

void CbctReconWidget::SLT_LoadPlanCT_mha() // m_spRecon -->m_spRefCT
{
  // typedef itk::ImageFileReader<ShortImageType> ReaderType;
  // ReaderType::Pointer reader = ReaderType::New();

  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_cbctrecon->m_strPathDirDefault,
      "Plan CT file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  if (!m_cbctrecon->LoadShortImageToUshort(fileName,
                                           m_cbctrecon->m_spRefCTImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  ImageCalculatorFilterType2::Pointer imageCalculatorFilter2 =
      ImageCalculatorFilterType2::New();
  imageCalculatorFilter2->SetImage(m_cbctrecon->m_spRefCTImg);
  imageCalculatorFilter2->Compute();

  auto minVal2 = static_cast<double>(imageCalculatorFilter2->GetMinimum());
  auto maxVal2 = static_cast<double>(imageCalculatorFilter2->GetMaximum());

  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  UShortImageType::SizeType imgDim =
      m_cbctrecon->m_spRefCTImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing =
      m_cbctrecon->m_spRefCTImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT, REGISTER_MANUAL_RIGID);
  m_cbctrecon->m_spCrntReconImg = m_cbctrecon->m_spRefCTImg;

  ui.lineEdit_Cur3DFileName->setText(fileName);
  m_cbctrecon->m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(initVal); // DrawRecon Imge is called
  ui.radioButton_graph_recon->setChecked(true);
}

void CbctReconWidget::SLT_ExportReconUSHORT() {
  if (m_cbctrecon->m_spCrntReconImg == nullptr) {
    std::cout << " no image to export" << std::endl;
    return;
  }

  QString strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "unsigned short meta image (*.mha)", nullptr,
      nullptr);
  if (strPath.length() <= 1) {
    return;
  }

  using WriterType = itk::ImageFileWriter<UShortImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(strPath.toLocal8Bit().constData());
  writer->SetUseCompression(true); // not exist in original code (rtkfdk)
  writer->SetInput(m_cbctrecon->m_spCrntReconImg);

  std::cout << "Writing is under progress...: "
            << strPath.toLocal8Bit().constData() << std::endl;
  writer->Update();
  std::cout << "Writing was successfully done" << std::endl;

  QString msgStr = QString("USHORT File Writing was successfully done");
  QMessageBox::information(this, "Procedure Done", msgStr);
}

// Function for independent projection his images
void CbctReconWidget::LoadRawHisImages() {

  QStringList files =
      QFileDialog::getOpenFileNames(this, "Select one or more files to open",
                                    m_cbctrecon->m_strPathDirDefault,
                                    "projection images (*.his,*.hnd,*.xim)");

  m_cbctrecon->m_iImgCnt = files.size();
  std::vector<std::string> fileVector;

  for (auto &cur_file : files) {
    fileVector.push_back(cur_file.toStdString());
  }

  if (m_cbctrecon->m_iImgCnt < 1) {
    return;
  }

  m_cbctrecon->ReleaseMemory();

  m_cbctrecon->m_arrYKImage.resize(m_cbctrecon->m_iImgCnt);
  using ReaderType = rtk::ProjectionsReader<FloatImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames(fileVector);
  reader->UpdateOutputInformation();
  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(reader->GetOutput());
  castFilter->Update();

  int width =
      castFilter->GetOutput()->GetLargestPossibleRegion().GetSize()[0]; // width
  int height = castFilter->GetOutput()
                   ->GetLargestPossibleRegion()
                   .GetSize()[1]; // height
  int sizePix = width * height;
  int sizeBuf = sizePix * sizeof(FloatImageType::PixelType);
  int bytesPerPix = qRound(sizeBuf / static_cast<double>(sizePix));

  size_t index = 0;
  for (auto &it : m_cbctrecon->m_arrYKImage) {
    const QString &strFile = files.at(index);

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

  m_cbctrecon->m_multiplyFactor = 1.0;

  ui.spinBoxImgIdx->setMinimum(0);
  ui.spinBoxImgIdx->setMaximum(m_cbctrecon->m_iImgCnt - 1);
  ui.spinBoxImgIdx->setValue(0);

  m_cbctrecon->SetMaxAndMinValueOfProjectionImage();
  SLT_InitializeGraphLim();

  SLT_DrawRawImages(); // Change FileName as well.. read spinbox value and draw
                       // image
}
// IO_END

void CbctReconWidget::SLT_SaveCurrentSetting() {
  if (!SaveCurrentSetting(m_cbctrecon->m_strPathDefaultConfigFile)) {
    std::cout << "Error! in SaveCurrentSetting" << std::endl;
    return;
  }
}

bool CbctReconWidget::SaveCurrentSetting(QString &strPathConfigFile) {
  QFileInfo fInfo(strPathConfigFile);
  if (!fInfo.exists()) {
    std::cout << "Config file not exist. will be created now" << std::endl;
  } else {
    std::cout << "Config file is found. it will be overwritten now"
              << std::endl;
  }
  auto &m_pDlgRegistration = m_dlgRegistration;

  std::ofstream fout;
  fout.open(strPathConfigFile.toLocal8Bit().constData());

  QString strRefmAs = ui.lineEdit_RefmAs->text();
  QString PostFOV_R = ui.lineEdit_PostFOV_R->text();
  QString PostTablePosY = ui.lineEdit_PostTablePosY->text();

  QString strBkFillCT = m_pDlgRegistration->ui.lineEditBkFillCT->text();
  QString strBkDetectCT = m_pDlgRegistration->ui.lineEditBkDetectCT->text();
  QString strBubFillCT = m_pDlgRegistration->ui.lineEditBubFillCT->text();

  QString strBkFillCBCT = m_pDlgRegistration->ui.lineEditBkFillCBCT->text();
  QString strBkDetectCBCT =
      m_pDlgRegistration->ui.lineEditBubDetectCBCT->text();
  QString strBubFillCBCT = m_pDlgRegistration->ui.lineEditBubFillCBCT->text();

  QString strCropContourName =
      m_pDlgRegistration->ui.lineEditCropContourName->text();

  QString strFOVPos = m_pDlgRegistration->ui.lineEditFOVPos->text();

  QString strArgument1 = m_pDlgRegistration->ui.lineEditArgument1->text();
  QString strArgument2 = m_pDlgRegistration->ui.lineEditArgument2->text();
  QString strArgument3 = m_pDlgRegistration->ui.lineEditArgument3->text();

  bool bExportFwd = ui.checkBox_ExportFwd->isChecked();
  bool bExportScat = ui.checkBox_ExportScat->isChecked();
  bool bExportCor = ui.checkBox_ExportCor->isChecked();
  bool bExportVolDICOM = ui.checkBox_ExportVolDICOM->isChecked();
  bool bCouchShiftAddToMacro = ui.checkBox_CouchShiftAddToMacro->isChecked();

  // From Registration GUI
  bool bCropBkgroundCT =
      m_pDlgRegistration->ui.checkBoxCropBkgroundCT->isChecked();
  bool bCropBkgroundCBCT =
      m_pDlgRegistration->ui.checkBoxCropBkgroundCBCT->isChecked();

  bool bFillBubbleCT = m_pDlgRegistration->ui.checkBoxFillBubbleCT->isChecked();
  bool bFillBubbleCBCT =
      m_pDlgRegistration->ui.checkBoxFillBubbleCBCT->isChecked();

  bool bUseROIForRigid =
      m_pDlgRegistration->ui.checkBoxUseROIForRigid->isChecked();
  bool bUseROIForDIR = m_pDlgRegistration->ui.checkBoxUseROIForDIR->isChecked();

  bool bRadioButton_mse = m_pDlgRegistration->ui.radioButton_mse->isChecked();
  bool bRadioButton_mi = m_pDlgRegistration->ui.radioButton_mi->isChecked();

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

bool CbctReconWidget::LoadCurrentSetting(QString &strPathConfigFile) {
  QFileInfo fInfo(strPathConfigFile);

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
  auto &m_pDlgRegistration = m_dlgRegistration;

  // int cnt = 0;
  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH); // Read out header
    QString tmpStr = QString(&str[0]);
    QStringList strList = tmpStr.split("\t"); // tab

    QString strHeader, strContent;
    bool bFlagContent = false;
    if (strList.count() == 2) {
      strHeader = strList.at(0);
      strContent = strList.at(1);

      bFlagContent = strContent.toInt() == 1;

      if (strHeader == "strRefmAs") {
        ui.lineEdit_RefmAs->setText(strContent);
      } else if (strHeader == "PostFOV_R") {
        ui.lineEdit_PostFOV_R->setText(strContent);
      } else if (strHeader == "PostTablePosY") {
        ui.lineEdit_PostTablePosY->setText(strContent);

      } else if (strHeader == "strBkFillCT") {
        m_pDlgRegistration->ui.lineEditBkFillCT->setText(strContent);
      } else if (strHeader == "strBkDetectCT") {
        m_pDlgRegistration->ui.lineEditBkDetectCT->setText(strContent);
      } else if (strHeader == "strBubFillCT") {
        m_pDlgRegistration->ui.lineEditBubFillCT->setText(strContent);

      } else if (strHeader == "strBkFillCBCT") {
        m_pDlgRegistration->ui.lineEditBkFillCBCT->setText(strContent);
      } else if (strHeader == "strBkDetectCBCT") {
        m_pDlgRegistration->ui.lineEditBubDetectCBCT->setText(strContent);
      } else if (strHeader == "strBubFillCBCT") {
        m_pDlgRegistration->ui.lineEditBubFillCBCT->setText(strContent);
      }

      if (strHeader == "strCropContourName") {
        m_pDlgRegistration->ui.lineEditCropContourName->setText(strContent);
      } else if (strHeader == "strFOVPos") {
        m_pDlgRegistration->ui.lineEditFOVPos->setText(strContent);

      } else if (strHeader == "strArgument1") {
        m_pDlgRegistration->ui.lineEditArgument1->setText(strContent);
      } else if (strHeader == "strArgument2") {
        m_pDlgRegistration->ui.lineEditArgument2->setText(strContent);
      } else if (strHeader == "strArgument3") {
        m_pDlgRegistration->ui.lineEditArgument3->setText(strContent);

      } else if (strHeader == "bExportFwd") {
        ui.checkBox_ExportFwd->setChecked(bFlagContent);
      } else if (strHeader == "bExportScat") {
        ui.checkBox_ExportScat->setChecked(bFlagContent);
      } else if (strHeader == "bExportCor") {
        ui.checkBox_ExportCor->setChecked(bFlagContent);

      } else if (strHeader == "bExportVolDICOM") {
        ui.checkBox_ExportVolDICOM->setChecked(bFlagContent);
      } else if (strHeader == "bCouchShiftAddToMacro") {
        ui.checkBox_CouchShiftAddToMacro->setChecked(bFlagContent);

      } else if (strHeader == "bCropBkgroundCT") {
        m_pDlgRegistration->ui.checkBoxCropBkgroundCT->setChecked(bFlagContent);
      } else if (strHeader == "bCropBkgroundCBCT") {
        m_pDlgRegistration->ui.checkBoxCropBkgroundCBCT->setChecked(
            bFlagContent);

      } else if (strHeader == "bFillBubbleCT") {
        m_pDlgRegistration->ui.checkBoxFillBubbleCT->setChecked(bFlagContent);
      } else if (strHeader == "bFillBubbleCBCT") {
        m_pDlgRegistration->ui.checkBoxFillBubbleCBCT->setChecked(bFlagContent);

      } else if (strHeader == "bUseROIForRigid") {
        m_pDlgRegistration->ui.checkBoxUseROIForRigid->setChecked(bFlagContent);
      } else if (strHeader == "bUseROIForDIR") {
        m_pDlgRegistration->ui.checkBoxUseROIForDIR->setChecked(bFlagContent);

      } else if (strHeader == "radioButton_mse") {
        m_pDlgRegistration->ui.radioButton_mse->setChecked(bFlagContent);
      } else if (strHeader == "radioButton_mi") {
        m_pDlgRegistration->ui.radioButton_mi->setChecked(bFlagContent);
      }
    }
  }
  fin.close();

  return true;
}
