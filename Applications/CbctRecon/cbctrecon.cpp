#include "cbctrecon.h"
#define USE_AVX false
#if USE_AVX
#include <immintrin.h>
#endif

#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif

#include <cstdio>

#ifdef USE_OPENMP
#include <omp.h>
#endif

#if ITK_VERSION_MAJOR >= 4
#include "gdcmUIDGenerator.h"
#else
#include "gdcm/src/gdcmFile.h"
#include "gdcm/src/gdcmUtil.h"
#endif

// related includes for forward projection
//#include "itkMedianImageFilter.h"
#if USE_CUDA
#include "rtkConjugateGradientConeBeamReconstructionFilter.h"
#include "rtkCudaForwardProjectionImageFilter.h"
#include "rtkSARTConeBeamReconstructionFilter.h"
#endif

// Local
#ifdef USE_OPENCL_RTK
#include "rtkOpenCLFDKConeBeamReconstructionFilter.h" // local
#endif

#include "DlgExternalCommand.h"
#include "DlgRegistration.h"
#include "OpenCLFFTFilter.h"

#define round(dTemp) (long(dTemp + (dTemp > 0 ? .5 : -.5)))


#include "io.hxx"
#include "compute.hxx"
#include "fdk.hxx"

// using namespace std;

// extern "C" void YKCudaMedianWrapFloat(float* CPUinput, float* CPUoutput,
// const int width, const int height, const int fWidth, const int fHeight);

struct TIFIFD {
  unsigned short TagID;
  unsigned short DataType;
  int DataCnt;
  int DataOrOffset;
};

struct RATIONAL {
  long a;
  long b;
};

double vectorMean(const std::vector<double> &vDouble);
double vectorSum(const std::vector<double> &vDouble);

CbctRecon::CbctRecon(QWidget *parent, Qt::WindowFlags flags)
    : QMainWindow(parent, flags) {
  ui.setupUi(this);

  m_dspYKReconImage = std::make_unique<YK16GrayImage>();
  m_dspYKImgProj = std::make_unique<YK16GrayImage>();

  // Badpixmap;
  m_pImgOffset =
      std::make_unique<YK16GrayImage>(DEFAULT_ELEKTA_PROJ_WIDTH, DEFAULT_ELEKTA_PROJ_HEIGHT);
  m_pImgGain =
      std::make_unique<YK16GrayImage>(DEFAULT_ELEKTA_PROJ_WIDTH, DEFAULT_ELEKTA_PROJ_HEIGHT);
  // Prepare Raw image

  m_bScanDirectionCW = true;

  // Disable cuda & opencl as defaults
  ui.radioButton_UseCPU->setChecked(true);
  ui.radioButton_UseCUDA->setDisabled(true);
  ui.radioButton_UseOpenCL->setDisabled(true);

  ui.pushButton_DoRecon->setDisabled(true);

  m_iTmpIdx = 60;

  m_fProjImgValueMax = 0.0; // value of float image
  m_fProjImgValueMin = 0.0;

  m_arrYKBufProj.clear();
  m_iCntSelectedProj = 0;

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
  m_pTableModel = nullptr;
  m_pDlgRegistration = std::make_unique<DlgRegistration>(this);
  // m_pDlgHistogram = new DlgHistogram(this);
  m_pDlgExternalCommand = std::make_unique<DlgExternalCommand>(this);

  m_structures = std::make_unique<StructureSet>();

  m_iFixedOffset_ScatterMap = 10000; // fixed! allows negative value of scatter
  // m_iFixedOffset_ScatterMap = 0;//fixed! allows negative value of scatter
  m_fResampleF = 1.0;
  m_fProjSpacingX = 0.4; // DEFAULT, will be updated during Load Proj selected
  m_fProjSpacingY = 0.4;

  // 20141017 QTIMER for sync
  m_Timer = new QTimer(this);
  connect(m_Timer, SIGNAL(timeout()), this, SLOT(SLT_TimerEvent()));
  m_busyTimer = false;

  m_strPathDirDefault = QDir::currentPath();
  std::cout << "Current Default Dir: "
            << m_strPathDirDefault.toLocal8Bit().constData() << std::endl;
  QString tmp_folder("tmp");
  init_DlgRegistration(tmp_folder); // to Setup plastimatch folder. this is
                                    // useful if registration will be only done

  // shell test
  // QString strCurFolder = "H:\\lib\\rtk\\NightlyBUILD64\\bin\\Release";
  // QString strCommand = QString("explorer %1").arg(strCurFolder);
  // std::cout << strCommand.toLocal8Bit().constData() << std::endl;
  //::system(strCommand.toLocal8Bit().constData());
  //::system("rtkfdk");

  //	if
  //(QProcess::execute(QString("H:\\lib\\rtk\\NightlyBUILD64\\bin\\Release\\rtkfdk"))
  //< 0) 	qDebug() << "Failed to run";
  m_strPathDirDefault =
      R"(D:\Program_data\01_20140827_CBCT_All\04_patient_phan_pelvis_M\IMAGES\img_1.3.46.423632.135786.1409186054.9_M20mAs6440)";

  QString strPathCurAppDir = QDir::currentPath(); // should be same as .exe file
  std::cout << "Current app path= "
            << strPathCurAppDir.toLocal8Bit().constData() << std::endl;
  m_strPathDefaultConfigFile = strPathCurAppDir + "/" + "DefaultConfig.cfg";

  if (!LoadCurrentSetting(m_strPathDefaultConfigFile)) // Update GUI
  {
    std::cout << "DefaultConfig.cfg is not found in the application folder. A "
                 "new one will be created"
              << std::endl;
    if (!SaveCurrentSetting(m_strPathDefaultConfigFile)) {
      std::cout << "Error in SaveCurrentSetting" << std::endl;
    }
  }
  m_bMacroContinue = true;
}

CbctRecon::~CbctRecon() {
  ReleaseMemory();
  delete m_Timer;
}

void CbctRecon::ReleaseMemory() {
  if (!m_arrYKImage.empty()) {
    // m_iImgCnt = 0;
    m_arrYKImage.clear();
  }

  if (m_pTableModel != nullptr) {
    m_pTableModel = nullptr;
  }

  if (!m_arrYKBufProj.empty()) {
    m_arrYKBufProj.clear();
    m_iCntSelectedProj = 0;
  }
}
void CbctRecon::SLT_LoadRawImages() {
  LoadRawHisImages();
}

// Hexa name ->decimal name

void CbctRecon::RenameFromHexToDecimal(QStringList &filenameList) {
  int size = filenameList.size();

  QString crntFilePath;

  for (int i = 0; i < size; i++) {
    crntFilePath = filenameList.at(i);
    QFileInfo fileInfo = QFileInfo(crntFilePath);
    QDir dir = fileInfo.absoluteDir();
    QString fileBase = fileInfo.baseName();
    QString newBaseName = HexStr2IntStr(fileBase);
    QString extStr = fileInfo.completeSuffix();

    QString newFileName = newBaseName.append(".").append(extStr);
    QString newPath = dir.absolutePath() + "/" + newFileName;

    // extract former part
    QFile::rename(crntFilePath, newPath);
  }
  // Extract
}

void CbctRecon::SLT_DrawRawImages() {
  int crntIdx = ui.spinBoxImgIdx->value();

  if (crntIdx >= m_iImgCnt) {
    return;
  }

  int windowMin = ui.sliderRawMin->value();
  int windowMax = ui.sliderRawMax->value();

  QFileInfo tmpInfo = QFileInfo(m_arrYKImage[crntIdx].m_strFilePath);
  ui.lineEditFileName->setText(tmpInfo.fileName());

  int width = m_arrYKImage[crntIdx].m_iWidth;
  int height = m_arrYKImage[crntIdx].m_iHeight;
  m_dspYKImgProj->CreateImage(width, height, 0);
  m_dspYKImgProj->CopyFromBuffer(m_arrYKImage[crntIdx].m_pData, width, height);

  m_dspYKImgProj->FillPixMapMinMax(windowMin, windowMax);
  ui.labelImageRaw->SetBaseImage(m_dspYKImgProj.get());
  ui.labelImageRaw->update();
}

void CbctRecon::SLT_DrawProjImages() {
  if (m_dspYKImgProj == nullptr) {
    return;
  }

  if (m_iImgCnt > 0) {
    SLT_DrawRawImages();
    //		SLT_DrawGraph();
    SLT_UpdateTable();
    return;
  }
  // Using slice iterator,
  // 1) Find the slice requested
  // 2) get dimension to create 2DYK16Image
  // 3) copy slice region to YK16 Image --> Cating: float to USHORT

  if (m_spProjImg3DFloat == nullptr) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<FloatImageType> it(
      m_spProjImg3DFloat, m_spProjImg3DFloat->GetRequestedRegion());

  FloatImageType::SizeType imgSize =
      m_spProjImg3DFloat->GetRequestedRegion().GetSize(); // 1016x1016 x z

  int width = imgSize[0];
  int height = imgSize[1];
  m_dspYKImgProj->CreateImage(width, height, 0);

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  int iNumSlice = 0;
  int iNumWidth = 0;
  int iNumHeight = 0;

  double realValGap = m_fProjImgValueMax - m_fProjImgValueMin;
  m_multiplyFactor = 0.0;

  if (realValGap > 0.0) {
    m_multiplyFactor = 65535.0 / realValGap;
  }

  int iReqSlice = ui.spinBoxImgIdx->value();

  while (!it.IsAtEnd()) {
    if (iNumSlice == iReqSlice) {
      iNumHeight = 0;
      while (!it.IsAtEndOfSlice()) {
        iNumWidth = 0;
        while (!it.IsAtEndOfLine()) {
          // double tmpVal = it.Get()*multiplyFactor;
          double tmpVal = it.Get();

          m_dspYKImgProj->m_pData[iNumWidth + width * iNumHeight] =
              static_cast<unsigned short>((tmpVal - m_fProjImgValueMin) *
                                          m_multiplyFactor);
          // it.Set() doesn't exist in the Const Iterator
          ++it;
          iNumWidth++;
        }
        it.NextLine();
        iNumHeight++;
      }
      // break;
    }
    it.NextSlice();
    iNumSlice++;
  }

  m_dspYKImgProj->FillPixMapMinMax(ui.sliderRawMin->value(),
                                   ui.sliderRawMax->value());

  ui.labelImageRaw->SetBaseImage(m_dspYKImgProj.get());
  ui.labelImageRaw->update();

  SLT_UpdateTable();
}

QString CbctRecon::HexStr2IntStr(QString &strHex) {
  int len = strHex.length();

  int tmpDecimal = 0;
  // int cnt = 0;
  int i = 0;

  for (i = len - 1; i >= 0; i--) {
    int tmpNum = 0;

    if (strHex.at(i) == 'a' || strHex.at(i) == 'A') {
      tmpNum = 10;
    } else if (strHex.at(i) == 'b' || strHex.at(i) == 'B') {
      tmpNum = 11;
    } else if (strHex.at(i) == 'c' || strHex.at(i) == 'C') {
      tmpNum = 12;
    } else if (strHex.at(i) == 'd' || strHex.at(i) == 'D') {
      tmpNum = 13;
    } else if (strHex.at(i) == 'e' || strHex.at(i) == 'E') {
      tmpNum = 14;
    } else if (strHex.at(i) == 'f' || strHex.at(i) == 'F') {
      tmpNum = 15;

    } else {
      QString tmpStr;
      tmpStr = strHex.at(i);
      tmpNum = tmpStr.toInt();
    }
    tmpDecimal = tmpDecimal + tmpNum * static_cast<int>(pow(16.0, len - 1 - i));
  }

  QString intStr = QString("%1").arg(tmpDecimal);

  return intStr;
  // return tmpDecimal;
  // m_str_10.Format("%d", tmpDecimal);
}

void CbctRecon::SLT_FileNameHex2Dec() {
  QStringList files = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open", m_strPathDirDefault,
      "projection images (*.his)");

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
    RenameFromHexToDecimal(files);
  }
}

void CbctRecon::SLT_MakeElektaXML() {
  // Define IMAGE.DBF path
  QString filePath_ImageDBF = QFileDialog::getOpenFileName(
      this, "SelectIMAGE.DBF file", m_strPathDirDefault,
      "Elekta DB file (*.dbf)", nullptr, nullptr);

  if (filePath_ImageDBF.length() < 2) {
    return;
  }

  QString filePath_FrameDBF = QFileDialog::getOpenFileName(
      this, "Select FRAME.DBF file", m_strPathDirDefault,
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

  QString genFilePath =
      MakeElektaXML(filePath_ImageDBF, filePath_FrameDBF, DICOM_UID);
  std::cout << "Generated ElektaXML path: "
            << genFilePath.toLocal8Bit().constData() << std::endl;
}


void CbctRecon::SLT_OpenOffsetFile() {
  // QString strPath = QFileDialog::getOpenFileNames(this,"Select one or more
  // files to open","/home","Images (*.raw)");
  QString strPath =
      QFileDialog::getOpenFileName(this, "Select a single file to open",
                                   m_strPathDirDefault, "raw image (*.raw)");

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_offsetPath->setText(strPath);
  if (hisIsUsed) {
    m_pImgOffset->LoadRawImage(strPath.toLocal8Bit().constData(),
                               DEFAULT_ELEKTA_PROJ_WIDTH,
                               DEFAULT_ELEKTA_PROJ_HEIGHT);
  } else {
    m_pImgOffset->LoadRawImage(strPath.toLocal8Bit().constData(),
                               DEFAULT_VARIAN_PROJ_WIDTH,
                               DEFAULT_VARIAN_PROJ_HEIGHT);
  }
}

void CbctRecon::SLT_OpenGainFile() {
  QString strPath =
      QFileDialog::getOpenFileName(this, "Select a single file to open",
                                   m_strPathDirDefault, "raw image (*.raw)");

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_gainPath->setText(strPath);
  m_pImgGain->LoadRawImage(strPath.toLocal8Bit().constData(),
                           DEFAULT_ELEKTA_PROJ_WIDTH,
                           DEFAULT_ELEKTA_PROJ_HEIGHT);
}

void CbctRecon::SLT_OpenBadpixelFile() {
  QString strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open", m_strPathDirDefault,
      "bad pixel map file (*.pmf)");

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_badpixelPath->setText(strPath);
  LoadBadPixelMap(strPath.toLocal8Bit().constData());
  // m_pImgGain->LoadRawImage(strPath.toLocal8Bit(),IMG_WIDTH, IMG_HEIGHT);
}

QString CbctRecon::CorrectSingleFile(const char *filePath) {
  // Load raw file
  YK16GrayImage rawImg(DEFAULT_ELEKTA_PROJ_WIDTH, DEFAULT_ELEKTA_PROJ_HEIGHT);
  rawImg.LoadRawImage(filePath, DEFAULT_ELEKTA_PROJ_WIDTH,
                      DEFAULT_ELEKTA_PROJ_HEIGHT);

  YK16GrayImage corrImg(DEFAULT_ELEKTA_PROJ_WIDTH, DEFAULT_ELEKTA_PROJ_HEIGHT);

  bool bDarkCorrApply = ui.checkBox_offsetOn->isChecked();
  bool bGainCorrApply = ui.checkBox_gainOn->isChecked();
  bool bDefectMapApply = ui.checkBox_badpixelOn->isChecked();

  // m_pParent->m_pCurrImageRaw->m_pData[i]

  if (!bDarkCorrApply && !bGainCorrApply) {
    for (int i = 0; i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT;
         i++) {
      corrImg.m_pData[i] = rawImg.m_pData[i]; // raw image
    }
  } else if (bDarkCorrApply && !bGainCorrApply) {
    if (m_pImgOffset->IsEmpty()) {
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg.m_pData[i] = rawImg.m_pData[i]; // raw image
      }
    } else {
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        if (rawImg.m_pData[i] > m_pImgOffset->m_pData[i]) {
          corrImg.m_pData[i] = rawImg.m_pData[i] - m_pImgOffset->m_pData[i];
        } else {
          corrImg.m_pData[i] = 0;
        }
      }
    }
  } else if (!bDarkCorrApply && bGainCorrApply) {
    if (m_pImgGain->IsEmpty()) {
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg.m_pData[i] = rawImg.m_pData[i]; // raw image
      }
    } else {
      // get a mean value for m_pGainImage
      double sum = 0.0;
      double MeanVal = 0.0;
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        sum = sum + m_pImgGain->m_pData[i];
      }
      MeanVal = sum / static_cast<double>(DEFAULT_ELEKTA_PROJ_WIDTH *
                                          DEFAULT_ELEKTA_PROJ_HEIGHT);

      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        if (m_pImgGain->m_pData[i] == 0) {
          corrImg.m_pData[i] = rawImg.m_pData[i];
        } else {
          corrImg.m_pData[i] = static_cast<unsigned short>(
              static_cast<double>(rawImg.m_pData[i]) /
              static_cast<double>(m_pImgGain->m_pData[i]) * MeanVal);
        }
      }
    }

    /*double CalibF = 0.0;

    for (int i = 0; i < xSize * ySize; i++) {
    CalibF =
    m_pParent->m_dlgControl->lineEditSingleCalibFactor->text().toDouble();
    pImageCorr[i] = (unsigned short)(pImageCorr[i] * CalibF);
    }*/
  }

  else if (bDarkCorrApply && bGainCorrApply) {
    bool bRawImage = false;
    if (m_pImgOffset->IsEmpty()) {
      bRawImage = true;
    }
    if (m_pImgGain->IsEmpty()) {
      bRawImage = true;
    }

    if (bRawImage) {
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg.m_pData[i] = rawImg.m_pData[i]; // raw image
      }
    } else // if not raw image
    {
      // get a mean value for m_pGainImage
      double sum = 0.0;
      double MeanVal = 0.0;
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        sum = sum + (m_pImgGain->m_pData[i] - m_pImgOffset->m_pData[i]);
      }
      MeanVal = sum / static_cast<double>(DEFAULT_ELEKTA_PROJ_WIDTH *
                                          DEFAULT_ELEKTA_PROJ_HEIGHT);

      double denom = 0.0;
      int iDenomLessZero = 0;
      int iDenomLessZero_RawIsGreaterThanDark = 0;
      int iDenomLessZero_RawIsSmallerThanDark = 0;
      int iDenomOK_RawValueMinus = 0;
      int iValOutOfRange = 0;

      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        denom = static_cast<double>(m_pImgGain->m_pData[i] -
                                    m_pImgOffset->m_pData[i]);

        if (denom <= 0) {
          iDenomLessZero++;

          if (rawImg.m_pData[i] > m_pImgOffset->m_pData[i]) {
            corrImg.m_pData[i] = rawImg.m_pData[i] - m_pImgOffset->m_pData[i];
            iDenomLessZero_RawIsGreaterThanDark++;
          } else {
            corrImg.m_pData[i] = 0;
            iDenomLessZero_RawIsSmallerThanDark++;
          }
        } else {
          double tmpVal = 0.0;
          tmpVal =
              (rawImg.m_pData[i] - m_pImgOffset->m_pData[i]) / denom * MeanVal;

          if (tmpVal < 0) {
            corrImg.m_pData[i] = 0;
            iDenomOK_RawValueMinus++;
          } else {
            if (tmpVal > 65535) { // 16bit max value
              iValOutOfRange++;
            }

            corrImg.m_pData[i] = static_cast<unsigned short>(tmpVal);
          }
        }
      } // end of for
    }   // end if not bRawImage

    // Apply Single Calib. Factor
    /*double CalibF = 0.0;

    for (int i = 0; i < xSize * ySize; i++) {
    CalibF =
    m_pParent->m_dlgControl->lineEditSingleCalibFactor->text().toDouble();
    pImageCorr[i] = (unsigned short)(pImageCorr[i] * CalibF);
    }*/

  } // else if (m_bDarkCorrApply && m_bGainCorrApply)

  /*******************Customized Thresholding after Gain Corr
   * ****************/
  // unsigned short customThreVal = ui.lineEdit_customThre->text().toDouble();
  // //13000  bool enableCustomThre = ui.checkBox_customThre->isChecked();

  /*if (enableCustomThre)
  {
  for (int i = 0; i < IMG_WIDTH* IMG_HEIGHT; i++)
  {
  if (corrImg.m_pData[i] >= customThreVal)
  corrImg.m_pData[i] = customThreVal;
  }
  }*/

  if (bDefectMapApply && !m_vPixelReplMap.empty()) // pixel replacement
  {
    BadPixReplacement(&corrImg);
  }

  // file naming & export
  // file end-fix

  // filePath
  // QString exportName = filePath;
  // corrImg.SaveDataAsRaw();
  QString endFix = "_CORR";

  QFileInfo srcFileInfo = QFileInfo(filePath);
  QDir dir = srcFileInfo.absoluteDir();
  QString baseName = srcFileInfo.baseName();
  QString extName = srcFileInfo.completeSuffix();

  QString newFileName = baseName.append(endFix).append(".").append(extName);
  QString newPath = dir.absolutePath() + "/" + newFileName;

  corrImg.SaveDataAsRaw(newPath.toLocal8Bit().constData());

  return newPath;
  // corrImg.ReleaseBuffer();
}

void CbctRecon::CorrectSingleFile(YK16GrayImage *pYKRawImg) {
  if (pYKRawImg == nullptr) {
    return;
  }

  // YK16GrayImage rawImg(IMG_WIDTH, IMG_HEIGHT);
  // rawImg.LoadRawImage(filePath, IMG_WIDTH, IMG_HEIGHT);

  YK16GrayImage corrImg(DEFAULT_ELEKTA_PROJ_WIDTH, DEFAULT_ELEKTA_PROJ_HEIGHT);

  bool bDarkCorrApply = ui.checkBox_offsetOn->isChecked();
  bool bGainCorrApply = ui.checkBox_gainOn->isChecked();
  bool bDefectMapApply = ui.checkBox_badpixelOn->isChecked();

  // m_pParent->m_pCurrImageRaw->m_pData[i]

  if (!bDarkCorrApply && !bGainCorrApply) {
    for (int i = 0; i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT;
         i++) {
      corrImg.m_pData[i] = pYKRawImg->m_pData[i]; // raw image
    }
  } else if (bDarkCorrApply && !bGainCorrApply) {
    if (m_pImgOffset->IsEmpty()) {
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg.m_pData[i] = pYKRawImg->m_pData[i]; // raw image
      }
    } else {
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        if (pYKRawImg->m_pData[i] > m_pImgOffset->m_pData[i]) {
          corrImg.m_pData[i] = pYKRawImg->m_pData[i] - m_pImgOffset->m_pData[i];
        } else {
          corrImg.m_pData[i] = 0;
        }
      }
    }
  } else if (!bDarkCorrApply && bGainCorrApply) {
    if (m_pImgGain->IsEmpty()) {
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg.m_pData[i] = pYKRawImg->m_pData[i]; // raw image
      }
    } else {
      // get a mean value for m_pGainImage
      double sum = 0.0;
      double MeanVal = 0.0;
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        sum = sum + m_pImgGain->m_pData[i];
      }
      MeanVal = sum / static_cast<double>(DEFAULT_ELEKTA_PROJ_WIDTH *
                                          DEFAULT_ELEKTA_PROJ_HEIGHT);

      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        if (m_pImgGain->m_pData[i] == 0) {
          corrImg.m_pData[i] = pYKRawImg->m_pData[i];
        } else {
          corrImg.m_pData[i] = static_cast<unsigned short>(
              static_cast<double>(pYKRawImg->m_pData[i]) /
              static_cast<double>(m_pImgGain->m_pData[i]) * MeanVal);
        }
      }
    }

    /*double CalibF = 0.0;

    for (int i = 0; i < xSize * ySize; i++) {
    CalibF =
    m_pParent->m_dlgControl->lineEditSingleCalibFactor->text().toDouble();
    pImageCorr[i] = (unsigned short)(pImageCorr[i] * CalibF);
    }*/
  }

  else if (bDarkCorrApply && bGainCorrApply) {
    bool bRawImage = false;
    if (m_pImgOffset->IsEmpty()) {
      bRawImage = true;
    }
    if (m_pImgGain->IsEmpty()) {
      bRawImage = true;
    }

    if (bRawImage) {
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg.m_pData[i] = pYKRawImg->m_pData[i]; // raw image
      }
    } else // if not raw image
    {
      // get a mean value for m_pGainImage
      double sum = 0.0;
      double MeanVal = 0.0;
      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        sum = sum + (m_pImgGain->m_pData[i] - m_pImgOffset->m_pData[i]);
      }
      MeanVal = sum / static_cast<double>(DEFAULT_ELEKTA_PROJ_WIDTH *
                                          DEFAULT_ELEKTA_PROJ_HEIGHT);

      double denom = 0.0;
      int iDenomLessZero = 0;
      int iDenomLessZero_RawIsGreaterThanDark = 0;
      int iDenomLessZero_RawIsSmallerThanDark = 0;
      int iDenomOK_RawValueMinus = 0;
      int iValOutOfRange = 0;

      for (int i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        denom = static_cast<double>(m_pImgGain->m_pData[i] -
                                    m_pImgOffset->m_pData[i]);

        if (denom <= 0) {
          iDenomLessZero++;

          if (pYKRawImg->m_pData[i] > m_pImgOffset->m_pData[i]) {
            corrImg.m_pData[i] =
                pYKRawImg->m_pData[i] - m_pImgOffset->m_pData[i];
            iDenomLessZero_RawIsGreaterThanDark++;
          } else {
            corrImg.m_pData[i] = 0;
            iDenomLessZero_RawIsSmallerThanDark++;
          }
        } else {
          double tmpVal = 0.0;
          tmpVal = (pYKRawImg->m_pData[i] - m_pImgOffset->m_pData[i]) / denom *
                   MeanVal;

          if (tmpVal < 0) {
            corrImg.m_pData[i] = 0;
            iDenomOK_RawValueMinus++;
          } else {
            if (tmpVal > 65535) { // 16bit max value
              iValOutOfRange++;
            }

            corrImg.m_pData[i] = static_cast<unsigned short>(tmpVal);
          }
        }
      } // end of for
    }   // end if not bRawImage

    // Apply Single Calib. Factor

  } // else if (m_bDarkCorrApply && m_bGainCorrApply)

  if (bDefectMapApply && !m_vPixelReplMap.empty()) // pixel replacement
  {
    BadPixReplacement(&corrImg);
  }

  // Replace old buffer with new one.
  pYKRawImg->CopyFromBuffer(corrImg.m_pData, corrImg.m_iWidth,
                            corrImg.m_iHeight);
}

void CbctRecon::LoadBadPixelMap(const char *filePath) {
  m_vPixelReplMap.clear();

  std::ifstream fin;
  fin.open(filePath);

  if (fin.fail()) {
    return;
  }

  char str[MAX_LINE_LENGTH];
  // memset(str, 0, MAX_LINE_LENGTH);

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString tmpStr = QString(&str[0]);

    if (tmpStr.contains("#ORIGINAL_X")) {
      break;
    }
  }

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString tmpStr = QString(&str[0]);

    QStringList strList = tmpStr.split("	");

    if (strList.size() == 4) {
      BADPIXELMAP tmpData{};
      tmpData.BadPixX = strList.at(0).toInt();
      tmpData.BadPixY = strList.at(1).toInt();
      tmpData.ReplPixX = strList.at(2).toInt();
      tmpData.ReplPixY = strList.at(3).toInt();
      m_vPixelReplMap.push_back(tmpData);
    }
  }
  fin.close();
}

void CbctRecon::BadPixReplacement(YK16GrayImage *targetImg) {
  if (m_vPixelReplMap.empty()) {
    return;
  }

  int oriIdx, replIdx;

  for (auto &it : m_vPixelReplMap) {
    BADPIXELMAP tmpData = it;
    oriIdx = tmpData.BadPixY * DEFAULT_ELEKTA_PROJ_WIDTH + tmpData.BadPixX;
    replIdx = tmpData.ReplPixY * DEFAULT_ELEKTA_PROJ_WIDTH + tmpData.ReplPixX;
    targetImg->m_pData[oriIdx] = targetImg->m_pData[replIdx];
  }
}

void CbctRecon::SLT_ApplyCalibration() {
  if (m_iImgCnt < 1) {
    return;
  }

  for (int i = 0; i < m_iImgCnt; i++) {
    CorrectSingleFile(&m_arrYKImage[i]); // pixel value will be changed
  }
  SLT_DrawRawImages();
}

void CbctRecon::SLT_DrawReconImage() {
  if (m_dspYKReconImage == nullptr) {
    return;
  }

  if (m_spCrntReconImg == nullptr) {
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
  duplicator->SetInputImage(m_spCrntReconImg);
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

  double originZ = m_spCrntReconImg->GetOrigin()[2];
  double spacingZ = m_spCrntReconImg->GetSpacing()[2];
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
  m_dspYKReconImage = YK16GrayImage::CopyItkImage2YKImage(
      pCrnt2D, std::move(m_dspYKReconImage)); // dimension should be same automatically.

  // m_dspYKReconImage->SaveDataAsRaw("D:\\RawFile.raw"); //410 410 OK

  PostApplyFOVDispParam();
  // SLT_UpdatePostProcDispObj();

  if (ui.checkBox_PostDispObjOn->isChecked()) {
    m_dspYKReconImage->m_bDrawFOVCircle = true;
    m_dspYKReconImage->m_bDrawTableLine = true;
  }

  else {
    m_dspYKReconImage->m_bDrawFOVCircle = false;
    m_dspYKReconImage->m_bDrawTableLine = false;
  }

  m_dspYKReconImage->FillPixMapMinMax(ui.sliderReconImgMin->value(),
                                      ui.sliderReconImgMax->value());
  ui.labelReconImage->SetBaseImage(m_dspYKReconImage.get());
  ui.labelReconImage->update();

  // SLT_DrawGraph();
  SLT_UpdateTable();
}

void CbctRecon::SLT_OpenElektaGeomFile() {
  QString strPath = QFileDialog::getOpenFileName(
      this, "Select a single file to open", m_strPathDirDefault,
      "Geometry file (*.xml)");

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_ElektaGeomPath->setText(strPath);
}

void CbctRecon::SLT_SetOutputPath() {
  QString strPath = QFileDialog::getSaveFileName(
      this, "File path to save", "D:\\", "meta 3D image data (*.mha)", nullptr,
      nullptr); // Filename don't need to exist

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_OutputFilePath->setText(strPath);
}


void CbctRecon::SLT_DoReconstruction() {
  if (ui.radioButton_UseCUDA->isChecked()) {
    CudaDoReconstructionFDK(REGISTER_RAW_CBCT);
  } else if (ui.radioButton_UseOpenCL->isChecked()) {
    OpenCLDoReconstructionFDK(REGISTER_RAW_CBCT);
  } else {
    DoReconstructionFDK(REGISTER_RAW_CBCT);
  }

  m_pDlgRegistration->UpdateListOfComboBox(0); // combo selection signalis
                                               // called
  m_pDlgRegistration->UpdateListOfComboBox(1);
  // m_pDlgRegistration->SelectComboExternal(0, REGISTER_RAW_CBCT); // will call
  // fixedImageSelected  m_pDlgRegistration->SelectComboExternal(1,
  // REGISTER_RAW_CBCT );

  // After first reconstruction, set Median size to 0 0 1 for scatter corrected
  // solution
  /* ui.lineEdit_PostMedSizeX->setText(QString("%1").arg(0.0));
  ui.lineEdit_PostMedSizeY->setText(QString("%1").arg(0.0));
  ui.lineEdit_PostMedSizeZ->setText(QString("%1").arg(1.0));*/
}

void CbctRecon::SLT_InitializeGraphLim() {
  // Set Max Min at graph
  if (ui.radioButton_graph_proj->isChecked()) {
    if (m_iImgCnt > 0) // if indep raw his images are loaded
    {
      int horLen = m_dspYKImgProj->m_iWidth;
      // int verLen = m_dspYKImgProj->m_iHeight;

      // set edit maxium min
      QString strXMin = QString("%1").arg(horLen);
      ui.lineEditXMin->setText("0");
      ui.lineEditXMax->setText(strXMin);

      QString strYMin, strYMax;
      strYMin = QString("%1").arg(m_fProjImgValueMin, 0, 'f', 1);
      strYMax = QString("%1").arg(m_fProjImgValueMax, 0, 'f', 1);

      ui.lineEditYMin->setText(strYMin);
      ui.lineEditYMax->setText(strYMax);
    }

    if (m_spProjImg3DFloat == nullptr) {
      return;
    }

    int horLen = m_spProjImg3DFloat->GetBufferedRegion().GetSize()[0];
    // int verLen = m_spProjImg3DFloat->GetBufferedRegion().GetSize()[1];

    // set edit maxium min
    QString strXMin = QString("%1").arg(horLen);
    ui.lineEditXMin->setText("0");
    ui.lineEditXMax->setText(strXMin);

    QString strYMin, strYMax;
    strYMin = QString("%1").arg(m_fProjImgValueMin, 0, 'f', 1);
    strYMax = QString("%1").arg(m_fProjImgValueMax, 0, 'f', 1);

    ui.lineEditYMin->setText(strYMin);
    ui.lineEditYMax->setText(strYMax);
  } else if (ui.radioButton_graph_recon->isChecked()) {
    if (m_spCrntReconImg == nullptr) {
      return;
    }

    int horLen = m_spCrntReconImg->GetBufferedRegion().GetSize()[0];
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

//
// void CbctRecon::SLT_GetProjectionProfile()
//{
//		m_dspYKReconImage->CreateImage(410,410,0);
//		ui.spinBoxReconImgSliceNo->setMinimum(0);
//		ui.spinBoxReconImgSliceNo->setMaximum(119);
//
//		m_iTmpIdx++;
//		ui.spinBoxReconImgSliceNo->setValue(m_iTmpIdx);
//
//
//		SLT_DrawReconImage();
//
//}
//
// void CbctRecon::SLT_GetReconImgProfile()
//{
//
//}

void CbctRecon::SLT_CopyTableToClipBoard() {
  qApp->clipboard()->clear();

  QStringList list;

  int rowCnt = m_pTableModel->rowCount();
  int columnCnt = m_pTableModel->columnCount();

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
      QStandardItem *item = m_pTableModel->item(j, i);
      list << item->text();
    }
    list << "\n";
  }

  qApp->clipboard()->setText(list.join("\t"));
}

void CbctRecon::SetProjDir(QString &strProjPath) {
  m_strPathGeomXML = "";
  m_strPathDirDefault = strProjPath;
  ui.lineEdit_HisDirPath->setText(strProjPath);

  UShortImageType::Pointer spNull;

  m_spCrntReconImg = spNull;     // fixed image // ID: RawCBCT
  m_spRawReconImg = spNull;      // just added --> when file is loaded
  m_spScatCorrReconImg = spNull; // just added --> after scatter correction

  FindAllRelevantPaths(strProjPath);
  init_DlgRegistration(m_strDCMUID);
}

void CbctRecon::SLT_SetHisDir() // Initialize all image buffer
{
  // Initializing..

  // Set folder --> then use RTK HIS Reader
  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.length() <= 1) {
    return;
  }

  SetProjDir(dirPath);

  m_vSelectedFileNames.clear();

  std::cout << "Push Load button to load projection images" << std::endl;

  // if (m_strPathGeomXML.length() > 1) //if geometry file was found
  //{
  //  SLT_LoadSelectedProjFiles();
  //}
  // else
  //{
  //  std::cout << "Geometry file is not ready. Find XML geometry file manually"
  //  << std::endl; ui.lineEdit_ElektaGeomPath->setText(QString(""));
  //}
}

QString getBowtiePath(QWidget *parent, const QDir &calDir) {
  return QFileDialog::getOpenFileName(
      parent, "Find air(+bowtie) filter image for subtraction",
      calDir.absolutePath(), "Projection (*.xim)", nullptr, nullptr);
}

std::tuple<bool, bool> CbctRecon::probeUser(const QString &guessDir) {

  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open CT DICOM Directory"), guessDir,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  bool dcm_success = false;
  if (!(dirPath.length() <= 1)) {
    Dcmtk_rt_study drs(dirPath.toLocal8Bit().constData());
    drs.load_directory(); // parse_directory();
    Plm_image plmImg;
    plmImg.set(drs.get_image());

    if (plmImg.have_image()) {
      // if (plmImg.load_native(dirPath.toLocal8Bit().constData())) {
      auto planCT_ss = drs.get_rtss(); // dies at end of scope...
      if (planCT_ss.get() != nullptr) {
        // ... so I copy to my own modern-C++ implementation
        m_structures->set_planCT_ss(planCT_ss.get());
      }

      ShortImageType::Pointer spShortImg = plmImg.itk_short();

      // Figure out whether this is NKI
      using ImageCalculatorFilterType =
          itk::MinimumMaximumImageCalculator<ShortImageType>;
      ImageCalculatorFilterType::Pointer imageCalculatorFilter =
          ImageCalculatorFilterType::New();
      // imageCalculatorFilter->SetImage(spShortImg);
      // imageCalculatorFilter->Compute();

      // double minVal0 = (double)(imageCalculatorFilter->GetMinimum());
      // double maxVal0 = (double)(imageCalculatorFilter->GetMaximum());

      // Thresholding
      using ThresholdImageFilterType =
          itk::ThresholdImageFilter<ShortImageType>;
      ThresholdImageFilterType::Pointer thresholdFilter =
          ThresholdImageFilterType::New();

      thresholdFilter->SetInput(spShortImg);
      thresholdFilter->ThresholdOutside(-1024, 3072); //--> 0 ~ 4095
      thresholdFilter->SetOutsideValue(-1024);
      thresholdFilter->Update();

      imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
      imageCalculatorFilter->Compute();

      auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
      auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

      std::cout << "Current Min and Max Values are	" << minVal << "	"
                << maxVal << std::endl;

      auto outputMinVal = static_cast<USHORT_PixelType>(minVal + 1024);
      auto outputMaxVal = static_cast<USHORT_PixelType>(maxVal + 1024);

      using RescaleFilterType =
          itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>;
      RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
      spRescaleFilter->SetInput(thresholdFilter->GetOutput());
      spRescaleFilter->SetOutputMinimum(outputMinVal);
      spRescaleFilter->SetOutputMaximum(outputMaxVal);
      spRescaleFilter->Update();

      // m_spRawReconImg = spRescaleFilter->GetOutput();
      m_spRefCTImg = spRescaleFilter->GetOutput();

      // UpdateReconImage(m_spRefCTImg, QString("DICOM reference image"));

      RegisterImgDuplication(REGISTER_REF_CT, REGISTER_MANUAL_RIGID);
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


void CbctRecon::SLT_LoadSelectedProjFiles() // main loading fuction for
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

  m_iImgCnt = 0; // should be reset
  m_iCntSelectedProj = 0;
  ReleaseMemory(); // only reset mem for indepent projection images

  std::string regexp;
  if (hisIsUsed) {
    regexp = ".[0-9a-fA-F].his";
  } else if (ximIsUsed) {
    regexp = "Proj_.*.xim";
  } else {
    regexp = "Proj_.*.hnd";
  }
  // char * regexp = ".*.his";
  std::vector<std::string> names;
  itk::RegularExpressionSeriesFileNames::Pointer regexpnames =
      itk::RegularExpressionSeriesFileNames::New();
  regexpnames->SetDirectory(dirPath.toLocal8Bit().constData());
  // regexpnames->SetNumericSort(false);
  regexpnames->SetNumericSort(true); // doesn't work with hexadecimal. and
                                     // [true/false] doesn't mean ascending or
                                     // descending
  regexpnames->SetRegularExpression(regexp);
  regexpnames->SetSubMatch(
      1); // SetSubMatch(0) led to sorting from last digit of the name
  // regexpnames->SetSubMatch(0); //SetSubMatch(0) led to sorting from last
  // digit of the name
  TRY_AND_EXIT_ON_ITK_EXCEPTION(names = regexpnames->GetFileNames());

  if (!IsFileNameOrderCorrect(names) && !ximIsUsed) {
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

  m_vExcludeProjIdx.clear();

  QString bowtiePath;

  if (!geomFileInfo.exists()) {
    std::cout << "Critical Error! geometry file is not existing. Please retry."
              << std::endl;
    return;
  }

  if (geomFileInfo.fileName() ==
      "_Frames.xml") // this is XVI XML. 
  {
    std::cout << "XVI Geometry File was found. This will be temporarily used:"
              << geomPath.toLocal8Bit().constData() << std::endl;
    LoadXVIGeometryFile(
        geomPath.toLocal8Bit().constData()); // will generate m_spFullGeometry
  } else if (geomFileInfo.fileName() ==
             "ProjectionInfo.xml") // this is OBI XML. 
  {
    std::cout
        << "Varian XML Geometry File was found. This will be temporarily used:"
        << geomPath.toLocal8Bit().constData() << std::endl;
    // ADDED BY AGRAVGAARD :::::::: Thief'd from RTK/applications/
    // rtkvarianobigemetry.cxx
    rtk::VarianObiGeometryReader::Pointer reader;
    reader = rtk::VarianObiGeometryReader::New();
    reader->SetXMLFileName(geomPath.toLocal8Bit().constData());
    reader->SetProjectionsFileNames(regexpnames->GetFileNames());
    TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->UpdateOutputData());
    // Write
    rtk::ThreeDCircularProjectionGeometryXMLFileWriter::Pointer xmlWriter =
        rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
    xmlWriter->SetFilename("RTKgeometry.xml");
    xmlWriter->SetObject(reader->GetGeometry());
    TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile());
    std::cout << "RTK standard Geometry XML File was created:"
              << "RTKgeometry.xml" << std::endl;
    LoadRTKGeometryFile("RTKgeometry.xml");
    // ::::::::::::::::::::::::::::LoadXMLGeometryFile(geomPath.toLocal8Bit().constData());
    // //will generate m_spFullGeometry
  } else if (geomFileInfo.fileName() ==
             "Scan.xml") // this is XIM XML. 
  {
    std::cout << "Varian Xim XML Geometry File was found. This will be "
                 "temporarily used:"
              << geomPath.toLocal8Bit().constData() << std::endl;
    // ADDED BY AGRAVGAARD :::::::: Thief'd from RTK/applications/
    // rtkvarianobigemetry.cxx
    rtk::VarianProBeamGeometryReader::Pointer reader;
    reader = rtk::VarianProBeamGeometryReader::New();
    reader->SetXMLFileName(geomPath.toLocal8Bit().constData());
    reader->SetProjectionsFileNames(regexpnames->GetFileNames());
    TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->UpdateOutputData());
    // Write
    rtk::ThreeDCircularProjectionGeometryXMLFileWriter::Pointer xmlWriter =
        rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
    xmlWriter->SetFilename("RTKgeometry.xml");
    xmlWriter->SetObject(reader->GetGeometry());
    TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile());
    std::cout << "RTK standard Geometry XML File was created:"
              << "RTKgeometry.xml" << std::endl;
    LoadRTKGeometryFile("RTKgeometry.xml");
    std::cout << "Done!";
    // ::::::::::::::::::::::::::::LoadXMLGeometryFile(geomPath.toLocal8Bit().constData());
    // //will generate m_spFullGeometry
  } else {
    std::cout << "RTK standard Geometry XML File was found:"
              << geomPath.toLocal8Bit().constData() << std::endl;
    LoadRTKGeometryFile(
        geomPath.toLocal8Bit().constData()); // will generate m_spFullGeometry
  }

  int iFullGeoDataSize = m_spFullGeometry->GetGantryAngles().size();
  if (iFullGeoDataSize < 1) {
    std::cout << "Not enough projection image (should be > 0)" << std::endl;
    return;
  }

  if (iFullGeoDataSize != fullCnt) {
    if (!ximIsUsed) {
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

  // 3) Seletively load projection file

  double meanGap = vectorMean(m_spFullGeometry->GetAngularGaps(
                       m_spFullGeometry->GetSourceAngles())) /
                   itk::Math::pi * 180.0;
  std::cout << "AngularGaps Mean (deg): " << meanGap << std::endl;
  double angleSum = vectorSum(m_spFullGeometry->GetAngularGaps(
                        m_spFullGeometry->GetSourceAngles())) /
                    itk::Math::pi * 180.0; // total rot. range from first angle
  std::cout << "AngularGaps Sum (deg): " << angleSum << std::endl;

  double gantryAngleInterval =
      ui.lineEdit_ManualProjAngleGap->text().toDouble();

  // if (ui.Radio_KeepOriginalAngles->isChecked())
  if (ui.Radio_ManualProjAngleGap->isChecked()) {
    // bManualGap = true;
    // std::cout << "Input angle gap in deg: " ;
    // cin >> gantryAngleInterval;

    if (gantryAngleInterval < meanGap) {
      std::cout << "Angle gap size is too small. Terminating the app"
                << std::endl;
      return;
      // bManualGap = false;
    }
  }

  ///////////////////////////////////Exclude outlier projection files

  std::vector<int> vExcludeIdx;
  if (!hisIsUsed) {
    vExcludeIdx.push_back(2);
    std::cout << "[Excluding-files-function] has exluded the second "
                 "projection, as it is often faulty in varian projections."
              << std::endl;
  } else {
    std::cout << "[Excluding-files-function] has been omitted. To reactivaite "
                 "it, please edit SLT_LoadSelectedProjFiles"
              << std::endl;
  }

  ///////////////////////////////////Exclude outlier projection files
  std::vector<int> vSelectedIdx;
  std::vector<int> vSelectedIdx_final;
  std::vector<int>::iterator itIdx;

  // double gantryAngleInterval =
  // ui.lineEdit_ManualProjAngleGap->text().toDouble();

  if (ui.Radio_ManualProjAngleGap->isChecked()) {
    // Select indices for recon
    // Generate norminal gantry values from the first angle
    double firstAngle = m_spFullGeometry->GetGantryAngles().at(0);
    double lastAngle =
        m_spFullGeometry->GetGantryAngles().at(iFullGeoDataSize - 1);

    std::vector<double> vNormAngles;

    int multiSize = round(angleSum / gantryAngleInterval) + 2;

    // CW only (179.xx -> 181.xx -> 359.xx --> 1.xx --> 179.xx), CCW should be
    // checked later
    for (int i = 0; i < multiSize; i++) {
      double curAngle = 0.0;

      if (m_bScanDirectionCW) {
        curAngle = firstAngle + i * gantryAngleInterval;
        if (curAngle >= 360.0) {
          curAngle = curAngle - 360.0;
        }
      } else {
        curAngle = firstAngle - i * gantryAngleInterval;
        if (curAngle < 0.0) {
          curAngle = curAngle + 360.0;
        }
      }
      // Don't add last gantry angle if their intervals are too small.

      // Last data will be added at the last part
      if (i > multiSize - 5) // last parts of the data
      {
        if (m_bScanDirectionCW) {
          if (curAngle <=
              lastAngle - gantryAngleInterval / 2.0) // from 5 latest indices,
          {
            vNormAngles.push_back(curAngle);
          }
        } else {
          if (curAngle >=
              lastAngle - gantryAngleInterval / 2.0) // from 5 latest indices,
          {
            vNormAngles.push_back(curAngle);
          }
        }
        // gantryAngleInterval/2.0 is given to remove "very near" value to the
        // last value
      } else {
        vNormAngles.push_back(curAngle);
      }
    }
    vNormAngles.push_back(lastAngle);

    for (double vNormAngle : vNormAngles) {
      std::cout << "Nominal proj. angle: ";
      std::cout << vNormAngle << std::endl;
    }

    // Collect appropriate indices
    GetSelectedIndices(m_spFullGeometry->GetGantryAngles(), vNormAngles,
                       vSelectedIdx, m_bScanDirectionCW, vExcludeIdx);

    for (itIdx = vSelectedIdx.begin(); itIdx != vSelectedIdx.end(); itIdx++) {
      std::cout << "Index: " << *itIdx << "     "
                << "GantryAngle: "
                << m_spFullGeometry->GetGantryAngles().at(*itIdx) << std::endl;
    }
  } else // not manual
  {
    for (int i = 0; i < iFullGeoDataSize; i++) {
      if (std::find(vExcludeIdx.begin(), vExcludeIdx.end(), i) ==
          vExcludeIdx.end()) { // if i is not included in vExcludeIdx
        vSelectedIdx.push_back(i);
      }
    }
  }
  // Another exlusion for kV off images

  vSelectedIdx_final.clear();

  // std::vector<int>::iterator itExclude;
  // for (itExclude = m_vExcludeProjIdx.begin(); itExclude !=
  // m_vExcludeProjIdx.end(); ++itExclude)
  //{
  //    int idx = (*itExclude);
  //    //if (std::find(m_vExcludeProjIdx.begin(), m_vExcludeProjIdx.end(),
  //    curIdx) == m_vExcludeProjIdx.end()) // if i is not included in
  //    vExcludeIdx
  //    //    vSelectedIdx_final.push_back(curIdx);

  //    std::cout << "Exclude " << idx << std::endl;
  //}

  std::vector<int>::iterator itFinal;
  int curIdx = 0;
  for (itFinal = vSelectedIdx.begin(); itFinal != vSelectedIdx.end();
       ++itFinal) {
    curIdx = (*itFinal);
    if (std::find(m_vExcludeProjIdx.begin(), m_vExcludeProjIdx.end(), curIdx) ==
        m_vExcludeProjIdx.end()) { // if i is not included in vExcludeIdx
      vSelectedIdx_final.push_back(curIdx);
    }
  }

  // Regenerate geometry object
  m_spCustomGeometry = GeometryType::New();

  // for (itIdx =vSelectedIdx.begin() ; itIdx != vSelectedIdx.end() ; itIdx++ )
  // #pragma omp parallel for private(itIdx) shared(m_spCustomGeometry)
  // schedule(static)
  for (itIdx = vSelectedIdx_final.begin(); itIdx != vSelectedIdx_final.end();
       itIdx++) {
    // 9 parameters are required
    const double curSID =
        m_spFullGeometry->GetSourceToIsocenterDistances().at(*itIdx);
    const double curSDD =
        m_spFullGeometry->GetSourceToDetectorDistances().at(*itIdx);
    double curGantryAngle = m_spFullGeometry->GetGantryAngles().at(*itIdx);
    const double kVAng =
        curGantryAngle * 180.0 *
        itk::Math::one_over_pi; // 360 / 2 = 180 radians to degrees
    double MVAng = kVAng - (hisIsUsed ? 0.0 : 90.0);
    if (MVAng < 0.0) {
      MVAng = MVAng + 360.0;
    }
    curGantryAngle = MVAng;

    const double curProjOffsetX =
        m_spFullGeometry->GetProjectionOffsetsX().at(*itIdx);
    const double curProjOffsetY =
        m_spFullGeometry->GetProjectionOffsetsY().at(*itIdx);

    const double curOutOfPlaneAngles =
        m_spFullGeometry->GetOutOfPlaneAngles().at(*itIdx);
    const double curInPlaneAngles =
        m_spFullGeometry->GetInPlaneAngles().at(*itIdx);

    const double curSrcOffsetX =
        m_spFullGeometry->GetSourceOffsetsX().at(*itIdx);
    const double curSrcOffsetY =
        m_spFullGeometry->GetSourceOffsetsY().at(*itIdx);

    m_spCustomGeometry->AddProjection(
        curSID, curSDD, curGantryAngle, curProjOffsetX,
        curProjOffsetY, // Flexmap, For Xim, these are 0
        curOutOfPlaneAngles,
        curInPlaneAngles,              // In elekta and varian, these are 0
        curSrcOffsetX, curSrcOffsetY); // In elekta and varian, these are 0
  }

  std::cout << "Total proj count: " << vSelectedIdx.size() << std::endl;
  std::cout << "Excluded proj count: " << m_vExcludeProjIdx.size() << std::endl;
  std::cout << "Final proj count: " << vSelectedIdx_final.size() << std::endl;

  // Regenerate fileNames and geometry object based on the selected indices.

  if (!m_vSelectedFileNames.empty()) {
    m_vSelectedFileNames.clear();
  }

  std::ofstream fout;
  fout.open("D:/DebugFileNames.txt");

  for (itIdx = vSelectedIdx_final.begin(); itIdx != vSelectedIdx_final.end();
       itIdx++) {
    std::string curStr = names.at(*itIdx);
    m_vSelectedFileNames.push_back(curStr);
    fout << curStr.c_str() << std::endl;
  }

  fout.close();

  ////YKTEMP
  ////std::cout << vSelectedIdx.size() << " is the number of selected index" <<
  /// std::endl; /std::cout << m_spCustomGeometry->GetGantryAngles().size() << "
  /// is the number of m_spCustomGeometry" << std::endl;
  // rtk::ThreeDCircularProjectionGeometryXMLFileWriter::Pointer xmlWriter =
  //    rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
  // xmlWriter->SetFilename("D:/FewProjGeom.xml");
  // xmlWriter->SetObject(m_spCustomGeometry);
  // TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile());
  ////Copy selected his files to a different folder

  // int countFiles = m_vSelectedFileNames.size();
  // for (int i = 0; i < countFiles; i++)
  //{
  //    QFileInfo fInfo(m_vSelectedFileNames.at(i).c_str());
  //    QString strDir = "D:/FewProjDir";
  //    QString strNewFilePath = strDir + "/" + fInfo.fileName();
  //    QFile::copy(fInfo.absoluteFilePath(), strNewFilePath);
  //}
  // std::cout << countFiles << " files were copied." << std::endl;
  ////YKTEMP

  // Reads the cone beam projections
  using ReaderType = rtk::ProjectionsReader<FloatImageType>;
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileNames(m_vSelectedFileNames);
  // TRY_AND_EXIT_ON_ITK_EXCEPTION(
  // std::thread calc_thread(read_projections, reader);
  std::thread calc_thread([&reader]() { reader->Update(); });
  // calc_thread.detach();

  std::cout << "Reader detached from main thread" << std::endl;
  // After reading the whole file,
  // HIS header should be saved

  using FilterReaderType =
      rtk::ProjectionsReader<FloatImage2DType>; // This one does a bit more than
                                                // we need, but mainly
                                                // statistical filters
  FilterReaderType::Pointer bowtiereader =
      FilterReaderType::New(); // we use is because we need the projections to
                               // be in the same unit (order of magnitude)

  QDir guessDir(geomFileInfo.absolutePath() + QString("/../"));
  //  Insta Recon, Dcm read
  std::tuple<bool, bool> answers;

  if (ximIsUsed) {
    QDir calDir(geomFileInfo.absolutePath() + QString("/Calibrations/"));
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
  } else {
    answers = probeUser(guessDir.absolutePath()); // ^^^
  }

  m_iCntSelectedProj = m_vSelectedFileNames.size();
  if (hisIsUsed) {
    std::cout << "Copying the HIS info to buffer." << std::endl;
    m_arrYKBufProj.resize(m_iCntSelectedProj);
    auto it_selected = m_vSelectedFileNames.begin();
    for (auto it = m_arrYKBufProj.begin();
         it != m_arrYKBufProj.end() &&
         it_selected != m_vSelectedFileNames.end();
         ++it, ++it_selected) {
      it->m_strFilePath = it_selected->c_str();
      it->CopyHisHeader(it_selected->c_str());
    }
  }

  calc_thread.join();

  std::cout << "Reader re-attached to main thread" << std::endl;

  if (bowtiePath.length() > 1 && ximIsUsed) {
    OpenCL_subtract3Dfrom2DbySlice_InPlace(
        static_cast<cl_float *>(reader->GetOutput()->GetBufferPointer()),
        static_cast<cl_float *>(bowtiereader->GetOutput()->GetBufferPointer()),
        reader->GetOutput()->GetLargestPossibleRegion().GetSize(),
        bowtiereader->GetOutput()->GetLargestPossibleRegion().GetSize());
  }

  std::cout << "ProjectionReader Get Spacing : "
            << reader->GetOutput()->GetSpacing() << std::endl;

  m_fProjSpacingX = reader->GetOutput()->GetSpacing()[0];
  m_fProjSpacingY = reader->GetOutput()->GetSpacing()[1];

  double originalMax = -1.0;
  double originalMin = -1.0;
  //                   µ/rho (cm^2/g)  * rho (g/cm^3) * path (cm)
  const double theoreticalMin =
      0.1541 * 1.225e-3 *
      sqrt(pow(m_spCustomGeometry->GetSourceToDetectorDistances()[0], 2) +
           pow(m_spCustomGeometry->GetSourceToDetectorDistances()[1],
               2) + // mm->cm
           pow(m_spCustomGeometry->GetSourceToDetectorDistances()[2], 2)) *
      0.1;

  const double correctionValue = GetMaxAndMinValueOfProjectionImage(
      originalMax, originalMin, reader->GetOutput()); // , theoreticalMin);
  std::cout << "Reader Max, Min=" << originalMax << "	" << originalMin
            << std::endl;

  if (correctionValue != 10.0 &&
      correctionValue != 20.0) { // 10 and 20 are error codes
    if ((originalMax - originalMin) > (log(65535.0f) - theoreticalMin)) {
      OpenCL_AddConst_MulConst_InPlace(
          static_cast<cl_float *>(reader->GetOutput()->GetBufferPointer()),
          reader->GetOutput()->GetLargestPossibleRegion().GetSize(),
          static_cast<cl_float>(-correctionValue), // 2048th lowest value
                                                   // (avoiding outliers simply)
          static_cast<cl_float>((log(65535.0f) - theoreticalMin) /
                                (originalMax - originalMin)));
    } else {
      OpenCL_AddConst_InPlace(
          static_cast<cl_float *>(reader->GetOutput()->GetBufferPointer()),
          reader->GetOutput()->GetLargestPossibleRegion().GetSize(),
          static_cast<cl_float>(-correctionValue)); // 2048th lowest value
                                                    // (avoiding outliers
                                                    // simply)
    }
    // Reset min max:
    originalMax = -1.0;
    originalMin = -1.0;
    if (GetMaxAndMinValueOfProjectionImage(originalMax, originalMin,
                                           reader->GetOutput()) != 20.0) {
      std::cout << "Reader Max, Min=" << originalMax << "	" << originalMin
                << std::endl;
    }
  }
  m_spProjImg3DFloat = reader->GetOutput(); // 1024 1024, line integ image

  m_fResampleF = ui.lineEdit_DownResolFactor->text().toDouble(); // 0.5

  if (m_fResampleF > 1 && m_fResampleF <= 0) {
    std::cout << "wrong resample factor. reset to 1.0" << std::endl;
    ui.lineEdit_DownResolFactor->setText("1.0");
    m_fResampleF = 1.0;
  }

  std::thread save_thread(saveImageAsMHA<FloatImageType>, m_spProjImg3DFloat);

  if (m_fResampleF != 1.0) {
    ResampleItkImage(m_spProjImg3DFloat, m_spProjImg3DFloat,
                     m_fResampleF); // was! BROKEN AF for .his where input size
                                    // != 1024 (tested with 1016) -> outputs
                                    // offset -inputoffset/refactor^2 and 4
                                    // pixels too few in x and y
  }

  if (!hisIsUsed && !ximIsUsed) { // -> hnd
    std::cout << "Fitted bowtie-filter correction ongoing..." << std::endl;
    SLT_DoBowtieCorrection();
  }

  ConvertLineInt2Intensity(m_spProjImg3DFloat, m_spProjImgRaw3D,
                           65535); // if X not 1024 == input size: out_offset =
                                   // in_offset + (1024*res_f -
                                   // X*res_f)*out_spacing     <- will still
                                   // break down at fw_projection

  FloatImageType::PointType originPt = m_spProjImg3DFloat->GetOrigin();
  FloatImageType::SizeType FloatImgSize =
      m_spProjImg3DFloat->GetBufferedRegion().GetSize();
  FloatImageType::SpacingType FloatImgSpacing =
      m_spProjImg3DFloat->GetSpacing();

  std::cout << "YKDEBUG: Origin" << originPt[0] << ", " << originPt[1] << ", "
            << originPt[2] << std::endl;
  std::cout << "YKDEBUG: Size" << FloatImgSize[0] << ", " << FloatImgSize[1]
            << ", " << FloatImgSize[2] << std::endl;
  std::cout << "YKDEBUG: Spacing" << FloatImgSpacing[0] << ", "
            << FloatImgSpacing[1] << ", " << FloatImgSpacing[2] << std::endl;

  std::cout << "Raw3DProj dimension "
            << m_spProjImgRaw3D->GetRequestedRegion().GetSize() << std::endl;

  std::cout << "Projection reading succeeded." << m_vSelectedFileNames.size()
            << " files were read" << std::endl;

  ui.pushButton_DoRecon->setEnabled(true);

  ui.spinBoxImgIdx->setMinimum(0);
  ui.spinBoxImgIdx->setMaximum(m_vSelectedFileNames.size() - 1);
  ui.spinBoxImgIdx->setValue(0); // it doesn't call Draw Event .. don't know
                                 // why.

  SetMaxAndMinValueOfProjectionImage(); // update min max projection image

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
  // Make sure the projections are saved before going out of scope.
  save_thread.join();
}

void CbctRecon::GetSelectedIndices(const std::vector<double> &vFullAngles,
                                   std::vector<double> &vNormAngles,
                                   std::vector<int> &vTargetIdx, bool bCW,
                                   std::vector<int> &vExcludingIdx) {
  // projection time. Begins with 179.xxx (CW)
  int latest_Idx = 0;

  double curVal = 0.0;
  double nextVal = 0.0;

  int sizeNom = vNormAngles.size();
  int sizeFull = vFullAngles.size();

  for (int i = 0; i < sizeNom; i++) {
    double tmpNominalValue = vNormAngles.at(i);

    for (int j = latest_Idx + 1; j < sizeFull - 1; j++) {
      int enExcludingMode = 0; // 0: safe,1: right is outlier, 2: left is
                               // outlier, 3: both are outlier

      // 1) Left point is outlier
      if (std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j) !=
              vExcludingIdx.end() &&
          std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j + 1) ==
              vExcludingIdx.end()) {
        enExcludingMode = 2;
      }
      // 2) Right point is outlier
      else if (std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j) ==
                   vExcludingIdx.end() &&
               std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j + 1) !=
                   vExcludingIdx.end()) {
        enExcludingMode = 1;
      }
      // 2) No outlier
      else if (std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j) ==
                   vExcludingIdx.end() &&
               std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j + 1) ==
                   vExcludingIdx.end()) {
        enExcludingMode = 0;
      }
      // 3) Both are outliers
      else if (std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j) !=
                   vExcludingIdx.end() &&
               std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j + 1) !=
                   vExcludingIdx.end()) {
        enExcludingMode = 3;
      }

      curVal = vFullAngles.at(j);
      nextVal = vFullAngles.at(j + 1);

      if (bCW) {
        // for full gantry angle value of 359.0 - 1.0 interface in CW
        if (curVal >
            nextVal + 0.2) // e.g.359)  - 0.5,, 0.2-->minimum angle diff.
        {
          if (tmpNominalValue < 100) {
            curVal = curVal - 360.0;
          } else if (tmpNominalValue > 260) {
            nextVal = nextVal + 360.0;
          }
        }
        if (tmpNominalValue >= curVal && tmpNominalValue <= nextVal) {

          // Add filtering
          // if j is among the excluding index list (e.g. outlier), just pass
          // it.

          double diffCur = fabs(tmpNominalValue - curVal);
          double diffNext = fabs(tmpNominalValue - nextVal);

          if (diffCur <= diffNext || enExcludingMode == 0 ||
              enExcludingMode == 1) {
            latest_Idx = j;
            vTargetIdx.push_back(latest_Idx);
          } else if (j != sizeFull - 2 &&
                     (enExcludingMode == 0 || enExcludingMode == 2)) {
            latest_Idx = j + 1;
            vTargetIdx.push_back(latest_Idx);
          } else {
            latest_Idx = j;
            // Skip to pushback
          }

          break;
        }
      } else // in CCW case
      {
        // for full gantry angle value of 1.0 - 359.0 interface in CCW
        if (curVal <
            nextVal + 0.01) // e.g.359)  - 0.5,, 0.2-->minimum angle diff.
        {
          if (tmpNominalValue < 100) { // for redundant check
            nextVal = nextVal - 360.0;
          } else if (tmpNominalValue > 260) { // for redundant check
            curVal = curVal + 360.0;
          }
        }

        // in CCW, next value should be smaller than curVal
        if (tmpNominalValue >= nextVal && tmpNominalValue <= curVal) {
          double diffCur = fabs(tmpNominalValue - curVal);
          double diffNext = fabs(tmpNominalValue - nextVal);

          if (diffCur <= diffNext || enExcludingMode == 0 ||
              enExcludingMode == 1) {
            latest_Idx = j;
            vTargetIdx.push_back(latest_Idx);
          } else if (j != sizeFull - 2 &&
                     (enExcludingMode == 0 || enExcludingMode == 2)) {
            latest_Idx = j + 1;
            vTargetIdx.push_back(latest_Idx);
          } else {
            latest_Idx = j;
            // Skip to pushback
          }
          break;
        }
      }
    }
  }
  // vTargetIdx.push_back(sizeFull-1); //omit the last image --> should be same
  // as first...
}

void CbctRecon::GetExcludeIndexByNames(
    const QString &outlierListPath, std::vector<std::string> &vProjFileFullPath,
    std::vector<int> &vExcludeIdx) {
  std::ifstream fin;
  fin.open(outlierListPath.toLocal8Bit().constData(), std::ios::in);
  if (static_cast<int>(fin.fail()) == 1) {
    return;
  }

  char str[MAX_LINE_LENGTH];

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    if (strlen(&str[0]) < 1) {
      continue;
    }

    std::cout << "Outlier file: " << &str[0] << std::endl;

    int curIdx = 0;
    for (auto &it : vProjFileFullPath) {
      if (strstr(it.c_str(), static_cast<const char *>(str)) != nullptr) {
        std::cout << "Detected in the list. Index = " << curIdx << std::endl;
        vExcludeIdx.push_back(curIdx);
        break;
      }
      curIdx++;
    }
  }
  fin.close();
}

void CbctRecon::SetMaxAndMinValueOfProjectionImage() // should be called
                                                     // whenever proj image is
                                                     // changed
{
  if (m_iImgCnt > 0) {
    m_fProjImgValueMax = 65535;
    m_fProjImgValueMin = 0;
    return;
  }

  if (m_spProjImg3DFloat == nullptr) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<FloatImageType> it(
      m_spProjImg3DFloat, m_spProjImg3DFloat->GetRequestedRegion());

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  m_fProjImgValueMin = 65535.0;
  m_fProjImgValueMax = -9999.0;

  while (!it.IsAtEnd()) {
    while (!it.IsAtEndOfSlice()) {
      while (!it.IsAtEndOfLine()) {
        double tmpVal = it.Get();

        if (m_fProjImgValueMax < tmpVal) {
          m_fProjImgValueMax = tmpVal;
        }

        if (m_fProjImgValueMin > tmpVal) {
          m_fProjImgValueMin = tmpVal;
        }

        ++it;
      }
      it.NextLine();
    }
    it.NextSlice();
  }
}

double CbctRecon::GetMaxAndMinValueOfProjectionImage(
    double &fProjImgValueMax, double &fProjImgValueMin,
    const FloatImageType::Pointer &projImage) //, const double theoreticalMin)
{
  if (projImage == nullptr) {
    fProjImgValueMax = -1.0;
    fProjImgValueMin = -1.0;
    return 20.0;
  }

  using MinMaxCalcType = itk::MinimumMaximumImageCalculator<FloatImageType>;
  MinMaxCalcType::Pointer MinMaxFilter = MinMaxCalcType::New();
  MinMaxFilter->SetImage(projImage);
  clock_t begin = std::clock();
  MinMaxFilter->Compute();
  fProjImgValueMin = MinMaxFilter->GetMinimum();
  fProjImgValueMax = MinMaxFilter->GetMaximum();
  clock_t end_time = std::clock();

  std::cout << "Min: " << fProjImgValueMin
            << " max: " << fProjImgValueMax << std::endl;

  return fProjImgValueMin;
}

void CbctRecon::SLT_DataProbeProj() {
  double dspWidth = ui.labelImageRaw->width();
  double dspHeight = ui.labelImageRaw->height();
  int dataWidth = 0;
  int dataHeight = 0;
  int dataX = 0;
  int dataY = 0;
  int dataZ = 0;
  double fProbeValue = 0.0;

  if (m_iImgCnt > 0) // there is indep loaded projection files
  {
    dataWidth = m_dspYKImgProj->m_iWidth;
    dataHeight = m_dspYKImgProj->m_iHeight;

    dataX = qRound(ui.labelImageRaw->x / dspWidth * dataWidth);
    dataY = qRound(ui.labelImageRaw->y / dspHeight * dataHeight);
    dataZ = ui.spinBoxImgIdx->value();
    fProbeValue = m_dspYKImgProj->m_pData[dataWidth * dataY + dataX];
  } else {
    if (m_spProjImg3DFloat == nullptr) {
      return;
    }

    dataWidth =
        static_cast<int>(m_spProjImg3DFloat->GetBufferedRegion().GetSize()[0]);
    dataHeight =
        static_cast<int>(m_spProjImg3DFloat->GetBufferedRegion().GetSize()[1]);

    // int crntIdx = ui.spinBoxImgIdx->value();
    // These are displayed data (just index data)
    dataX = qRound(ui.labelImageRaw->x / dspWidth * dataWidth);
    dataY = qRound(ui.labelImageRaw->y / dspHeight * dataHeight);
    dataZ = ui.spinBoxImgIdx->value();

    // fProbeValue = m_dspYKImgProj->m_pData[dataWidth*dataY +
    // dataX]/m_multiplyFactor;
    fProbeValue =
        m_dspYKImgProj->m_pData[dataWidth * dataY + dataX] / m_multiplyFactor +
        m_fProjImgValueMin;
  }

  QString dspText;
  dspText = QString("(%1, %2, %3): %4")
                .arg(dataX)
                .arg(dataY)
                .arg(dataZ)
                .arg(fProbeValue, 0, 'f', 2);
  ui.lineEdit_DataProbe_Proj->setText(dspText);
}

void CbctRecon::SLT_DataProbeRecon() {
  if (m_spCrntReconImg == nullptr) {
    return;
  }

  double dspWidth = ui.labelReconImage->width();
  double dspHeight = ui.labelReconImage->height();

  auto dataWidth =
      static_cast<int>(m_spCrntReconImg->GetBufferedRegion().GetSize()[0]);
  auto dataHeight =
      static_cast<int>(m_spCrntReconImg->GetBufferedRegion().GetSize()[1]);

  // int crntIdx = ui.spinBoxImgIdx->value();
  // These are displayed data (just index data)

  int dataX = qRound(ui.labelReconImage->x / dspWidth * dataWidth);
  int dataY = qRound(ui.labelReconImage->y / dspHeight * dataHeight);
  int dataZ = ui.spinBoxReconImgSliceNo->value();

  unsigned short iProbeValue =
      m_dspYKReconImage->m_pData[dataWidth * dataY + dataX];
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


void CbctRecon::SLT_DrawGraph() // based on profile
{
  if (m_pTableModel == nullptr) {
    return;
  }

  // Draw only horizontal, center

  QVector<double> vAxisX; // can be rows or columns
  QVector<double> vAxisY;

  // QStandardItemModel 	m_pTableModel.item()
  int dataLen = m_pTableModel->rowCount();

  if (dataLen < 1) {
    return;
  }

  // std::cout << "check graph 1" << std::endl;
  ui.customPlot->clearGraphs();

  double minX = 9999.0;
  double maxX = -1.0;

  for (int i = 0; i < dataLen; i++) {
    QStandardItem *tableItem1 = m_pTableModel->item(i, 0);
    QStandardItem *tableItem2 = m_pTableModel->item(i, 1);
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

void CbctRecon::SLT_UpdateTable() {
  if (m_spCrntReconImg == nullptr) {
    ui.radioButton_graph_proj->setChecked(true);
  }

  // std::cout << "check 1" << std::endl;
  YK16GrayImage *pYKImg = nullptr;
  double fMultiPlyFactor = 1.0;
  double fMinValue = 0.0;

  if (ui.radioButton_graph_proj->isChecked()) {
    pYKImg = m_dspYKImgProj.get(); // you may look, but no touching!

    if (m_iImgCnt > 0) { // if indep image
      fMultiPlyFactor = 1.0;
    } else {
      fMultiPlyFactor = m_multiplyFactor;
      fMinValue = m_fProjImgValueMin;
    }
  } else {
    pYKImg = m_dspYKReconImage.get();
    fMultiPlyFactor = 1.0;
    fMinValue = 0.0;
  }
  if (pYKImg == nullptr) {
    return;
  }

  // std::cout << "check 2" << std::endl;

  if (m_pTableModel != nullptr) {
    delete m_pTableModel;
    m_pTableModel = nullptr;
  }
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
  m_pTableModel =
      new QStandardItemModel(rowSize, columnSize, this); // 2 Rows and 3 Columns

  // for (int i = 0 ; i<columnSize ; i++)
  //{
  // QFileInfo tmpInfo = QFileInfo(m_arrYKImage[i].m_strFilePath);
  // m_pTableModel->setHorizontalHeaderItem(0, new
  // QStandardItem(QString("Index"))); m_pTableModel->setHorizontalHeaderItem(0,
  // new QStandardItem(QString("Profile")));
  m_pTableModel->setHorizontalHeaderItem(
      0, new QStandardItem(QString("Position(mm)")));
  m_pTableModel->setHorizontalHeaderItem(1,
                                         new QStandardItem(QString("Value")));
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
    if (m_spCrntReconImg != nullptr) {
      UShortImageType::PointType tmpOrigin = m_spCrntReconImg->GetOrigin();
      UShortImageType::SpacingType tmpSpacing = m_spCrntReconImg->GetSpacing();
      originX = tmpOrigin[0];
      originY = tmpOrigin[1];
      spacingX = tmpSpacing[0];
      spacingY = tmpSpacing[1];
    }
  }

  // std::cout << "check 6" << std::endl;

  QVector<double> vPos;
  if (ui.radioButton_Profile_Hor->isChecked()) {
    for (int i = 0; i < rowSize; i++) {
      vPos.push_back(originX + i * spacingX);
    }
  } else {
    for (int i = 0; i < rowSize; i++) {
      vPos.push_back(originY + i * spacingY);
    }
  }

  QVector<double> vProfile;
  if (ui.radioButton_Profile_Hor->isChecked()) {
    pYKImg->GetProfileData(vProfile, DIRECTION_HOR);
  } else {
    pYKImg->GetProfileData(vProfile, DIRECTION_VER);
  }

  // int i = fixedY;
  for (int i = 0; i < rowSize; i++) {
    auto tmpVal1 = static_cast<qreal>(vPos[i]);
    m_pTableModel->setItem(i, 0, new QStandardItem(QString("%1").arg(tmpVal1)));

    qreal tmpVal2 =
        static_cast<qreal>(vProfile[i]) / fMultiPlyFactor + fMinValue;
    m_pTableModel->setItem(i, 1, new QStandardItem(QString("%1").arg(tmpVal2)));
  }

  ui.tableViewReconImgProfile->setModel(m_pTableModel); // also for proj

  // std::cout << "check 7" << std::endl;
  SLT_DrawGraph();
}

bool CbctRecon::IsFileNameOrderCorrect(std::vector<std::string> &vFileNames) {
  // regardless of whether number or hexa codes,
  // we can convert it from number to hexa number and compare the order

  int size = vFileNames.size();

  if (size < 2) {
    return false;
  }

  std::vector<int> arrNum(size);

  QString crntFilePath;

  int index = 0;
  while (index < size) {
    crntFilePath = vFileNames.at(index++).c_str();
    QFileInfo fileInfo = QFileInfo(crntFilePath);
    QDir dir = fileInfo.absoluteDir();
    QString file_basename = fileInfo.baseName();
    QString newBaseName = HexStr2IntStr(file_basename);
    arrNum.push_back(newBaseName.toInt());
  }

  /*bool bOrderOK = true;
  for (int i = 0; i < size - 1; i++)
  {
      if (arrNum[i] >= arrNum[i + 1]) {
          bOrderOK = false;
              }
  }*/

  auto index_of_nonascending =
      std::adjacent_find(arrNum.begin(), arrNum.end(), std::greater<>());

  return (index_of_nonascending == arrNum.end());
}

// Mouse Left Click
void CbctRecon::SLT_CalculateROI_Recon() {
  if (m_dspYKReconImage == nullptr) {
    return;
  }

  if (m_spCrntReconImg == nullptr) {
    return;
  }

  double dspWidth = ui.labelReconImage->width();
  double dspHeight = ui.labelReconImage->height();

  int dataWidth = m_dspYKReconImage->m_iWidth;
  int dataHeight = m_dspYKReconImage->m_iHeight;
  if (dataWidth * dataHeight == 0) {
    return;
  }

  // int crntIdx = ui.spinBoxImgIdx->value();
  // These are displayed data (just index data)

  int dataX = qRound(ui.labelReconImage->x / dspWidth * dataWidth);
  int dataY = qRound(ui.labelReconImage->y / dspHeight * dataHeight);
  int dataZ = ui.spinBoxReconImgSliceNo->value();

  auto originX = static_cast<double>(m_spCrntReconImg->GetOrigin()[0]);
  auto originY = static_cast<double>(m_spCrntReconImg->GetOrigin()[1]);
  auto originZ = static_cast<double>(m_spCrntReconImg->GetOrigin()[2]);

  auto spacingX = static_cast<double>(m_spCrntReconImg->GetSpacing()[0]);
  auto spacingY = static_cast<double>(m_spCrntReconImg->GetSpacing()[1]);
  auto spacingZ = static_cast<double>(m_spCrntReconImg->GetSpacing()[2]);

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

  m_dspYKReconImage->SetProfileProbePos(dataX, dataY);
  if (ui.radioButton_Profile_Hor->isChecked()) {
    m_dspYKReconImage->m_bDrawProfileX = true;
    m_dspYKReconImage->m_bDrawProfileY = false;
  } else {
    m_dspYKReconImage->m_bDrawProfileX = false;
    m_dspYKReconImage->m_bDrawProfileY = true;
  }

  // m_dspYKReconImage value itself
  int ROI_size = ui.lineEdit_ROI_size->text().toInt();
  if (ROI_size < 0) {
    return;
  }

  if (ROI_size > 0) {
    m_dspYKReconImage->setROI(
        qRound(dataX - ROI_size / 2.0), qRound(dataY - ROI_size / 2.0),
        qRound(dataX + ROI_size / 2.0), qRound(dataY + ROI_size / 2.0));
    m_dspYKReconImage->CalcImageInfo_ROI();
    m_dspYKReconImage->DrawROIOn(true);

    // m_dspYKImgProj->m_pData[iNumWidth + width*iNumHeight] = (unsigned
    // short)((tmpVal- m_fProjImgValueMin)*m_multiplyFactor);

    QString strMean;
    strMean = QString("%1").arg(m_dspYKReconImage->m_fPixelMean_ROI, 0, 'f', 2);
    QString strSD;
    strSD = QString("%1").arg(m_dspYKReconImage->m_fPixelSD_ROI, 0, 'f', 2);
    ui.lineEdit_ROI_mean->setText(strMean);
    ui.lineEdit_ROI_SD->setText(strSD);
  } else {
    m_dspYKReconImage->DrawROIOn(false);
  }

  SLT_DrawReconImage();
}

void CbctRecon::SLT_CalculateROI_Proj() {
  if (m_dspYKImgProj == nullptr) {
    return;
  }

  double dspWidth = ui.labelImageRaw->width();
  double dspHeight = ui.labelImageRaw->height();

  int dataWidth = m_dspYKImgProj->m_iWidth;
  int dataHeight = m_dspYKImgProj->m_iHeight;
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

  m_dspYKImgProj->SetProfileProbePos(dataX, dataY);

  if (ui.radioButton_Profile_Hor->isChecked()) {
    m_dspYKImgProj->m_bDrawProfileX = true;
    m_dspYKImgProj->m_bDrawProfileY = false;
  } else {
    m_dspYKImgProj->m_bDrawProfileX = false;
    m_dspYKImgProj->m_bDrawProfileY = true;
  }

  // m_dspYKReconImage value itself
  int ROI_size = ui.lineEdit_ROI_size->text().toInt();
  if (ROI_size < 0) {
    return;
  }

  if (ROI_size > 0) {
    m_dspYKImgProj->setROI(
        qRound(dataX - ROI_size / 2.0), qRound(dataY - ROI_size / 2.0),
        qRound(dataX + ROI_size / 2.0), qRound(dataY + ROI_size / 2.0));
    m_dspYKImgProj->CalcImageInfo_ROI();
    m_dspYKImgProj->DrawROIOn(true);
    QString strMean;
    // strMean.sprintf("%5.1f", m_dspYKImgProj->m_fPixelMean_ROI);
    strMean = QString("%1").arg(
        (m_dspYKImgProj->m_fPixelMean_ROI / m_multiplyFactor) +
            m_fProjImgValueMin,
        0, 'f', 2);

    QString strSD;
    strSD = QString("%1").arg(m_dspYKImgProj->m_fPixelSD_ROI / m_multiplyFactor,
                              0, 'f', 2);
    ui.lineEdit_ROI_mean->setText(strMean);
    ui.lineEdit_ROI_SD->setText(strSD);
  } else {
    m_dspYKImgProj->DrawROIOn(false);
  }

  SLT_DrawProjImages();
}



void CbctRecon::SLT_GoForcedProbePos() // when forced probe button was clicked
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
    if (m_spProjImg3DFloat == nullptr) {
      return;
    }

    originX = m_spProjImg3DFloat->GetOrigin()[0];
    originY = m_spProjImg3DFloat->GetOrigin()[1];
    originZ = m_spProjImg3DFloat->GetOrigin()[2];

    spacingX = m_spProjImg3DFloat->GetSpacing()[0];
    spacingY = m_spProjImg3DFloat->GetSpacing()[1];
    spacingZ = m_spProjImg3DFloat->GetSpacing()[2];

    sliceIdx = qRound((fForcedProbePosZ - originZ) / spacingZ);

    if (sliceIdx < 0 || sliceIdx >= m_iImgCnt) {
      return;
    }

    ui.spinBoxImgIdx->setValue(sliceIdx); // Draw function is called

    dspWidth = ui.labelImageRaw->width();
    dspHeight = ui.labelImageRaw->height();

    dataWidth = m_dspYKImgProj->m_iWidth;
    dataHeight = m_dspYKImgProj->m_iHeight;

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
    if (m_spCrntReconImg == nullptr) {
      return;
    }

    originX = m_spCrntReconImg->GetOrigin()[0];
    originY = m_spCrntReconImg->GetOrigin()[1];
    originZ = m_spCrntReconImg->GetOrigin()[2];

    spacingX = m_spCrntReconImg->GetSpacing()[0];
    spacingY = m_spCrntReconImg->GetSpacing()[1];
    spacingZ = m_spCrntReconImg->GetSpacing()[2];

    sliceIdx = qRound((fForcedProbePosZ - originZ) / spacingZ);

    if (sliceIdx < 0 ||
        sliceIdx >= static_cast<int>(
                        m_spCrntReconImg->GetBufferedRegion().GetSize()[2])) {
      return;
    }

    ui.spinBoxReconImgSliceNo->setValue(sliceIdx); // Draw function is called

    dspWidth = ui.labelReconImage->width();
    dspHeight = ui.labelReconImage->height();

    dataWidth = m_dspYKReconImage->m_iWidth;
    dataHeight = m_dspYKReconImage->m_iHeight;

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

void CbctRecon::PostApplyFOVDispParam() {
  if (m_dspYKReconImage == nullptr) {
    return;
  }

  float physPosX = ui.lineEdit_PostFOV_X->text().toFloat();
  float physPosY = ui.lineEdit_PostFOV_Y->text().toFloat();

  float physRadius = ui.lineEdit_PostFOV_R->text().toFloat();
  float physTablePosY = ui.lineEdit_PostTablePosY->text().toFloat();

  UShortImageType::PointType origin = m_spCrntReconImg->GetOrigin();
  UShortImageType::SpacingType spacing = m_spCrntReconImg->GetSpacing();
  // UShortImageType::SizeType size =
  // m_spCrntReconImg->GetBufferedRegion().GetSize();

  int pixPosX = qRound((physPosX - static_cast<double>(origin[0])) /
                       static_cast<double>(spacing[0]));
  int pixPosY = qRound((physPosY - static_cast<double>(origin[1])) /
                       static_cast<double>(spacing[1]));

  int pixRadius = qRound(physRadius / static_cast<double>(spacing[0]));

  // int pixWidth = qRound((qreal) size[0]);
  // int pixHeight = qRound((qreal) size[1]);

  int pixTableY = qRound((physTablePosY - static_cast<double>(origin[1])) /
                         static_cast<double>(spacing[1]));

  if (pixPosX >= 0 && pixPosY < m_dspYKReconImage->m_iWidth && pixPosY >= 0 &&
      pixPosY < m_dspYKReconImage->m_iHeight && pixRadius > 0 &&
      pixRadius < m_dspYKReconImage->m_iWidth && pixTableY >= 0 &&
      pixTableY < m_dspYKReconImage->m_iHeight) {
    m_dspYKReconImage->m_ptFOVCenter.setX(pixPosX); // data pos
    m_dspYKReconImage->m_ptFOVCenter.setY(pixPosY);

    m_dspYKReconImage->m_iFOVRadius = pixRadius;
    m_dspYKReconImage->m_iTableTopPos = pixTableY;
  }
}

void CbctRecon::SLT_PostApplyFOVDispParam() {
  // PostApplyFOVDispParam();
  SLT_DrawReconImage();
}

void CbctRecon::CropSupInf(UShortImageType::Pointer &sp_Img,
                           float physPosInfCut, float physPosSupCut) {
  if (sp_Img == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position
  UShortImageType::PointType origin = sp_Img->GetOrigin();
  UShortImageType::SpacingType spacing = sp_Img->GetSpacing();
  UShortImageType::SizeType size = sp_Img->GetBufferedRegion().GetSize();

  std::cout << "Old Origin" << origin << std::endl;
  std::cout << "Old spacing" << spacing << std::endl;
  std::cout << "Old size" << size << std::endl;

  UShortImageType::SizeType sizeLower{}, sizeUpper{}; // not index. this is
                                                      // width of pixels that
                                                      // will be taken away.
  sizeLower[0] = 0;
  sizeLower[1] = 0;
  sizeLower[2] = 0;
  /*indexUpper[0] = size[0] - 1;
  indexUpper[1] = size[1] - 1;
  indexUpper[2] = size[2] - 1;*/
  sizeUpper[0] = 0;
  sizeUpper[1] = 0;
  sizeUpper[2] = 0;

  double minPosSI = origin[2];
  double maxPosSI = origin[2] + (size[2] - 1) * spacing[2];

  if (minPosSI >= physPosInfCut) {
    physPosInfCut = minPosSI;
  }
  if (maxPosSI <= physPosSupCut) {
    physPosSupCut = maxPosSI;
  }

  if (physPosSupCut <= physPosInfCut) {
    return;
  }

  ////calc index
  sizeLower[2] = qRound((physPosInfCut - minPosSI) / spacing[2]);
  sizeUpper[2] = qRound((maxPosSI - physPosSupCut) / spacing[2]);
  //
  using CropImageFilterType =
      itk::CropImageFilter<UShortImageType, UShortImageType>;
  CropImageFilterType::Pointer CropFilter = CropImageFilterType::New();

  CropFilter->SetInput(sp_Img);
  CropFilter->SetLowerBoundaryCropSize(sizeLower);
  CropFilter->SetUpperBoundaryCropSize(sizeUpper);

  CropFilter->Update();

  if (sp_Img == m_spRawReconImg) {
    sp_Img = CropFilter->GetOutput();
    m_spRawReconImg = sp_Img;
  }

  if (sp_Img == m_spRefCTImg) {
    sp_Img = CropFilter->GetOutput();
    m_spRefCTImg = sp_Img;
    m_spManualRigidCT = sp_Img;
  }

  UShortImageType::PointType origin_new = sp_Img->GetOrigin();
  UShortImageType::SpacingType spacing_new = sp_Img->GetSpacing();
  UShortImageType::SizeType size_new = sp_Img->GetBufferedRegion().GetSize();

  // origin_new[2] = physPosInfCut;
  // sp_Img->SetOrigin(origin_new);

  std::cout << "New Origin" << origin_new << std::endl;
  std::cout << "New spacing" << spacing_new << std::endl;
  std::cout << "New size" << size_new << std::endl;

  std::cout << "LowPos[mm, index] = " << physPosInfCut << ", " << sizeLower[2]
            << std::endl;
  std::cout << "UpperPos[mm, index] = " << physPosSupCut << ", " << sizeUpper[2]
            << std::endl;
  std::cout << "Cropping SI has been successfully done." << std::endl;

  // Result: same image after cropping
  /*
      sizeDiff[0] = CropFilter->GetOutput()->GetBufferedRegion().GetSize()[0] -
     pYKImageROI->m_iWidth; sizeDiff[1] =
     CropFilter->GetOutput()->GetBufferedRegion().GetSize()[1] -
     pYKImageROI->m_iHeight;

      if (sizeDiff[0] != 0 || sizeDiff[1] != 0)
      {
      std::cout << "Cross-correlation error! template size is not matching ROI
     image even after cropping" << std::endl; return;
      }

      rescaleFilter->SetInput(CropFilter->GetOutput());*/
}

// mm
void CbctRecon::CropFOV3D(UShortImageType::Pointer &sp_Img, float physPosX,
                          float physPosY, float physRadius,
                          float physTablePosY) {
  if (sp_Img == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position
  UShortImageType::PointType origin = sp_Img->GetOrigin();
  UShortImageType::SpacingType spacing = sp_Img->GetSpacing();
  // UShortImageType::SizeType size = sp_Img->GetBufferedRegion().GetSize();

  // itk::ImageSliceConstIteratorWithIndex<FloatImageType> it (m_spReconImg,
  // m_spReconImg->GetRequestedRegion());
  itk::ImageSliceIteratorWithIndex<UShortImageType> it(
      sp_Img, sp_Img->GetBufferedRegion());

  // ImageSliceConstIteratorWithIndex<ImageType> it( image,
  // image->GetRequestedRegion() );
  // UShortImageType::SizeType imgSize = sp_Img->GetBufferedRegion().GetSize();
  // //1016x1016 x z

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

        if (pow(crntPhysX - physPosX, 2.0) + pow(crntPhysY - physPosY, 2.0) >=
            pow(physRadius, 2.0)) {
          //(*it) = (unsigned short)0; //air value
          it.Set(0);
        }

        if (crntPhysY >= physTablePosY) {
          it.Set(0);
        }
        ++it;
        iPosX++;
      }
      it.NextLine();
      iPosY++;
    }
    it.NextSlice();
    iNumSlice++;
  }
}

void CbctRecon::SLT_DoPostProcessing() {
  if (m_spCrntReconImg == nullptr) {
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

  CropFOV3D(m_spCrntReconImg, physPosX, physPosY, physRadius, physTablePosY);

  SLT_DrawReconImage();
}

void CbctRecon::SLT_PostProcCropInv() {
  if (m_spCrntReconImg == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position

  double physPosX = ui.lineEdit_PostFOV_X->text().toDouble();
  double physPosY = ui.lineEdit_PostFOV_Y->text().toDouble();

  double physRadius = ui.lineEdit_PostFOV_R->text().toDouble();
  // double physTablePosY = ui.lineEdit_PostTablePosY->text().toDouble();

  UShortImageType::PointType origin = m_spCrntReconImg->GetOrigin();
  UShortImageType::SpacingType spacing = m_spCrntReconImg->GetSpacing();
  // UShortImageType::SizeType size =
  // m_spCrntReconImg->GetBufferedRegion().GetSize();

  // itk::ImageSliceConstIteratorWithIndex<FloatImageType> it (m_spReconImg,
  // m_spReconImg->GetRequestedRegion());
  itk::ImageSliceIteratorWithIndex<UShortImageType> it(
      m_spCrntReconImg, m_spCrntReconImg->GetRequestedRegion());

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


void CbctRecon::SLT_ExportALL_DCM_and_SHORT_HU_and_calc_WEPL() {

  if (m_spCrntReconImg == nullptr) {
    return;
  }

  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), m_strPathDirDefault,
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
      strPatientID = m_strDCMUID;
    }
    // strPatientID = m_strDCMUID + "_" + strEndFix;
  } else {
    strPatientID = m_strDCMUID;
  }

  if (strPatientID.isEmpty()) {
    return;
  }

  for (int i = 0; i < m_pDlgRegistration->ui.comboBoxImgFixed->count(); i++) {
    m_pDlgRegistration->ui.comboBoxImgFixed->setCurrentIndex(i);
    QString strDirName = m_pDlgRegistration->ui.comboBoxImgFixed->currentText();
    bool tmpResult = crntDir.mkdir(strDirName); // what if the directory exists?
    if (!tmpResult) {
      std::cout
          << "DICOM dir seems to exist already. Files will be overwritten."
          << std::endl;
    }

    QString strSavingFolder = dirPath + "/" + strDirName;
    QString strFullName = strLastName + ", " + strFirstName;
    m_pDlgRegistration->LoadImgFromComboBox(0, strDirName);

    SaveUSHORTAsSHORT_DICOM(m_pDlgRegistration->m_spFixed, strPatientID,
                            strFullName, strSavingFolder);
    QString mhaFileName = strSavingFolder + "/" + strDirName + ".mha";
    ExportReconSHORT_HU(m_pDlgRegistration->m_spFixed, mhaFileName);
  }
  SLT_GeneratePOIData();
  QString angle_end_one("1");
  ui.lineEdit_AngEnd->setText(angle_end_one);
  SLT_ExportAngularWEPL_byFile();
}

void CbctRecon::SLT_ExportReconSHORT_HU() {
  QString strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "signed short meta image (*.mha)", nullptr,
      nullptr);
  if (strPath.length() <= 1) {
    return;
  }
  ExportReconSHORT_HU(m_spCrntReconImg, strPath);

}

void CbctRecon::CopyDictionary(itk::MetaDataDictionary &fromDict,
                               itk::MetaDataDictionary &toDict) {
  using DictionaryType = itk::MetaDataDictionary;

  DictionaryType::ConstIterator itr = fromDict.Begin();
  DictionaryType::ConstIterator end = fromDict.End();
  using MetaDataStringType = itk::MetaDataObject<std::string>;

  while (itr != end) {
    itk::MetaDataObjectBase::Pointer entry = itr->second;

    MetaDataStringType::Pointer entryvalue =
        dynamic_cast<MetaDataStringType *>(entry.GetPointer());
    if (entryvalue != nullptr) {
      std::string tagkey = itr->first;
      std::string tagvalue = entryvalue->GetMetaDataObjectValue();
      itk::EncapsulateMetaData<std::string>(toDict, tagkey, tagvalue);
    }
    ++itr;
  }
}

inline float BeamHardModel(const double val, const double a, const double b,
                           const double c, const double d) {
  return a * pow(val, 3.0) + b * pow(val, 2.0) + c * val + d;
}

inline double HndBeamHardModel(const double val) {
  // a * x^3 + b * x^2 + c * x + d
  return 6.0e-08 * pow(val, 3.0) - 1.0e-08 * pow(val, 2.0) - 5.0e-07 * val +
         8.0e-01;
}

inline double HisBeamHardModel(const double val) {
  // a * x^3 + b * x^2 + c * x + d
  return 9.321e-05 * pow(val, 3.0) - 2.609e-03 * pow(val, 2.0) +
         3.374e-02 * val + 9.691e-01;
}

inline double
XimBeamHardModel(const double val) { // a * x^3 + b * x^2 + c * x + d
  return 6.0e-08 * pow(val, 3.0) + 9.0e-5 * pow(val, 2.0) + 1.0e-2 * val + 0.8;
}

void hndBeamHardening(float *pBuffer, const int nPix) {
#pragma omp parallel for
  for (int i = 0; i < nPix; i++) {
    pBuffer[i] = static_cast<float>(
        ((pBuffer[i] < 1.189) ? pBuffer[i]
                              : pBuffer[i] * HndBeamHardModel(pBuffer[i])) +
        1.47);
  }
}
void hisBeamHardening(float *pBuffer, const int nPix) {
#pragma omp parallel for
  for (int i = 0; i < nPix; i++) {
    pBuffer[i] = static_cast<float>(
        ((pBuffer[i] < 1.189) ? pBuffer[i]
                              : pBuffer[i] * HisBeamHardModel(pBuffer[i])));
  }
}
void ximBeamHardening(float *pBuffer, const int nPix) {
#pragma omp parallel for
  for (int i = 0; i < nPix; i++) {
    pBuffer[i] = static_cast<float>(
        ((pBuffer[i] < 1.189) ? pBuffer[i]
                              : pBuffer[i] * XimBeamHardModel(pBuffer[i])) -
        1.47);
  }
}
void CbctRecon::DoBeamHardeningCorrection() {
  if (m_spProjImg3DFloat == nullptr) {
    return;
  }

  // FloatImageType m_spProjImg3D: float image

  // double crntVal = 0.0;
  // double corrF = 0.0;

  // If Hnd: // REMEMBER to change in above inlined functions, here is only for
  // debug !!
  double poly3_a = 6.0e-08;
  double poly3_b = -1.0e-08;
  double poly3_c = -5.0e-07;
  double poly3_d = 8.0e-01;
  double poly3_e = 1.47;

  // double corrVal = 0.0;
  if (hisIsUsed) {
    poly3_a = 9.321e-05;
    poly3_b = -2.609e-03;
    poly3_c = 3.374e-02;
    poly3_d = 9.691e-01;
    poly3_e = 0.0;
  } else if (ximIsUsed) {
    poly3_a = 6.0e-8;
    poly3_b = 9.0e-5;
    poly3_c = 1.0e-2;
    poly3_d = 0.8;
    poly3_e = -1.47;
  }
  // Shortening factor 0.9 is applied
  /*double poly3_a = 11.24e-05;
  double poly3_b = -29.28e-04;
  double poly3_c = 35.48e-03;
  double poly3_d = 9.701e-01;*/

  // Shortening factor 0.7 is applied
  /*double poly3_a = 1.660e-04;
  double poly3_b = -3.699e-03;
  double poly3_c = 3.923e-02;
  double poly3_d = 9.727e-01;*/

  std::cout << "Beam hardening corrF poly curve:" << poly3_a << "	"
            << poly3_b << "	" << poly3_c << "	" << poly3_d << "  "
            << poly3_e << std::endl;

  float *pImgBuffer = m_spProjImg3DFloat->GetBufferPointer();
  const itk::Size<3U> pImgSize =
      m_spProjImg3DFloat->GetLargestPossibleRegion().GetSize();
  const unsigned int nPix = pImgSize[0] * pImgSize[1] * pImgSize[2];

  if (hisIsUsed) { // we have now helped the compiler to allow aggressive
                   // optimization
    hisBeamHardening(pImgBuffer, nPix);
  } else if (ximIsUsed) {
    ximBeamHardening(pImgBuffer, nPix);
  } else {
    hndBeamHardening(pImgBuffer, nPix);
  }
}

void CbctRecon::SLT_DoBHC() {

  std::cout << "Beam hardening correction is under progress.." << std::endl;
  DoBeamHardeningCorrection(); // only for m_spProjImg3D
  SetMaxAndMinValueOfProjectionImage();

  SLT_DrawProjImages();
}

template <typename T> int sgn(T val) { return (T(0) < val) - (val < T(0)); }
double heaviside(double x) { return .5 * sgn(x) + 0.5; }

double fullFan_subFunction(double a, double b, double c, double d, double x) {
  return c - sqrt(abs(pow(a, 2) - pow(x * d - b, 2))) *
                 heaviside(x * d - b + a) * heaviside(-(x * d - b - a));
}

double fullFan_Function(double a, double b, double c, double d, double e,
                        double x) {
  return (fullFan_subFunction(a, b, c, d, x - 3. * e) +
          fullFan_subFunction(a, b, c, d, x - 2. * e) +
          fullFan_subFunction(a, b, c, d, x - e) +
          fullFan_subFunction(a, b, c, d, x) +
          fullFan_subFunction(a, b, c, d, x + e) +
          fullFan_subFunction(a, b, c, d, x + 2. * e) +
          fullFan_subFunction(a, b, c, d, x + 3. * e)) *
         .142857; // = 1/7
}

void CbctRecon::SLT_DoBowtieCorrection() {
  if (m_spProjImg3DFloat == nullptr) {
    return;
  }

  if (ximIsUsed || hisIsUsed) {
    std::cout
        << "Bow tie filtering should not be used for His data or Xim data!!"
        << std::endl;
    return;
  }

  QStringList strList = ui.comboBox_fBTcor->currentText().split(';');

  if (strList.length() != 4 && !ui.checkBox_Fullfan->isChecked()) {
    std::cout << "Wrong number of arguments!" << std::endl
              << "Must be a;b;c;d -> d. / (1 + exp(-b.*(x - a))) + c"
              << std::endl;
    return;
  }
  if (strList.length() != 5 && ui.checkBox_Fullfan->isChecked()) {
    std::cout << "Wrong number of arguments!" << std::endl
              << "Must be a;b;c;d;e -> c - sqrt(abs(a^2-((x+/-e)*d-b)^2)) * "
                 "heaviside((x+/-e)*d-b+a) * heaviside(-((x+/-e)*d-b-a))"
              << std::endl;
    return;
  }
  double poly3_a = strList.at(0).toDouble(); // 264.6; //comboBox_fBTcor
  double poly3_b = strList.at(1).toDouble(); // 0.06258;
  double poly3_c = strList.at(2).toDouble(); // 2.502;
  double poly3_d = strList.at(3).toDouble(); // 1.455;
  double poly3_e;
  if (ui.checkBox_Fullfan->isChecked()) {
    poly3_e = strList.at(4).toDouble();
  }

  FloatImageType::SizeType imgSize =
      m_spProjImg3DFloat->GetLargestPossibleRegion().GetSize();
  if (ui.checkBox_Fullfan->isChecked()) {
    std::cout << "Bow-tie correction curve:" << poly3_d << " / ( 1 + exp(-"
              << poly3_b << " * (x - " << poly3_a << "))) + " << poly3_c
              << " + " << poly3_d << " / ( 1 + exp(" << poly3_b << " * (x - "
              << poly3_e << ")))" << std::endl;
  } else {
    std::cout << "Bow-tie correction curve:" << poly3_d << " / ( 1 + exp(-"
              << poly3_b << " * (x - " << poly3_a << "))) + " << poly3_c
              << std::endl;
  }

  using iteratorType = itk::ImageRegionIteratorWithIndex<FloatImageType>;
  iteratorType it(m_spProjImg3DFloat, m_spProjImg3DFloat->GetRequestedRegion());

  double crntVal = 0.0;
  double corrF = 0.0;
  double x_idx;
  it.GoToBegin();
  if (ui.checkBox_Fullfan->isChecked()) {
    while (!it.IsAtEnd()) {
      crntVal = static_cast<double>(it.Get()); // (65535 / exp(it.Get())); //raw
                                               // mu_t = ln(65535/I) <-> I =
                                               // 65535 / exp(mu_t)
      x_idx = static_cast<double>(it.GetIndex()[0]) *
              (512.0 / imgSize[0]); // 512 from current fit -> conversion to be
                                    // consistent with downResFactor
      // if (crntVal > (poly3_c - poly3_a)), negative values are fine, don't
      // worry
      corrF =
          fullFan_Function(poly3_a, poly3_b, poly3_c, poly3_d, poly3_e, x_idx);
      it.Set(static_cast<float>(crntVal -
                                corrF)); // (log(65535 / (crntVal - corrF))));
      ++it;
    }
  } else {
    while (!it.IsAtEnd()) {
      crntVal = static_cast<double>(it.Get()); // (65535 / exp(it.Get())); //raw
                                               // mu_t = ln(65535/I) <-> I =
                                               // 65535 / exp(mu_t)
      x_idx = static_cast<double>(it.GetIndex()[0]) *
              (512.0 / imgSize[0]); // 512 from current fit -> conversion to be
                                    // consistent with downResFactor
      // if (crntVal > poly3_c), negative values are fine, don't worry
      corrF = poly3_d / (1.0 + exp(-poly3_b * (x_idx - poly3_a))) + poly3_c;
      it.Set(static_cast<float>(crntVal -
                                corrF)); // (log(65535 / (crntVal - corrF))));
      ++it;
    }
  }
  SetMaxAndMinValueOfProjectionImage();
  SLT_DrawProjImages();
  std::cout << "Bow-tie correction done." << std::endl;
}

void CbctRecon::SLT_ViewRegistration() // default showing function
{
  m_pDlgRegistration->UpdateListOfComboBox(0); // combo selection signalis
                                               // called
  m_pDlgRegistration->UpdateListOfComboBox(1);
  m_pDlgRegistration->show();
}

void CbctRecon::SLT_ViewHistogram() // default showing function
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

void CbctRecon::Draw2DFrom3DDouble(UShortImageType::Pointer &spFixedImg,
                                   UShortImageType::Pointer &spMovingImg,
                                   enPLANE enPlane, double pos,
                                   YK16GrayImage &YKFixed,
                                   YK16GrayImage &YKMoving) {
  if ((spFixedImg == nullptr) || (spMovingImg == nullptr)) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(
      spFixedImg, spFixedImg->GetRequestedRegion());

  UShortImageType::SizeType imgSize =
      spFixedImg->GetRequestedRegion().GetSize(); // 1016x1016 x z
  // UShortImageType::SizeType imgSizeBuf =
  // spFixedImg->GetBufferedRegion().GetSize(); //1016x1016 x z
  // UShortImageType::SizeType imgSizeLargest =
  // spFixedImg->GetLargestPossibleRegion().GetSize(); //1016x1016 x z

  UShortImageType::PointType imgOrigin = spFixedImg->GetOrigin();
  UShortImageType::SpacingType imgSpacing = spFixedImg->GetSpacing();

  int width = 0;
  int height = 0;
  int iReqSlice = 0;
  int iCntSlice = 0;

  // For moving image
  using ResampleFilterType =
      itk::ResampleImageFilter<UShortImageType, UShortImageType>;
  ResampleFilterType::Pointer filter = ResampleFilterType::New();

  filter->SetInput(spMovingImg);

  using TransformType = itk::AffineTransform<double, 3>;
  TransformType::Pointer transform = TransformType::New();
  filter->SetTransform(transform);

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<UShortImageType, double>;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  filter->SetInterpolator(interpolator);
  filter->SetDefaultPixelValue(0);

  UShortImageType::DirectionType direction;
  direction.SetIdentity();
  filter->SetOutputDirection(direction);

  UShortImageType::SpacingType movingSpacing = imgSpacing;
  UShortImageType::PointType movingOrigin = imgOrigin;
  UShortImageType::SizeType movingSize = imgSize;

  switch (enPlane) {
  case PLANE_AXIAL:
    width = imgSize[0];
    height = imgSize[1];
    iCntSlice = imgSize[2];
    iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?

    movingSpacing[2] = 1.0;
    movingOrigin[2] = pos;
    movingSize[2] = 1;
    // Resample Here! make corresponding 2D image for Moving image
    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);

    break;
  case PLANE_FRONTAL:
    width = imgSize[0];
    height = imgSize[2];
    iCntSlice = imgSize[1];
    iReqSlice = qRound((pos - imgOrigin[1]) / imgSpacing[1]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(2); // y?

    movingSpacing[1] = 1.0;
    movingOrigin[1] = pos;
    movingSize[1] = 1;

    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[2]);
    break;
  case PLANE_SAGITTAL:
    width = imgSize[1];
    height = imgSize[2];
    iCntSlice = imgSize[0];
    iReqSlice = qRound((pos - imgOrigin[0]) / imgSpacing[0]);
    it.SetFirstDirection(1);  // x?
    it.SetSecondDirection(2); // y?

    movingSpacing[0] = 1.0;
    movingOrigin[0] = pos;
    movingSize[0] = 1;

    YKFixed.SetSpacing(imgSpacing[1], imgSpacing[2]);
    break;

  default:
    std::cout << "default should not passed by" << std::endl;
    width = imgSize[0];
    height = imgSize[1];
    iCntSlice = imgSize[2];
    iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?
    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);
    break;
  }

  filter->SetOutputSpacing(movingSpacing);
  filter->SetOutputOrigin(movingOrigin);
  filter->SetSize(movingSize);
  filter->Update();

  YKFixed.CreateImage(width, height, 0);
  // std::cout << "Before MovingImg Creation " << std::endl;

  YKMoving.CreateImage(width, height, 0); // exactly same dimension

  itk::ImageRegionConstIterator<UShortImageType> itMoving(
      filter->GetOutput(), filter->GetOutput()->GetBufferedRegion());
  int cnt = 0;

  // this simple code will cause flip of the image in frontal and sagittal image
  for (itMoving.GoToBegin(); !itMoving.IsAtEnd(); ++itMoving) {
    YKMoving.m_pData[cnt] = itMoving.Get();
    cnt++;
  }
  if (enPlane != PLANE_AXIAL) {
    YKMoving.EditImage_Flip();
  }

  // std::cout << "tot pixel no: " << cnt << std::endl;

  // YK16GrayImage::CopyItkImage2YKImage(filter->GetOutput(), &YKMoving);

  // std::cout << "After MovingImg Creation " << std::endl;

  it.GoToBegin();

  int iNumSlice = 0;
  int iNumWidth = 0;
  int iNumHeight = 0;

  if (iReqSlice < 0 || iReqSlice >= iCntSlice) {
    return;
  }

  while (!it.IsAtEnd()) {
    if (iNumSlice == iReqSlice) {
      iNumHeight = 0;

      while (!it.IsAtEndOfSlice()) {
        iNumWidth = 0;
        while (!it.IsAtEndOfLine()) {
          // double tmpVal = it.Get()*multiplyFactor;
          // ShortImageType::PixelType fixedImgVal = it.Get();
          UShortImageType::PixelType fixedImgVal = it.Get();
          UShortImageType::IndexType pixelIdxFixed{};
          UShortImageType::PointType pixelPhysPt;
          pixelIdxFixed = it.GetIndex();

          // spFixedImg->TransformIndexToPhysicalPoint (pixelIdxFixed,
          // pixelPhysPt);
          // spMovingImg->TransformPhysicalPointToIndex(pixelPhysPt,pixelIdxMoving);

          // Fill YKMoving image
          // calculate the position of this iterator

          // double movingImgVal = spMovingImg->GetPixel()

          // double movingImgVal = (double)spMovingImg->GetPixel( pixelIdxMoving
          // );  unsigned short movingImgVal = 0;

          if (enPlane == PLANE_AXIAL) {
            YKFixed.m_pData[iNumWidth + width * iNumHeight] = fixedImgVal;
            // YKMoving.m_pData[iNumWidth + width*iNumHeight] = movingImgVal;
          } else {
            YKFixed.m_pData[iNumWidth + width * (height - iNumHeight - 1)] =
                fixedImgVal;
            // YKMoving.m_pData[iNumWidth + width*(height - iNumHeight-1)] =
            // movingImgVal;
          }

          ++it;
          iNumWidth++;
        }
        it.NextLine();
        iNumHeight++;
      }
      break;
    }
    it.NextSlice();
    iNumSlice++;
  }

  //    std::cout << "YK Images were filled" << std::endl;
}

// Actually just an overload of the above function
void CbctRecon::Draw2DFrom3DDouble(UShortImageType::Pointer &spFixedImg,
                                   UShortImageType::Pointer &spMovingImg,
                                   enPLANE enPlane, double pos,
                                   AG17RGBAImage &YKFixed,
                                   AG17RGBAImage &YKMoving) {
  if ((spFixedImg == nullptr) || (spMovingImg == nullptr)) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(
      spFixedImg, spFixedImg->GetRequestedRegion());

  UShortImageType::SizeType imgSize =
      spFixedImg->GetRequestedRegion().GetSize(); // 1016x1016 x z
  // UShortImageType::SizeType imgSizeBuf =
  // spFixedImg->GetBufferedRegion().GetSize(); //1016x1016 x z
  // UShortImageType::SizeType imgSizeLargest =
  // spFixedImg->GetLargestPossibleRegion().GetSize(); //1016x1016 x z

  UShortImageType::PointType imgOrigin = spFixedImg->GetOrigin();
  UShortImageType::SpacingType imgSpacing = spFixedImg->GetSpacing();

  int width = 0;
  int height = 0;
  int iReqSlice = 0;
  int iCntSlice = 0;

  // For moving image
  using ResampleFilterType =
      itk::ResampleImageFilter<UShortImageType, UShortImageType>;
  ResampleFilterType::Pointer filter = ResampleFilterType::New();

  filter->SetInput(spMovingImg);

  using TransformType = itk::AffineTransform<double, 3>;
  TransformType::Pointer transform = TransformType::New();
  filter->SetTransform(transform);

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<UShortImageType, double>;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  filter->SetInterpolator(interpolator);
  filter->SetDefaultPixelValue(0);

  // const double outputSpacing[2] = { 1.0, 1.0 };
  // const double outputOrigin[2] = { 0.0, 0.0 };

  UShortImageType::DirectionType direction;
  direction.SetIdentity();
  filter->SetOutputDirection(direction);

  // ResampledImgType2D::SizeType outSize;

  UShortImageType::SpacingType movingSpacing = imgSpacing;
  UShortImageType::PointType movingOrigin = imgOrigin;
  UShortImageType::SizeType movingSize = imgSize;

  switch (enPlane) {
  case PLANE_AXIAL:
    width = imgSize[0];
    height = imgSize[1];
    iCntSlice = imgSize[2];
    iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?

    movingSpacing[2] = 1.0;
    movingOrigin[2] = pos;
    movingSize[2] = 1;
    // Resample Here! make corresponding 2D image for Moving image
    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);

    break;
  case PLANE_FRONTAL:
    width = imgSize[0];
    height = imgSize[2];
    iCntSlice = imgSize[1];
    iReqSlice = qRound((pos - imgOrigin[1]) / imgSpacing[1]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(2); // y?

    movingSpacing[1] = 1.0;
    movingOrigin[1] = pos;
    movingSize[1] = 1;

    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[2]);
    break;
  case PLANE_SAGITTAL:
    width = imgSize[1];
    height = imgSize[2];
    iCntSlice = imgSize[0];
    iReqSlice = qRound((pos - imgOrigin[0]) / imgSpacing[0]);
    it.SetFirstDirection(1);  // x?
    it.SetSecondDirection(2); // y?

    movingSpacing[0] = 1.0;
    movingOrigin[0] = pos;
    movingSize[0] = 1;

    YKFixed.SetSpacing(imgSpacing[1], imgSpacing[2]);
    break;

  default:
    std::cout << "default should not passed by" << std::endl;
    width = imgSize[0];
    height = imgSize[1];
    iCntSlice = imgSize[2];
    iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?
    YKFixed.SetSpacing(imgSpacing[0], imgSpacing[1]);
    break;
  }

  filter->SetOutputSpacing(movingSpacing);
  filter->SetOutputOrigin(movingOrigin);
  filter->SetSize(movingSize);
  filter->Update();

  YKFixed.CreateImage(width, height, 0);
  // std::cout << "Before MovingImg Creation " << std::endl;

  YKMoving.CreateImage(width, height, 0); // exactly same dimension

  itk::ImageRegionConstIterator<UShortImageType> itMoving(
      filter->GetOutput(), filter->GetOutput()->GetBufferedRegion());
  int cnt = 0;

  // this simple code will cause flip of the image in frontal and sagittal image
  for (itMoving.GoToBegin(); !itMoving.IsAtEnd(); ++itMoving) {
    YKMoving.m_pData[cnt] = itMoving.Get();
    cnt++;
  }
  if (enPlane != PLANE_AXIAL) {
    YKMoving.EditImage_Flip();
  }

  it.GoToBegin();

  int iNumSlice = 0;
  int iNumWidth = 0;
  int iNumHeight = 0;

  if (iReqSlice < 0 || iReqSlice >= iCntSlice) {
    return;
  }

  while (!it.IsAtEnd()) {
    if (iNumSlice == iReqSlice) {
      iNumHeight = 0;

      while (!it.IsAtEndOfSlice()) {
        iNumWidth = 0;
        while (!it.IsAtEndOfLine()) {
          UShortImageType::PixelType fixedImgVal = it.Get();

          if (enPlane == PLANE_AXIAL) {
            YKFixed.m_pData[iNumWidth + width * iNumHeight] = fixedImgVal;
          } else {
            YKFixed.m_pData[iNumWidth + width * (height - iNumHeight - 1)] =
                fixedImgVal;
          }

          ++it;
          iNumWidth++;
        }
        it.NextLine();
        iNumHeight++;
      }
      break;
    }
    it.NextSlice();
    iNumSlice++;
  }

  //    std::cout << "YK Images were filled" << std::endl;
}

void CbctRecon::Draw2DFrom3D(UShortImageType::Pointer &pImg, enPLANE direction,
                             double pos, YK16GrayImage &Output2D) {
  if (pImg == nullptr) {
    return;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it(
      pImg, pImg->GetRequestedRegion());

  UShortImageType::SizeType imgSize =
      pImg->GetRequestedRegion().GetSize(); // 1016x1016 x z
  // UShortImageType::SizeType imgSizeBuf = pImg->GetBufferedRegion().GetSize();
  // //1016x1016 x z UShortImageType::SizeType imgSizeLargest =
  // pImg->GetLargestPossibleRegion().GetSize(); //1016x1016 x z

  UShortImageType::PointType imgOrigin = pImg->GetOrigin();
  UShortImageType::SpacingType imgSpacing = pImg->GetSpacing();

  int width = 0;
  int height = 0;
  int iReqSlice = 0;
  int iCntSlice = 0;

  switch (direction) {
  case PLANE_AXIAL:
    width = imgSize[0];
    height = imgSize[1];
    iCntSlice = imgSize[2];
    iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?
    Output2D.SetSpacing(imgSpacing[0], imgSpacing[1]);
    break;
  case PLANE_FRONTAL:
    width = imgSize[0];
    height = imgSize[2];
    iCntSlice = imgSize[1];
    iReqSlice = qRound((pos - imgOrigin[1]) / imgSpacing[1]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(2); // y?
    Output2D.SetSpacing(imgSpacing[0], imgSpacing[2]);
    break;
  case PLANE_SAGITTAL:
    width = imgSize[1];
    height = imgSize[2];
    iCntSlice = imgSize[0];
    iReqSlice = qRound((pos - imgOrigin[0]) / imgSpacing[0]);
    it.SetFirstDirection(1);  // x?
    it.SetSecondDirection(2); // y?
    Output2D.SetSpacing(imgSpacing[1], imgSpacing[2]);
    break;
  default:
    std::cout << "default should not be passed by" << std::endl;
    width = imgSize[0];
    height = imgSize[1];
    iCntSlice = imgSize[2];
    iReqSlice = qRound((pos - imgOrigin[2]) / imgSpacing[2]);
    it.SetFirstDirection(0);  // x?
    it.SetSecondDirection(1); // y?
    Output2D.SetSpacing(imgSpacing[0], imgSpacing[1]);
    break;
  }

  Output2D.CreateImage(width, height, 0);
  it.GoToBegin();

  int iNumSlice = 0;
  int iNumWidth = 0;
  int iNumHeight = 0;

  if (iReqSlice < 0 || iReqSlice >= iCntSlice) {
    return;
  }

  while (!it.IsAtEnd()) {
    if (iNumSlice == iReqSlice) {
      iNumHeight = 0;
      while (!it.IsAtEndOfSlice()) {
        iNumWidth = 0;
        while (!it.IsAtEndOfLine()) {
          // double tmpVal = it.Get()*multiplyFactor;
          double tmpVal = it.Get();

          if (direction == PLANE_AXIAL) {
            Output2D.m_pData[iNumWidth + width * iNumHeight] = tmpVal;
          } else {
            Output2D.m_pData[iNumWidth + width * (height - iNumHeight - 1)] =
                tmpVal;
          }

          ++it;
          iNumWidth++;
        }
        it.NextLine();
        iNumHeight++;
      }
      break;
    }
    it.NextSlice();
    iNumSlice++;
  }
}

void CbctRecon::RegisterImgDuplication(enREGI_IMAGES src,
                                       enREGI_IMAGES target) {

  UShortImageType::Pointer tmpSrc;

  switch (src) {
  case REGISTER_REF_CT:
    tmpSrc = m_spRefCTImg;
    break;
  default:
    std::cerr << "You are using a non-valid target!" << std::endl;
    return;
  }

  if (tmpSrc == nullptr) {
    std::cout << "src image is empty" << std::endl;
    return;
  }

  // Duplication for registration. Starting point is manual Rigid CT image
  using DuplicatorType = itk::ImageDuplicator<UShortImageType>;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(tmpSrc);
  duplicator->Update();

  switch (target) {
  case REGISTER_MANUAL_RIGID:
    m_spManualRigidCT = duplicator->GetOutput();
    break;
  default:
    std::cerr << "You are using a non-valid target!" << std::endl;
    return;
  }
  // Duplication for : End
}
void CbctRecon::FindAllRelevantPaths(
    const QString &pathProjHisDir) // called following SLT_SetHisDir
{
  // in case of eletka, img_UID
  // QString aa;
  // std::cout<< "ddd " << aa.toLocal8Bit().constData() << std::endl;

  m_strDCMUID = "";
  m_strPathPatientDir = "";
  m_strPatientDirName = "";
  m_strPathFRAME_DBF = "";
  m_strPathIMAGE_DBF = "";
  m_strPathGeomXML = "";
  m_strPathPlanCTDir = "";
  m_strPathRS = "";
  m_strPathRS_CBCT = "";
  m_strPathElektaINI = "";
  m_strPathElektaINIXVI2 = "";

  m_strPathPlan = "";

  m_strPathIMAGES = "";

  QDir curHisDir(pathProjHisDir);
  QDir movingDir(pathProjHisDir);
  hisIsUsed = true;

  if (!curHisDir.dirName().contains("img_", Qt::CaseSensitive) &&
      !curHisDir.dirName().contains("fwd_", Qt::CaseSensitive) &&
      !curHisDir.dirName().contains("sca_", Qt::CaseSensitive) &&
      !curHisDir.dirName().contains("cor_", Qt::CaseSensitive)) {
    if (curHisDir.dirName().contains("Scan0", Qt::CaseSensitive)) {
      std::cout << "XML set by guessing: Scan0/../ProjectionInfo.xml"
                << std::endl;
      m_strPathGeomXML =
          curHisDir.absolutePath() + "/../" + "ProjectionInfo.xml";
      ui.lineEdit_ElektaGeomPath->setText(m_strPathGeomXML);
      std::cout << "Patient DIR set to: Scan0/../../" << std::endl;
      m_strPathGeomXML = curHisDir.absolutePath() + "/../";
      ximIsUsed = false;
      hisIsUsed = false;
      return;
    }
    if (curHisDir.absolutePath().contains("Acquisitions", Qt::CaseSensitive)) {
      std::cout << "XML set by guessing: Acquisitions/../Scan.xml" << std::endl;
      m_strPathGeomXML = curHisDir.absolutePath() + "/../../" + "Scan.xml";
      ui.lineEdit_ElektaGeomPath->setText(m_strPathGeomXML);
      std::cout << "Patient DIR set to: Acquisitions/../../../" << std::endl;
      m_strPathGeomXML = curHisDir.absolutePath() + "/../../";
      ximIsUsed = true;
      hisIsUsed = false;
      return;
    }

    std::cout << "Projection folder should have format [img_UID]" << std::endl;
    std::cout << "XML file cannot be made" << std::endl;
    return;
  }

  QString tmpStr = curHisDir.dirName();
  QStringList strListDir = tmpStr.split("_");
  m_strDCMUID = strListDir.at(1);

  // m_strDCMUID = curHisDir.dirName().right(curHisDir.dirName().length() - 4);

  if (!movingDir.cdUp()) // projDir ==> IMAGES
  {
    std::cout << "no upper dir" << std::endl;
    return;
  }
  QDir tmpDir_IMAGES(movingDir.absolutePath());
  m_strPathIMAGES = tmpDir_IMAGES.absolutePath();

  if (!movingDir.cdUp()) // IMAGES ==> patient_402-02-78
  {
    std::cout << "no upper dir" << std::endl;
    return;
  }
  QDir tmpDir_PatientFolder(movingDir.absolutePath());

  if (!movingDir
           .cdUp()) // patient_402-02-78 ==> Data folder where DBF files are.
  {
    std::cout << "no upper dir" << std::endl;
    return;
  }

  m_strPatientDirName = tmpDir_PatientFolder.dirName();
  m_strPathPatientDir = tmpDir_PatientFolder.absolutePath();

  QDir tmpDir_RootFolder(movingDir.absolutePath()); // root folder

  if (tmpDir_RootFolder.absolutePath().length() > 1) {
    m_strPathDirDefault = tmpDir_RootFolder.absolutePath();
  }

  // option 1: already made rtk xml file
  QString tmpPathRTKGeometry = tmpDir_RootFolder.absolutePath() + "/" +
                               "ElektaGeom_" + m_strDCMUID + ".xml";
  QFileInfo rtkGeomInfo(tmpPathRTKGeometry);

  // option 2
  QString pathXVIGeometryXML = curHisDir.absolutePath() + "/" + "_Frames.xml";
  QFileInfo xviGeomInfo(pathXVIGeometryXML);

  if (rtkGeomInfo
          .exists()) // The best option: rtk geometry file is already existing
  {
    std::cout << "RTK XLM file is found" << std::endl;
    m_strPathGeomXML = tmpPathRTKGeometry;
  } else if (xviGeomInfo.exists()) // 2nd option:_Frames.xml already exists in
                                   // each projection folder for > XVI5.0.2
  {
    std::cout << "XVI XLM file is found" << std::endl;
    // YKdebug: Later, it
    m_strPathGeomXML = pathXVIGeometryXML;

    // Later, it should be

    // XVIXMLReader.SetFile(pathXVIGeometryXML)
    // XVIXMLReader.GenerateData()

    // rtk::ThreeDCircularProjectionGeometryXMLFileWriter::Pointer xmlWriter
    // xmlWriter->SetFilename( tmpPathRTKGeometry ) //as option 1
    // xmlWriter->SetObject(XVIXMLReader->GetGeometry());
    // TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile());
    // m_strPathGeomXML = tmpPathRTKGeometry;
  } else {
    QString tmpStrPath1 = m_strPathPatientDir;
    QString tmpStrPath2 = m_strPathPatientDir;

    // 1st priority: DBF files saved in each "patient" folder --> just in case
    // data are collected separately
    QFileInfo fInfo_FrameDBF = QFileInfo(tmpStrPath1.append("/FRAME.DBF"));
    QFileInfo fInfo_ImageDBF = QFileInfo(tmpStrPath2.append("/IMAGE.DBF"));

    if (!fInfo_FrameDBF.exists() || !fInfo_ImageDBF.exists()) {
      std::cout << "No found in the patient folder. DBF files can be saved in "
                   "each individual patient as well. Continues to search them "
                   "again in root folder(standard)"
                << std::endl;

      fInfo_FrameDBF =
          QFileInfo(tmpDir_RootFolder.absolutePath().append("/FRAME.DBF"));
      fInfo_ImageDBF =
          QFileInfo(tmpDir_RootFolder.absolutePath().append("/IMAGE.DBF"));

      if (!fInfo_FrameDBF.exists() || !fInfo_ImageDBF.exists()) {
        std::cout << "DBF files were not found" << std::endl;
        std::cout << "XML file cannot be made" << std::endl;
        return;
      }
    } else {
      std::cout << "DBF files are found in the individual patient directory."
                << std::endl;
    }
    m_strPathFRAME_DBF = fInfo_FrameDBF.absoluteFilePath();
    m_strPathIMAGE_DBF = fInfo_ImageDBF.absoluteFilePath();

    m_strPathGeomXML = "";
    m_strPathGeomXML = MakeElektaXML(
        m_strPathIMAGE_DBF, m_strPathFRAME_DBF,
        m_strDCMUID); // if DBF files exist but UID is not found, it will crash

    if (m_strPathGeomXML.length() < 1) {
      std::cout << "No releated data in DBF file" << std::endl;
      return;
    }
  }

  // std::cout << "Root folder: " <<
  // tmpDir_RootFolder.absolutePath().toLocal8Bit().constData() << std::endl;
  // std::cout << "Root folder2: " <<
  // tmpDir_RootFolder.absolutePath().append("\\FRAME.DBF").toLocal8Bit().constData()
  // << std::endl;  std::cout <<
  // fInfo_FrameDBF.absoluteFilePath().toLocal8Bit().constData() << std::endl;

  // GenerateXMLFunc.

  // Search for the geometry XML file: naming convention:
  // ElektaGeom_DICOMUID.xml  after Generation of the XML from DBF files
  //  movingDir = tmpDir_IMAGES;
  // tmpDir_PatientFolder;

  int enDirStructure_Type = 0;
  // 0: in patient DIR --> 3 folders(CT_SET, DICOM_PLAN, IMAGES)
  // 1: // Patient DIR ==> IMAGES --> CT_SET / DICOM_PLAN
  // 2: NO CT image

  QDir tmpDIR_CTSET =
      QDir(tmpDir_PatientFolder.absolutePath().append("/CT_SET"));

  if (tmpDIR_CTSET.exists()) {
    enDirStructure_Type = 0;
  } else {
    QString tmpStrPathCTSET = m_strPathIMAGES;
    tmpDIR_CTSET = QDir(tmpStrPathCTSET.append("/CT_SET"));

    if (tmpDIR_CTSET.exists()) {
      enDirStructure_Type = 1;
    } else {
      enDirStructure_Type = 2;
    }
  }

  // QString strPathCTSet = m_strPathIMAGES.append("/CT_SET");

  // switch (enDirStructure_Type)
  // {
  // case 0:
  //  movingDir = tmpDir_IMAGES;
  // break;
  // case 1:
  // movingDir = tmpDir_IMAGES;
  // break;
  // case 2:
  // std::cout << "No CT DICOM folder exist. Proceeding w/o CT" << std::endl;
  // break;
  // }
  //

  QFileInfoList listDir;
  if (enDirStructure_Type != 2) {
    listDir = tmpDIR_CTSET.entryInfoList(QDir::Dirs, QDir::Name);
    if (listDir.size() <= 2) // only /. and /.. exist
    {
      std::cout << "No CT DICOM folder exist. Proceeding w/o CT" << std::endl;
    } else {
      m_strPathPlanCTDir =
          listDir.at(2).absoluteFilePath(); // . , .. , real DICOM folder
    }
  }

  // for (int i = 0 ; i < listDir.size() ; i++)
  // {
  ////std::cout << listDir.at(i).absolutePath().toLocal8Bit().constData() <<
  /// std::endl; //this returns Dir, not itself
  // QString tmpPath = listDir.at(i).absoluteFilePath();
  // }

  QFileInfoList listFile = tmpDIR_CTSET.entryInfoList(
      QDir::Files, QDir::Name); // search for DICOM RS file

  if (listFile.empty()) {
    std::cout << "No CT DICOM RS file exist. proceeding w/o RS" << std::endl;
    // return;
  } else {
    for (const auto &i : listFile) {
      if (i.suffix().contains("DCM", Qt::CaseInsensitive)) {
        m_strPathRS = i.absoluteFilePath();
        break;
      }
    }
  }

  QDir tmpDIR_DCM_Plan =
      QDir(tmpDir_PatientFolder.absolutePath().append("/DICOM_PLAN"));

  if (tmpDIR_DCM_Plan.exists()) {
    enDirStructure_Type = 0;
  } else {
    QString tmpStrPathCTSET = m_strPathIMAGES;
    tmpDIR_DCM_Plan = QDir(tmpStrPathCTSET.append("/DICOM_PLAN"));

    if (tmpDIR_DCM_Plan.exists()) {
      enDirStructure_Type = 1;
    } else {
      enDirStructure_Type = 2;
    }
  }

  QFileInfoList listFileDCMPlan;
  if (enDirStructure_Type != 2) {
    listFileDCMPlan = tmpDIR_DCM_Plan.entryInfoList(QDir::Files, QDir::Name);
    if (listFileDCMPlan.empty()) // should be /. and /.. and one dcm file
    {
      std::cout << "No DCM plan file exists. Proceeding w/o dicom plan"
                << std::endl;
    } else if (listFileDCMPlan.size() > 1) {
      std::cout << "Warning! More than one DCM plan file exists. First DCM "
                   "plan file will be used"
                << std::endl;
    } else {
    }

    for (const auto &i : listFileDCMPlan) {
      if (i.suffix().contains("DCM", Qt::CaseInsensitive)) {
        m_strPathPlan = i.absoluteFilePath();
        break; // get fisrt one only
      }
    }
  }

  QDir movingDirCBCTRS;

  if (enDirStructure_Type == 0) {
    movingDirCBCTRS = tmpDir_PatientFolder;
  } else if (enDirStructure_Type == 1) {
    movingDirCBCTRS = tmpDir_IMAGES;
  }

  if (!movingDirCBCTRS.cd("CBCT_RS")) {
    std::cout << "no CBCT_RS dir exists. Proceed with out CBCT RS image"
              << std::endl;
  } else {
    QFileInfoList listFile2 =
        movingDirCBCTRS.entryInfoList(QDir::Files, QDir::Name);

    if (listFile2.empty()) {
      std::cout << "No CBCT DICOM RS file exist. proceeding w/o RS"
                << std::endl;
      // return;
      m_strPathRS_CBCT = "";
    } else {
      for (const auto &i : listFile2) {
        if (i.suffix().contains("DCM", Qt::CaseInsensitive)) {
          m_strPathRS_CBCT = i.absoluteFilePath();
          break;
        }
      }
    }
  }
  ui.lineEdit_PathCBCTSkinPath->setText(m_strPathRS_CBCT);

  QString strPathAcqParamDir = pathProjHisDir + "/Reconstruction";
  QDir tmpAcqParamDir = QDir(strPathAcqParamDir);

  if (tmpAcqParamDir.exists()) {
    QFileInfoList listFileAcqParam = tmpAcqParamDir.entryInfoList(
        QDir::Files, QDir::Name); // search for DICOM RS file

    int iMinNameLength = 9999;

    int iMaxNameLength = 0;
    int iCnt_INIXVI = 0;

    QString strPathINIXVI_long;
    for (const auto &i : listFileAcqParam) {
      // suffix:*.tar.gz ==> gz only
      if (i.suffix().contains("INI", Qt::CaseInsensitive)) {
        QString tmpPath = i.absoluteFilePath();

        if (tmpPath.length() < iMinNameLength) {
          iMinNameLength = tmpPath.length();
          m_strPathElektaINI = tmpPath;
        }
      }

      QString StrSuffix = i.completeSuffix();

      if (StrSuffix.contains("INI.XVI", Qt::CaseInsensitive)) {
        iCnt_INIXVI++;

        QString tmpPath2 = i.absoluteFilePath();

        if (tmpPath2.length() > iMaxNameLength) {
          iMaxNameLength = tmpPath2.length();
          strPathINIXVI_long = tmpPath2;
        }
      }
    }

    if (iCnt_INIXVI == 2) {
      m_strPathElektaINIXVI2 = strPathINIXVI_long;
    }
  }

  std::cout << "m_strDCMUID: " << m_strDCMUID.toLocal8Bit().constData()
            << std::endl;
  std::cout << "m_strPathPatientDir: "
            << m_strPathPatientDir.toLocal8Bit().constData() << std::endl;
  std::cout << "m_strPatientDirName: "
            << m_strPatientDirName.toLocal8Bit().constData() << std::endl;
  std::cout << "m_strPathFRAME_DBF: "
            << m_strPathFRAME_DBF.toLocal8Bit().constData() << std::endl;
  std::cout << "m_strPathIMAGE_DBF: "
            << m_strPathIMAGE_DBF.toLocal8Bit().constData() << std::endl;
  std::cout << "m_strPathGeomXML: "
            << m_strPathGeomXML.toLocal8Bit().constData() << std::endl;
  std::cout << "m_strPathPlanCTDir: "
            << m_strPathPlanCTDir.toLocal8Bit().constData() << std::endl;
  std::cout << "m_strPathRS: " << m_strPathRS.toLocal8Bit().constData()
            << std::endl;
  std::cout << "m_strPathRS_CBCT: "
            << m_strPathRS_CBCT.toLocal8Bit().constData() << std::endl;
  std::cout << "m_strPathPlan: " << m_strPathPlan.toLocal8Bit().constData()
            << std::endl;
  std::cout << "m_strPathElektaINI: "
            << m_strPathElektaINI.toLocal8Bit().constData() << std::endl;
  std::cout << "m_strPathElektaINIXVI2: "
            << m_strPathElektaINIXVI2.toLocal8Bit().constData() << std::endl;

  float kVp = 0.0;
  float mA = 0.0;
  float ms = 0.0;
  GetXrayParamFromINI(m_strPathElektaINI, kVp, mA, ms);

  if (kVp * mA * ms != 0) {
    // update GUI
    std::cout << "Updating current mAs setting from INI file: "
              << "kVp= " << kVp << ", mA= " << mA << ", ms= " << ms
              << std::endl;
    ui.lineEdit_CurmAs->setText(QString("%1, %2").arg(mA).arg(ms));
  }

  VEC3D couch_trans = {-999, -999,
                       -999}; // mm. In the text file, these values are in cm.
  VEC3D couch_rot = {-999, -999,
                     -999}; // mm. In the text file, these values are in cm.

  bool res =
      GetCouchShiftFromINIXVI(m_strPathElektaINIXVI2, &couch_trans, &couch_rot);

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

  // YKTEMP: delete if the UI is changed
  ui.lineEdit_ElektaGeomPath->setText(m_strPathGeomXML);
}

void CbctRecon::init_DlgRegistration(
    QString &strDCM_UID) // init dlgRegistrations
{
  m_pDlgRegistration->initDlgRegistration(
      strDCM_UID); // NULLing all temporary spImage
}

// output spProjCT3D => intensity value, not line integral
void CbctRecon::ForwardProjection(UShortImageType::Pointer &spVolImg3D,
                                  GeometryType::Pointer &spGeometry,
                                  UShortImageType::Pointer &spProjCT3D,
                                  bool bSave) {
  if (spVolImg3D == nullptr) {
    std::cout << "ERROR! No 3D-CT file. Load 3D CT file first" << std::endl;
    return;
  }

  if (m_iCntSelectedProj < 1 && bSave) {
    std::cout << "Error! No projection image is loaded" << std::endl;
    return;
  }

  if (spGeometry->GetGantryAngles().empty()) {
    std::cout << "No geometry!" << std::endl;
    return;
  }

#ifndef USE_CUDA
  CPU_ForwardProjection(spVolImg3D, spGeometry, spProjCT3D,
                        bSave); // final moving image
#else
  if (ui.radioButton_UseCUDA->isChecked()) {
    CUDA_ForwardProjection(spVolImg3D, spGeometry, spProjCT3D,
                           bSave); // final moving image
  } else {
    CPU_ForwardProjection(spVolImg3D, spGeometry, spProjCT3D,
                          bSave); // final moving image
  }
#endif // !USE_CUDA
}

#if USE_CUDA
// output spProjCT3D => intensity value, not line integral
void CbctRecon::CUDA_ForwardProjection(UShortImageType::Pointer &spVolImg3D,
                                       GeometryType::Pointer &spGeometry,
                                       UShortImageType::Pointer &spProjCT3D,
                                       bool bSave) {

  // m_spProjCTImg --> spProjCT3D

  FloatImageType::Pointer spResultProjImageFloat;
  // Euler Transformation for RTK's weird orientation

  // int iNumOfProjections = 0;

  {
    // 0) CT image Transformation
    UShortImageType::SizeType size_original =
        spVolImg3D->GetLargestPossibleRegion().GetSize();
    UShortImageType::SpacingType spacing_original = spVolImg3D->GetSpacing();

    // Same image type from original image -3D & float
    UShortImageType::IndexType start_trans{};
    start_trans[0] = 0;
    start_trans[1] = 0;
    start_trans[2] = 0;

    UShortImageType::SizeType size_trans{};
    size_trans[0] = size_original[1]; // X //512
    size_trans[1] = size_original[2]; // Y  //512
    size_trans[2] = size_original[0]; // Z //300

    std::cout << " size_trans" << size_trans << std::endl;

    UShortImageType::SpacingType spacing_trans;
    spacing_trans[0] = spacing_original[1];
    spacing_trans[1] = spacing_original[2];
    spacing_trans[2] = spacing_original[0];

    std::cout << " spacing_trans" << spacing_trans << std::endl;

    UShortImageType::PointType Origin_trans;
    Origin_trans[0] = -0.5 * size_trans[0] * spacing_trans[0];
    Origin_trans[1] = -0.5 * size_trans[1] * spacing_trans[1];
    Origin_trans[2] = -0.5 * size_trans[2] * spacing_trans[2];

    UShortImageType::RegionType region_trans;
    region_trans.SetSize(size_trans);
    region_trans.SetIndex(start_trans);

    using FilterType = itk::FlipImageFilter<UShortImageType>;
    FilterType::Pointer flipFilter = FilterType::New();
    using FlipAxesArrayType = FilterType::FlipAxesArrayType;

    FlipAxesArrayType arrFlipAxes;
    arrFlipAxes[0] = true;
    arrFlipAxes[1] = false;
    arrFlipAxes[2] = false;

    flipFilter->SetFlipAxes(arrFlipAxes);
    flipFilter->SetInput(spVolImg3D); // plan CT, USHORT image

    using TransformType = itk::Euler3DTransform<double>;
    TransformType::Pointer transform = TransformType::New();

    TransformType::ParametersType param;
    param.SetSize(6);
    param.put(0, itk::Math::pi / -2.0); // rot X // 0.5 = PI/2
    param.put(1, 0);                    // rot Y
    param.put(2, itk::Math::pi / 2.0);  // rot Z
    param.put(3, 0.0);                  // Trans X mm
    param.put(4, 0.0);                  // Trans Y mm
    param.put(5, 0.0);                  // Trans Z mm

    TransformType::ParametersType fixedParam(3); // rotation center
    fixedParam.put(0, 0);
    fixedParam.put(1, 0);
    fixedParam.put(2, 0);

    transform->SetParameters(param);
    transform->SetFixedParameters(fixedParam); // Center of the Transform

    std::cout << "Transform matrix:"
              << "\n"
              << transform->GetMatrix() << std::endl;

    using ResampleFilterType =
        itk::ResampleImageFilter<UShortImageType, UShortImageType>;
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    resampler->SetInput(flipFilter->GetOutput());
    resampler->SetSize(size_trans);
    resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
    resampler->SetOutputSpacing(spacing_trans); // 1 1 1
    resampler->SetOutputDirection(
        flipFilter->GetOutput()->GetDirection()); // image normal?
    resampler->SetTransform(transform);

    using CastFilterType =
        itk::CastImageFilter<UShortImageType,
                             FloatImageType>; // Maybe not inplace filter
    CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(resampler->GetOutput());

    // Default value
    double calibF_A = 1.0;
    double calibF_B = 0.0;

    std::cout << "Temporary forcing CT# applied for tissue" << std::endl;

    std::cout << "CBCT calibration Factor(Recommended: 1, 0): A = " << calibF_A
              << "  B= " << calibF_B << std::endl;
    using MultiplyImageFilterType =
        itk::MultiplyImageFilter<FloatImageType, FloatImageType,
                                 FloatImageType>;
    MultiplyImageFilterType::Pointer multiplyImageFilter =
        MultiplyImageFilterType::New();
    multiplyImageFilter->SetInput(castFilter->GetOutput());
    multiplyImageFilter->SetConstant(calibF_A / 65535.0);

    using AddImageFilterType =
        itk::AddImageFilter<FloatImageType, FloatImageType, CUDAFloatImageType>;
    AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
    addImageFilter->SetInput1(multiplyImageFilter->GetOutput());
    double addingVal = calibF_B / 65535.0;
    addImageFilter->SetConstant2(addingVal);
    addImageFilter->Update(); // will generate map of real_mu (att.coeff)

    CUDAFloatImageType::Pointer spCTImg_mu;
    spCTImg_mu = addImageFilter->GetOutput();

    // 2) Prepare empty projection images //Should be corresonponding to raw
    // projection images

    // Create a stack of empty projection images
    using ConstantImageSourceType =
        rtk::ConstantImageSource<CUDAFloatImageType>; // Output: FLoat image =
                                                      // may be mu_t =
                                                      // log(I_0/I)
    ConstantImageSourceType::Pointer constantImageSource =
        ConstantImageSourceType::New();

    ConstantImageSourceType::SizeType size{};
    ConstantImageSourceType::SpacingType spacing;
    ConstantImageSourceType::PointType origin;

    std::cout << "Setting-up vacant projection image data" << std::endl;

    // a) size
    // std::cout << "chk1" << std::endl;
    size[0] = m_spProjImg3DFloat->GetBufferedRegion()
                  .GetSize()[0]; // qRound((double)DEFAULT_W*m_fResampleF);
    size[1] = m_spProjImg3DFloat->GetBufferedRegion()
                  .GetSize()[1]; // qRound((double)DEFAULT_H*m_fResampleF);
    size[2] = spGeometry->GetGantryAngles().size();
    // iNumOfProjections = size[2];

    // b) spacing
    spacing[0] = m_fProjSpacingX / m_fResampleF; // typical HIS file
    spacing[1] = m_fProjSpacingY / m_fResampleF;
    spacing[2] = 1.0;

    // c) Origin: can center be the image center? or should be related to the CT
    // image???
    origin[0] = spacing[0] * (size[0] - 1) * -0.5;
    origin[1] = spacing[1] * (size[1] - 1) * -0.5;
    origin[2] = 0.0;

    constantImageSource->SetOrigin(origin);
    constantImageSource->SetSpacing(spacing);

    FloatImageType::DirectionType imageDirection;
    imageDirection.SetIdentity(); // no effect
    constantImageSource->SetDirection(imageDirection);
    constantImageSource->SetSize(size);
    constantImageSource->SetConstant(1.0);
    constantImageSource->UpdateOutputInformation();
    std::cout << "Canvas for projection image is ready to write" << std::endl;

    // 4) Prepare CT image to be projected

    // Create forward projection image filter
    rtk::CudaForwardProjectionImageFilter<CUDAFloatImageType,
                                          CUDAFloatImageType>::Pointer
        CudaForwardProjection; // Float to Float
    CudaForwardProjection =
        rtk::CudaForwardProjectionImageFilter<CUDAFloatImageType,
                                              CUDAFloatImageType>::New();

    itk::TimeProbe projProbe;
    std::cout << "Forward projection is now ongoing" << std::endl;

    CudaForwardProjection->SetInput(
        constantImageSource
            ->GetOutput()); // Canvas. projection image will be saved here.
    CudaForwardProjection->SetInput(1, spCTImg_mu); // reference plan CT image
    CudaForwardProjection->SetGeometry(spGeometry);
    projProbe.Start();
    TRY_AND_EXIT_ON_ITK_EXCEPTION(CudaForwardProjection->Update());
    projProbe.Stop();
    spResultProjImageFloat = CudaForwardProjection->GetOutput();
    std::cout << "Forward projection done in:	" << projProbe.GetMean() << ' '
              << projProbe.GetUnit() << '.' << std::endl;
  } // release all the memory

  // From Float to USHORT and line integral to intensity

  spProjCT3D = UShortImageType::New(); // later
  UShortImageType::SizeType projCT_size =
      spResultProjImageFloat->GetLargestPossibleRegion()
          .GetSize(); // 1024 1024 350
  UShortImageType::IndexType projCT_idxStart =
      spResultProjImageFloat->GetLargestPossibleRegion().GetIndex(); // 0 0 0
  UShortImageType::SpacingType projCT_spacing =
      spResultProjImageFloat->GetSpacing(); // 0.4 0.4 1.0
  UShortImageType::PointType projCT_origin =
      spResultProjImageFloat->GetOrigin(); //-204.6 -204.6 -174.5

  // Copy informations from spResultProjImageFloat
  FloatImageType::RegionType projCT_region;
  projCT_region.SetSize(projCT_size);
  projCT_region.SetIndex(projCT_idxStart);

  spProjCT3D->SetRegions(projCT_region);
  spProjCT3D->SetSpacing(projCT_spacing);
  spProjCT3D->SetOrigin(projCT_origin);

  spProjCT3D->Allocate();
  spProjCT3D->FillBuffer(0);

  // Calculation process
  itk::ImageRegionConstIterator<FloatImageType> itSrc(
      spResultProjImageFloat, spResultProjImageFloat->GetRequestedRegion());
  itk::ImageRegionIterator<UShortImageType> itTarg(
      spProjCT3D, spProjCT3D->GetRequestedRegion()); // writing

  itSrc.GoToBegin();
  itTarg.GoToBegin();

  float fProjVal = 0.0;
  double tmpConvVal = 0.0;

  // Convert line integral to intensity value (I0/I = exp(mu_t)) --> I =
  // I0/exp(mu_t)
  while (!itSrc.IsAtEnd() && !itTarg.IsAtEnd()) {
    fProjVal = itSrc.Get();                 // mu_t //63.5 --> 6.35
    tmpConvVal = (65535.0 / exp(fProjVal)); // physically true

    if (tmpConvVal <= 0.0) {
      itTarg.Set(0);
    } else if (tmpConvVal > 65535.0) {
      itTarg.Set(65535);
    } else {
      itTarg.Set(static_cast<unsigned short>(tmpConvVal));
    }

    ++itSrc;
    ++itTarg;
  }

  // spProjCT3D: USHORT IMAGE of intensity. Not inverted (physical intensity)

  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout << " Patient DIR Path: "
              << m_strPathPatientDir.toLocal8Bit().constData() << std::endl;

    bool manuallySelectedDir = false; // <- just to make sure I don't break
                                      // usecases of the older version.
    if (m_strPathPatientDir.isEmpty()) {
      std::cout << "File save error!: No patient DIR name" << std::endl;

      m_strPathPatientDir = QFileDialog::getExistingDirectory(
          this, tr("Open Directory"), ".",
          QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
      if (m_strPathPatientDir.length() <= 1) {
        return;
      }
      manuallySelectedDir = true;
    }

    // Get current folder
    QString subdir_images("IMAGES");
    QString strCrntDir =
        m_strPathPatientDir + "/" + subdir_images; // current Proj folder

    // Make a sub directory
    QDir crntDir(strCrntDir);

    if (!crntDir.exists()) {
      if (manuallySelectedDir) {
        QDir current_dir(m_strPathPatientDir);
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

    QString fwdDirName = "fwd_" + m_strDCMUID;

    bool tmpResult = crntDir.mkdir(fwdDirName); // what if the directory exists?

    if (!tmpResult) {
      std::cout << "FwdProj directory seems to exist already. Files will be "
                   "overwritten."
                << std::endl;
    }

    QString strSavingFolder = strCrntDir + "/" + fwdDirName;
    SaveProjImageAsHIS(spProjCT3D, m_arrYKBufProj, strSavingFolder,
                       m_fResampleF);
  }
}
#endif

// output spProjCT3D => intensity value, not line integral
void CbctRecon::CPU_ForwardProjection(UShortImageType::Pointer &spVolImg3D,
                                      GeometryType::Pointer &spGeometry,
                                      UShortImageType::Pointer &spProjCT3D,
                                      bool bSave) {
  // m_spProjCTImg --> spProjCT3D

  FloatImageType::Pointer spResultProjImageFloat;
  // Euler Transformation for RTK's weird orientation

  // int iNumOfProjections = 0;

  {
    // 0) CT image Transformation
    UShortImageType::SizeType size_original =
        spVolImg3D->GetLargestPossibleRegion().GetSize();
    UShortImageType::SpacingType spacing_original = spVolImg3D->GetSpacing();

    // Same image type from original image -3D & float
    UShortImageType::IndexType start_trans{};
    start_trans[0] = 0;
    start_trans[1] = 0;
    start_trans[2] = 0;

    UShortImageType::SizeType size_trans{};
    size_trans[0] = size_original[1]; // X //512
    size_trans[1] = size_original[2]; // Y  //512
    size_trans[2] = size_original[0]; // Z //300

    std::cout << " size_trans" << size_trans << std::endl;

    UShortImageType::SpacingType spacing_trans;
    spacing_trans[0] = spacing_original[1];
    spacing_trans[1] = spacing_original[2];
    spacing_trans[2] = spacing_original[0];

    std::cout << " spacing_trans" << spacing_trans << std::endl;

    UShortImageType::PointType Origin_trans;
    Origin_trans[0] = -0.5 * size_trans[0] * spacing_trans[0];
    Origin_trans[1] = -0.5 * size_trans[1] * spacing_trans[1];
    Origin_trans[2] = -0.5 * size_trans[2] * spacing_trans[2];

    UShortImageType::RegionType region_trans;
    region_trans.SetSize(size_trans);
    region_trans.SetIndex(start_trans);

    using FilterType = itk::FlipImageFilter<UShortImageType>;
    FilterType::Pointer flipFilter = FilterType::New();
    using FlipAxesArrayType = FilterType::FlipAxesArrayType;

    FlipAxesArrayType arrFlipAxes;
    arrFlipAxes[0] = true;
    arrFlipAxes[1] = false;
    arrFlipAxes[2] = false;

    flipFilter->SetFlipAxes(arrFlipAxes);
    flipFilter->SetInput(spVolImg3D); // plan CT, USHORT image

    using TransformType = itk::Euler3DTransform<double>;
    TransformType::Pointer transform = TransformType::New();

    TransformType::ParametersType param;
    param.SetSize(6);
    param.put(0, itk::Math::pi / -2.0); // rot X // 0.5 = PI/2
    param.put(1, 0);                    // rot Y
    param.put(2, itk::Math::pi / 2.0);  // rot Z
    param.put(3, 0.0);                  // Trans X mm
    param.put(4, 0.0);                  // Trans Y mm
    param.put(5, 0.0);                  // Trans Z mm

    TransformType::ParametersType fixedParam(3); // rotation center
    fixedParam.put(0, 0);
    fixedParam.put(1, 0);
    fixedParam.put(2, 0);

    transform->SetParameters(param);
    transform->SetFixedParameters(fixedParam); // Center of the Transform

    std::cout << "Transform matrix:"
              << "\n"
              << transform->GetMatrix() << std::endl;

    using ResampleFilterType =
        itk::ResampleImageFilter<UShortImageType, UShortImageType>;
    ResampleFilterType::Pointer resampler = ResampleFilterType::New();

    resampler->SetInput(flipFilter->GetOutput());
    resampler->SetSize(size_trans);
    resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
    resampler->SetOutputSpacing(spacing_trans); // 1 1 1
    resampler->SetOutputDirection(
        flipFilter->GetOutput()->GetDirection()); // image normal?
    resampler->SetTransform(transform);

    using CastFilterType =
        itk::CastImageFilter<UShortImageType,
                             FloatImageType>; // Maybe not inplace filter
    CastFilterType::Pointer castFilter = CastFilterType::New();
    castFilter->SetInput(resampler->GetOutput());

    /*
    //Default value
    double calibF_A = 1.0;
    double calibF_B = 0.0;

    std::cout << "Temporary forcing CT# applied for tissue" << std::endl;

    std::cout << "CBCT calibration Factor(Recommended: 1, 0): A = " << calibF_A
    << "  B= " << calibF_B << std::endl; typedef
    itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>
    MultiplyImageFilterType; MultiplyImageFilterType::Pointer
    multiplyImageFilter = MultiplyImageFilterType::New();
    multiplyImageFilter->SetInput(castFilter->GetOutput());
    multiplyImageFilter->SetConstant(calibF_A / 65535.0);

    typedef itk::AddImageFilter <FloatImageType, FloatImageType, FloatImageType>
    AddImageFilterType; AddImageFilterType::Pointer addImageFilter =
    AddImageFilterType::New();
    addImageFilter->SetInput1(multiplyImageFilter->GetOutput());
    double addingVal = calibF_B / 65535.0;
    addImageFilter->SetConstant2(addingVal);
    addImageFilter->Update(); //will generate map of real_mu (att.coeff)
    */
    FloatImageType::Pointer spCTImg_mu;
    spCTImg_mu = castFilter->GetOutput(); // addImageFilter->GetOutput();

    // 2) Prepare empty projection images //Should be corresonponding to raw
    // projection images

    // Create a stack of empty projection images
    using ConstantImageSourceType =
        rtk::ConstantImageSource<FloatImageType>; // Output: FLoat image = may
                                                  // be mu_t = log(I_0/I)
    ConstantImageSourceType::Pointer constantImageSource =
        ConstantImageSourceType::New();

    ConstantImageSourceType::SizeType size{};
    ConstantImageSourceType::SpacingType spacing;
    ConstantImageSourceType::PointType origin;

    std::cout << "Setting-up vacant projection image data" << std::endl;

    // a) size
    // std::cout << "chk1" << std::endl;
    size[0] = m_spProjImg3DFloat->GetBufferedRegion()
                  .GetSize()[0]; // qRound((double)DEFAULT_W*m_fResampleF);
    size[1] = m_spProjImg3DFloat->GetBufferedRegion()
                  .GetSize()[1]; // qRound((double)DEFAULT_H*m_fResampleF);
    size[2] = spGeometry->GetGantryAngles().size();
    // iNumOfProjections = size[2];

    // b) spacing
    spacing[0] = m_fProjSpacingX / m_fResampleF; // typical HIS file
    spacing[1] = m_fProjSpacingY / m_fResampleF;
    spacing[2] = 1.0;

    // c) Origin: can center be the image center? or should be related to the CT
    // image???
    origin[0] = spacing[0] * (size[0] - 1) * -0.5;
    origin[1] = spacing[1] * (size[1] - 1) * -0.5;
    origin[2] = 0.0;

    constantImageSource->SetOrigin(origin);
    constantImageSource->SetSpacing(spacing);

    FloatImageType::DirectionType imageDirection;
    imageDirection.SetIdentity(); // no effect
    constantImageSource->SetDirection(imageDirection);
    constantImageSource->SetSize(size);
    constantImageSource->SetConstant(1.0);
    constantImageSource->UpdateOutputInformation();
    std::cout << "Canvas for projection image is ready to write" << std::endl;

    // 4) Prepare CT image to be projected
    // int fwdMethod = en_Joseph; //later, it will be coming from the GUI
    // std::cout << "projection algorithm (0:Joseph, 1: CUDA, 2:RayCast ): " <<
    // fwdMethod << std::endl;

    // Create forward projection image filter
    rtk::ForwardProjectionImageFilter<FloatImageType, FloatImageType>::Pointer
        forwardProjection; // Float to Float

    forwardProjection =
        rtk::JosephForwardProjectionImageFilter<FloatImageType,
                                                FloatImageType>::New();
    // forwardProjection =
    // rtk::RayCastInterpolatorForwardProjectionImageFilter<FloatImageType,
    // FloatImageType>::New();

    itk::TimeProbe projProbe;
    std::cout << "Forward projection is now ongoing" << std::endl;
    forwardProjection->SetInput(
        constantImageSource
            ->GetOutput()); // Canvas. projection image will be saved here.
    forwardProjection->SetInput(1, spCTImg_mu); // reference plan CT image
    forwardProjection->SetGeometry(spGeometry);
    projProbe.Start();
    TRY_AND_EXIT_ON_ITK_EXCEPTION(forwardProjection->Update());
    projProbe.Stop();
    spResultProjImageFloat = forwardProjection->GetOutput();
    std::cout << "Forward projection done in:	" << projProbe.GetMean() << ' '
              << projProbe.GetUnit() << '.' << std::endl;
  } // release all the memory

  // From Float to USHORT and line integral to intensity

  spProjCT3D = UShortImageType::New(); // later
  UShortImageType::SizeType projCT_size =
      spResultProjImageFloat->GetLargestPossibleRegion()
          .GetSize(); // 1024 1024 350
  UShortImageType::IndexType projCT_idxStart =
      spResultProjImageFloat->GetLargestPossibleRegion().GetIndex(); // 0 0 0
  UShortImageType::SpacingType projCT_spacing =
      spResultProjImageFloat->GetSpacing(); // 0.4 0.4 1.0
  UShortImageType::PointType projCT_origin =
      spResultProjImageFloat->GetOrigin(); //-204.6 -204.6 -174.5

  // Copy informations from spResultProjImageFloat
  FloatImageType::RegionType projCT_region;
  projCT_region.SetSize(projCT_size);
  projCT_region.SetIndex(projCT_idxStart);

  spProjCT3D->SetRegions(projCT_region);
  spProjCT3D->SetSpacing(projCT_spacing);
  spProjCT3D->SetOrigin(projCT_origin);

  spProjCT3D->Allocate();
  spProjCT3D->FillBuffer(0);

  // Calculation process
  itk::ImageRegionConstIterator<FloatImageType> itSrc(
      spResultProjImageFloat, spResultProjImageFloat->GetRequestedRegion());
  itk::ImageRegionIterator<UShortImageType> itTarg(
      spProjCT3D, spProjCT3D->GetRequestedRegion()); // writing

  itSrc.GoToBegin();
  itTarg.GoToBegin();

  float fProjVal = 0.0;
  double tmpConvVal = 0.0;

  // Convert line integral to intensity value (I0/I = exp(mu_t)) --> I =
  // I0/exp(mu_t)
  while (!itSrc.IsAtEnd() && !itTarg.IsAtEnd()) {
    fProjVal = itSrc.Get();                 // mu_t //63.5 --> 6.35
    tmpConvVal = (65535.0 / exp(fProjVal)); // physically true

    if (tmpConvVal <= 0.0) {
      itTarg.Set(0);
    } else if (tmpConvVal > 65535.0) {
      itTarg.Set(65535);
    } else {
      itTarg.Set(static_cast<unsigned short>(tmpConvVal));
    }

    ++itSrc;
    ++itTarg;
  }

  // spProjCT3D: USHORT IMAGE of intensity. Not inverted (physical intensity)

  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout << " Patient DIR Path: "
              << m_strPathPatientDir.toLocal8Bit().constData() << std::endl;

    if (m_strPathPatientDir.isEmpty()) {
      std::cout << "File save error!: No patient DIR name" << std::endl;
      return;
    }

    // Get current folder
    QString strCrntDir =
        m_strPathPatientDir + "/" + "IMAGES"; // current Proj folder

    // Make a sub directory
    QDir crntDir(strCrntDir);

    if (!crntDir.exists()) {
      std::cout << "File save error: The specified folder does not exist."
                << std::endl;
      return;
    }

    QString fwdDirName = "fwd_" + m_strDCMUID;

    bool tmpResult = crntDir.mkdir(fwdDirName); // what if the directory exists?

    if (!tmpResult) {
      std::cout << "FwdProj directory seems to exist already. Files will be "
                   "overwritten."
                << std::endl;
    }

    QString strSavingFolder = strCrntDir + "/" + fwdDirName;
    SaveProjImageAsHIS(spProjCT3D, m_arrYKBufProj, strSavingFolder,
                       m_fResampleF);
  }
}

void CbctRecon::SaveProjImageAsHIS(UShortImageType::Pointer &spProj3D,
                                   std::vector<YK16GrayImage> arrYKImage,
                                   QString &strSavingFolder, double resampleF) {
  std::cout << "Starting Saving files" << std::endl;
  FILE *fd = nullptr;

  UShortImageType::Pointer targetImg3D;
  double restoreResampleF = 1.0 / resampleF;

  if (resampleF != 1.0) {
    std::cout << "restore the  resampled image by applying a factor of "
              << restoreResampleF << std::endl;
    ResampleItkImage(spProj3D, targetImg3D, restoreResampleF);
  } else {
    targetImg3D = spProj3D;
  }

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it_FwdProj(
      targetImg3D, targetImg3D->GetRequestedRegion());

  it_FwdProj.SetFirstDirection(0);
  it_FwdProj.SetSecondDirection(1);
  it_FwdProj.GoToBegin();

  for (auto it = arrYKImage.begin();
       it != arrYKImage.end() && !it_FwdProj.IsAtEnd(); ++it) {
    QFileInfo crntFileInfo(it->m_strFilePath);

    QString crntFileName = crntFileInfo.fileName();
    QString crntPath = strSavingFolder + "/" + crntFileName;

#ifdef WIN32
    if (fopen_s(&fd, crntPath.toLocal8Bit().constData(), "wb") == 0) {
      std::cerr << "Could not open file: " << crntPath.toLocal8Bit().constData()
                << " for writing!" << std::endl;
      return;
    }
#else
    fd = fopen(crntPath.toLocal8Bit().constData(), "wb");
#endif
    fwrite(it->m_pElektaHisHeader, 100, 1,
           fd); // this buffer only include header info

    // int imgSize = m_arrYKBufProj[i].m_iWidth * m_arrYKImage[i].m_iHeight;

    // Search matching slice using slice iterator for m_spProjCTImg
    while (!it_FwdProj.IsAtEndOfSlice()) {
      while (!it_FwdProj.IsAtEndOfLine()) {
        auto tmpVal = static_cast<unsigned short>(it_FwdProj.Get());
        tmpVal = 65535 - tmpVal; // inverse is done here

        fwrite(&tmpVal, 2, 1, fd);
        ++it_FwdProj;
      }
      it_FwdProj.NextLine();
    }
    fclose(fd);

    it_FwdProj.NextSlice();
    // std::cout << "Now saving " << i+1 << " th file: " <<
    // crntFileName.toLocal8Bit().constData() << std::endl;
  }

  std::cout << "Saving completed" << std::endl;
}

void CbctRecon::SLT_DoScatterCorrection_APRIORI() {

  bool bExportProj_Fwd = ui.checkBox_ExportFwd->isChecked();
  bool bExportProj_Scat = ui.checkBox_ExportScat->isChecked();
  bool bExportProj_Cor = ui.checkBox_ExportCor->isChecked();

  // ForwardProjection(m_spRefCTImg, m_spCustomGeometry, m_spProjImgCT3D,
  // false); //final moving image
  if (m_pDlgRegistration->m_spMoving != nullptr) {
    ForwardProjection(m_pDlgRegistration->m_spMoving, m_spCustomGeometry,
                      m_spProjImgCT3D, bExportProj_Fwd); // final moving image
  } else if (m_spRefCTImg != nullptr) {
    std::cout << "No Moving image in Registration is found. Ref CT image will "
                 "be used instead"
              << std::endl;
    ForwardProjection(m_spRefCTImg, m_spCustomGeometry, m_spProjImgCT3D,
                      bExportProj_Fwd); // final moving image
  } else {
    std::cout << "Error!: No ref image for forward projection is found."
              << std::endl;
    return;
  }

  // YKTEMP
  std::cout << "ProjImgCT Size = "
            << m_spProjImgCT3D->GetBufferedRegion().GetSize()[0] << ", "
            << m_spProjImgCT3D->GetBufferedRegion().GetSize()[1] << ", "
            << m_spProjImgCT3D->GetBufferedRegion().GetSize()[2] << std::endl;
  std::cout << "ProjImgCT origin = " << m_spProjImgCT3D->GetOrigin()[0] << ", "
            << m_spProjImgCT3D->GetOrigin()[1] << ", "
            << m_spProjImgCT3D->GetOrigin()[2] << std::endl;
  std::cout << "ProjImgCT spacing = " << m_spProjImgCT3D->GetSpacing()[0]
            << ", " << m_spProjImgCT3D->GetSpacing()[1] << ", "
            << m_spProjImgCT3D->GetSpacing()[2] << std::endl;

  // double scaResam = ui.lineEdit_scaResam->text().toDouble();
  double scaMedian = ui.lineEdit_scaMedian->text().toDouble();
  double scaGaussian = ui.lineEdit_scaGaussian->text().toDouble();

  std::cout << "Generating scatter map is ongoing..." << std::endl;

  std::cout << "To account for the mAs values, the intensity scale factor of "
            << GetRawIntensityScaleFactor()
            << "will be multiplied during scatter correction to avoid negative "
               "scatter"
            << std::endl;

  GenScatterMap_PriorCT(m_spProjImgRaw3D, m_spProjImgCT3D, m_spProjImgScat3D,
                        scaMedian, scaGaussian, m_iFixedOffset_ScatterMap,
                        bExportProj_Scat); // void GenScatterMap2D_PriorCT()
  m_spProjImgCT3D->Initialize();           // memory saving

  std::cout << "Scatter correction is in progress..." << std::endl;

  int postScatMedianSize = ui.lineEdit_scaPostMedian->text().toInt();
  ScatterCorr_PrioriCT(m_spProjImgRaw3D, m_spProjImgScat3D, m_spProjImgCorr3D,
                       m_iFixedOffset_ScatterMap, postScatMedianSize,
                       bExportProj_Cor);
  m_spProjImgScat3D->Initialize(); // memory saving

  std::cout << "AfterCorrectionMacro is ongoing..." << std::endl;
  AfterScatCorrectionMacro();
  std::cout << "FINISHED!Scatter correction: CBCT DICOM files are saved"
            << std::endl;
}

void CbctRecon::CalculateIntensityScaleFactorFromMeans(
    UShortImageType::Pointer &spProjRaw3D,
    UShortImageType::Pointer &spProjCT3D) {
  using StatFilterType = itk::StatisticsImageFilter<UShortImageType>;
  double rawMean = 0.0;
  double ctMean = 0.0;
#pragma omp parallel sections
  {
#pragma omp section
      {StatFilterType::Pointer statFilter = StatFilterType::New();
  statFilter->SetInput(spProjRaw3D);
  statFilter->Update();

  rawMean = statFilter->GetMean();
}
#pragma omp section
{
  StatFilterType::Pointer statFilter = StatFilterType::New();
  statFilter->SetInput(spProjCT3D);
  statFilter->Update();

  ctMean = statFilter->GetMean();
}
}
;
const double scaling = ctMean / rawMean;

ui.lineEdit_CurmAs->setText(QString("%1,20").arg((64 * 40 / 20) / scaling));

ui.lineEdit_RefmAs->setText(QString("64,40"));
}

int divisible_by_235_const(const int size) {
  int input_size = size;
  bool ok = true;
  while (ok) {
    if (input_size % 2 == 0) {
      // ok so far
      input_size /= 2; // compiler should optimize so division and modulo is
                       // calculated simultaneously
    } else if (input_size % 3 == 0) {
      // ok so far
      input_size /= 3;
    } else if (input_size % 5 == 0) {
      // ok so far
      input_size /= 5;
    } else {
      ok = false;
    }
  }
  return input_size;
}

int get_padding(int input_size) {
  int cur_padding = 0;
  while (true) {
    // Padding necessary
    if (divisible_by_235_const(input_size + cur_padding) != 1) {
      // While is broken if divisible, else more padding needed.
      cur_padding += 2;
    } else {
      break;
    }
  }
  return cur_padding;
}

#ifdef LOWPASS_FFT
// template<typename T>
FloatImage2DType::Pointer LowPassFFT(FloatImage2DType::Pointer &input,
                                     const double sigmaValue) {
  const unsigned int Dimension = input->GetImageDimension();
  using RealImageType = FloatImage2DType;

  // Some FFT filter implementations, like VNL's, need the image size to be a
  // multiple of small prime numbers.
  using PadFilterType = itk::WrapPadImageFilter<RealImageType, RealImageType>;
  PadFilterType::Pointer padFilter = PadFilterType::New();
  padFilter->SetInput(input);
  PadFilterType::SizeType padding{};
  auto input_size = input->GetLargestPossibleRegion().GetSize();

  for (unsigned int dim = 0; dim < Dimension; dim++) {
    padding[dim] = get_padding(input_size[dim]);
    // Even though the size is usually 512 or 384 which gives padding 0,
    // it would not really speed up the code much to have the checks for these
    // sizes.
  }
  padFilter->SetPadUpperBound(padding);

  using ForwardFFTFilterType = itk::ForwardFFTImageFilter<RealImageType>;
  using ComplexImageType = ForwardFFTFilterType::OutputImageType;
  ForwardFFTFilterType::Pointer forwardFFTFilter = ForwardFFTFilterType::New();
  forwardFFTFilter->SetInput(padFilter->GetOutput());
  forwardFFTFilter->UpdateOutputInformation();

  // A Gaussian is used here to create a low-pass filter.
  using GaussianSourceType = itk::GaussianImageSource<RealImageType>;
  GaussianSourceType::Pointer gaussianSource = GaussianSourceType::New();
  gaussianSource->SetNormalized(false);
  gaussianSource->SetScale(1.0);
  ComplexImageType::ConstPointer transformedInput =
      forwardFFTFilter->GetOutput();
  const ComplexImageType::RegionType inputRegion(
      transformedInput->GetLargestPossibleRegion());
  const ComplexImageType::SizeType inputSize = inputRegion.GetSize();
  const ComplexImageType::SpacingType inputSpacing =
      transformedInput->GetSpacing();
  const ComplexImageType::PointType inputOrigin = transformedInput->GetOrigin();
  const ComplexImageType::DirectionType inputDirection =
      transformedInput->GetDirection();
  gaussianSource->SetSize(inputSize);
  gaussianSource->SetSpacing(inputSpacing);
  gaussianSource->SetOrigin(inputOrigin);
  gaussianSource->SetDirection(inputDirection);
  GaussianSourceType::ArrayType sigma;
  GaussianSourceType::PointType mean;
  sigma.Fill(sigmaValue);
  for (unsigned int ii = 0; ii < Dimension; ++ii) {
    const double halfLength = inputSize[ii] * inputSpacing[ii] / 2.0;
    sigma[ii] *= halfLength;
    mean[ii] = inputOrigin[ii] + halfLength;
  }
  mean = inputDirection * mean;
  gaussianSource->SetSigma(sigma);
  gaussianSource->SetMean(mean);

  using FFTShiftFilterType =
      itk::FFTShiftImageFilter<RealImageType, RealImageType>;
  FFTShiftFilterType::Pointer fftShiftFilter = FFTShiftFilterType::New();
  fftShiftFilter->SetInput(gaussianSource->GetOutput());

  using MultiplyFilterType =
      itk::MultiplyImageFilter<ComplexImageType, RealImageType,
                               ComplexImageType>;
  MultiplyFilterType::Pointer multiplyFilter = MultiplyFilterType::New();
  multiplyFilter->SetInput1(forwardFFTFilter->GetOutput());
  multiplyFilter->SetInput2(fftShiftFilter->GetOutput());

  using InverseFilterType =
      itk::InverseFFTImageFilter<ComplexImageType, RealImageType>;
  InverseFilterType::Pointer inverseFFTFilter = InverseFilterType::New();
  inverseFFTFilter->SetInput(multiplyFilter->GetOutput());
  inverseFFTFilter->Update();
  return inverseFFTFilter->GetOutput();
}
#endif

// spProjRaw3D: raw intensity value (0-65535), spProjCT3D: raw intensity value
// (0-65535)
void CbctRecon::GenScatterMap_PriorCT(UShortImageType::Pointer &spProjRaw3D,
                                      UShortImageType::Pointer &spProjCT3D,
                                      UShortImageType::Pointer &spProjScat3D,
                                      double /*medianRadius*/,
                                      double gaussianSigma,
                                      int nonNegativeScatOffset, bool bSave) {
  // Scatter map: should be 2D to use 2D median, Gaussian filters
  if (m_iCntSelectedProj < 1) {
    std::cout << "error: no count of proj image" << std::endl;
    return;
  }

  if ((spProjRaw3D == nullptr) || (spProjCT3D == nullptr)) {
    std::cout << "error: proj image 3D is not ready" << std::endl;
    return;
  }

  using SizeType = UShortImageType::SizeType;
  SizeType size1 = spProjRaw3D->GetRequestedRegion().GetSize();
  SizeType size2 = spProjCT3D->GetRequestedRegion().GetSize();

  std::cout << "Raw3DProj Size= " << size1 << std::endl;
  std::cout << "spProjCT Size= " << size2 << std::endl;

  bool bHighResolMacro = false; // raw imag= 1024, scattermap = 512
  if (size1[0] != size2[0] || size1[1] != size2[1] || size1[2] != size2[2]) {
    std::cout << "Raw and CT projection dimension are not matching. under the "
                 "high resolution macro?"
              << std::endl;
    // Why 2.0 shouldn't it be 1/downresfactor ?
    if (size1[0] ==
            static_cast<SizeType::SizeValueType>(qRound(size2[0] * 2.0)) &&
        size1[1] ==
            static_cast<SizeType::SizeValueType>(qRound(size2[1] * 2.0)) &&
        size1[2] == size2[2]) {
      bHighResolMacro = true;
    } else {
      return;
    }
  }

  UShortImageType::Pointer spTmpProjRaw3D;

  if (bHighResolMacro) {
    std::cout << "bHighResolMacro is unexpectedly on" << std::endl;
    ResampleItkImage(spProjRaw3D, spTmpProjRaw3D, 0.5);
  } else {
    spTmpProjRaw3D = spProjRaw3D;
  }

  AllocateByRef(spTmpProjRaw3D, spProjScat3D);
  // AllocateByRef(spProjCT3D, spProjScat3D);

  // std::cout << "Scat3D size = " <<
  // spProjScat3D->GetRequestedRegion().GetSize() << std::endl;

  UShortImageType::SizeType imgSize =
      spTmpProjRaw3D->GetRequestedRegion().GetSize();
  // UShortImageType::SizeType imgSize =
  // spProjCT3D->GetRequestedRegion().GetSize();

  // UShortImageType::SizeType imgSize =
  // spSrcImg3D->GetBufferedRegion().GetSize();  Create spProjScat3D with same
  // dimension of the spProjRaw3D
  int iSizeZ = imgSize[2];

  // std::cout << "resample factor " << resF2D << std::endl;
  CalculateIntensityScaleFactorFromMeans(spProjRaw3D, spProjCT3D);
  double mAs_correctionFactor = GetRawIntensityScaleFactor();
  for (int i = 0; i < iSizeZ; i++) {
    FloatImage2DType::Pointer spImg2DRaw;
    FloatImage2DType::Pointer spImg2DPrim;
    FloatImage2DType::Pointer spImg2DScat;

    Get2DFrom3D(spTmpProjRaw3D, spImg2DRaw, i,
                PLANE_AXIAL); // simple conversion between ushort 3D to float 2D
                              // (using casting, not log): input/output: 0-65535
    Get2DFrom3D(spProjCT3D, spImg2DPrim, i, PLANE_AXIAL);

    // Dimension should be matched
    AllocateByRef(spImg2DRaw, spImg2DScat);

    itk::ImageRegionConstIteratorWithIndex<FloatImage2DType> it_Src1(
        spImg2DRaw, spImg2DRaw->GetRequestedRegion());
    itk::ImageRegionConstIteratorWithIndex<FloatImage2DType> it_Src2(
        spImg2DPrim, spImg2DPrim->GetRequestedRegion());
    itk::ImageRegionIteratorWithIndex<FloatImage2DType> it_Tar(
        spImg2DScat, spImg2DScat->GetRequestedRegion());

    // int cnt1 = 0; int cnt2 = 0;
    for (it_Src1.GoToBegin(), it_Src2.GoToBegin(), it_Tar.GoToBegin();
         !it_Src1.IsAtEnd() && !it_Src2.IsAtEnd() && !it_Tar.IsAtEnd();
         ++it_Src1, ++it_Src2, ++it_Tar) {
      float intensityValScat =
          it_Src1.Get() * mAs_correctionFactor -
          it_Src2.Get(); // raw intensity * mAs_CF - primary intensity
      it_Tar.Set(intensityValScat); // float 	  //allow minus value
    }

#ifdef LOWPASS_FFT
    spImg2DScat = LowPassFFT(spImg2DScat, gaussianSigma);
#else
    // ResampleItkImage2D(spImg2DScat, spImg2DScat, resF2D);
    using MedianFilterType =
        itk::MedianImageFilter<FloatImage2DType, FloatImage2DType>;

    MedianFilterType::Pointer medianFilterX = MedianFilterType::New();
    MedianFilterType::InputSizeType radiusX{};
    radiusX[0] = medianRadius;
    radiusX[1] = 0;
    medianFilterX->SetRadius(radiusX);
    medianFilterX->SetInput(spImg2DScat);
    // medianFilterX->Update();
    // spImg2DScat = medianFilterX->GetOutput();

    MedianFilterType::Pointer medianFilterY = MedianFilterType::New();
    MedianFilterType::InputSizeType radiusY{};
    radiusY[0] = 0;
    radiusY[1] = medianRadius;
    medianFilterY->SetRadius(radiusY);
    medianFilterY->SetInput(medianFilterX->GetOutput());
    medianFilterY->Update();
    spImg2DScat = medianFilterY->GetOutput();

    using SmoothingFilterType =
        itk::SmoothingRecursiveGaussianImageFilter<FloatImage2DType,
                                                   FloatImage2DType>;
    SmoothingFilterType::Pointer gaussianFilter = SmoothingFilterType::New();
    // gaussianFilter->SetInput(medianFilter->GetOutput());
    gaussianFilter->SetInput(spImg2DScat);
    if (hisIsUsed) {
      gaussianFilter->SetSigma(
          gaussianSigma); // filter specific setting for 512x 512 image
    } else {
      SmoothingFilterType::SigmaArrayType gaussianSigmaArray;
      gaussianSigmaArray[0] = gaussianSigma;
      gaussianSigmaArray[1] = .75 * gaussianSigma;
      gaussianFilter->SetSigmaArray(
          gaussianSigmaArray); // filter specific setting for 512x384 (varian/2)
    }
      // gaussianFilter->Update();
      // spImg2DScat = gaussianFilter->GetOutput();
#endif

    using AddImageFilterType =
        itk::AddImageFilter<FloatImage2DType, FloatImage2DType,
                            FloatImage2DType>;
    AddImageFilterType::Pointer addFilter = AddImageFilterType::New();
#ifdef LOWPASS_FFT
    addFilter->SetInput1(spImg2DScat);
#else
    addFilter->SetInput1(gaussianFilter->GetOutput());
#endif
    addFilter->SetConstant2(static_cast<float>(nonNegativeScatOffset));
    addFilter->Update();
    spImg2DScat = addFilter->GetOutput(); // even after the offset applied, -
                                          // value is still possible

    // float to unsigned short
    Set2DTo3D(spImg2DScat, spProjScat3D, i,
              PLANE_AXIAL); // input/Output: 0-65535 intensity valuesno mu_t to
                            // intensity converion is involved

    int unit = qRound(iSizeZ / 10.0);
    if (i % unit == 0) {
      std::cout << "Generating scatter map: "
                << (i / static_cast<double>(unit)) * 10.0 << " % is done"
                << std::endl;
    }
  } // end of for

  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout << "Patient DIR Path: "
              << m_strPathPatientDir.toLocal8Bit().constData() << std::endl;

    if (m_strPathPatientDir.isEmpty() && hisIsUsed) {
      std::cout << "File save error!: No patient DIR name" << std::endl;
      return;
    }

    // Get current folder
    QString strCrntDir;
    if (hisIsUsed) {
      strCrntDir = m_strPathPatientDir + "/" + "IMAGES"; // current Proj folder
    } else {
      strCrntDir = m_pDlgRegistration->m_strPathPlastimatch;
    }

    // Make a sub directory
    QDir crntDir(strCrntDir);

    if (!crntDir.exists()) {
      std::cout << "File save error: The specified folder does not exist."
                << std::endl;
      return;
    }

    QString scatDirName = "sca_" + m_strDCMUID;

    bool tmpResult = crntDir.mkdir(scatDirName); // what if the directory
                                                 // exists?

    if (!tmpResult) {
      std::cout << "Scatter map directory seems to exist already. Files will "
                   "be overwritten."
                << std::endl;
    }

    QString strSavingFolder = strCrntDir + "/" + scatDirName;
    if (hisIsUsed) {
      SaveProjImageAsHIS(spProjScat3D, m_arrYKBufProj, strSavingFolder,
                         m_fResampleF);
    } else {
      using imagewritertype = itk::ImageFileWriter<UShortImageType>;
      imagewritertype::Pointer imagewriter = imagewritertype::New();
      imagewriter->SetInput(spProjScat3D);
      imagewriter->SetFileName(
          QString(strSavingFolder + "/scatter.mha").toStdString());
      imagewriter->Update();
    }
  }
}

void CbctRecon::ScatterCorr_PrioriCT(UShortImageType::Pointer &spProjRaw3D,
                                     UShortImageType::Pointer &spProjScat3D,
                                     UShortImageType::Pointer &m_spProjCorr3D,
                                     int nonNegativeScatOffset, int postMedian,
                                     bool bSave) {
  // Scatter map: should be 2D to use 2D median, Gaussian filters
  if (m_iCntSelectedProj < 1) {
    std::cout << "error: no count of proj image" << std::endl;
    return;
  }

  if ((spProjRaw3D == nullptr) || (spProjScat3D == nullptr)) {
    std::cout << "Error: proj image 3D is not ready" << std::endl;
    return;
  }

  UShortImageType::SizeType size1 = spProjRaw3D->GetRequestedRegion().GetSize();
  UShortImageType::SizeType size2 =
      spProjScat3D->GetRequestedRegion().GetSize();

  std::cout << "Raw3DProj Size= " << size1 << std::endl;
  std::cout << "spProjScat3D Size= " << size2 << std::endl;

  bool bHighResolMacro = false;

  if (size1[0] != size2[0] || size1[1] != size2[1] || size1[2] != size2[2]) {
    std::cout << "Raw and scatter projection dimension are not matching. under "
                 "the high resolution macro?"
              << std::endl;

    if (static_cast<int>(size1[0]) == qRound(size2[0] * 2.0) &&
        static_cast<int>(size1[1]) == qRound(size2[1] * 2.0) &&
        size1[2] == size2[2]) {
      bHighResolMacro = true;
    } else {
      return;
    }
  }

  UShortImageType::Pointer spTmpProjScat3D;

  if (bHighResolMacro) {
    ResampleItkImage(spProjScat3D, spTmpProjScat3D, 2.0);
  } else {
    spTmpProjScat3D = spProjScat3D;
  }

  AllocateByRef(spProjRaw3D, m_spProjCorr3D);

  UShortImageType::SizeType imgSize =
      spProjRaw3D->GetRequestedRegion().GetSize();

  // UShortImageType::SizeType imgSize =
  // spSrcImg3D->GetBufferedRegion().GetSize();  Create spProjScat3D with same
  // dimension of the spProjRaw3D
  int iSizeZ = imgSize[2];

  // std::cout << "resample factor " << resF2D << std::endl;.

  double mAs_correctionFactor = GetRawIntensityScaleFactor();
  for (int i = 0; i < iSizeZ; i++) {
    FloatImage2DType::Pointer spImg2DRaw;
    FloatImage2DType::Pointer spImg2DScat;
    FloatImage2DType::Pointer spImg2DCorr;

    Get2DFrom3D(spProjRaw3D, spImg2DRaw, i, PLANE_AXIAL);
    Get2DFrom3D(spTmpProjScat3D, spImg2DScat, i, PLANE_AXIAL);

    // Dimension should be matched
    AllocateByRef(spImg2DRaw, spImg2DCorr);

    itk::ImageRegionConstIteratorWithIndex<FloatImage2DType> it_Src1(
        spImg2DRaw, spImg2DRaw->GetRequestedRegion());
    itk::ImageRegionConstIteratorWithIndex<FloatImage2DType> it_Src2(
        spImg2DScat, spImg2DScat->GetRequestedRegion());
    itk::ImageRegionIteratorWithIndex<FloatImage2DType> it_Tar(
        spImg2DCorr, spImg2DCorr->GetRequestedRegion());

    for (it_Src1.GoToBegin(), it_Src2.GoToBegin(), it_Tar.GoToBegin();
         !it_Src1.IsAtEnd() && !it_Src2.IsAtEnd() && !it_Tar.IsAtEnd();
         ++it_Src1, ++it_Src2, ++it_Tar) {
      float rawVal = it_Src1.Get() * mAs_correctionFactor;
      float scatVal = it_Src2.Get() - nonNegativeScatOffset;
      float corrVal = rawVal - scatVal;

      if (corrVal < 1.0f) {
        corrVal = 1.0f; // Overflow control
      }
      if (corrVal >
          65534.0f) { // 65535 -->(inversion) --> 0 --> LOg (65536 / 0) = ERROR!
        corrVal = 65534.0f;
      }

      it_Tar.Set(corrVal); // float // later, add customSPR
      // corrVal = (float)(rawVal - customSPR*scatterVal);
    }
    // Post Median filtering

    if (bHighResolMacro) {
      postMedian = postMedian * 2.0;
    }

    if (postMedian >= 2) // YK2015
    {
      using MedianFilterType =
          itk::MedianImageFilter<FloatImage2DType, FloatImage2DType>;
      MedianFilterType::Pointer medianFilter = MedianFilterType::New();
      MedianFilterType::InputSizeType radius{};

      radius[0] = qRound(postMedian / 2.0);
      radius[1] = qRound(postMedian / 2.0);

      /*	if (ui.radioButton_UseCUDA->isChecked())
              {
              int wndX = radius[0] * 2 + 1;
              int wndY = radius[1] * 2 + 1;

              cudaMedianFilter2DITK(spImg2DCorr, wndX, wndY);
              }
              else
              {*/
      medianFilter->SetRadius(radius);
      medianFilter->SetInput(spImg2DCorr);
      medianFilter->Update();
      spImg2DCorr = medianFilter->GetOutput();
      //}
    }

    Set2DTo3D(spImg2DCorr, m_spProjCorr3D, i, PLANE_AXIAL); // float2D to USHORT

    int unit = qRound(iSizeZ / 10.0);
    if (i % unit == 0) {
      std::cout << "Applying scatter correction: "
                << (i / static_cast<double>(unit)) * 10.0 << " % is done"
                << std::endl;
    }
  }

  if (bSave) {
    // Saving part: save as his file in sub-folder of raw image
    std::cout << "Files are being saved" << std::endl;
    std::cout << "Patient DIR Path: "
              << m_strPathPatientDir.toLocal8Bit().constData() << std::endl;

    if (m_strPathPatientDir.isEmpty()) {
      std::cout << "File save error!: No patient DIR name" << std::endl;
      return;
    }

    // Get current folder
    QString strCrntDir =
        m_strPathPatientDir + "/" + "IMAGES"; // current Proj folder

    // Make a sub directory
    QDir crntDir(strCrntDir);

    if (!crntDir.exists()) {
      std::cout << "File save error: The specified folder does not exist."
                << std::endl;
      return;
    }

    QString scatDirName = "cor_" + m_strDCMUID;

    bool tmpResult = crntDir.mkdir(scatDirName); // what if the directory
                                                 // exists?

    if (!tmpResult) {
      std::cout << "Corrected projection directory seems to exist already. "
                   "Files will be overwritten."
                << std::endl;
    }

    QString strSavingFolder = strCrntDir + "/" + scatDirName;

    if (!bHighResolMacro) {
      SaveProjImageAsHIS(m_spProjCorr3D, m_arrYKBufProj, strSavingFolder,
                         m_fResampleF);
    } else {
      SaveProjImageAsHIS(m_spProjCorr3D, m_arrYKBufProj, strSavingFolder, 1.0);
    }
  }
  // spProjScat3D->Initialize(); //memory release
}

// spSrcImg3D: usually projImage in USHORT type
void CbctRecon::Get2DFrom3D(UShortImageType::Pointer &spSrcImg3D,
                            FloatImage2DType::Pointer &spTargetImg2D, int idx,
                            enPLANE iDirection) {
  if (spSrcImg3D == nullptr) {
    return;
  }

  int idxHor, idxVer, idxZ;

  switch (iDirection) {
  case PLANE_AXIAL:
    idxHor = 0;
    idxVer = 1;
    idxZ = 2;
    break;
  case PLANE_FRONTAL:
    idxHor = 0;
    idxVer = 2;
    idxZ = 1;
    break;
  case PLANE_SAGITTAL:
    idxHor = 1;
    idxVer = 2;
    idxZ = 0;
    break;
  }

  // Create 2D target image based on geometry of 3D
  UShortImageType::SizeType imgDim = spSrcImg3D->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing = spSrcImg3D->GetSpacing();
  // UShortImageType::PointType origin = spSrcImg3D->GetOrigin();

  // int width = imgDim[idxHor];
  // int height = imgDim[idxVer];
  int zSize = imgDim[idxZ];
  // std::cout << "Get2DFrom3D zSize = " << zSize << std::endl;

  if (idx < 0 || idx >= zSize) {
    std::cout << "Error! idx is out of the range" << std::endl;
    return;
  }

  FloatImage2DType::IndexType idxStart{};
  idxStart[0] = 0;
  idxStart[1] = 0;

  FloatImage2DType::SizeType size2D{};
  size2D[0] = imgDim[idxHor];
  size2D[1] = imgDim[idxVer];

  FloatImage2DType::SpacingType spacing2D;
  spacing2D[0] = spacing[idxHor];
  spacing2D[1] = spacing[idxVer];

  FloatImage2DType::PointType origin2D;
  //  origin2D[0] = origin[idxHor];
  //  origin2D[1] = origin[idxVer];
  origin2D[0] = size2D[0] * spacing2D[0] / -2.0;
  origin2D[1] = size2D[1] * spacing2D[1] / -2.0;

  FloatImage2DType::RegionType region;
  region.SetSize(size2D);
  region.SetIndex(idxStart);

  // spTargetImg2D is supposed to be empty.
  if (spTargetImg2D != nullptr) {
    std::cout
        << "something is here in target image. is it gonna be overwritten?"
        << std::endl;
  }

  spTargetImg2D = FloatImage2DType::New();
  spTargetImg2D->SetRegions(region);
  spTargetImg2D->SetSpacing(spacing2D);
  spTargetImg2D->SetOrigin(origin2D);

  spTargetImg2D->Allocate();
  spTargetImg2D->FillBuffer(0);

  // std::cout << "src size = " << spSrcImg3D->GetRequestedRegion().GetSize() <<
  // " " << std::endl;  std::cout << "target image size = " <<
  // spTargetImg2D->GetRequestedRegion().GetSize() << " " << std::endl;

  itk::ImageSliceConstIteratorWithIndex<UShortImageType> it_3D(
      spSrcImg3D, spSrcImg3D->GetRequestedRegion());
  // itk::ImageRegionIteratorWithIndex<FloatImageType2D> it_2D (spTargetImg2D,
  // spTargetImg2D->GetRequestedRegion());
  itk::ImageRegionIterator<FloatImage2DType> it_2D(
      spTargetImg2D, spTargetImg2D->GetRequestedRegion());

  it_3D.SetFirstDirection(idxHor);
  it_3D.SetSecondDirection(idxVer);

  it_3D.GoToBegin();
  it_2D.GoToBegin();

  for (int i = 0; i < zSize && !it_3D.IsAtEnd(); i++) {
    /*QFileInfo crntFileInfo(arrYKImage[i].m_strFilePath);
    QString crntFileName = crntFileInfo.fileName();
    QString crntPath = strSavingFolder + "/" + crntFileName;*/
    // Search matching slice using slice iterator for m_spProjCTImg
    // std::cout << "Get2DFrom3D: Slide= " << i  << " ";

    if (i == idx) {
      while (!it_3D.IsAtEndOfSlice()) // Error here why?
      {
        while (!it_3D.IsAtEndOfLine()) {
          auto tmpVal = static_cast<float>(
              it_3D.Get()); // in proj image case, this is intensity
          it_2D.Set(tmpVal);
          ++it_2D;
          ++it_3D;
        } // while2
        it_3D.NextLine();
      } // while1
      break;
    } // end if
    it_3D.NextSlice();
  } // end of for

  // std::cout << "cnt = " << cnt << " TotCnt " << cntTot << std::endl;
  /*YK16GrayImage tmpYK;
  tmpYK.UpdateFromItkImageFloat(spTargetImg2D);
  QString str = QString("D:\\testYK\\InsideFunc_%1.raw").arg(idx);
  tmpYK.SaveDataAsRaw(str.toLocal8Bit().constData());*/
}

void CbctRecon::Set2DTo3D(FloatImage2DType::Pointer &spSrcImg2D,
                          UShortImageType::Pointer &spTargetImg3D, int idx,
                          enPLANE iDirection) {
  if ((spSrcImg2D == nullptr) ||
      (spTargetImg3D == nullptr)) { // Target image should be also ready.
    return;
  }

  int idxHor, idxVer, idxZ;

  switch (iDirection) {
  case PLANE_AXIAL:
    idxHor = 0;
    idxVer = 1;
    idxZ = 2;
    break;
  case PLANE_FRONTAL:
    idxHor = 0;
    idxVer = 2;
    idxZ = 1;
    break;
  case PLANE_SAGITTAL:
    idxHor = 1;
    idxVer = 2;
    idxZ = 0;
    break;
  }

  FloatImage2DType::SizeType imgDim2D =
      spSrcImg2D->GetBufferedRegion().GetSize();
  // FloatImage2DType::SpacingType spacing2D = spSrcImg2D->GetSpacing();
  // FloatImage2DType::PointType origin2D = spSrcImg2D->GetOrigin();

  UShortImageType::SizeType imgDim3D =
      spTargetImg3D->GetBufferedRegion().GetSize();
  // UShortImageType::SpacingType spacing3D = spTargetImg3D->GetSpacing();
  // UShortImageType::PointType origin3D = spTargetImg3D->GetOrigin();

  // Filtering
  if (imgDim2D[0] != imgDim3D[idxHor] || imgDim2D[1] != imgDim3D[idxVer] ||
      idx < 0 || idx >= static_cast<int>(imgDim3D[idxZ])) {
    std::cout << "Error: image dimensions is not matching" << std::endl;
    std::cout << "2D= " << imgDim2D << std::endl;
    std::cout << "3D= " << imgDim3D << std::endl;
    return;
  }
  /*int width = imgDim[idxHor];
  int height  = imgDim[idxVer];*/

  // itk::ImageRegionConstIteratorWithIndex<FloatImageType2D> it_2D (spSrcImg2D,
  // spSrcImg2D->GetRequestedRegion());
  itk::ImageRegionConstIterator<FloatImage2DType> it_2D(
      spSrcImg2D, spSrcImg2D->GetRequestedRegion());
  itk::ImageSliceIteratorWithIndex<UShortImageType> it_3D(
      spTargetImg3D, spTargetImg3D->GetRequestedRegion());

  it_3D.SetFirstDirection(idxHor);
  it_3D.SetSecondDirection(idxVer);
  it_3D.GoToBegin();

  int zSize = imgDim3D[idxZ];

  it_2D.GoToBegin();

  float fVal2D = 0.0;
  unsigned short outputVal = 0;

  for (int i = 0; i < zSize && !it_3D.IsAtEnd(); i++) {
    /*QFileInfo crntFileInfo(arrYKImage[i].m_strFilePath);
    QString crntFileName = crntFileInfo.fileName();
    QString crntPath = strSavingFolder + "/" + crntFileName;*/
    // Search matching slice using slice iterator for m_spProjCTImg
    if (i == idx) {
      while (!it_3D.IsAtEndOfSlice()) {
        while (!it_3D.IsAtEndOfLine()) {
          fVal2D = it_2D.Get();

          if (fVal2D < 0.0) {
            outputVal = 0U;
          } else if (fVal2D > 65535.0) {
            outputVal = 65535U;
          } else {
            outputVal = static_cast<unsigned short>(qRound(fVal2D));
          }

          it_3D.Set(outputVal);
          // float tmpVal = (float)(it_3D.Get()); //in proj image case, this is
          // intensity  it_2D.Set(tmpVal);
          ++it_2D;
          ++it_3D;
        } // while2
        it_3D.NextLine();
      } // while1
      break;
    }
    //
    it_3D.NextSlice();
  } // end of for
}

// void CbctRecon::Get2DFrom3D( FloatImageType::Pointer& spSrcImg3D,
// FloatImageType2D::Pointer& spTargetImg2D, enPLANE iDirection)
//{
//
//}

// From line integral to raw intensity
// bkIntensity is usually 65535
void CbctRecon::ConvertLineInt2Intensity(
    FloatImageType::Pointer &spProjLineInt3D,
    UShortImageType::Pointer &spProjIntensity3D, int bkIntensity) {
  if (spProjLineInt3D == nullptr) {
    return;
  }
  // FloatImageType::IMageRegionIteratorWithIndex

  AllocateByRef(spProjLineInt3D, spProjIntensity3D);

  itk::ImageRegionConstIteratorWithIndex<FloatImageType> it_Src(
      spProjLineInt3D, spProjLineInt3D->GetRequestedRegion());
  itk::ImageRegionIteratorWithIndex<UShortImageType> it_Tar(
      spProjIntensity3D, spProjIntensity3D->GetRequestedRegion());

  for (it_Src.GoToBegin(), it_Tar.GoToBegin();
       !it_Src.IsAtEnd() && !it_Tar.IsAtEnd(); ++it_Src, ++it_Tar) {
    float intensityVal = exp(static_cast<double>(it_Src.Get()) * -1.0) *
                         static_cast<double>(bkIntensity);

    if (intensityVal <= 1.0) {
      intensityVal = 1.0;
    }
    if (intensityVal >= 65534) {
      intensityVal = 65534.0;
    }

    it_Tar.Set(static_cast<unsigned short>(intensityVal));
  }
}

void CbctRecon::ConvertIntensity2LineInt(
    UShortImageType::Pointer &spProjIntensity3D,
    FloatImageType::Pointer &spProjLineInt3D, int bkIntensity) {
  if (spProjIntensity3D == nullptr) {
    return;
  }
  // FloatImageType::IMageRegionIteratorWithIndex

  AllocateByRef(spProjIntensity3D, spProjLineInt3D);

  itk::ImageRegionConstIteratorWithIndex<UShortImageType> it_Src(
      spProjIntensity3D, spProjIntensity3D->GetRequestedRegion());
  itk::ImageRegionIteratorWithIndex<FloatImageType> it_Tar(
      spProjLineInt3D, spProjLineInt3D->GetRequestedRegion());

  const auto background_intensity = static_cast<double>(bkIntensity);

  for (it_Src.GoToBegin(), it_Tar.GoToBegin();
       !it_Src.IsAtEnd() && !it_Tar.IsAtEnd(); ++it_Src, ++it_Tar) {
    // mu = ln(I_0/I) OR mu = ln(I/I0)
    float mu_t_val =
        log(background_intensity / static_cast<double>(it_Src.Get()));
    it_Tar.Set(mu_t_val);
  }
}

void CbctRecon::AllocateByRef(FloatImageType::Pointer &spRefImg3D,
                              FloatImageType::Pointer &spTarImg3D) {
  FloatImageType::SizeType sizeSrc = spRefImg3D->GetBufferedRegion().GetSize();
  FloatImageType::IndexType startSrc =
      spRefImg3D->GetBufferedRegion().GetIndex();

  FloatImageType::SpacingType spacingSrc = spRefImg3D->GetSpacing();
  FloatImageType::PointType originSrc = spRefImg3D->GetOrigin();

  FloatImageType::RegionType region;
  region.SetSize(sizeSrc);
  region.SetIndex(startSrc);

  spTarImg3D = FloatImageType::New();

  spTarImg3D->SetRegions(region);
  spTarImg3D->SetSpacing(spacingSrc);
  spTarImg3D->SetOrigin(originSrc);

  spTarImg3D->Allocate();
  spTarImg3D->FillBuffer(0);
}

void CbctRecon::AllocateByRef(UShortImageType::Pointer &spRefImg3D,
                              UShortImageType::Pointer &spTarImg3D) {
  UShortImageType::SizeType sizeSrc = spRefImg3D->GetBufferedRegion().GetSize();
  UShortImageType::IndexType startSrc =
      spRefImg3D->GetBufferedRegion().GetIndex();

  UShortImageType::SpacingType spacingSrc = spRefImg3D->GetSpacing();
  UShortImageType::PointType originSrc = spRefImg3D->GetOrigin();

  UShortImageType::RegionType region;
  region.SetSize(sizeSrc);
  region.SetIndex(startSrc);

  spTarImg3D = UShortImageType::New();

  spTarImg3D->SetRegions(region);
  spTarImg3D->SetSpacing(spacingSrc);
  spTarImg3D->SetOrigin(originSrc);

  spTarImg3D->Allocate();
  spTarImg3D->FillBuffer(0);
}

void CbctRecon::AllocateByRef(FloatImage2DType::Pointer &spRefImg2D,
                              FloatImage2DType::Pointer &spTarImg2D) {
  FloatImage2DType::SizeType sizeSrc =
      spRefImg2D->GetBufferedRegion().GetSize();
  FloatImage2DType::IndexType startSrc =
      spRefImg2D->GetBufferedRegion().GetIndex();

  FloatImage2DType::SpacingType spacingSrc = spRefImg2D->GetSpacing();
  FloatImage2DType::PointType originSrc = spRefImg2D->GetOrigin();

  FloatImage2DType::RegionType region;
  region.SetSize(sizeSrc);
  region.SetIndex(startSrc);

  spTarImg2D = FloatImage2DType::New();

  spTarImg2D->SetRegions(region);
  spTarImg2D->SetSpacing(spacingSrc);
  spTarImg2D->SetOrigin(originSrc);

  spTarImg2D->Allocate();
  spTarImg2D->FillBuffer(0);
}

void CbctRecon::AllocateByRef(UShortImageType::Pointer &spRefImg3D,
                              FloatImageType::Pointer &spTarImg3D) {
  UShortImageType::SizeType sizeSrc = spRefImg3D->GetBufferedRegion().GetSize();
  UShortImageType::IndexType startSrc =
      spRefImg3D->GetBufferedRegion().GetIndex();
  UShortImageType::SpacingType spacingSrc = spRefImg3D->GetSpacing();
  UShortImageType::PointType originSrc = spRefImg3D->GetOrigin();

  FloatImageType::SizeType size = FloatImageType::SizeType(sizeSrc);
  FloatImageType::IndexType start = FloatImageType::IndexType(startSrc);
  FloatImageType::SpacingType spacing = FloatImageType::SpacingType(spacingSrc);
  FloatImageType::PointType origin = FloatImageType::PointType(originSrc);

  FloatImageType::RegionType region;
  region.SetSize(size);
  region.SetIndex(start);

  spTarImg3D = FloatImageType::New();

  spTarImg3D->SetRegions(region);
  spTarImg3D->SetSpacing(spacing);
  spTarImg3D->SetOrigin(origin);

  spTarImg3D->Allocate();
  spTarImg3D->FillBuffer(0);
}

void CbctRecon::AllocateByRef(FloatImageType::Pointer &spRefImg3D,
                              UShortImageType::Pointer &spTarImg3D) {
  FloatImageType::SizeType sizeSrc = spRefImg3D->GetBufferedRegion().GetSize();
  FloatImageType::IndexType startSrc =
      spRefImg3D->GetBufferedRegion().GetIndex();
  FloatImageType::SpacingType spacingSrc = spRefImg3D->GetSpacing();
  FloatImageType::PointType originSrc = spRefImg3D->GetOrigin();

  UShortImageType::SizeType size = UShortImageType::SizeType(sizeSrc);
  UShortImageType::IndexType start = UShortImageType::IndexType(startSrc);
  UShortImageType::SpacingType spacing =
      UShortImageType::SpacingType(spacingSrc);
  UShortImageType::PointType origin = UShortImageType::PointType(originSrc);

  UShortImageType::RegionType region;

  region.SetSize(size);
  region.SetIndex(start);

  spTarImg3D = UShortImageType::New();

  spTarImg3D->SetRegions(region);
  spTarImg3D->SetSpacing(spacing);
  spTarImg3D->SetOrigin(origin);

  spTarImg3D->Allocate();
  spTarImg3D->FillBuffer(0);
}
// it works! new memory will be allocated for spTarImg
void CbctRecon::ResampleItkImage(FloatImageType::Pointer &spSrcImg,
                                 FloatImageType::Pointer &spTarImg,
                                 double resFactor) {
  if (spSrcImg == nullptr) {
    return;
  }

  //  std::cout << "original Origin: " << spSrcImg2D->GetOrigin() << std::endl;
  using ResampleImageFilterType =
      itk::ResampleImageFilter<FloatImageType, FloatImageType, float>;
  ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();

  resample->SetOutputDirection(spSrcImg->GetDirection());

  using TransformType = itk::AffineTransform<float, 3>;
  TransformType::Pointer transform = TransformType::New();

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<FloatImageType, float>;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  resample->SetInterpolator(interpolator);
  if ((hisIsUsed && DEFAULT_ELEKTA_PROJ_HEIGHT ==
                        spSrcImg->GetBufferedRegion().GetSize()[1]) ||
      (!hisIsUsed && DEFAULT_VARIAN_PROJ_HEIGHT ==
                         spSrcImg->GetBufferedRegion().GetSize()[1])) {
    resample->SetDefaultPixelValue(50);
  } else {
    resample->SetDefaultPixelValue(0);
  }

  FloatImageType::SizeType inputSize =
      spSrcImg->GetLargestPossibleRegion().GetSize();
  FloatImageType::SizeType outputSize{};
  outputSize[0] = qRound(inputSize[0] * resFactor);
  outputSize[1] = qRound(inputSize[1] * resFactor);
  outputSize[2] = inputSize[2];
  resample->SetSize(outputSize);

  FloatImageType::SpacingType outputSpacing;
  outputSpacing[0] =
      spSrcImg->GetSpacing()[0] *
      (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
  outputSpacing[1] =
      spSrcImg->GetSpacing()[1] *
      (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
  outputSpacing[2] = spSrcImg->GetSpacing()[2];
  resample->SetOutputSpacing(outputSpacing);

  FloatImageType::PointType outputOrigin = spSrcImg->GetOrigin(); // Float image
  resample->SetOutputOrigin(outputOrigin);

  resample->SetInput(spSrcImg);
  transform->SetIdentity();
  resample->SetTransform(transform);

  resample->Update();

  // resample->GetOutput()->SetOrigin(prevOrigin);
  spTarImg = resample->GetOutput(); // is it copied? or replaced?
}

void CbctRecon::ResampleItkImage(UShortImageType::Pointer &spSrcImg,
                                 UShortImageType::Pointer &spTarImg,
                                 double resFactor) {
  if (spSrcImg == nullptr) {
    return;
  }

  //  std::cout << "original Origin: " << spSrcImg2D->GetOrigin() << std::endl;
  using ResampleImageFilterType =
      itk::ResampleImageFilter<UShortImageType, UShortImageType, float>;
  ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();

  resample->SetOutputDirection(spSrcImg->GetDirection());

  using TransformType = itk::AffineTransform<float, 3>;
  TransformType::Pointer transform = TransformType::New();

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<UShortImageType, float>;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  resample->SetInterpolator(interpolator);

  resample->SetDefaultPixelValue(50);

  UShortImageType::SizeType inputSize =
      spSrcImg->GetLargestPossibleRegion().GetSize();
  UShortImageType::SizeType outputSize{};
  outputSize[0] = qRound(inputSize[0] * resFactor);
  outputSize[1] = qRound(inputSize[1] * resFactor);
  outputSize[2] = inputSize[2];
  resample->SetSize(outputSize);

  UShortImageType::SpacingType outputSpacing;
  outputSpacing[0] =
      spSrcImg->GetSpacing()[0] *
      (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
  outputSpacing[1] =
      spSrcImg->GetSpacing()[1] *
      (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
  outputSpacing[2] = spSrcImg->GetSpacing()[2];
  resample->SetOutputSpacing(outputSpacing);

  UShortImageType::PointType outputOrigin = spSrcImg->GetOrigin(); // Float
                                                                   // image
  resample->SetOutputOrigin(outputOrigin);

  resample->SetInput(spSrcImg);
  transform->SetIdentity();
  resample->SetTransform(transform);

  resample->Update();

  // resample->GetOutput()->SetOrigin(prevOrigin);
  spTarImg = resample->GetOutput(); // is it copied? or replaced?
}

void CbctRecon::ResampleItkImage2D(FloatImage2DType::Pointer &spSrcImg2D,
                                   FloatImage2DType::Pointer &spTarImg2D,
                                   double resFactor) {
  if (spSrcImg2D == nullptr) {
    std::cout << "ERROR! SrcImage is empty" << std::endl;
    return;
  }

  //  std::cout << "original Origin: " << spSrcImg2D->GetOrigin() << std::endl;
  using ResampleImageFilterType =
      itk::ResampleImageFilter<FloatImage2DType, FloatImage2DType, float>;
  ResampleImageFilterType::Pointer resample = ResampleImageFilterType::New();

  resample->SetOutputDirection(spSrcImg2D->GetDirection());

  // outputSpacing[2] = spSrcImg2D->GetSpacing()[2];

  // std::cout << "Output spacing: " << outputSpacing << std::endl;

  // FloatImageType2D::Pointer input = spSrcImg2D;
  // FloatImageType2D::PointType prevOrigin = input->GetOrigin(); //-204.6 -
  // 204.6  0

  // input->SetOrigin(outputOrigin);

  /* std::cout << "outputSize " << outputSize << std::endl;
   std::cout << "OutputSpacing " << outputSpacing << std::endl;
   std::cout << "OutputOrigin " << outputOrigin << std::endl;*/

  // Resample the image
  // typedef itk::IdentityTransform<float, 2> TransformType;

  using TransformType = itk::AffineTransform<float, 2>;
  TransformType::Pointer transform = TransformType::New();

  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<FloatImage2DType, float>;
  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  resample->SetInterpolator(interpolator);

  resample->SetDefaultPixelValue(50);

  FloatImage2DType::SizeType inputSize =
      spSrcImg2D->GetLargestPossibleRegion().GetSize();
  FloatImage2DType::SizeType outputSize{};
  outputSize[0] = qRound(inputSize[0] * resFactor);
  outputSize[1] = qRound(inputSize[1] * resFactor);
  resample->SetSize(outputSize);

  FloatImage2DType::SpacingType outputSpacing;
  outputSpacing[0] =
      spSrcImg2D->GetSpacing()[0] *
      (static_cast<double>(inputSize[0]) / static_cast<double>(outputSize[0]));
  outputSpacing[1] =
      spSrcImg2D->GetSpacing()[1] *
      (static_cast<double>(inputSize[1]) / static_cast<double>(outputSize[1]));
  resample->SetOutputSpacing(outputSpacing);

  FloatImage2DType::PointType outputOrigin =
      spSrcImg2D->GetOrigin(); // Float image
  resample->SetOutputOrigin(outputOrigin);

  resample->SetInput(spSrcImg2D);
  transform->SetIdentity();
  resample->SetTransform(transform);

  resample->Update();

  // resample->GetOutput()->SetOrigin(prevOrigin);
  spTarImg2D = resample->GetOutput(); // is it copied? or replaced?

  // std::cout << "resampled Origin: " << spTarImg2D->GetOrigin() << std::endl;
}

void CbctRecon::AfterScatCorrectionMacro() {
  // Original projection file can be replaced by the corrected one
  // Current projection map (float) used for the reconstruction is:
  // m_spProjImg3DFloat and this is resampled one
  ConvertIntensity2LineInt(m_spProjImgCorr3D, m_spProjImg3DFloat, 65535U);

  int iSizeZ = m_spProjImg3DFloat->GetRequestedRegion().GetSize()[2];

  // Update UI
  ui.pushButton_DoRecon->setEnabled(true);
  ui.spinBoxImgIdx->setMinimum(0);
  ui.spinBoxImgIdx->setMaximum(iSizeZ - 1);
  ui.spinBoxImgIdx->setValue(0);
  SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();
  SLT_DrawProjImages(); // Update Table is called

  // return;

  // Do reconstruction + //Update GUI

  // Regardeless of previous setting, The Truncation should not be applied!

  // Truncation is invalidated inside the function
  if (ui.radioButton_UseCUDA->isChecked()) {
    CudaDoReconstructionFDK(REGISTER_COR_CBCT);
  } else {
    DoReconstructionFDK(REGISTER_COR_CBCT);
  }
  // Skin removal (using CT contour w/ big margin)

  std::cout
      << "Post  FDK reconstruction is done. Moving on to post skin removal"
      << std::endl;

  m_pDlgRegistration->PostSkinRemovingCBCT(m_spRawReconImg);
  m_pDlgRegistration->PostSkinRemovingCBCT(m_spScatCorrReconImg);

  // 20151208 Removal of high intensity skin mask
  // Main issue: raw CBCT projection includes mask, deformed CT doesn't include
  // mask. In case of weight loss, mask signal is independent from skin contour,
  // but deformed CT cannot have that signal.  Therefore, after the subtraction
  // (CBCTcor projections), there is always a big peak. DIR quality doesn't
  // matter because it cannot 'create' mask signal anyway.  Assumption: near the
  // skin contour, this kind of discrepancy is not expected.
  // m_pDlgRegistration->ThermoMaskRemovingCBCT(m_spRawReconImg,
  // m_spScatCorrReconImg, threshold_HU);

  m_pDlgRegistration->UpdateListOfComboBox(0); // combo selection signalis
                                               // called
  m_pDlgRegistration->UpdateListOfComboBox(1);
  m_pDlgRegistration->SelectComboExternal(
      0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  m_pDlgRegistration->SelectComboExternal(1, REGISTER_COR_CBCT);

  m_pDlgRegistration
      ->SLT_DoLowerMaskIntensity(); // it will check the check button.
  std::cout << "Updating ReconImage..";
  QString updated_text = QString("Scatter corrected CBCT");
  UpdateReconImage(m_spScatCorrReconImg, updated_text); // main GUI update

  // Save Image as DICOM
  if (ui.checkBox_ExportVolDICOM->isChecked()) {
    // Get current folder
    QString strCrntDir = m_strPathPatientDir + "/" + "IMAGES" + "/" + "cor_" +
                         m_strDCMUID; // current Proj folder
    QDir crntDir(strCrntDir);
    QString SubDirName = "Reconstruction";
    bool tmpResult = crntDir.mkdir(SubDirName); // what if the directory exists?
    if (!tmpResult) {
      std::cout
          << "DICOM dir seems to exist already. Files will be overwritten."
          << std::endl;
    }
    QString strSavingFolder = strCrntDir + "/" + SubDirName;
    QString updated_text_ct = QString("PriorCT_ScatterCorr");
    SaveUSHORTAsSHORT_DICOM(m_spScatCorrReconImg, m_strDCMUID, updated_text_ct,
                            strSavingFolder);
    // Export as DICOM (using plastimatch) folder?
  }
  std::cout << "Exiting AfterScatCorrectionMacro.";
}

// called whenver recon 3D image for display changes.
void CbctRecon::UpdateReconImage(UShortImageType::Pointer &spNewImg,
                                 QString &fileName) {
  m_spCrntReconImg = spNewImg;

  UShortImageType::PointType origin_new = m_spCrntReconImg->GetOrigin();
  UShortImageType::SpacingType spacing_new = m_spCrntReconImg->GetSpacing();
  UShortImageType::SizeType size_new =
      m_spCrntReconImg->GetBufferedRegion().GetSize();

  std::cout << "New Origin" << origin_new << std::endl;
  std::cout << "New spacing" << spacing_new << std::endl;
  std::cout << "New size" << size_new << std::endl;

  ui.lineEdit_Cur3DFileName->setText(fileName);

  UShortImageType::SizeType size =
      m_spCrntReconImg->GetRequestedRegion().GetSize();

  m_dspYKReconImage->CreateImage(size[0], size[1], 0);

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

void CbctRecon::SLT_TempAudit() {
  if (m_spRawReconImg != nullptr) {
    std::cout << "m_spRawReconImg " << m_spRawReconImg << std::endl;
  }

  if (m_spRefCTImg != nullptr) {
    std::cout << "m_spRefCTImg " << m_spRefCTImg << std::endl;
  }

  if (m_spCrntReconImg != nullptr) {
    std::cout << "m_spCrntReconImg " << m_spCrntReconImg << std::endl;
  }
}

void CbctRecon::SaveUSHORTAsSHORT_DICOM(UShortImageType::Pointer &spImg,
                                        QString &strPatientID,
                                        QString &strPatientName,
                                        QString &strPathTargetDir) {
  if (spImg == nullptr) {
    return;
  }

  ShortImageType::Pointer spShortImg;
  ConvertUshort2Short(spImg, spShortImg);

  Plm_image plm_img(spShortImg);

  QString endFix = "_DCM";

  QString newDirPath = strPathTargetDir + "/" + strPatientID + "_DCM";

  QDir dirNew(newDirPath);
  if (!dirNew.exists()) {
    dirNew.mkdir(".");
  }

  Rt_study_metadata rsm;
  rsm.set_patient_id(strPatientID.toLocal8Bit().constData());
  rsm.set_patient_name(strPatientName.toLocal8Bit().constData());

  plm_img.save_short_dicom(newDirPath.toLocal8Bit().constData(), &rsm);
}

void CbctRecon::ConvertUshort2Short(UShortImageType::Pointer &spImgUshort,
                                    ShortImageType::Pointer &spImgShort) {
  // std::cout << "Before filter spImgUshort" << spImgUshort << std::endl;

  using ThresholdImageFilterType = itk::ThresholdImageFilter<UShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilter =
      ThresholdImageFilterType::New();
  thresholdFilter->SetInput(spImgUshort);
  thresholdFilter->ThresholdOutside(0, 4096); //--> 0 ~ 4095
  thresholdFilter->SetOutsideValue(0);
  thresholdFilter->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<UShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
  imageCalculatorFilter->Compute();
  auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());
  // std::cout <<"Min and Max Values are	" << minVal << "	" <<
  // maxVal
  // << std::endl; //should be 0 and 4096 std::cout <<"Min and Max Values are
  // " << minVal << "	" << maxVal << std::endl;

  // Min value is always 3024 --> outside the FOV
  auto outputMinVal = static_cast<SHORT_PixelType>(minVal - 1024);
  auto outputMaxVal = static_cast<SHORT_PixelType>(maxVal - 1024);

  // USHORT_PixelType outputMinVal = (USHORT_PixelType)(minVal - minVal);
  // USHORT_PixelType outputMaxVal = (USHORT_PixelType) (maxVal - minVal);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<UShortImageType, ShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();

  spImgShort = spRescaleFilter->GetOutput();

  // std::cout << "After filter spImgUshort" << spImgUshort << std::endl;
  // std::cout << "After filter spImgShort" << spImgShort << std::endl;
}

void CbctRecon::SLT_LoadPlanCT_USHORT() {
  // typedef itk::ImageFileWriter<FloatImageType> WriterType;
  using ReaderType = itk::ImageFileReader<UShortImageType>;
  ReaderType::Pointer reader = ReaderType::New();

  QString fileName =
      QFileDialog::getOpenFileName(this, "Open Image", m_strPathDirDefault,
                                   "Projection file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  reader->SetFileName(fileName.toLocal8Bit().constData());
  reader->Update();

  m_spRefCTImg = reader->GetOutput();
  QString ref_ct = QString("RefCT");
  UpdateReconImage(m_spRefCTImg, ref_ct);

  RegisterImgDuplication(REGISTER_REF_CT, REGISTER_MANUAL_RIGID);
}

void CbctRecon::SLT_CalcAndSaveAngularWEPL() // single point
{
  std::vector<WEPLData> vOutputWEPL;

  double fAngleGap = ui.lineEdit_WEPL_AngRes->text().toDouble();
  double fAngleStart = ui.lineEdit_AngStart->text().toDouble();
  double fAngleEnd = ui.lineEdit_AngEnd->text().toDouble();

  VEC3D curPOI{};
  curPOI.x = ui.lineEdit_ForcedProbePosX->text().toDouble(); // in mm
  curPOI.y = ui.lineEdit_ForcedProbePosY->text().toDouble();
  curPOI.z = ui.lineEdit_ForcedProbePosZ->text().toDouble();

  GetAngularWEPL_SinglePoint(m_spCrntReconImg, fAngleGap, fAngleStart,
                             fAngleEnd, curPOI, 0, vOutputWEPL, true);
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
void CbctRecon::SLT_DoScatterCorrectionUniform() {
  if (m_spProjImg3DFloat == nullptr) {
    return;
  }

  UShortImageType::Pointer spIntensityRaw;
  ConvertLineInt2Intensity(m_spProjImg3DFloat, spIntensityRaw, 65535);

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

  ConvertIntensity2LineInt(spIntensityUniformCorr, m_spProjImg3DFloat, 65535);

  // ConvertLineInt2Intensity(m_spProjImg3DFloat, m_spProjImgRaw3D, 65535);

  m_spProjImgRaw3D = spIntensityUniformCorr;

  ui.spinBoxImgIdx->setValue(0);
  SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();

  this->SLT_DrawProjImages(); // Update Table is called

  std::cout << "FINISHED!: Uniform Scatter correction for raw projection "
               "images (Boallaard method) is completed. Proceed to "
               "reconstruction"
            << std::endl;
}

void CbctRecon::SLT_FileExportShortDICOM_CurrentImg() {

  if (m_spCrntReconImg == nullptr) {
    return;
  }

  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), m_strPathDirDefault,
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
      strPatientID = m_strDCMUID;
    }
    // strPatientID = m_strDCMUID + "_" + strEndFix;
  } else {
    strPatientID = m_strDCMUID;
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
  SaveUSHORTAsSHORT_DICOM(m_spCrntReconImg, strPatientID, strFullName,
                          strSavingFolder);
}

void CbctRecon::SLT_AddConstHUToCurImg() {
  if (m_spCrntReconImg == nullptr) {
    return;
  }
  int addingVal = ui.lineEdit_AddConstHU->text().toInt();
  AddConstHU(m_spCrntReconImg, addingVal);
  QString updated_text = QString("Added%1").arg(addingVal);
  UpdateReconImage(m_spCrntReconImg, updated_text);
}

double CbctRecon::GetRawIntensityScaleFactor() {
  // GetRawIntensity Scale Factor
  double rawIntensityScaleF = 1.0;

  double fRef_mAs = 0.0;
  double fCur_mAs = 0.0;
  QString strRef_mAs = ui.lineEdit_RefmAs->text();
  QStringList listmAsRef = strRef_mAs.split(",");
  if (listmAsRef.length() == 2) {
    fRef_mAs = listmAsRef.at(0).toDouble() * listmAsRef.at(1).toDouble();
  }
  QString strCur_mAs = ui.lineEdit_CurmAs->text();
  QStringList listmAsCur = strCur_mAs.split(",");
  if (listmAsCur.length() == 2) {
    fCur_mAs = listmAsCur.at(0).toDouble() * listmAsCur.at(1).toDouble();
  }
  if (fRef_mAs * fCur_mAs != 0) {
    rawIntensityScaleF = fRef_mAs / fCur_mAs;
  }

  return rawIntensityScaleF;
  // if 64 40 ref, 40 40 cur --> scaleF = 1.6
  // raw intensity X scaleF ==> raw intensity increased --> this avoids negative
  // scatter map
}

void CbctRecon::SLT_SetCBCTSkinRSPath() {
  QString strPath =
      QFileDialog::getOpenFileName(this, "Open RS file", m_strPathDirDefault,
                                   "DICOM RS (*.dcm)", nullptr, nullptr);

  if (strPath.length() <= 1) {
    return;
  }

  ui.lineEdit_PathCBCTSkinPath->setText(strPath);
}
void CbctRecon::SLT_CropSkinUsingThreshold() {
  std::cout << "Overwriting of values below threshold to air ";
  if (m_spCrntReconImg == nullptr) {
    return;
  }

  using threshFilterType =
      itk::BinaryThresholdImageFilter<UShortImageType, UShortImageType>;
  threshFilterType::Pointer threshFilter = threshFilterType::New();
  threshFilter->SetInput(m_spCrntReconImg);

  threshFilter->SetOutsideValue(0);
  threshFilter->SetInsideValue(1);
  threshFilter->SetLowerThreshold(ui.lineEdit_Threshold->text().toInt());
  threshFilter->Update();
  UShortImageType::Pointer spCrntImgMask = threshFilter->GetOutput();
  using iteratorType = itk::ImageRegionIteratorWithIndex<UShortImageType>;
  iteratorType it(spCrntImgMask, spCrntImgMask->GetRequestedRegion());
  UShortImageType::SizeType imgDims =
      spCrntImgMask->GetBufferedRegion().GetSize();

  it.GoToBegin();
  while (!it.IsAtEnd()) {
    int z_idx = it.GetIndex()[2];
    if (z_idx == static_cast<int>(imgDims[2] - 10) || z_idx == 10) {
      it.Set(1.0f);
    }
    ++it;
  }

  using HoleFillingFilterType = itk::BinaryFillholeImageFilter<UShortImageType>;
  HoleFillingFilterType::Pointer HoleFillingFilter =
      HoleFillingFilterType::New();
  HoleFillingFilter->SetForegroundValue(1);
  HoleFillingFilter->SetFullyConnected(false);
  std::cout << "Threshold filtering.. ";
  HoleFillingFilter->SetInput(spCrntImgMask);

  using StructElementType =
      itk::BinaryBallStructuringElement<USHORT_PixelType, 3>;
  using ErodeFilterType =
      itk::BinaryErodeImageFilter<UShortImageType, UShortImageType,
                                  StructElementType>;
  using DilateFilterType =
      itk::BinaryDilateImageFilter<UShortImageType, UShortImageType,
                                   StructElementType>;
  ErodeFilterType::Pointer binaryErode = ErodeFilterType::New();
  binaryErode->SetErodeValue(1);
  StructElementType erodeStructElement;
  erodeStructElement.SetRadius(ui.lineEdit_ErodeRadius->text().toInt());
  erodeStructElement.CreateStructuringElement();
  binaryErode->SetKernel(erodeStructElement);
  std::cout << "filling holes.. ";
  binaryErode->SetInput(HoleFillingFilter->GetOutput());

  DilateFilterType::Pointer binaryDilate = DilateFilterType::New();
  binaryDilate->SetDilateValue(1);
  StructElementType dilateStructElement;
  dilateStructElement.SetRadius(ui.lineEdit_DilateRadius->text().toInt());
  dilateStructElement.CreateStructuringElement();
  binaryDilate->SetKernel(dilateStructElement);
  std::cout << "eroding dirt.. ";
  binaryDilate->SetInput(binaryErode->GetOutput());
  std::cout << "Skin mask is being created..." << std::endl;

  using MaskFilterType =
      itk::MaskImageFilter<UShortImageType, UShortImageType, UShortImageType>;
  MaskFilterType::Pointer MaskFilter = MaskFilterType::New();
  MaskFilter->SetMaskingValue(0);
  std::cout << "Dilating.. ";
  MaskFilter->SetMaskImage(binaryDilate->GetOutput());
  QString update_text = QString("Thresh-based skin cropped image");
  if (m_spCrntReconImg == m_spRawReconImg) {
    MaskFilter->SetInput(m_spRawReconImg);
    MaskFilter->Update();
    m_spRawReconImg = MaskFilter->GetOutput();
    UpdateReconImage(m_spRawReconImg, update_text);
  } else if (m_spCrntReconImg == m_spRefCTImg) {
    MaskFilter->SetInput(m_spRefCTImg);
    MaskFilter->Update();
    m_spRefCTImg = MaskFilter->GetOutput();
    UpdateReconImage(m_spRefCTImg, update_text);
  } else if (m_spCrntReconImg == m_spScatCorrReconImg) {
    MaskFilter->SetInput(m_spScatCorrReconImg);
    MaskFilter->Update();
    m_spScatCorrReconImg = MaskFilter->GetOutput();
    UpdateReconImage(m_spScatCorrReconImg, update_text);
  }
}

void CbctRecon::SLT_CropSkinUsingRS() {
  QString strPathRS = ui.lineEdit_PathCBCTSkinPath->text();
  if (strPathRS.length() < 1) {
    return;
  }

  double croppingMargin = ui.lineEdit_SkinMargin->text().toDouble();
  QString update_text = QString("RS-based skin cropped image");
  if (m_spCrntReconImg == m_spRawReconImg) {
    m_pDlgRegistration->CropSkinUsingRS(m_spRawReconImg, strPathRS,
                                        croppingMargin);
    UpdateReconImage(m_spRawReconImg, update_text);
  } else if (m_spCrntReconImg == m_spRefCTImg) {
    m_pDlgRegistration->CropSkinUsingRS(m_spRefCTImg, strPathRS,
                                        croppingMargin);
    UpdateReconImage(m_spRefCTImg, update_text);
  } else if (m_spCrntReconImg == m_spScatCorrReconImg) {
    m_pDlgRegistration->CropSkinUsingRS(m_spScatCorrReconImg, strPathRS,
                                        croppingMargin);
    UpdateReconImage(m_spScatCorrReconImg, update_text);
  }
}

// Below version is optimized for many points and much faster
void CbctRecon::ExportAngularWEPL_byFile(QString &strPathOutput) {
  if (strPathOutput.length() < 1) {
    return;
  }

  double fAngleGap = ui.lineEdit_WEPL_AngRes->text().toDouble();
  double fAngleStart = ui.lineEdit_AngStart->text().toDouble();
  double fAngleEnd = ui.lineEdit_AngEnd->text().toDouble();

  if (m_vPOI_DCM.empty()) {
    std::cout << "No POI data is prepared. Load them first" << std::endl;
    return;
  }

  if (m_spRawReconImg == nullptr) {
    std::cout << "Error: no Raw Recon image is found" << std::endl;
    return;
  }

  if (m_spScatCorrReconImg == nullptr) {
    std::cout << "Warning: no ScatCorrReconImg is found" << std::endl;
  }

  if (m_spManualRigidCT == nullptr) {
    std::cout << "Warning: no ManualRigidCT is found" << std::endl;
  }
  if (m_spAutoRigidCT == nullptr) {
    std::cout << "Warning: no AutoRigidCT is found" << std::endl;
  }
  if (m_spDeformedCT_Final == nullptr) {
    std::cout << "Warning: no DeformedCT is found" << std::endl;
  }

  std::vector<WEPLData> vOutputWEPL_manual;
  std::vector<WEPLData> vOutputWEPL_auto_rigid;
  std::vector<WEPLData> vOutputWEPL_deform;
  std::vector<WEPLData> vOutputWEPL_rawCBCT;
  std::vector<WEPLData> vOutputWEPL_corCBCT;
#pragma omp parallel sections
  {
#pragma omp section
    {
      GetAngularWEPL_window(m_spRawReconImg, fAngleGap, fAngleStart, fAngleEnd,
                            vOutputWEPL_rawCBCT, true); // mandatory
      std::cout << "Done: (RAW)";
    }
#pragma omp section
    {
      if (m_spScatCorrReconImg != nullptr) {
        try {
          GetAngularWEPL_window(m_spScatCorrReconImg, fAngleGap, fAngleStart,
                                fAngleEnd, vOutputWEPL_corCBCT, true);
          std::cout << " (COR)";
        } catch (std::exception &e) {
          std::cout << " (COR) failed!!: e=" << e.what() << std::endl;
        }
      }
    }
#pragma omp section
    {
      if (m_spManualRigidCT != nullptr) {
        try {
          GetAngularWEPL_window(m_spManualRigidCT, fAngleGap, fAngleStart,
                                fAngleEnd, vOutputWEPL_manual, true);
          std::cout << " (MAN)";
        } catch (std::exception &e) {
          std::cout << " (MAN) failed!!: e=" << e.what() << std::endl;
        }
      }
    }
#pragma omp section
    {
      if (m_spAutoRigidCT != nullptr) {
        try {
          GetAngularWEPL_window(m_spAutoRigidCT, fAngleGap, fAngleStart,
                                fAngleEnd, vOutputWEPL_auto_rigid, true);
          std::cout << " (AUT)";
        } catch (std::exception &e) {
          std::cout << " (AUT) failed!!: e=" << e.what() << std::endl;
        }
      }
    }
#pragma omp section
    {
      if (m_spDeformedCT_Final != nullptr) {
        try {
          GetAngularWEPL_window(m_spDeformedCT_Final, fAngleGap, fAngleStart,
                                fAngleEnd, vOutputWEPL_deform, true);
          std::cout << " (DEF)";
        } catch (std::exception &e) {
          std::cout << " (DEF) failed!!: e=" << e.what() << std::endl;
        }
      }
    }
  }
  std::cout << std::endl;

  std::cout << "Saving results...";

  std::ofstream fout;
  fout.open(strPathOutput.toLocal8Bit().constData());

  int cntWEPL = vOutputWEPL_rawCBCT.size();

  fout << "Point Index"
       << "\t"
       << "Gantry Angle"
       << "\t"
       << "Sample Number"
       << "\t"
       << "RawCBCT"
       << "\t";

  if ((m_spScatCorrReconImg != nullptr) &&
      static_cast<int>(vOutputWEPL_corCBCT.size()) == cntWEPL) {
    fout << "CorrCBCT"
         << "\t";
  }
  if ((m_spManualRigidCT != nullptr) &&
      static_cast<int>(vOutputWEPL_manual.size()) == cntWEPL) {
    fout << "ManualRigidCT"
         << "\t";
  }
  if ((m_spAutoRigidCT != nullptr) &&
      static_cast<int>(vOutputWEPL_auto_rigid.size()) == cntWEPL) {
    fout << "AutoRigidCT"
         << "\t";
  }
  if ((m_spDeformedCT_Final != nullptr) &&
      static_cast<int>(vOutputWEPL_deform.size()) == cntWEPL) {
    fout << "DeformedCT"
         << "\t";
  }
  fout << std::endl;

  for (int i = 0; i < cntWEPL; i++) {
    fout << vOutputWEPL_rawCBCT.at(i).ptIndex << "\t"
         << vOutputWEPL_rawCBCT.at(i).fGanAngle << "\t" << i << "\t"
         << vOutputWEPL_rawCBCT.at(i).fWEPL << "\t";

    if ((m_spScatCorrReconImg != nullptr) &&
        static_cast<int>(vOutputWEPL_corCBCT.size()) == cntWEPL) {
      fout << vOutputWEPL_corCBCT.at(i).fWEPL << "\t";
    }
    if ((m_spManualRigidCT != nullptr) &&
        static_cast<int>(vOutputWEPL_manual.size()) == cntWEPL) {
      fout << vOutputWEPL_manual.at(i).fWEPL << "\t";
    }
    if ((m_spAutoRigidCT != nullptr) &&
        static_cast<int>(vOutputWEPL_auto_rigid.size()) == cntWEPL) {
      fout << vOutputWEPL_auto_rigid.at(i).fWEPL << "\t";
    }
    if ((m_spDeformedCT_Final != nullptr) &&
        static_cast<int>(vOutputWEPL_deform.size()) == cntWEPL) {
      fout << vOutputWEPL_deform.at(i).fWEPL << "\t";
    }

    fout << std::endl;
  }
  fout.close();
  std::cout << "done!" << std::endl;
}

void CbctRecon::SLT_ExportAngularWEPL_byFile() {
  // export arrWEPL
  QString filePath = QFileDialog::getSaveFileName(
      this, "Save data", m_strPathDirDefault, "txt image file (*.txt)", nullptr,
      nullptr); // Filename don't need to exist

  if (filePath.length() < 1) {
    return;
  }

  ExportAngularWEPL_byFile(filePath);
}

void CbctRecon::GetAngularWEPL_window(UShortImageType::Pointer &spUshortImage,
                                      float fAngleGap, float fAngleStart,
                                      float fAngleEnd,
                                      std::vector<WEPLData> &vOutputWEPLData,
                                      bool bAppend) {
  if (spUshortImage == nullptr) {
    return;
  }

  if (fAngleGap <= 0) {
    return;
  }

  if (!bAppend) {
    vOutputWEPLData.clear();
  }

  FloatImageType::Pointer wepl_image = ConvertUshort2WeplFloat(spUshortImage);

  const double fullAngle = fAngleEnd - fAngleStart;
  const int sizeAngles = qRound(fullAngle / fAngleGap);

  const std::array<double, 3> pixel_size = {{spUshortImage->GetSpacing()[0],
                                             spUshortImage->GetSpacing()[1],
                                             spUshortImage->GetSpacing()[2]}};
  const std::array<size_t, 3> cubedim = {
      {spUshortImage->GetLargestPossibleRegion().GetSize()[0],
       spUshortImage->GetLargestPossibleRegion().GetSize()[1],
       spUshortImage->GetLargestPossibleRegion().GetSize()[2]}};

  const double couch = 0.0;

  for (int i = 0; i < sizeAngles; i++) {
    const double gantry = fAngleStart + i * fAngleGap;
    const std::array<double, 3> basis = get_basis_from_angles(gantry, couch);
    size_t loop_idx = 0;
    // int z_slice = -5000; // for progress and debug
    for (auto poi_it = m_vPOI_DCM.begin(); poi_it < m_vPOI_DCM.end();
         poi_it++, loop_idx++) {
      UShortImageType::PointType cur_point;
      cur_point[0] = poi_it->x;
      cur_point[1] = poi_it->y;
      cur_point[2] = poi_it->z;

      UShortImageType::IndexType cur_idx{};
      spUshortImage->TransformPhysicalPointToIndex(cur_point, cur_idx);
      /* // Turn on when debugging
      if (cur_idx[2] != z_slice) {
              z_slice = cur_idx[2];
              std::cout << z_slice << ", " << std::endl;
      }
      */
      const std::array<size_t, 3> point_id = {
          {static_cast<size_t>(cur_idx[0]), static_cast<size_t>(cur_idx[1]),
           static_cast<size_t>(cur_idx[2])}};

      WEPLData wepl_data{};
      wepl_data.fWEPL =
          WEPL_from_point(point_id, basis, pixel_size, cubedim, wepl_image);
      wepl_data.ptIndex = loop_idx;
      wepl_data.fGanAngle = gantry;

      vOutputWEPLData.push_back(wepl_data);
    }
  }
};

void CbctRecon::GetAngularWEPL_SinglePoint(
    UShortImageType::Pointer &spUshortImage, float fAngleGap, float fAngleStart,
    float fAngleEnd, VEC3D calcPt, int curPtIdx,
    std::vector<WEPLData> &vOutputWEPLData, bool bAppend) {
  if (spUshortImage == nullptr) {
    return;
  }

  if (fAngleGap <= 0) {
    return;
  }

  FloatImageType::Pointer wepl_image = ConvertUshort2WeplFloat(spUshortImage);

  const double fullAngle = fAngleEnd - fAngleStart;
  const auto sizeAngles = static_cast<size_t>(qRound(fullAngle / fAngleGap));

  UShortImageType::PointType calc_point;
  calc_point[0] = calcPt.x;
  calc_point[1] = calcPt.y;
  calc_point[2] = calcPt.z;
  UShortImageType::IndexType calc_idx{};

  if (!spUshortImage->TransformPhysicalPointToIndex(calc_point, calc_idx)) {
    std::cerr << "Point was outside image!" << std::endl;
  }

  // 1) Generate parms according to the angle e.g 360 parms
  const std::array<size_t, 3> isoTarget = {{static_cast<size_t>(calc_idx[0]),
                                            static_cast<size_t>(calc_idx[1]),
                                            static_cast<size_t>(calc_idx[2])}};

  std::cout << "Target: ( " << isoTarget[0] << ", " << isoTarget[1] << ", "
            << isoTarget[2] << " ), sizeAngles: " << sizeAngles << std::endl;

  const std::string stdout_file = "WEPL_stdout.txt";

  std::ofstream ofs(stdout_file); // Open stdout_file for writing
  if (ofs.is_open()) {
    std::cerr
        << "couldn't open file: " << stdout_file << " for writing!\n"
        << "Are you running this app from a folder without write permissions?"
        << std::endl;
    return;
  }
  const std::array<double, 3> pixel_size = {{wepl_image->GetSpacing()[0],
                                             wepl_image->GetSpacing()[1],
                                             wepl_image->GetSpacing()[2]}};

  FloatImageType::RegionType region = wepl_image->GetLargestPossibleRegion();
  const std::array<size_t, 3> cube_dim = {
      {region.GetSize()[0], region.GetSize()[1], region.GetSize()[2]}};

  std::vector<WEPLData> stArrWEPL(sizeAngles);

  size_t i = 0;
  for (auto it : stArrWEPL) {
    // YKTEMP Should be updated according to recent update of plastimatch
    const float curAngle = fAngleStart + i * fAngleGap;

    const std::array<double, 3> basis = get_basis_from_angles(curAngle, 0.0);

    ofs << std::fixed << std::setprecision(3) << curAngle << ", [" << basis[0]
        << ", " << basis[1] << ", " << basis[2] << "]: ";

    it.fWEPL = WEPL_from_point(isoTarget, basis, pixel_size, cube_dim,
                               wepl_image); // get_rgdepth
    it.fGanAngle = curAngle;
    it.ptIndex = curPtIdx;
    ofs << std::fixed << std::setprecision(5) << it.fWEPL << "\n";
    i++;
  }

  ofs.close();

  if (!bAppend) {
    vOutputWEPLData.clear();
  }
  vOutputWEPLData.insert(vOutputWEPLData.end(), stArrWEPL.begin(),
                         stArrWEPL.end());
  /*
for (int i = 0; i < sizeAngles; i++)
{
  vOutputWEPLData.push_back(stArrWEPL[i]);
}*/

  // delete[] stArrWEPL;
}

void CbctRecon::GetAngularWEPL_MultiPoint(
    UShortImageType::Pointer &spUshortImage, float fAngleGap, float fAngleStart,
    float fAngleEnd, std::vector<WEPLData> &vOutputWEPLData, bool bAppend) {
  if (spUshortImage == nullptr) {
    return;
  }

  if (fAngleGap <= 0) {
    return;
  }

  ShortImageType::Pointer spShortImg;
  ConvertUshort2Short(spUshortImage, spShortImg);

  using CastFilterType = itk::CastImageFilter<ShortImageType, FloatImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(spShortImg);
  castFilter->Update();

  // Plm_image::Pointer ct_vol = Plm_image::New (spShortImg);
  Plm_image::Pointer ct_vol = Plm_image::New(
      castFilter->GetOutput()); // Oh ohh.. not actually a pointer now is it?

  double fullAngle = fAngleEnd - fAngleStart;
  int sizeAngles = qRound(fullAngle / fAngleGap);

  // 1) Generate parms according to the angle e.g 360 parms
  double isoTarget[3];
  // isoTarget[0] = ui.lineEdit_ForcedProbePosX->text().toDouble(); //in mm
  // isoTarget[1] = ui.lineEdit_ForcedProbePosY->text().toDouble();
  // isoTarget[2] = ui.lineEdit_ForcedProbePosZ->text().toDouble();

  float srcDistance = 2200.0; // in mm, 2.2 m
  float srcProton[3] = {0.0, 0.0, 0.0};

  float ap_distance = 1900.0; // 300 mm from the target //offset from the source
                              // float ap_distance = 1000.0; //

  float ap_spacing[2] = {1.0, 1.0}; // resolution
  const plm_long ap_dim[2] = {1, 1};
  float ap_center[2] = {1, 1};

  float ray_step = 1.0; // mm

  double curAngle = 0.0;
  VEC3D curPOI{};

  size_t sizePOI = m_vPOI_DCM.size();
  std::vector<WEPLData> stArrWEPL;
  stArrWEPL.reserve(sizePOI * sizeAngles);

  // std::cout << "Everything is OK!? Let's calculate WEPL...\n" << std::endl;

  const std::string stdout_file = "WEPL_stdout.txt";

  std::ofstream ofs(stdout_file); // Open stdout_file for writing
  if (ofs.is_open()) {
    std::cerr
        << "couldn't open file: " << stdout_file << " for writing!\n"
        << "Are you running this app from a folder without write permissions?"
        << std::endl;
    return;
  }

#pragma omp parallel for
  for (int p = 0; p < static_cast<int>(sizePOI); p++) {
    curPOI = m_vPOI_DCM.at(p);
    isoTarget[0] = curPOI.x; // Sorry, might do better later..
    isoTarget[1] = curPOI.y;
    isoTarget[2] = curPOI.z;

    itk::Point<double, 3U> itk_isoTarget;
    itk_isoTarget[0] = curPOI.x;
    itk_isoTarget[1] = curPOI.y;
    itk_isoTarget[2] = curPOI.z;

    UCharImageType::RegionType region = spShortImg->GetLargestPossibleRegion();

    UCharImageType::Pointer target_bool = UCharImageType::New();
    target_bool->SetRegions(region);
    target_bool->SetDirection(spShortImg->GetDirection());
    target_bool->SetSpacing(spShortImg->GetSpacing());
    target_bool->SetOrigin(spShortImg->GetOrigin());
    target_bool->Allocate();
    UCharImageType::IndexType pixelIndex{};
    target_bool->TransformPhysicalPointToIndex(itk_isoTarget, pixelIndex);
    target_bool->SetPixel(pixelIndex, 255);
    // GetAngularWEPL_SinglePoint(m_spCrntReconImg, fAngleGap, curPOI, i,
    // vOutputWEPL, true);  if (m_spRawReconImg)
    // GetAngularWEPL_MultiPoint(m_spRawReconImg, fAngleGap, fAngleStart,
    // fAngleEnd, curPOI, i, vOutputWEPL_rawCBCT, true);//mandatory

    ofs << "POI idx: " << p
        << ", calculation for all given angles startng...\n";

    // for (int i = qRound(fAngleStart / fAngleGap); i < (sizeAngles +
    // qRound(fAngleStart / fAngleGap)); i++)
    for (int i = 0; i < sizeAngles; i++) {
      // YKTEMP Should be updated according to recent update of plastimatch
      curAngle = fAngleStart + i * fAngleGap;

      srcProton[0] =
          isoTarget[0] + (srcDistance * sin(curAngle * M_PI / 180.0));
      srcProton[1] =
          isoTarget[1] - (srcDistance * cos(curAngle * M_PI / 180.0));
      srcProton[2] = isoTarget[2];

      // std::cout << "\nDefine scene.." << std::endl;
      auto scene = std::make_shared<Rt_plan>();
      // std::cout << "New beam.." << std::endl;
      Rt_beam *newBeam = scene->append_beam(); // . to ->
      scene->set_patient(ct_vol);              // . to ->
      newBeam->get_aperture()->set_distance(ap_distance);
      newBeam->get_aperture()->set_distance(ap_distance);
      newBeam->get_aperture()->set_spacing(&ap_spacing[0]);
      newBeam->get_aperture()->set_dim(&ap_dim[0]);
      newBeam->get_aperture()->set_center(&ap_center[0]);

      newBeam->set_step_length(ray_step);
      newBeam->set_isocenter_position(&isoTarget[0]);
      newBeam->set_source_position(&srcProton[0]);

      // std::cout << "\nPrep beam.." << std::endl;
      // scene.prepare_beam_for_calc(newBeam);

      scene->set_target(target_bool);
      // std::string flavor("b");
      // newBeam->set_flavor(flavor); // b for beta to get fastest dose calc
      // scene->compute_dose(newBeam);

      Plm_image::Pointer ct_hu = Plm_image::New();
      ct_hu->set_volume(scene->get_patient_volume());

      using Float_pair_list = std::list<std::pair<float, float>>;
      Float_pair_list lookup;
      lookup.emplace_back(NLMIN(float), 0.0f);
      lookup.emplace_back(-1000.0f, 0.00106f);
      lookup.emplace_back(0.0f, 1.0f);
      lookup.emplace_back(41.46f, 1.048674f);
      lookup.emplace_back(NLMAX(float), 0.005011f);

      Volume::Pointer psp = volume_adjust(ct_hu->get_volume(), lookup);
      Plm_image::Pointer patient_psp = Plm_image::New(psp);

      newBeam->prepare_for_calc(ct_hu, patient_psp, newBeam->get_target());

      // scene->create_patient_psp(); // hu to stopping power
      // std::string flavor("b");
      // newBeam->set_flavor(flavor); // b for beta to get fastest dose calc
      // scene->compute_dose(newBeam); // ONLY way to call prepare_for_calc from
      // the "outside" (without ugly tricks)
      // newBeam->prepare_for_calc(ct_vol, newBeam->get_ct_psp(),
      // newBeam->get_target());
      // wed_ct_compute in wed_main

      Rpl_volume *rpl_vol = newBeam->rsp_accum_vol; // rpl_vol;

      // rpl_vol->compute_proj_wed_volume()

      Proj_volume *proj_vol = rpl_vol->get_proj_volume();
      // float *proj_wed_vol_img = (float*) rpl_vol->proj_wed_vol->img;
      // std::cout << "\nA lot of const doubles..";

      const double *src = proj_vol->get_src();
      const double *iso = proj_vol->get_iso();
      const double sid_length =
          proj_vol->get_proj_matrix()->sid; // distance from source to aperture
      double src_iso_vec[3];
      vec3_sub3(&src_iso_vec[0], src, iso);
      const double src_iso_distance = vec3_len(&src_iso_vec[0]);
      const double ap_iso_distance = src_iso_distance - sid_length;
      // std::cout << " front clipping plane..";
      double base_rg_dist =
          ap_iso_distance - rpl_vol->get_front_clipping_plane();

      const double base_dist =
          proj_vol->get_proj_matrix()->sid; // distance from source to aperture

      // const plm_long *ires = proj_vol->get_image_dim();

      // int ap_ij[2]; //ray index of rvol
      // plm_long ap_idx = 0;  //ray number always 0 here

      Ray_data *ray_data;
      double ray_ap[3];     // std::vector from src to ray intersection with ap
                            // plane
      double ray_ap_length; // length of std::vector from src to ray
                            // intersection with ap plane
      double rglength; // length that we insert into get_rgdepth for each ray
                       // std::cout << " get ray data.. " << std::endl;
      rpl_vol->compute_ray_data(); // computes p2 <- called in compute dose for
                                   // some flavors of beam
      ray_data = rpl_vol->get_ray_data();

      Ray_data *ray_data_single = &ray_data[0]; // ap_idx];

      /* Coordinate of ray intersection with aperture plane */
      double *ap_xyz = &ray_data_single->p2[0];
      vec3_sub3(&ray_ap[0], ap_xyz, src);
      ray_ap_length = vec3_len(&ray_ap[0]);
      rglength = base_rg_dist * (ray_ap_length / base_dist);

      plm_long ap_idx_default[2] = {0, 0};
      WEPLData curWEPLData{};
      curWEPLData.fWEPL =
          rpl_vol->get_value(ap_idx_default, rglength); // rgdepth
      curWEPLData.fGanAngle = curAngle;
      curWEPLData.ptIndex = p;

      stArrWEPL.push_back(curWEPLData);
      // This should be in the Rt_beam destructor, but it isn't
      delete newBeam->rsp_accum_vol;
      delete newBeam->hu_samp_vol;
      delete newBeam->sigma_vol;
      delete newBeam->rpl_vol_lg;
      delete newBeam->rpl_vol_samp_lg;
      delete newBeam->sigma_vol_lg;
      delete newBeam->dose_rv; // rpl_dose_vol;

      delete newBeam;
      // delete scene; not necessary for shared_ptr
    }
  }

  // std::cout << "Return results for saving.." << std::endl;

  if (!bAppend) {
    vOutputWEPLData.clear();
  }
  vOutputWEPLData.insert(vOutputWEPLData.end(), stArrWEPL.begin(),
                         stArrWEPL.end());
  /*
  for (int p = 0; p < sizePOI; p++)
  {
          for (auto i = int(fAngleStart / fAngleGap); i < sizeAngles; i++)
          {
                  vOutputWEPLData.push_back(stArrWEPL[i - qRound(fAngleStart /
  fAngleGap) + p * sizeAngles]);
          }
  }

  delete[] stArrWEPL;
  */
  ct_vol->free();
  ofs.close();
}

void CbctRecon::SLT_GeneratePOIData() // it fills m_vPOI_DCM
{
  if (!m_vPOI_DCM.empty()) {
    m_vPOI_DCM.clear();
  }
  /*
  VEC3D imgDims{};
  imgDims.x = ui.lineEdit_outImgDim_LR->text().toInt();
  imgDims.y = ui.lineEdit_outImgDim_AP->text().toInt();
  imgDims.z = ui.lineEdit_outImgDim_SI->text().toInt();
  VEC3D imgSpac{};
  imgSpac.x = ui.lineEdit_outImgSp_LR->text().toDouble();
  imgSpac.y = ui.lineEdit_outImgSp_AP->text().toDouble();
  imgSpac.z = ui.lineEdit_outImgSp_SI->text().toDouble();
  */
  UShortImageType::SizeType imgSize =
      m_spCrntReconImg->GetLargestPossibleRegion().GetSize();
  VEC3D imgDims = {static_cast<double>(imgSize[0]),
                   static_cast<double>(imgSize[1]),
                   static_cast<double>(imgSize[2])};
  VEC3D imgSpac{m_spCrntReconImg->GetSpacing()[0],
                m_spCrntReconImg->GetSpacing()[1],
                m_spCrntReconImg->GetSpacing()[2]};

  if (ui.checkBox_AP->isChecked()) {
    for (size_t k = 2; k < (imgDims.z - 2); k++) {
      for (size_t i = 2; i < (imgDims.x - 2); i++) {
        VEC3D fPOI = {i * imgSpac.x - ((imgSpac.x * imgDims.x) / 2.),
                      ui.lineEdit_PostTablePosY->text().toDouble(),
                      k * imgSpac.z - ((imgSpac.z * imgDims.z) / 2.)};
        m_vPOI_DCM.push_back(fPOI);
      }
    }
  } else {
    for (size_t k = 2; k < (imgDims.z - 2); k++) {
      for (size_t j = 2; j < (imgDims.y - 2); j++) {
        VEC3D fPOI = {2, j * imgSpac.y - ((imgSpac.y * imgDims.y) / 2.),
                      k * imgSpac.z - ((imgSpac.z * imgDims.z) / 2.)};
        m_vPOI_DCM.push_back(fPOI);
      }
    }
  }
  std::cout << "POI data generated! last value: [" << m_vPOI_DCM.back().x
            << ", " << m_vPOI_DCM.back().y << ", " << m_vPOI_DCM.back().z << "]"
            << std::endl;
}

void CbctRecon::SLT_LoadPOIData() // it fills m_vPOI_DCM
{
  if (!m_vPOI_DCM.empty()) {
    m_vPOI_DCM.clear();
  }

  QString filePath =
      QFileDialog::getOpenFileName(this, "POI data file", m_strPathDirDefault,
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

    m_vPOI_DCM.push_back(fPOI);
  }
  for (int i = 0; i < static_cast<int>(m_vPOI_DCM.size()); i++) {
    std::cout << "Data " << i << "	" << m_vPOI_DCM.at(i).x << ", "
              << m_vPOI_DCM.at(i).y << ", " << m_vPOI_DCM.at(i).z << std::endl;
  }
  std::cout << "POI data has been loaded. " << m_vPOI_DCM.size()
            << " data points are read" << std::endl;
  fin.close();
}

void CbctRecon::SLT_StartSyncFromSharedMem() {
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

void CbctRecon::SLT_StopSyncFromSharedMem() {
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

void CbctRecon::SLT_TimerEvent() {

#ifndef _WIN32
  std::cerr << "Function not implemented for non-windows systems!!"
            << std::endl;
  return;
#else
  if (m_busyTimer) {
    return;
  }

  if (m_arrYKImage.empty() || m_iImgCnt != 1) {
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
    m_arrYKImage.at(0).m_pData[i] =
        ((charBuf[idxA] << 8) | charBuf[idxB]); // little endian
  }

  ui.spinBoxImgIdx->setValue(0);
  SLT_DrawRawImages();

  CloseHandle(handle);

  m_busyTimer = false;
#endif
}

void CbctRecon::SLTM_ViewExternalCommand() { m_pDlgExternalCommand->show(); }

void CbctRecon::LoadExternalFloatImage(QString &strPath, bool bConversion) {
  using ReaderType = itk::ImageFileReader<FloatImageType>;
  ReaderType::Pointer reader = ReaderType::New();

  // QString filePath = strPath;

  reader->SetFileName(strPath.toLocal8Bit().constData());
  reader->Update();

  FloatImageType::Pointer spCrntImg = reader->GetOutput();

  // Float image
  std::cout << "Float image has been loaded" << std::endl;

  if (bConversion) {
    TransformationRTK2IEC(spCrntImg);
  }

  using AbsImageFilterType =
      itk::AbsImageFilter<FloatImageType, FloatImageType>;
  AbsImageFilterType::Pointer absImgFilter = AbsImageFilterType::New();
  absImgFilter->SetInput(spCrntImg); // 20140206 modified it was a bug

  using MultiplyImageFilterType =
      itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  MultiplyImageFilterType::Pointer multiplyImageFilter =
      MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(absImgFilter->GetOutput());
  multiplyImageFilter->SetConstant(65536); // calculated already

  using CastFilterType = itk::CastImageFilter<FloatImageType, UShortImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(multiplyImageFilter->GetOutput());
  castFilter->Update();
  m_spRawReconImg = castFilter->GetOutput();

  QString strCrntFileName;
  QFileInfo outFileInfo(strPath);
  strCrntFileName = outFileInfo.fileName();

  UpdateReconImage(m_spRawReconImg, strCrntFileName);
}

void CbctRecon::TransformationRTK2IEC(FloatImageType::Pointer &spSrcTarg) {
  FloatImageType::SizeType sizeOutput =
      spSrcTarg->GetBufferedRegion().GetSize();
  FloatImageType::SpacingType spacingOutput = spSrcTarg->GetSpacing();

  // Transformation is applied
  std::cout << "Euler 3D Transformation: from RTK-procuded volume to standard "
               "DICOM coordinate"
            << std::endl;
  // Same image type from original image -3D & float
  FloatImageType::IndexType start_trans{};
  start_trans[0] = 0;
  start_trans[1] = 0;
  start_trans[2] = 0;

  FloatImageType::SizeType size_trans{};
  size_trans[0] = sizeOutput[0]; // X //410
  size_trans[1] = sizeOutput[2]; // Y  // 410
  size_trans[2] = sizeOutput[1]; // Z // 120?

  FloatImageType::SpacingType spacing_trans;
  spacing_trans[0] = spacingOutput[0];
  spacing_trans[1] = spacingOutput[2];
  spacing_trans[2] = spacingOutput[1];

  FloatImageType::PointType Origin_trans;
  Origin_trans[0] = -0.5 * size_trans[0] * spacing_trans[0];
  Origin_trans[1] = -0.5 * size_trans[1] * spacing_trans[1];
  Origin_trans[2] = -0.5 * size_trans[2] * spacing_trans[2];

  FloatImageType::RegionType region_trans;
  region_trans.SetSize(size_trans);
  region_trans.SetIndex(start_trans);

  /* 2) Prepare Target image */
  const FloatImageType::Pointer &targetImg = spSrcTarg;

  /* 3) Configure transform */
  using TransformType = itk::Euler3DTransform<double>;
  TransformType::Pointer transform = TransformType::New();

  TransformType::ParametersType param;
  param.SetSize(6);
  // MAXIMUM PARAM NUMBER: 6!!!
  param.put(0, 0.0);                  // rot X // 0.5 = PI/2
  param.put(1, itk::Math::pi / 2.0);  // rot Y
  param.put(2, itk::Math::pi / -2.0); // rot Z
  param.put(3, 0.0);                  // Trans X mm
  param.put(4, 0.0);                  // Trans Y mm
  param.put(5, 0.0);                  // Trans Z mm

  TransformType::ParametersType fixedParam(3); // rotation center
  fixedParam.put(0, 0);
  fixedParam.put(1, 0);
  fixedParam.put(2, 0);

  transform->SetParameters(param);
  transform->SetFixedParameters(fixedParam); // Center of the Transform

  std::cout << "Transform matrix:"
            << "	" << std::endl;
  std::cout << transform->GetMatrix() << std::endl;

  using ResampleFilterType =
      itk::ResampleImageFilter<FloatImageType, FloatImageType>;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();
  // FloatImageType::RegionType fixedImg_Region =
  // fixedImg->GetLargestPossibleRegion().GetSize();

  resampler->SetInput(targetImg);
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(targetImg->GetDirection()); // image normal?
  resampler->SetTransform(transform);

  // LR flip

  std::cout << "LR flip filter is being applied" << std::endl;

  using FilterType = itk::FlipImageFilter<FloatImageType>;

  FilterType::Pointer flipFilter = FilterType::New();
  using FlipAxesArrayType = FilterType::FlipAxesArrayType;

  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;

  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(resampler->GetOutput());

  flipFilter->Update();

  spSrcTarg = flipFilter->GetOutput();
}

void CbctRecon::SLTM_LoadDICOMdir() {
  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open Directory"), m_strPathDirDefault,
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.length() <= 1) {
    return;
  }

  Plm_image plmImg;
  plmImg.load_native(dirPath.toLocal8Bit().constData());
  ShortImageType::Pointer spShortImg = plmImg.itk_short();

  // Figure out whether this is NKI
  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  // imageCalculatorFilter->SetImage(spShortImg);
  // imageCalculatorFilter->Compute();

  // double minVal0 = (double)(imageCalculatorFilter->GetMinimum());
  // double maxVal0 = (double)(imageCalculatorFilter->GetMaximum());

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilter =
      ThresholdImageFilterType::New();

  thresholdFilter->SetInput(spShortImg);
  thresholdFilter->ThresholdOutside(-1024, 3072); //--> 0 ~ 4095
  thresholdFilter->SetOutsideValue(-1024);
  thresholdFilter->Update();

  imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
  imageCalculatorFilter->Compute();

  auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Current Min and Max Values are	" << minVal << "	"
            << maxVal << std::endl;

  auto outputMinVal = static_cast<USHORT_PixelType>(minVal + 1024);
  auto outputMaxVal = static_cast<USHORT_PixelType>(maxVal + 1024);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();

  // m_spRawReconImg = spRescaleFilter->GetOutput();
  m_spRefCTImg = spRescaleFilter->GetOutput();
  QString update_text = QString("DICOM reference image");
  UpdateReconImage(m_spRefCTImg, update_text);

  RegisterImgDuplication(REGISTER_REF_CT, REGISTER_MANUAL_RIGID);
}

void CbctRecon::SLTM_LoadRTKoutput() {
  QString filePath = QFileDialog::getOpenFileName(
      this, "Open Image", m_strPathDirDefault, "rtk output float image (*.mha)",
      nullptr, nullptr);
  LoadExternalFloatImage(filePath, true);
}
// Only can be used for m_spRawRecon
void CbctRecon::MedianFilterByGUI() {
  if (m_spCrntReconImg == nullptr) {
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
    using FilterType = itk::MedianImageFilter<UShortImageType, UShortImageType>;
    FilterType::Pointer medFilter = FilterType::New();

    // this is radius. 1 --> median window 3
    std::cout << "Post median(3D) filtering is under progress..Size(radius X Y "
                 "Z) is = "
              << indexRadius << std::endl;

    medFilter->SetRadius(indexRadius);
    medFilter->SetInput(m_spCrntReconImg);
    medFilter->Update();

    m_spCrntReconImg = medFilter->GetOutput();
    std::cout << "median filtering has been done" << std::endl;

    QString prevFileName = ui.lineEdit_Cur3DFileName->text();

    UpdateReconImage(m_spCrntReconImg, prevFileName.append("_med"));
  } else {
    std::cout << "Not valid median window" << std::endl;
  }
}

void CbctRecon::SLT_OutPathEdited() {
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

// Only can be used for m_spRawRecon
void CbctRecon::FileExportByGUI() // USHORT
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
    writer->SetInput(m_spRawReconImg);

    std::cout << "Writing the image to: "
              << outputFilePath.toLocal8Bit().constData() << std::endl;

    TRY_AND_EXIT_ON_ITK_EXCEPTION(writer->Update());

    std::cout << std::endl;
    std::cout << "The output image was successfully saved" << std::endl;
  }
}

void CbctRecon::SLT_MedianFilterDoNow() { MedianFilterByGUI(); }

void CbctRecon::SLT_Export2DDose_TIF() // 2D dose from current displayed image
                                       // of reconstruction
{
  if (m_dspYKReconImage == nullptr) {
    return;
  }

  if (m_spCrntReconImg == nullptr) {
    return;
  }

  QString strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "signed short meta image (*.tif)", nullptr,
      nullptr);
  if (strPath.length() <= 1) {
    return;
  }

  auto originLeft = static_cast<double>(m_spCrntReconImg->GetOrigin()[0]);
  auto originTop =
      static_cast<double>(m_spCrntReconImg->GetOrigin()[1]); // not sure...

  auto spacingX = static_cast<double>(m_spCrntReconImg->GetSpacing()[0]);
  auto spacingY =
      static_cast<double>(m_spCrntReconImg->GetSpacing()[1]); // not sure...

  if (!SaveDoseGrayImage(strPath.toLocal8Bit().constData(),
                         m_dspYKReconImage->m_iWidth,
                         m_dspYKReconImage->m_iHeight, spacingX, spacingY,
                         originLeft, originTop, m_dspYKReconImage->m_pData)) {
    std::cout << "Failed in save gray dose file" << std::endl;
  } else {
    std::cout << "image exported successfully." << std::endl;
  }
}

void CbctRecon::SLTM_Export2DDoseMapAsMHA() {
  if (m_dspYKReconImage == nullptr) {
    return;
  }

  if (m_spCrntReconImg == nullptr) {
    return;
  }

  QString strPath = QFileDialog::getSaveFileName(
      this, "Save Image", "", "itk compatible meta image (*.mha)", nullptr,
      nullptr);
  if (strPath.length() <= 1) {
    return;
  }

  auto originLeft = static_cast<double>(m_spCrntReconImg->GetOrigin()[0]);
  auto originTop =
      static_cast<double>(m_spCrntReconImg->GetOrigin()[1]); // not sure...

  auto spacingX = static_cast<double>(m_spCrntReconImg->GetSpacing()[0]);
  auto spacingY =
      static_cast<double>(m_spCrntReconImg->GetSpacing()[1]); // not sure...

  // Export float 2D image
  FloatImage2DType::Pointer doseImg2D = FloatImage2DType::New();
  FloatImage2DType::SizeType doseSize{};
  doseSize[0] = m_dspYKReconImage->m_iWidth;
  doseSize[1] = m_dspYKReconImage->m_iHeight;

  FloatImage2DType::IndexType doseStart{};
  doseStart[0] = 0;
  doseStart[1] = 0;

  FloatImage2DType::RegionType doseRegion;
  doseRegion.SetSize(doseSize);
  doseRegion.SetIndex(doseStart);

  FloatImage2DType::SpacingType doseSpacing;
  doseSpacing[0] = spacingX;
  doseSpacing[1] = spacingY;

  FloatImage2DType::PointType doseOrigin;
  doseOrigin[0] = originLeft;
  doseOrigin[1] = originTop;

  doseImg2D->SetRegions(doseRegion);
  doseImg2D->SetSpacing(doseSpacing);
  doseImg2D->SetOrigin(doseOrigin);

  doseImg2D->Allocate();
  doseImg2D->FillBuffer(0);

  double factor_ushort2float = 0.01; // cGy --> Gy

  itk::ImageRegionIterator<FloatImage2DType> it(
      doseImg2D, doseImg2D->GetLargestPossibleRegion());

  float pixel_val = 0.0f;
  int i = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    pixel_val = static_cast<double>(m_dspYKReconImage->m_pData[i]) *
                factor_ushort2float;
    it.Set(pixel_val);
    i++;
  }
  // YK201502
  using WriterType = itk::ImageFileWriter<FloatImage2DType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(strPath.toLocal8Bit().constData());
  writer->SetUseCompression(true);
  writer->SetInput(doseImg2D);
  writer->Update();

  std::cout << "File was exported successfully" << std::endl;
}

void CbctRecon::SLTM_ExportProjGeometryTXT() {
  // if (!m_spFullGeometry)
  //	return;

  if (m_spCustomGeometry ==
      nullptr) { // will be filled after Projection load button is pushed
    return;
  }

  QString strPath = QFileDialog::getSaveFileName(
      this, "Save text file", "", "text (*.txt)", nullptr, nullptr);

  if (strPath.length() <= 1) {
    return;
  }

  std::vector<double>::const_iterator itAng, itShiftX, itShiftY;

  int cntAngle = m_spCustomGeometry->GetGantryAngles().size();
  int cntShiftX = m_spCustomGeometry->GetProjectionOffsetsX().size();
  int cntShiftY = m_spCustomGeometry->GetProjectionOffsetsY().size();

  if (cntAngle <= 0) {
    std::cout << "Error! no angle std::vector is found" << std::endl;
    return;
  }

  if (cntAngle != cntShiftX || cntAngle != cntShiftY) {
    std::cout << "Error! Angle number and shift number are not matching."
              << std::endl;
    return;
  }

  itShiftX = m_spCustomGeometry->GetProjectionOffsetsX().begin();
  itShiftY = m_spCustomGeometry->GetProjectionOffsetsY().begin();

  std::ofstream fout;
  fout.open(strPath.toLocal8Bit().constData());

  fout << "MV_Gantry_Angle"
       << "	"
       << "PanelShiftX(mm)"
       << "	"
       << "PanelShiftY(mm)" << std::endl;

  for (itAng = m_spCustomGeometry->GetGantryAngles().begin();
       itAng != m_spCustomGeometry->GetGantryAngles().end(); itAng++) {
    fout << (*itAng) << "	" << (*itShiftX) << "	" << (*itShiftY)
         << std::endl;

    itShiftX++;
    itShiftY++;
  }

  fout.close();
}

void CbctRecon::LoadXVIGeometryFile(const char *filePath) {
  QString strFilePath = filePath;

  m_spFullGeometry = GeometryType::New();

  FLEXDATA flxData{};

  /* We'll parse the example.xml */
  auto *file = new QFile(strFilePath);
  /* If we can't open it, let's show an error message. */
  if (!file->open(QIODevice::ReadOnly | QIODevice::Text)) {
    return;
  }
  /* QXmlStreamReader takes any QIODevice. */
  QXmlStreamReader xml(file);
  // QList< QMap<QString, QString> > persons;
  /* We'll parse the XML until we reach end of it.*/

  m_vExcludeProjIdx.clear();
  int iIdx = 0;

  while (!xml.atEnd() && !xml.hasError()) {
    /* Read next element.*/
    QXmlStreamReader::TokenType token = xml.readNext();
    /* If token is just StartDocument, we'll go to next.*/
    if (token == QXmlStreamReader::StartDocument) {
      continue;
    }
    /* If token is StartElement, we'll see if we can read it.*/
    if (token == QXmlStreamReader::StartElement) {
      /* If it's named persons, we'll go to the next.*/
      if (xml.name() == "Frames") {
        continue;
      }
      /* If it's named person, we'll dig the information from there.*/
      if (xml.name() == "Frame") {
        flxData = XML_parseFrameForXVI5(xml);
        // m_vRetroFlexmap.push_back(flxData);

        if (flxData.fGanAngle < 0) {
          flxData.fGanAngle = flxData.fGanAngle + 360.0;
        }

        flxData.fPanelOffsetX = -flxData.fPanelOffsetX;
        flxData.fPanelOffsetY = -flxData.fPanelOffsetY;

        if (!flxData.bKV_On) {
          m_vExcludeProjIdx.push_back(iIdx);
        }

        /*               if (flxData.bKV_On)
                           m_vExcludeProjIdx.push_back(iIdx);*/

        ////Image qual test
        // flxData.fGanAngle = -flxData.fGanAngle;
        // if (flxData.fGanAngle < 0)
        //	flxData.fGanAngle = flxData.fGanAngle + 360.0;
        // flxData.fPanelOffsetX = -flxData.fPanelOffsetX;
        // flxData.fPanelOffsetY = -flxData.fPanelOffsetY;

        m_spFullGeometry->AddProjection(1000.0, 1536.0, flxData.fGanAngle,
                                        flxData.fPanelOffsetX,
                                        flxData.fPanelOffsetY, // Flexmap
                                        0.0, 0.0,  // In elekta, these are 0
                                        0.0, 0.0); // In elekta, these are 0

        iIdx++;
      }
    }
  }
  /* Error handling. */
  if (xml.hasError()) {
    QMessageBox::critical(this, "SLT_SetFlexmap", xml.errorString(),
                          QMessageBox::Ok);
  }
  /* Removes any device() or data from the reader
   * and resets its internal state to the initial state. */
  xml.clear();
}

FLEXDATA CbctRecon::XML_parseFrameForXVI5(QXmlStreamReader &xml) {
  FLEXDATA tmpResult{};
  tmpResult.fGanAngle = 0.0;
  tmpResult.fPanelOffsetX = 0.0;
  tmpResult.fPanelOffsetY = 0.0;
  tmpResult.bKV_On = true;
  tmpResult.bMV_On = false;

  /* Let's check that we're really getting a person. */
  if (xml.tokenType() != QXmlStreamReader::StartElement &&
      xml.name() == "Frame") {
    return tmpResult;
  }
  /* Let's get the attributes for person */
  // QXmlStreamAttributes attributes = xml.attributes();
  /* Let's check that person has id attribute. */
  // if (attributes.hasAttribute("id")) {
  //	/* We'll add it to the map. */
  //	person["id"] = attributes.value("id").toString();
  //}
  /* Next element... */
  xml.readNext();
  /*
   * We're going to loop over the things because the order might change.
   * We'll continue the loop until we hit an EndElement named person.
   */
  while (!(xml.tokenType() == QXmlStreamReader::EndElement &&
           xml.name() == "Frame")) {
    QStringRef tmpXmlName = (xml.name());
    QString strTmpXMLName = QString(tmpXmlName.toLocal8Bit().constData());
    // int tmpType = (int)(xml.tokenType());

    QString tmpStr;
    if (xml.tokenType() == QXmlStreamReader::StartElement) {
      /* We've found first name. */
      if (xml.name() == "Seq") {
        tmpStr = XML_GetSingleItemString(xml);
      }
      /* We've found surname. */
      else if (xml.name() == "DeltaMS") {
        tmpStr = XML_GetSingleItemString(xml);
      }
      /* We've found email. */
      else if (xml.name() == "HasPixelFactor") {
        tmpStr = XML_GetSingleItemString(xml);
      }
      /* We've found website. */
      else if (xml.name() == "PixelFactor") {
        tmpStr = XML_GetSingleItemString(xml);
      } else if (xml.name() == "GantryAngle") {
        tmpStr = XML_GetSingleItemString(xml);
        tmpResult.fGanAngle = tmpStr.toDouble();
      } else if (xml.name() == "Exposed") {
        tmpStr = XML_GetSingleItemString(xml);
        if (tmpStr == "True") {
          tmpResult.bKV_On = true;
        } else {
          tmpResult.bKV_On = false;
        }
      }
      /*else if (xml.name() == "Exposed") {
          tmpStr = XML_GetSingleItemString(xml);
          if (tmpStr == "True")
              tmpResult.bKV_On = true;
          else
              tmpResult.bKV_On = false;
      }*/
      else if (xml.name() == "MVOn") {
        tmpStr = XML_GetSingleItemString(xml);
        if (tmpStr == "True") {
          tmpResult.bMV_On = true;
        } else {
          tmpResult.bMV_On = false;
        }
      } else if (xml.name() == "UCentre") {
        tmpStr = XML_GetSingleItemString(xml);
        tmpResult.fPanelOffsetX = tmpStr.toDouble();
      } else if (xml.name() == "VCentre") {
        tmpStr = XML_GetSingleItemString(xml);
        tmpResult.fPanelOffsetY = tmpStr.toDouble();
      } else if (xml.name() == "Inactive") {
        tmpStr = XML_GetSingleItemString(xml);
      }
    }
    xml.readNext();
  }
  return tmpResult;
}

QString CbctRecon::XML_GetSingleItemString(QXmlStreamReader &xml) {
  QString strResult = "";
  /* We need a start element, like <foo> */
  if (xml.tokenType() != QXmlStreamReader::StartElement) {
    return strResult;
  }

  /* Let's read the name... */
  QString elementName = xml.name().toString();
  /* ...go to the next. */
  xml.readNext();
  /*
   * This elements needs to contain Characters so we know it's
   * actually data, if it's not we'll leave.
   */
  if (xml.tokenType() != QXmlStreamReader::Characters) {
    return strResult;
  }
  strResult = xml.text().toString();
  return strResult;
}

void CbctRecon::SLTM_ForwardProjection() {
  if (m_spRawReconImg == nullptr) {
    return;
  }

  GeometryType::Pointer crntGeometry = GeometryType::New();

  if (m_spCustomGeometry == nullptr) {
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

    ForwardProjection(m_spRawReconImg, crntGeometry, m_spProjImgRaw3D, false);
    // Save proj3D;

    // QString outputPath = "D:/ProjTemplate.mha";
    // QString outputPath = "D:/2D3DRegi/FwdProj_0.mha";
    QString outputPath = QFileDialog::getSaveFileName(
        this, "File path to save", m_strPathDirDefault,
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
    writer->SetInput(m_spProjImgRaw3D);
    writer->Update();

    return;
  }
  // if there is a geometry

  int cntProj = m_spCustomGeometry->GetGantryAngles().size();

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
    double curSID = m_spCustomGeometry->GetSourceToIsocenterDistances().at(i);
    double curSDD = m_spCustomGeometry->GetSourceToDetectorDistances().at(i);
    double curGantryAngle = m_spCustomGeometry->GetGantryAngles().at(i);
    double kVAng = curGantryAngle * 360. / (2. * itk::Math::pi);
    double MVAng = kVAng - (hisIsUsed ? 0.0 : 90.0);
    if (MVAng < 0.0) {
      MVAng = MVAng + 360.0;
    }
    curGantryAngle = MVAng;

    double curProjOffsetX = m_spCustomGeometry->GetProjectionOffsetsX().at(i);
    double curProjOffsetY = m_spCustomGeometry->GetProjectionOffsetsY().at(i);

    double curOutOfPlaneAngles =
        m_spCustomGeometry->GetOutOfPlaneAngles().at(i);
    double curInPlaneAngles = m_spCustomGeometry->GetInPlaneAngles().at(i);

    double curSrcOffsetX = m_spCustomGeometry->GetSourceOffsetsX().at(i);
    double curSrcOffsetY = m_spCustomGeometry->GetSourceOffsetsY().at(i);

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

  ForwardProjection(m_spRawReconImg, crntGeometry, m_spProjImgRaw3D, true);

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

void CbctRecon::SLTM_FineResolScatterCorrectrionMacro() {
  float curResampleF = m_fResampleF;
  ui.lineEdit_DownResolFactor->setText("1.0");
  SLT_LoadSelectedProjFiles();
  m_fResampleF = curResampleF;
  ui.lineEdit_DownResolFactor->setText(QString("%1").arg(m_fResampleF));

  // Scatter correction

  SLT_DoScatterCorrection_APRIORI();
}

void CbctRecon::SLTM_FullScatterCorrectionMacroAP() // single. should be called
                                                    // after HIS folder is
                                                    // defined
{
  if (m_strPathPatientDir.length() < 2) {
    return;
  }

  enREGI_IMAGES enRegImg = REGISTER_DEFORM_FINAL;
  bool bFullResolForFinalRecon = false;

  bool bIntensityShift = true;

  FullScatterCorrectionMacroSingle(m_strPathPatientDir, enRegImg,
                                   bFullResolForFinalRecon, false,
                                   bIntensityShift);
}

void CbctRecon::SLTM_BatchScatterCorrectionMacroAP() {
  // Scatter parameters
  QTime batchmodeTime = QTime::currentTime();

  // 1) Get img_ file lists
  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open IMAGES Directory"), m_strPathDirDefault,
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
      this, tr("Open Output Directory"), m_strPathDirDefault,
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

  /*QMessageBox msgBox;
  QString strMsg = "Full-resolution reconstruction after scatter generation?";
  msgBox.setText(strMsg);
  msgBox.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  int res = msgBox.exec();
  if (res == QMessageBox::Yes)
  {
  bFullResolForFinalRecon = true;
  }*/

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
      SetProjDir(curProjDirPath);
      FullScatterCorrectionMacroSingle(strOutDirPath, enRegImg,
                                       bFullResolForFinalRecon,
                                       bExportShortImages, bIntensityShift);
      cntHisDir++;
    }
  }

  /*QMessageBox msgBoxFinal;
  msgBoxFinal.setText("All Done!");
  msgBoxFinal.exec();*/

  float elapsedSec = batchmodeTime.elapsed() / 1000.0;

  std::cout << "Batch mode calculation is done! "
            << QString::number(elapsedSec, 'f', 2).toLocal8Bit().constData()
            << " seconds was spent for " << cntHisDir << " cases" << std::endl;
}

bool CbctRecon::FullScatterCorrectionMacroSingle(QString &outputDirPath,
                                                 enREGI_IMAGES enFwdRefImg,
                                                 bool bFullResolRecon,
                                                 bool bExportImages,
                                                 bool bCBCT_IntensityShift) {
  if (m_strDCMUID.length() < 1) {
    return false;
  }

  m_bMacroContinue = true;

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
    CropFOV3D(m_spRawReconImg, physPosX, physPosY, physRadius, physTablePosY);
  }

  SLT_ViewRegistration();

  m_pDlgRegistration->SLT_PreProcessCT();
  if (!m_bMacroContinue) {
    std::cout << "Stopped during MacroSingle due to error in PreProcessCT"
              << std::endl;
    return false;
  }

  QString strSuffix;
  switch (enFwdRefImg) {
  case REGISTER_MANUAL_RIGID:
    m_pDlgRegistration->SLT_ManualMoveByDCMPlan();
    m_pDlgRegistration
        ->SLT_ConfirmManualRegistration(); // skin cropping for
                                           // CBCT. only works when
                                           // CBCT_skin crop is on
    strSuffix = strSuffix + "_man";

    if (bCBCT_IntensityShift) {
      m_pDlgRegistration->SLT_IntensityNormCBCT();
    }

    break;
  case REGISTER_AUTO_RIGID:
    m_pDlgRegistration->SLT_ManualMoveByDCMPlan();
    m_pDlgRegistration->SLT_ConfirmManualRegistration(); // skin cropping

    if (bCBCT_IntensityShift) {
      m_pDlgRegistration->SLT_IntensityNormCBCT();
    }

    // OPtional
    if (bFOVCropping) {
      CropFOV3D(m_spManualRigidCT, physPosX, physPosY, physRadius,
                physTablePosY);
    }

    m_pDlgRegistration->SLT_DoRegistrationRigid();
    strSuffix = strSuffix + "_rigid";
    break;
  case REGISTER_DEFORM_FINAL:
    std::cout << "REGISTER_DEFORM_FINAL was chosen." << std::endl;
    m_pDlgRegistration->SLT_ManualMoveByDCMPlan();
    m_pDlgRegistration->SLT_ConfirmManualRegistration(); // skin cropping

    if (bCBCT_IntensityShift) {
      std::cout << "IntensityShift is underway" << std::endl;
      m_pDlgRegistration->SLT_IntensityNormCBCT();
    }

    if (bFOVCropping) {
      CropFOV3D(m_spManualRigidCT, physPosX, physPosY, physRadius,
                physTablePosY);
    }

    m_pDlgRegistration->SLT_DoRegistrationRigid();
    m_pDlgRegistration->SLT_DoRegistrationDeform();
    strSuffix = strSuffix + "_defrm";
    break;

  case REGISTER_DEFORM_SKIP_AUTORIGID:
    m_pDlgRegistration->SLT_ManualMoveByDCMPlan();
    m_pDlgRegistration->SLT_ConfirmManualRegistration(); // skin cropping
    // m_pDlgRegistration->SLT_DoRegistrationRigid();

    if (bCBCT_IntensityShift) {
      m_pDlgRegistration->SLT_IntensityNormCBCT();
    }

    if (bFOVCropping) {
      CropFOV3D(m_spManualRigidCT, physPosX, physPosY, physRadius,
                physTablePosY);
    }

    m_pDlgRegistration->SLT_DoRegistrationDeform();
    strSuffix = strSuffix + "_defrm_skipRigid";
    break;

  default:
    break;
  }

  if (bFullResolRecon) // if fullResol recon is on, load original proj files
                       // again
  {
    float curResampleF = m_fResampleF;
    ui.lineEdit_DownResolFactor->setText("1.0");
    SLT_LoadSelectedProjFiles();
    m_fResampleF = curResampleF;
    ui.lineEdit_DownResolFactor->setText(QString("%1").arg(m_fResampleF));

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
  QString outputPath_rawCBCT =
      outputDirPath + "/" + m_strDCMUID + strSuffix + "_rawCBCT.mha";
  QString outputPath_corrCBCT =
      outputDirPath + "/" + m_strDCMUID + strSuffix + "_corrCBCT.mha";
  QString outputPath_manCT =
      outputDirPath + "/" + m_strDCMUID + strSuffix + "_manCT.mha";
  QString outputPath_rigidCT =
      outputDirPath + "/" + m_strDCMUID + strSuffix + "_rigidCT.mha";
  QString outputPath_deformCT =
      outputDirPath + "/" + m_strDCMUID + strSuffix + "_deformCT.mha";

  if (bExportImages) {
    ExportReconSHORT_HU(m_spRawReconImg, outputPath_rawCBCT);
    ExportReconSHORT_HU(m_spScatCorrReconImg, outputPath_corrCBCT);
    ExportReconSHORT_HU(m_spManualRigidCT, outputPath_manCT);
    ExportReconSHORT_HU(m_spAutoRigidCT, outputPath_rigidCT);
    ExportReconSHORT_HU(m_spDeformedCT_Final, outputPath_deformCT);
  }

  // 2) Calculate batched WEPL points
  if (!m_vPOI_DCM.empty()) {
    QString outputTxtPath =
        outputDirPath + "/" + m_strDCMUID + strSuffix + "_WEPL.txt";
    ExportAngularWEPL_byFile(outputTxtPath);
  }
  return true;
}

bool CbctRecon::GetCouchShiftFromINIXVI(QString &strPathINIXVI, VEC3D *pTrans,
                                        VEC3D *pRot) {
  QFileInfo fInfo(strPathINIXVI);
  if (!fInfo.exists()) {
    return false;
  }

  std::ifstream fin;
  fin.open(strPathINIXVI.toLocal8Bit().constData());

  if (fin.fail()) {
    return false;
  }

  char str[MAX_LINE_LENGTH];

  float couch_Lat_cm = 0.0;
  float couch_Long_cm = 0.0;
  float couch_Vert_cm = 0.0;

  float couch_Pitch = 0.0;
  float couch_Yaw = 0.0;
  float couch_Roll = 0.0;

  bool bFound = false;
  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString tmpStr = QString(&str[0]);
    QStringList strListParam = tmpStr.split("=");

    QString tagName, strVal;

    if (strListParam.count() == 2) {
      tagName = strListParam.at(0);
      strVal = strListParam.at(1);
      tagName = tagName.trimmed();
      strVal = strVal.trimmed();

      if (tagName == "CouchShiftLat") {
        couch_Lat_cm = strVal.toFloat();
        bFound = true;
      } else if (tagName == "CouchShiftLong") {
        couch_Long_cm = strVal.toFloat();
      } else if (tagName == "CouchShiftHeight") {
        couch_Vert_cm = strVal.toFloat();
      } else if (tagName == "CouchPitch") {
        couch_Pitch = strVal.toFloat();
      } else if (tagName == "CouchRoll") {
        couch_Yaw = strVal.toFloat();
      } else if (tagName == "CouchYaw") {
        couch_Roll = strVal.toFloat();
      }
    }
  }
  fin.close();

  if (!bFound) {
    return false;
  }

  // Warning!! dicom convention!
  pTrans->x = couch_Lat_cm * 10.0; // sign should be checked
  // pTrans->y = couch_Vert_cm*10.0; //sign should be checked // IEC-->DICOM is
  // already accounted for..but sign!
  pTrans->y = couch_Vert_cm * (-10.0); // consistent with Tracking software
  pTrans->z = couch_Long_cm * 10.0;    // sign should be checked

  pRot->x = couch_Pitch;
  pRot->y = couch_Yaw;
  pRot->z = couch_Roll;
  // x,y,z: dicom
  return true;
}

bool CbctRecon::GetXrayParamFromINI(QString &strPathINI, float &kVp, float &mA,
                                    float &ms) {
  QFileInfo info = QFileInfo(strPathINI);

  kVp = 0.0;
  mA = 0.0;
  ms = 0.0;

  if (!info.exists()) {
    return false;
  }

  // TubeMA=64.0000
  // TubeKV = 120.0000
  // TubeKVLength = 40.0000
  std::ifstream fin;
  fin.open(strPathINI.toLocal8Bit().constData());

  if (fin.fail()) {
    return false;
  }

  char str[MAX_LINE_LENGTH];

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    QString tmpStr = QString(&str[0]);
    QStringList strListParam = tmpStr.split("=");

    QString tagName;
    QString strVal;

    if (strListParam.count() == 2) {
      tagName = strListParam.at(0);
      strVal = strListParam.at(1);
      tagName = tagName.trimmed();
      strVal = strVal.trimmed();

      if (tagName == "TubeMA") {
        mA = strVal.toFloat();
      }
      if (tagName == "TubeKVLength") {
        ms = strVal.toFloat();
      }
      if (tagName == "TubeKV") {
        kVp = strVal.toFloat();
      }
    }
  }
  fin.close();

  return !(kVp == 0 || mA == 0 || ms == 0);
}

void CbctRecon::GenerateCylinderMask(UShortImageType::Pointer &spImgCanvas,
                                     float fDcmPosX, float fDcmPosY,
                                     float fRadius) {
  if (spImgCanvas == nullptr) {
    return;
  }
  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position
  UShortImageType::PointType origin = spImgCanvas->GetOrigin();
  UShortImageType::SpacingType spacing = spImgCanvas->GetSpacing();
  // UShortImageType::SizeType size =
  // spImgCanvas->GetBufferedRegion().GetSize();

  // itk::ImageSliceConstIteratorWithIndex<FloatImageType> it (m_spReconImg,
  // m_spReconImg->GetRequestedRegion());
  itk::ImageSliceIteratorWithIndex<UShortImageType> it(
      spImgCanvas, spImgCanvas->GetRequestedRegion());

  // ImageSliceConstIteratorWithIndex<ImageType> it( image,
  // image->GetRequestedRegion() );  UShortImageType::SizeType imgSize =
  // spImgCanvas->GetRequestedRegion().GetSize(); //1016x1016 x z

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

        if (pow(crntPhysX - fDcmPosX, 2.0) + pow(crntPhysY - fDcmPosY, 2.0) >=
            pow(fRadius, 2.0)) {
          //(*it) = (unsigned short)0; //air value
          it.Set(0);
        } else {
          it.Set(1);
        }

        ++it;
        iPosX++;
      }
      it.NextLine();
      iPosY++;
    }
    it.NextSlice();
    iNumSlice++;
  }
}

float CbctRecon::GetMeanIntensity(UShortImageType::Pointer &spImg,
                                  float sphereR, float *sdIntensity) {
  if (spImg == nullptr) {
    return -1.0;
  }

  float meanIntensity = 0.0;

  // 1) region iterator, set 0 for all pixels outside the circle and below the
  // table top, based on physical position
  UShortImageType::PointType origin = spImg->GetOrigin();
  UShortImageType::SpacingType spacing = spImg->GetSpacing();
  // UShortImageType::SizeType size = spImg->GetBufferedRegion().GetSize();

  itk::ImageSliceIteratorWithIndex<UShortImageType> it(
      spImg, spImg->GetRequestedRegion());
  // UShortImageType::SizeType imgSize = spImg->GetRequestedRegion().GetSize();
  // //1016x1016 x z

  // int width = imgSize[0];
  // int height = imgSize[1];

  it.SetFirstDirection(0);  // x?
  it.SetSecondDirection(1); // y?
  it.GoToBegin();

  int iNumSlice = 0;
  int iPosX = 0;
  int iPosY = 0;

  double crntPhysX = 0.0;
  double crntPhysY = 0.0;
  double crntPhysZ = 0.0;

  double pixSum = 0.0;
  int iCnt = 0;

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
        crntPhysZ = iNumSlice * static_cast<double>(spacing[2]) +
                    static_cast<double>(origin[2]);

        if (pow(crntPhysX, 2.0) + pow(crntPhysY, 2.0) + pow(crntPhysZ, 2.0) <
            pow(sphereR, 2.0)) {
          pixSum = pixSum + static_cast<double>(it.Get());
          iCnt++;
        }
        ++it;
        iPosX++;
      }
      it.NextLine();
      iPosY++;
    }
    it.NextSlice();
    iNumSlice++;
  }

  if (iCnt > 0) {
    meanIntensity = pixSum / static_cast<double>(iCnt);
  } else {
    meanIntensity = -1.0;
  }

  if (sdIntensity == nullptr) {
    return meanIntensity;
  }

  double devSum = 0.0;
  it.GoToBegin();

  iNumSlice = 0;

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
        crntPhysZ = iNumSlice * static_cast<double>(spacing[2]) +
                    static_cast<double>(origin[2]);

        if (pow(crntPhysX, 2.0) + pow(crntPhysY, 2.0) + pow(crntPhysZ, 2.0) <
            pow(sphereR, 2.0)) {
          devSum = devSum +
                   pow((static_cast<double>(it.Get()) - meanIntensity), 2.0);
        }
        ++it;
        iPosX++;
      }
      it.NextLine();
      iPosY++;
    }
    it.NextSlice();
    iNumSlice++;
  }

  if (iCnt > 0) {
    *sdIntensity = sqrt(devSum / static_cast<double>(iCnt));
  } else {
    *sdIntensity = -1.0;
  }

  return meanIntensity;
}

void CbctRecon::AddConstHU(UShortImageType::Pointer &spImg, int HUval) {

  using iteratorType = itk::ImageRegionIteratorWithIndex<UShortImageType>;
  iteratorType it(spImg, spImg->GetRequestedRegion());

  it.GoToBegin();

  int crntVal = 0;
  int newVal = 0;

  while (!it.IsAtEnd()) {
    crntVal = static_cast<int>(it.Get());

    newVal = HUval + crntVal;

    if (newVal <= 0) {
      newVal = 0;
    }

    if (newVal >= 4095) {
      newVal = 4095;
    }

    it.Set(static_cast<unsigned short>(newVal));
    ++it;
  }
}

void CbctRecon::SLT_OpenPhaseData() {
  if (!m_vPhaseFloat.empty()) {
    m_vPhaseFloat.clear();
  }

  // Open file
  QString filePath =
      QFileDialog::getOpenFileName(this, "Open phase text", m_strPathDirDefault,
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
    m_vPhaseFloat.push_back(tmpPhase);
    phaseCnt++;
    phaseSum = phaseSum + tmpPhase;
  }
  fin.close();
  std::cout << "NumOfPhaseData[Float]= " << phaseCnt << "  Mean Phase value= "
            << phaseSum / static_cast<double>(phaseCnt) << std::endl;
}

void CbctRecon::SLT_Export4DCBCT() {
  if (m_spCustomGeometry == nullptr) {
    std::cout << "Error! no Geometry information loaded yet" << std::endl;
    return;
  }

  int NumOfGanAngle = m_spCustomGeometry->GetGantryAngles().size();
  int NumOfPhase = m_vPhaseFloat.size();

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

  QString strDirForXML = m_strPathDirDefault; // where xml file is located
  // QString strUID ;//P00102030P + m_strDCMUID
  QString strDirForProj = m_strPathIMAGES;

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
    if (!ResortCBCTProjection(vPhaseBinsSelected, strDirForXML, strDirForProj,
                              m_strDCMUID, m_vPhaseFloat, m_spCustomGeometry,
                              m_vSelectedFileNames)) {
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

bool CbctRecon::ResortCBCTProjection(std::vector<int> &vIntPhaseBinSelected,
                                     QString &strPathForXML,
                                     QString &strPathProjRoot, QString &strUID,
                                     std::vector<float> &vFloatPhaseFull,
                                     GeometryType::Pointer &spGeomFull,
                                     std::vector<std::string> &vProjPathsFull) {
  if (vIntPhaseBinSelected.empty()) {
    return false;
  }

  int NumOfPhaseFull = vFloatPhaseFull.size();
  int NumOfGeomFull = spGeomFull->GetGantryAngles().size();
  int NumOfProjFileFull = vProjPathsFull.size();

  if (NumOfPhaseFull != NumOfGeomFull || NumOfGeomFull != NumOfProjFileFull) {
    std::cout << "Num of data is not matching:"
              << " NumOfPhaseFull= " << NumOfPhaseFull
              << " NumOfGeomFull= " << NumOfGeomFull
              << " NumOfProjFileFull= " << NumOfProjFileFull << std::endl;
    return false;
  }

  // Check Root dir is set

  if (strUID.length() < 1) {
    return false;
  }

  QDir dirSaveXML(strPathForXML);
  QDir dirSaveProj(strPathProjRoot);

  if (!dirSaveXML.exists() || !dirSaveProj.exists()) {
    std::cout << "Error! Directories don't exist" << std::endl;
    return false;
  }
  // Generate a new UID
  QString strUID_Endfix;
  int iNumOfSelPhase = vIntPhaseBinSelected.size();

  strUID_Endfix = "P";
  QChar zero('0');
  for (int i = 0; i < iNumOfSelPhase; i++) {
    QString strNum;
    strNum = QString("%1").arg(vIntPhaseBinSelected.at(i), 2, 10, zero);
    strUID_Endfix = strUID_Endfix + strNum;
  }
  strUID_Endfix = strUID_Endfix + "P"; // UID...P00102030405060P
  QString strNewUID = strUID + strUID_Endfix;

  // Create a subDir
  QDir curProjRoot(strPathProjRoot);
  QString strSubDirName = "img_" + strNewUID;
  curProjRoot.mkdir(strSubDirName);

  QString strPathProj = strPathProjRoot + "/" + strSubDirName;

  QDir projDir(strPathProj);
  if (!projDir.exists()) {
    std::cout << "no Proj Dir exists" << std::endl;
    return false;
  }

  QDir xmlDir(strPathForXML);
  if (!xmlDir.exists()) {
    std::cout << "no XML Dir exists" << std::endl;
    return false;
  }

  // strPathProj
  // strPathForXML
  std::vector<int> vSelectedIdxTemp;
  std::vector<int> vSelectedIdxFin;

  for (int i = 0; i < iNumOfSelPhase; i++) {
    AppendInPhaseIndex(vIntPhaseBinSelected.at(i), vFloatPhaseFull,
                       vSelectedIdxTemp);
  }
  // Remove redandancy

  sort(vSelectedIdxTemp.begin(), vSelectedIdxTemp.end()); // hopefully,
                                                          // ascending
  std::cout << "sorting check" << std::endl;
  std::cout << "0 " << vSelectedIdxTemp.at(0) << std::endl;
  std::cout << "1 " << vSelectedIdxTemp.at(1) << std::endl;

  std::vector<int>::iterator it;

  int prevVal = -1;
  for (it = vSelectedIdxTemp.begin(); it != vSelectedIdxTemp.end(); ++it) {
    if ((*it) > prevVal) {
      vSelectedIdxFin.push_back(*it);
    }

    prevVal = (*it);
  }

  GeometryType::Pointer spSubGeometry = GeometryType::New();

  std::vector<int>::iterator itIdx;

  for (itIdx = vSelectedIdxFin.begin(); itIdx != vSelectedIdxFin.end();
       itIdx++) {
    std::cout << "cur Idx=" << (*itIdx) << std::endl;
    // 9 parameters are required
    double curSID = spGeomFull->GetSourceToIsocenterDistances().at(*itIdx);
    double curSDD = spGeomFull->GetSourceToDetectorDistances().at(*itIdx);
    double curGantryAngle = spGeomFull->GetGantryAngles().at(*itIdx);

    double curProjOffsetX = spGeomFull->GetProjectionOffsetsX().at(*itIdx);
    double curProjOffsetY = spGeomFull->GetProjectionOffsetsY().at(*itIdx);

    double curOutOfPlaneAngles = spGeomFull->GetOutOfPlaneAngles().at(*itIdx);
    double curInPlaneAngles = spGeomFull->GetInPlaneAngles().at(*itIdx);

    double curSrcOffsetX = spGeomFull->GetSourceOffsetsX().at(*itIdx);
    double curSrcOffsetY = spGeomFull->GetSourceOffsetsY().at(*itIdx);

    spSubGeometry->AddProjection(
        curSID, curSDD, curGantryAngle, curProjOffsetX,
        curProjOffsetY,                        // Flexmap
        curOutOfPlaneAngles, curInPlaneAngles, // In elekta, these are 0
        curSrcOffsetX, curSrcOffsetY);         // In elekta, these are 0
  }
  // Export spSubGeometry
  rtk::ThreeDCircularProjectionGeometryXMLFileWriter::Pointer xmlWriter =
      rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();

  QString geomFileName = "ElektaGeom_" + strNewUID + ".xml";
  QString geomFilePath = strPathForXML + "/" + geomFileName;

  xmlWriter->SetFilename(geomFilePath.toLocal8Bit().constData());
  xmlWriter->SetObject(spSubGeometry);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile());
  // Copy selected his files to a different folder

  for (itIdx = vSelectedIdxFin.begin(); itIdx != vSelectedIdxFin.end();
       itIdx++) {
    QString strPathProjOriginal = vProjPathsFull.at(*itIdx).c_str();
    // Copy this file to target dir

    QFileInfo fInfo(strPathProjOriginal);
    QString strPathProjNew = strPathProj + "/" + fInfo.fileName();
    QFile::copy(fInfo.absoluteFilePath(), strPathProjNew);
  }

  //    std::vector<float>& vFloatPhaseFull, GeometryType::Pointer& spGeomFull,
  //    std::vector<std::string>& vProjPathsFull
  std::cout << vSelectedIdxFin.size() << " files were copied." << std::endl;

  return true;
}

void CbctRecon::AppendInPhaseIndex(int iPhase,
                                   std::vector<float> &vFloatPhaseFull,
                                   std::vector<int> &vOutputIndex, int margin) {

  int iNumOfPhase = vFloatPhaseFull.size();

  int iCurPhase = 0;

  int startPhase1;
  int endPhase1;

  int startPhase2;
  int endPhase2;

  for (int i = 0; i < iNumOfPhase; i++) {
    iCurPhase = qRound(vFloatPhaseFull.at(i) * 100.0);
    // determine wether it is within the range

    if (iPhase < margin) // if 5 --> 0 ~ 10%, IF 4--> 99 ~ 09
    {
      startPhase2 = iPhase + 100 - margin;
      endPhase2 = 100;

      startPhase1 = 0;
      endPhase1 = iPhase + margin;
    } else {
      startPhase1 = iPhase - margin;
      endPhase1 = iPhase + margin;

      startPhase2 = 1; // reverse
      endPhase2 = 0;
    }

    if ((iCurPhase >= startPhase1 && iCurPhase <= endPhase1) ||
        (iCurPhase >= startPhase2 && iCurPhase <= endPhase2)) {
      vOutputIndex.push_back(i);
    }
  }
}

void CbctRecon::SLT_LoadCBCTcorrMHA() {
  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_strPathDirDefault, "Short image file (*.mha)",
      nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  LoadShort3DImage(fileName, REGISTER_COR_CBCT);

  // std::cout << m_spScatCorrReconImg->GetBufferedRegion().GetSize() <<
  // std::endl;

  SLT_DrawReconImage();
}

void CbctRecon::SLT_LoadCTrigidMHA() {
  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_strPathDirDefault, "Short image file (*.mha)",
      nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  LoadShort3DImage(fileName, REGISTER_AUTO_RIGID);

  SLT_DrawReconImage();
}

void CbctRecon::SLT_LoadCTdeformMHA() {
  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_strPathDirDefault, "Short image file (*.mha)",
      nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  LoadShort3DImage(fileName, REGISTER_DEFORM_FINAL);

  SLT_DrawReconImage();
}

void CbctRecon::LoadShort3DImage(QString &filePath, enREGI_IMAGES enTarget) {
  QFileInfo fInfo(filePath);
  if (!fInfo.exists()) {
    return;
  }

  UShortImageType::Pointer spImg;

  if (!LoadShortImageToUshort(filePath, spImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  switch (enTarget) {
  case REGISTER_RAW_CBCT:
    m_spRawReconImg = spImg;
    break;
  case REGISTER_COR_CBCT:
    m_spScatCorrReconImg = spImg;
    break;
  case REGISTER_MANUAL_RIGID:
    m_spManualRigidCT = spImg;
    break;
  case REGISTER_AUTO_RIGID:
    m_spAutoRigidCT = spImg;
    break;
  case REGISTER_DEFORM_FINAL:
    m_spDeformedCT_Final = spImg;
    break;
  default:
    m_spRawReconImg = spImg;
    break;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  ImageCalculatorFilterType2::Pointer imageCalculatorFilter2 =
      ImageCalculatorFilterType2::New();
  // imageCalculatorFilter2->SetImage(m_spReconImg);
  imageCalculatorFilter2->SetImage(spImg);
  imageCalculatorFilter2->Compute();

  auto minVal2 = static_cast<double>(imageCalculatorFilter2->GetMinimum());
  auto maxVal2 = static_cast<double>(imageCalculatorFilter2->GetMaximum());

  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  UShortImageType::SizeType imgDim = spImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing = spImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  m_spCrntReconImg = spImg;

  ui.lineEdit_Cur3DFileName->setText(filePath);
  m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(initVal); // DrawRecon Imge is called
  ui.radioButton_graph_recon->setChecked(true);

  m_pDlgRegistration->UpdateListOfComboBox(0); // combo selection signalis
                                               // called
  m_pDlgRegistration->UpdateListOfComboBox(1);
  m_pDlgRegistration->SelectComboExternal(
      0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  m_pDlgRegistration->SelectComboExternal(1, enTarget);
}

// trans: mm, dicom order
// COuch shift values: directlry come from the INI.XVI file only multiplied
// by 10.0
void CbctRecon::ImageTransformUsingCouchCorrection(
    UShortImageType::Pointer &spUshortInput,
    UShortImageType::Pointer &spUshortOutput, VEC3D couch_trans,
    VEC3D couch_rot) {
  // couch_trans, couch_rot--> as it is from the text file. only x 10.0 was
  // applied
  if (spUshortInput == nullptr) {
    return;
  }

  using FilterType = itk::ResampleImageFilter<UShortImageType, UShortImageType>;
  FilterType::Pointer filter = FilterType::New();

  using TransformType = itk::AffineTransform<double, 3>;
  TransformType::Pointer transform = TransformType::New();
  filter->SetTransform(transform);
  using InterpolatorType =
      itk::NearestNeighborInterpolateImageFunction<UShortImageType, double>;

  InterpolatorType::Pointer interpolator = InterpolatorType::New();
  filter->SetInterpolator(interpolator);

  filter->SetDefaultPixelValue(0);

  //  const double spacing[3] = { 1.0, 1.0, 1.0 };
  UShortImageType::SpacingType spacing = spUshortInput->GetSpacing();

  filter->SetOutputSpacing(spacing);

  UShortImageType::PointType origin = spUshortInput->GetOrigin();

  filter->SetOutputOrigin(origin);

  UShortImageType::DirectionType direction;
  direction.SetIdentity();
  filter->SetOutputDirection(direction);

  UShortImageType::SizeType size =
      spUshortInput->GetLargestPossibleRegion().GetSize();
  filter->SetSize(size);
  filter->SetInput(spUshortInput);

  // NOTE: In couch shift reading
  // pTrans->x = couch_Lat_cm*10.0; //sign should be checked
  // pTrans->y = couch_Vert_cm*10.0; //sign should be checked // IEC-->DICOM is
  // already accounted for..but sign!  pTrans->z = couch_Long_cm*10.0; //sign
  // should be checked  pRot->x = couch_Pitch;  pRot->y = couch_Yaw;  pRot->z =
  // couch_Roll;

  TransformType::OutputVectorType translation;
  translation[0] = -couch_trans.x; // X translation in millimeters
  // translation[1] = +couch_trans.y; //so far so good// This is because when
  // IEC->DICOM, sign was not changed during reading the text file
  translation[1] = -couch_trans.y; // Consistent with Tracking software
  translation[2] = -couch_trans.z; // empirically found

  TransformType::OutputVectorType rotation;
  rotation[0] = -couch_rot.x; // X translation in millimeters
  rotation[1] = -couch_rot.y;
  rotation[2] = -couch_rot.z;

  transform->Translate(
      translation); // original position - (couch shift value in DICOM)
  // transform->Rotate3D(rotation);
  filter->Update();

  spUshortOutput = filter->GetOutput();
  // std::cout << "affine transform is successfully done" << std::endl;
}

void CbctRecon::SLT_DoCouchCorrection() {
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

  ImageTransformUsingCouchCorrection(m_spRawReconImg, m_spRawReconImg,
                                     couchShiftTrans, couchShiftRot);
  ImageTransformUsingCouchCorrection(m_spScatCorrReconImg, m_spScatCorrReconImg,
                                     couchShiftTrans, couchShiftRot);
  ImageTransformUsingCouchCorrection(m_spDeformedCT_Final, m_spDeformedCT_Final,
                                     couchShiftTrans, couchShiftRot);
  ImageTransformUsingCouchCorrection(m_spAutoRigidCT, m_spAutoRigidCT,
                                     couchShiftTrans, couchShiftRot);

  m_pDlgRegistration->UpdateListOfComboBox(0); // combo selection signalis
                                               // called
  m_pDlgRegistration->UpdateListOfComboBox(1);
  m_pDlgRegistration->SelectComboExternal(
      0, REGISTER_RAW_CBCT); // will call fixedImageSelected
  m_pDlgRegistration->SelectComboExternal(1, REGISTER_COR_CBCT);

  m_spCrntReconImg = m_spScatCorrReconImg;
  SLT_DrawReconImage();

  std::cout << "Couch shift and rotation was successfully applied."
            << std::endl;
}

// Multiple mha files
void CbctRecon::SLTM_WELPCalcMultipleFiles() {
  // Singed short
  QStringList listFilePath = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open", m_strPathDirDefault,
      "signed short 3D images (*.mha)");

  int iCntFiles = listFilePath.count();
  if (iCntFiles < 1) {
    return;
  }

  int iCntPOI = m_vPOI_DCM.size();

  if (iCntPOI < 1) {
    std::cout << "There is no POI file loaded." << std::endl;
    SLT_LoadPOIData();
  }
  iCntPOI = m_vPOI_DCM.size();
  if (iCntPOI < 1) {
    std::cout << "Error! still no POI" << std::endl;
    return;
  }

  QString strPathOutText = QFileDialog::getSaveFileName(
      this, "File path to save", m_strPathDirDefault, "WEPL_value (*.txt)",
      nullptr, nullptr); // Filename don't need to exist
  if (strPathOutText.length() <= 1) {
    return;
  }

  std::vector<std::vector<WEPLData>> vArrOutputWEPL;
  vArrOutputWEPL.resize(iCntFiles);

  for (int i = 0; i < iCntFiles; i++) {
    GetWEPLDataFromSingleFile(listFilePath.at(i), m_vPOI_DCM,
                              vArrOutputWEPL.at(i));
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

void CbctRecon::GetWEPLDataFromSingleFile(const QString &filePath,
                                          std::vector<VEC3D> &vPOI,
                                          std::vector<WEPLData> &vOutputWEPL) {

  int iCntPOI = vPOI.size();

  if (iCntPOI < 1) {
    return;
  }

  float fAngleGap = 1.0;
  double fAngleStart = ui.lineEdit_AngStart->text().toDouble();
  double fAngleEnd = ui.lineEdit_AngEnd->text().toDouble();

  UShortImageType::Pointer spImg;

  QString strFilePath = filePath;
  if (!LoadShortImageToUshort(strFilePath, spImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
    return;
  }

  for (int i = 0; i < iCntPOI; i++) {
    VEC3D curPOI = vPOI.at(i);
    // append mode
    GetAngularWEPL_SinglePoint(spImg, fAngleGap, fAngleStart, fAngleEnd, curPOI,
                               i, vOutputWEPL, true); // mandatory
  }
}

void CbctRecon::SLTM_ScatterCorPerProjRef() // load text file
{
  if (m_strListPerProjRefVol.empty()) {
    std::cout << "Error! Ref Vol list is not ready yet. Load it first"
              << std::endl;
    return;
  }

  // Find the mask file
  // QString strPath_mskSkinCT_final;
  // QString strPath_mskSkinCT_autoRegi_exp =
  // m_pDlgRegistration->m_strPathPlastimatch + "/msk_skin_CT_autoRegi_exp.mha";
  // QFileInfo maskInfoAuto(strPath_mskSkinCT_autoRegi_exp);

  // QString strPath_mskSkinCT_manualRegi_exp =
  // m_pDlgRegistration->m_strPathPlastimatch + "/msk_skin_CT_manRegi_exp.mha";
  // QFileInfo maskInfoManual(strPath_mskSkinCT_manualRegi_exp);

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

  ////std::cout << "Plastimatch Path " <<
  /// m_strPathPlastimatch.toLocal8Bit().constData() << std::endl;

  // if (m_pDlgRegistration->m_strPathPlastimatch.length() < 1)
  //{
  //    std::cout << "NO plastimatch Dir was defined. CorrCBCT will not be saved
  //    automatically" << std::endl; return;
  //}
  // Forward proj

  // Make a canvas
  m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize();

  if (m_spProjImgRaw3D == nullptr) {
    std::cout << "ERRORRR! m_spProjImgRaw3D" << std::endl;
    return;
  }

  m_spProjImgCT3D = UShortImageType::New(); // later
  UShortImageType::SizeType projCT_size =
      m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize(); // 1024 1024 350
  UShortImageType::IndexType projCT_idxStart =
      m_spProjImgRaw3D->GetLargestPossibleRegion().GetIndex(); // 0 0 0
  UShortImageType::SpacingType projCT_spacing =
      m_spProjImgRaw3D->GetSpacing(); // 0.4 0.4 1.0
  UShortImageType::PointType projCT_origin =
      m_spProjImgRaw3D->GetOrigin(); //-204.6 -204.6 -174.5

  UShortImageType::RegionType projCT_region;
  projCT_region.SetSize(projCT_size);
  projCT_region.SetIndex(projCT_idxStart);

  m_spProjImgCT3D->SetRegions(projCT_region);
  m_spProjImgCT3D->SetSpacing(projCT_spacing);
  m_spProjImgCT3D->SetOrigin(projCT_origin);

  m_spProjImgCT3D->Allocate();
  m_spProjImgCT3D->FillBuffer(0);

  // YKTEMP
  std::cout << "ProjImgCT Size = "
            << m_spProjImgCT3D->GetBufferedRegion().GetSize()[0] << ", "
            << m_spProjImgCT3D->GetBufferedRegion().GetSize()[1] << ", "
            << m_spProjImgCT3D->GetBufferedRegion().GetSize()[2] << std::endl;
  std::cout << "ProjImgCT origin = " << m_spProjImgCT3D->GetOrigin()[0] << ", "
            << m_spProjImgCT3D->GetOrigin()[1] << ", "
            << m_spProjImgCT3D->GetOrigin()[2] << std::endl;
  std::cout << "ProjImgCT spacing = " << m_spProjImgCT3D->GetSpacing()[0]
            << ", " << m_spProjImgCT3D->GetSpacing()[1] << ", "
            << m_spProjImgCT3D->GetSpacing()[2] << std::endl;

  int iCntRefVol = m_strListPerProjRefVol.count();

  if (iCntRefVol < 1) {
    std::cout << "Error! no volume data for loading" << std::endl;
    return;
  }

  auto flexCnt = static_cast<int>(m_spCustomGeometry->GetGantryAngles().size());
  if (flexCnt != iCntRefVol) {
    std::cout << "Error! flex count doesn't match" << std::endl;
    return;
  }

  double curMVAngle = 0.0;
  double curPanelOffsetX = 0.0;
  double curPanelOffsetY = 0.0;

  for (int i = 0; i < iCntRefVol; i++) {
    // Load volume: Short image
    ShortImageType::Pointer spOutputShort_raw = ShortImageType::New();
    // ShortImageType::Pointer spOutputShort_threshold = ShortImageType::New();
    UShortImageType::Pointer spOutputUshort = UShortImageType::New();
    // UShortImageType::Pointer spOutputUshort_register =
    // UShortImageType::New();
    UShortImageType::Pointer spUshortRotated = UShortImageType::New();
    FloatImageType::Pointer spAttFloat = FloatImageType::New();

    QString strDirPath = m_strListPerProjRefVol.at(i);

    if (!LoadShortImageDirOrFile(strDirPath, spOutputShort_raw)) {
      std::cout << "Error! in " << i
                << " th image. File couldn't be found. Path= "
                << strDirPath.toLocal8Bit().constData() << std::endl;
      return;
    }

    ConvertShort2Ushort(spOutputShort_raw, spOutputUshort);

    RotateImgBeforeFwd(spOutputUshort,
                       spUshortRotated); // IEC to RTK w/ kVGantry
    ConvertUshort2AttFloat(spUshortRotated, spAttFloat);

    curMVAngle = m_spCustomGeometry->GetGantryAngles().at(i);
    curPanelOffsetX = m_spCustomGeometry->GetProjectionOffsetsX().at(i);
    curPanelOffsetY = m_spCustomGeometry->GetProjectionOffsetsY().at(i);

    SingleForwardProjection(spAttFloat, curMVAngle, curPanelOffsetX,
                            curPanelOffsetY, m_spProjImgCT3D, i);
    std::cout << "Proj: " << i << "/" << iCntRefVol << std::endl;
  }

  /* typedef itk::ImageFileWriter<UShortImageType> WriterType;
   WriterType::Pointer writer = WriterType::New();
   writer->SetFileName("D:/TmpProjCT3D.mha");
   writer->SetUseCompression(true);
   writer->SetInput(m_spProjImgCT3D);
   writer->Update();
   */
  double scaMedian = ui.lineEdit_scaMedian->text().toDouble();
  double scaGaussian = ui.lineEdit_scaGaussian->text().toDouble();

  std::cout << "Generating scatter map is ongoing..." << std::endl;

  std::cout << "To account for the mAs values, the intensity scale factor of "
            << GetRawIntensityScaleFactor()
            << "will be multiplied during scatter correction to avoid negative "
               "scatter"
            << std::endl;

  GenScatterMap_PriorCT(m_spProjImgRaw3D, m_spProjImgCT3D, m_spProjImgScat3D,
                        scaMedian, scaGaussian, m_iFixedOffset_ScatterMap,
                        false);  // void GenScatterMap2D_PriorCT()
  m_spProjImgCT3D->Initialize(); // memory saving

  std::cout << "Scatter correction is in progress..." << std::endl;

  int postScatMedianSize = ui.lineEdit_scaPostMedian->text().toInt();
  ScatterCorr_PrioriCT(m_spProjImgRaw3D, m_spProjImgScat3D, m_spProjImgCorr3D,
                       m_iFixedOffset_ScatterMap, postScatMedianSize, true);
  m_spProjImgScat3D->Initialize(); // memory saving

  std::cout << "AfterCorrectionMacro is ongoing..." << std::endl;
  AfterScatCorrectionMacro();
  std::cout << "FINISHED!Scatter correction: CBCT DICOM files are saved"
            << std::endl;

  ////1) Export current CBCT file
  // QString filePathCBCT = m_strPathPlastimatch + "/" + "CorrCBCT.mha";
  // //usually corrected one  QString filePathCBCT_noSkin = m_strPathPlastimatch
  // +
  // "/" + "CorrCBCT_final.mha"; //usually corrected one

  // typedef itk::ImageFileWriter<UShortImageType> writerType;
  // writerType::Pointer writer = writerType::New();
  // writer->SetFileName(filePathCBCT.toLocal8Bit().constData());
  // writer->SetUseCompression(true);
  // writer->SetInput(spCBCT);

  // std::cout << "Writing the CBCT file" << std::endl;
  // writer->Update();

  // QFileInfo CBCTInfo(filePathCBCT);
  // if (!CBCTInfo.exists())
  //{
  //    std::cout << "No CBCT file to read. Maybe prior writing failed" <<
  //    std::endl; return;
  //}

  ////ERROR HERE! delete the temporry folder.
  // std::cout << "Delete the temporary folder if it crashes" << std::endl;

  ////4) eliminate the air region (temporarily)
  ////Mask_parms parms_msk;
  ////DIMENSION SHOULD BE MATCHED!!!! BETWEEN raw CBCT and Mask files
  // Mask_operation mask_option = MASK_OPERATION_MASK;
  // QString input_fn = filePathCBCT.toLocal8Bit().constData();
  // QString mask_fn = strPath_mskSkinCT_final.toLocal8Bit().constData();
  // QString output_fn = filePathCBCT_noSkin.toLocal8Bit().constData();
  // float mask_value = 0.0; //unsigned short
  // plm_mask_main(mask_option, input_fn, mask_fn, output_fn, mask_value);
}

// refer to YKPRoc later
// void YKPROC::ForwardProjection(FloatImageType::Pointer& spVolImgFloat, float
// fMVGanAngle, float panelOffsetX, float panelOffsetY, ,
// UShortImageType::Pointer& spProj3D)  iSliceIdx == Proj index 0 - 364
void CbctRecon::SingleForwardProjection(FloatImageType::Pointer &spVolImgFloat,
                                        float fMVGanAngle, float panelOffsetX,
                                        float panelOffsetY,
                                        UShortImageType::Pointer &spProjImg3D,
                                        int iSliceIdx) {
  if (spVolImgFloat == nullptr) {
    return;
  }
  if (spProjImg3D == nullptr) {
    return;
  }

  // 2) Prepare empty projection images //Should be corresonponding to raw
  // projection images

  // Create a stack of empty projection images
  using ConstantImageSourceType =
      rtk::ConstantImageSource<FloatImageType>; // Output: FLoat image = may be
                                                // mu_t = log(I_0/I)
  ConstantImageSourceType::Pointer constantImageSource =
      ConstantImageSourceType::New();

  ConstantImageSourceType::SizeType size{};
  ConstantImageSourceType::SpacingType spacing;
  ConstantImageSourceType::PointType origin;

  // std::cout << "Setting-up vacant projection image data" << std::endl;

  // a) size
  // std::cout << "chk1" << std::endl;
  size[0] = spProjImg3D->GetBufferedRegion().GetSize()[0];
  size[1] = spProjImg3D->GetBufferedRegion().GetSize()[1];
  size[2] = 1;

  int totalProjSize = spProjImg3D->GetBufferedRegion().GetSize()[2];
  if (iSliceIdx >= totalProjSize) {
    std::cout << "Error! totalProjSize= " << totalProjSize
              << " iSliceIdx= " << iSliceIdx << std::endl;
  }

  // b) spacing
  spacing[0] = spProjImg3D->GetSpacing()[0];
  spacing[1] = spProjImg3D->GetSpacing()[1];
  spacing[2] = 1.0;

  // c) Origin: can center be the image center? or should be related to the CT
  // image???
  /*origin[0] = spacing[0] * (size[0] - 1) * -0.5;
  origin[1] = spacing[1] * (size[1] - 1) * -0.5;
  origin[2] = 0.0;*/

  origin[0] = spProjImg3D->GetOrigin()[0];
  origin[1] = spProjImg3D->GetOrigin()[1];
  origin[2] = 0.0;

  constantImageSource->SetOrigin(origin);
  constantImageSource->SetSpacing(spacing);

  FloatImageType::DirectionType imageDirection;
  imageDirection.SetIdentity(); // no effect
  constantImageSource->SetDirection(imageDirection);
  constantImageSource->SetSize(size);
  constantImageSource->SetConstant(1.0);
  constantImageSource->UpdateOutputInformation();

  // std::cout << "Canvas for projection image is ready to write" << std::endl;

  // 4) Prepare CT image to be projected
  int fwdMethod = en_CudaRayCast; // later, it will be coming from the GUI
  //    std::cout << "projection algorithm (0:Joseph, 1: CUDA, 2:RayCast ): " <<
  //    fwdMethod << std::endl;

  // Create forward projection image filter
#if USE_CUDA
  rtk::CudaForwardProjectionImageFilter<CUDAFloatImageType,
                                        CUDAFloatImageType>::Pointer
      CudaForwardProjection; // Float to Float
#endif                       // CUDA_FOUND

  rtk::ForwardProjectionImageFilter<FloatImageType, FloatImageType>::Pointer
      forwardProjection; // Float to Float

  switch (fwdMethod) {
  case (en_Joseph):
    forwardProjection =
        rtk::JosephForwardProjectionImageFilter<FloatImageType,
                                                FloatImageType>::New();
    break;
  case (en_CudaRayCast):
#if USE_CUDA
    CudaForwardProjection =
        rtk::CudaForwardProjectionImageFilter<CUDAFloatImageType,
                                              CUDAFloatImageType>::New();
#else
    std::cerr << "The program has not been compiled with cuda option"
              << std::endl;
    return;
#endif
    break;
//#if RTK_MINOR_VERSION < 4
//  case (en_RayCastInterpolator):
//    forwardProjection = rtk::RayCastInterpolatorForwardProjectionImageFilter<
//        FloatImageType, FloatImageType>::New();
//    break;
//#endif

  default:
    std::cerr << "Unhandled --method value." << std::endl;
    return;
  }

  GeometryType::Pointer spGeometry = GeometryType::New();

  // 9 parameters are required
  double curSAD = 1000.0; // SourceToIsocenterDistances
  double curSDD = 1536.0;
  double curGantryAngle = fMVGanAngle; // MV

  double curProjOffsetX = panelOffsetX;
  double curProjOffsetY = panelOffsetY;

  double curOutOfPlaneAngles = 0.0;
  double curInPlaneAngles = 0.0;

  double curSrcOffsetX = 0.0;
  double curSrcOffsetY = 0.0;

  spGeometry->AddProjection(
      curSAD, curSDD, curGantryAngle, curProjOffsetX, curProjOffsetY, // Flexmap
      curOutOfPlaneAngles, curInPlaneAngles, // In elekta, these are 0
      curSrcOffsetX, curSrcOffsetY);         // In elekta, these are 0

  itk::TimeProbe projProbe;
  // std::cout << "Forward projection is now ongoing" << std::endl;
  FloatImageType::Pointer resultFwdImg;
  if (fwdMethod == en_CudaRayCast) {
    forwardProjection->SetInput(
        constantImageSource
            ->GetOutput()); // Canvas. projection image will be saved here.
    forwardProjection->SetInput(1, spVolImgFloat); // reference plan CT image
    forwardProjection->SetGeometry(spGeometry);

    projProbe.Start();
    TRY_AND_EXIT_ON_ITK_EXCEPTION(forwardProjection->Update())
    projProbe.Stop();

    resultFwdImg = forwardProjection->GetOutput();
  } else {
    forwardProjection->SetInput(
        constantImageSource
            ->GetOutput()); // Canvas. projection image will be saved here.
    forwardProjection->SetInput(1, spVolImgFloat); // reference plan CT image
    forwardProjection->SetGeometry(spGeometry);

    projProbe.Start();
    TRY_AND_EXIT_ON_ITK_EXCEPTION(forwardProjection->Update())
    projProbe.Stop();

    resultFwdImg = forwardProjection->GetOutput();
  }
  std::cout << "Forward projection done by in method ID = " << fwdMethod
            << " in:	" << projProbe.GetMean() << ' ' << projProbe.GetUnit()
            << '.' << std::endl;

  // normalization or shift

  // typedef itk::MinimumMaximumImageCalculator<FloatImageType>
  // MinMaxCalculatorType;  MinMaxCalculatorType::Pointer spCalculator =
  // MinMaxCalculatorType::New();  spCalculator->SetImage(resultFwdImg);
  // spCalculator->Compute();

  // float minValAtt = spCalculator->GetMinimum();
  // float maxValAtt = spCalculator->GetMaximum();

  // float maxValProj = (65535.0 / exp(minValAtt));
  // float minValProj = (65535.0 / exp(maxValAtt)); //physically true

  // float valOffset = maxValProj - 65535.0; //not possible! always <=65535
  // if (valOffset < 0)
  //    valOffset = 0.0;

  // std::cout << "MaxValProj=" << maxValProj << " MInval=" << minValProj << "
  // ValOffset = " << valOffset << std::endl;

  itk::ImageRegionConstIterator<FloatImageType> itSrc(
      resultFwdImg, resultFwdImg->GetRequestedRegion()); // 2D

  float fProjVal = 0.0; // mu_t, the lower means the higher attn.
  double tmpConvVal = 0.0;

  // Convert line integral to intensity value (I0/I = exp(mu_t)) --> I =
  // I0/exp(mu_t)

  /*if (resultFwdImg->GetBufferedRegion().GetSize()[0] != pYKImage2D->m_iWidth)
      return;*/

  itSrc.GoToBegin();

  itk::ImageSliceIteratorWithIndex<UShortImageType> it_FwdProj3D(
      spProjImg3D, spProjImg3D->GetRequestedRegion());

  it_FwdProj3D.SetFirstDirection(0);
  it_FwdProj3D.SetSecondDirection(1);
  it_FwdProj3D.GoToBegin();

  int curSliceIdx = 0;

  while (!it_FwdProj3D.IsAtEnd()) {
    if (curSliceIdx == iSliceIdx) {
      // Search matching slice using slice iterator for m_spProjCTImg
      while (!it_FwdProj3D.IsAtEndOfSlice() && !itSrc.IsAtEnd()) {
        while (!it_FwdProj3D.IsAtEndOfLine() && !itSrc.IsAtEnd()) {
          fProjVal = itSrc.Get();                 // mu_t //63.5 --> 6.35
          tmpConvVal = (65535.0 / exp(fProjVal)); // intensity value

          unsigned short val = 0;
          if (tmpConvVal <= 0.0) {
            val = 0;
          } else if (tmpConvVal > 65535.0) {
            val = 65535;
          } else {
            val = static_cast<unsigned short>(tmpConvVal);
          }

          // unsigned short tmpVal = (unsigned short)(it_FwdProj3D.Get());
          // tmpVal = 65535 - tmpVal; //inverse is done here

          it_FwdProj3D.Set(val);

          ++it_FwdProj3D;
          ++itSrc;
        }
        it_FwdProj3D.NextLine();
      }
      it_FwdProj3D.NextSlice();
    }
    it_FwdProj3D.NextSlice();
    curSliceIdx++;
  }
  // Save this file

} // release all the memory

void CbctRecon::SLTM_LoadPerProjRefList() {
  QString filePath =
      QFileDialog::getOpenFileName(this, "PerProjVol list", m_strPathDirDefault,
                                   "File path list (*.txt)", nullptr, nullptr);

  if (filePath.length() < 1) {
    return;
  }

  std::ifstream fin;
  fin.open(filePath.toLocal8Bit().constData(), std::ios::in);
  if (fin.fail()) {
    return;
  }

  m_strListPerProjRefVol.clear();

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

    m_strListPerProjRefVol.push_back(strList.at(3));
  }

  std::cout << m_strListPerProjRefVol.count() << " image paths were found"
            << std::endl;

  fin.close();
}

// double CbctRecon::CropSkinUsingRS(UShortImageType::Pointer& spImgUshort,
// QString& strPathRS, double cropMargin )
//{
//   //if (m_pParent->m_strPathRS.isEmpty())
//	//return;
//  /* End of [1]Segment air region*/
//
//  //plastimatch convert --input E:\PlastimatchData\DicomEg\OLD\RS.dcm
//  --output-ss-img E:\PlastimatchData\DicomEg\OLD\ssimg_all2.mha
//  --output-ss-list E:\PlastimatchData\DicomEg\OLD\sslist_all.txt
//  --referenced-ct E:\PlastimatchData\DicomEg\OLD\CT
//
//  //plm_clp_parse (&parms, &parse_fn, &usage_fn, argc, argv, 1)
//  //plastimatch segment --input E:\PlastimatchData\DicomEg\OLD\CT --output-img
//  E:\PlastimatchData\DicomEg\OLD\msk_bubbles_oldCT.mha --lower-threshold -600
//
//  //do_command_warp(argc, argv);
//
//
//
//
//}

bool SaveDoseGrayImage(
    const char *filePath, int width, int height, double spacingX,
    double spacingY, double originLeft_mm, double originTop_mm,
    unsigned short *pData) // export dose array to a specified file (16bit TIF)
{
  // Global variables
  long m_iSubFileType = 0;
  short m_iWidth = width;
  short m_iHeight = height;
  short m_iBitsPerSample = 16;
  short m_iCompression = 1;
  short m_iPhotometric = 0;
  long m_iStripOffset = 1024;
  short m_iSamplePerPixel = 1;
  long m_iRowsPerStrip = height;
  long m_iStripByteCnts = qRound(width * height * 2.0);

  short m_iResolUnit = 2;
  short m_iPgNum = 0; // or 1?
  unsigned short m_iMinSampleVal = 0;
  unsigned short m_iMaxSampleVal = 65535U; // old: 255

  RATIONAL m_rXResol{}; // spacingX in dpi
  RATIONAL m_rYResol{}; // spacingY

  RATIONAL m_rXPos{};
  RATIONAL m_rYPos{};

  m_rXResol.b = 10000000;
  m_rYResol.b = 10000000;

  m_rXResol.a =
      static_cast<long>(qRound(1 / spacingX * 25.4 * m_rXResol.b)); // dpi
  m_rYResol.a = static_cast<long>(qRound(1 / spacingY * 25.4 * m_rYResol.b));

  int m_iNextOffset = 0;

  if (pData == nullptr) {
    return false;
  }

  // Set Center
  QPoint dataPt;
  dataPt.setX(qRound(m_iWidth / 2.0));
  dataPt.setY(qRound(m_iHeight / 2.0));

  m_rXPos.b = 10000000;
  m_rYPos.b = 10000000;

  // double fLeftPosMM = -dataPt.x()*spacingX;
  // double fTopPosMM = dataPt.y()*spacingY;
  double fLeftPosMM = originLeft_mm;
  double fTopPosMM = -originTop_mm;

  m_rXPos.a = static_cast<long>(qRound(fLeftPosMM / 25.4 * m_rXPos.b));
  m_rYPos.a = static_cast<long>(qRound(fTopPosMM / 25.4 * m_rYPos.b));
  // Set Center

  FILE *fd = nullptr;

#ifdef WIN32
  if (fopen_s(&fd, filePath, "wb") == 0) {
    std::cerr << "Could not open file: " << filePath << " for writing!"
              << std::endl;
    return false;
  }
#else
  fd = fopen(filePath, "wb");
  if (fd == nullptr) {
    std::cerr << "Could not open file: " << filePath << " for writing!"
              << std::endl;
    return false;
  }
#endif

  long MarkerUpper;
  long MarkerLower;

  MarkerUpper = 0x002A4949;
  MarkerLower = 0x00000008;

  fwrite(&MarkerUpper, sizeof(long), 1, fd); // 4
  fwrite(&MarkerLower, sizeof(long), 1, fd); // 8

  // int IFDSize = GetValidIFDCnt();
  int IFDSize = 18;

  fwrite(&IFDSize, sizeof(unsigned short), 1, fd); // 10

  // auto* IFDArr = new TIFIFD[IFDSize];
  std::vector<TIFIFD> IFDarr;
  IFDarr.reserve(IFDSize);

  int offsetX;
  int offsetY = 0;

  /*int idx = 0;
  int TagID = 0;
  int dataType = 0;
  int DataCnt = 0;
  int dataVal = 0;
      */
  if (m_iSubFileType >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 254;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iSubFileType;
    IFDarr.push_back(tififd_tmp);
  }

  if (m_iWidth >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 256;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iWidth;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iHeight >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 257;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iHeight;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iBitsPerSample >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 258;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iBitsPerSample;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iCompression >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 259;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iCompression;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iPhotometric >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 262;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iPhotometric;
    IFDarr.push_back(tififd_tmp); // 1로 강제 지정
    // dataVal = 0; //0으로 강제 지정
    // dataVal이 초기값이면 insert 안함
  }
  if (m_iStripOffset >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 273;
    tififd_tmp.DataType = 4;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iStripOffset;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iSamplePerPixel >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 277;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iSamplePerPixel;
    IFDarr.push_back(tififd_tmp);
    // 1로강제지정
    // dataVal = 1;
    // dataVal이 초기값이면 insert 안함
  }
  if (m_iRowsPerStrip >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 278;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iRowsPerStrip;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iStripByteCnts >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 279;
    tififd_tmp.DataType = 4;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iStripByteCnts;
    IFDarr.push_back(tififd_tmp);
    /*if (m_iSamplePerPixel == 1)
    dataVal = m_iStripByteCnts;
    else if (m_iSamplePerPixel == 3)
    dataVal = (int)(m_iStripByteCnts/3.0);
    */
    // dataVal이 초기값이면 insert 안함
  }
  if (m_rXResol.a != 0) {
    offsetX = 8 + 2 + (12 * IFDSize) + 4;

    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 282;
    tififd_tmp.DataType = 5;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = offsetX; // maximum
    IFDarr.push_back(tififd_tmp);      // dataVal이 초기값이면 insert 안함
  }
  if (m_rYResol.a != 0) {
    offsetY = 8 + 2 + (12 * IFDSize) + 4 + 8;

    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 283;
    tififd_tmp.DataType = 5;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = offsetY;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }

  // IFDSize 단위 데이터 몇개인지 나타냄
  // 20111226추가 //center를 표시
  if (m_rXPos.a != 0) {
    offsetX = 8 + 2 + (12 * IFDSize) + 4 + 8 + 8;

    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 286;
    tififd_tmp.DataType = 5;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = offsetX; // maximum
    IFDarr.push_back(tififd_tmp);      // dataVal이 초기값이면 insert 안함
  }
  if (m_rYPos.a != 0) {
    offsetY = 8 + 2 + (12 * IFDSize) + 4 + 8 + 8 + 8;

    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 287;
    tififd_tmp.DataType = 5;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = offsetY;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }

  //////
  if (m_iMinSampleVal >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 280;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iMinSampleVal;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iMaxSampleVal >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 281;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iMaxSampleVal;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iResolUnit >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 296;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 1;
    tififd_tmp.DataOrOffset = m_iResolUnit;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  if (m_iPgNum >= 0) {
    TIFIFD tififd_tmp{};
    tififd_tmp.TagID = 297;
    tififd_tmp.DataType = 3;
    tififd_tmp.DataCnt = 2;
    tififd_tmp.DataOrOffset = m_iPgNum;
    IFDarr.push_back(tififd_tmp); // dataVal이 초기값이면 insert 안함
  }
  /*
for (int i = 0; i < IFDSize; i++)
{
  fwrite(&IFDArr[i], sizeof(TIFIFD), 1, fd);
}
  */
  for (auto &&it : IFDarr) {
    fwrite(&it, sizeof(TIFIFD), 1, fd);
  }
  fwrite(&m_iNextOffset, 4, 1, fd);

  fwrite(&m_rXResol, 8, 1, fd);
  fwrite(&m_rYResol, 8, 1, fd);

  fwrite(&m_rXPos, 8, 10, fd);
  fwrite(&m_rYPos, 8, 1, fd);

  int iDummySize = 0;
  iDummySize = 1024 - (offsetY + 8);

  // char tmpDummy [802]; // 1024 -222

  auto *tmpDummy = new char[iDummySize];
  memset(tmpDummy, 0, iDummySize);
  fwrite(tmpDummy, iDummySize, 1, fd); //`까지 0으로 채움

  delete[] tmpDummy;
  // delete[] IFDArr;

  int imgSize = m_iWidth * m_iHeight;
  // fwrite(m_pImage, imgSize, 1, fd);

  //쓰기용 버퍼 생성
  // unsigned short* writeBuf = new unsigned short[imgSize];

  // for (int i = 0 ; i<imgSize ; i++)
  //{
  //	//fread(&m_pImage[i], 2, 1, fd);
  //	if (pData[i] < 0)
  //		writeBuf[i] = 0;
  //	else if  (pData[i] > 65535)
  //		writeBuf[i] = 65535;
  //	else
  //		writeBuf[i] = pData[i];  //gray 이미지를 건드는 것이다!
  //}

  for (int i = 0; i < imgSize; i++) {
    fwrite(&pData[i], 2, 1, fd);
  }
  fclose(fd);

  return true;
}

double vectorMean(const std::vector<double> &vDouble) {
  int size = vDouble.size();
  if (size <= 0) {
    return 0.0;
  }

  double sum = 0.0;
  std::vector<double>::const_iterator it;

  for (it = vDouble.begin(); it != vDouble.end(); it++) {
    sum = sum + (*it);
  }

  return sum / static_cast<double>(size);
}

double vectorSum(const std::vector<double> &vDouble) {
  int size = vDouble.size();
  if (size <= 0) {
    return 0.0;
  }

  double sum = 0.0;
  std::vector<double>::const_iterator it;

  for (it = vDouble.begin(); it != vDouble.end(); it++) {
    sum = sum + (*it);
  }

  return sum;
}

// Projection image Median filtering using CUDA

// Dir or File
bool CbctRecon::LoadShortImageDirOrFile(
    QString &strPathDir, ShortImageType::Pointer &spOutputShortImg) {
  QFileInfo fInfo(strPathDir);
  if (!fInfo.exists()) {
    return false;
  }

  Plm_image plmImg;
  plmImg.load_native(strPathDir.toLocal8Bit().constData());
  ShortImageType::Pointer spShortImg = plmImg.itk_short();

  // Figure out whether this is NKI
  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(spShortImg);
  imageCalculatorFilter->Compute();

  /* double minVal0 = (double)(imageCalculatorFilter->GetMinimum());
  double maxVal0 = (double)(imageCalculatorFilter->GetMaximum());*/

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;

  ThresholdImageFilterType::Pointer thresholdFilterAbove =
      ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(spShortImg);
  thresholdFilterAbove->ThresholdAbove(3072);
  thresholdFilterAbove->SetOutsideValue(3072);

  ThresholdImageFilterType::Pointer thresholdFilterBelow =
      ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(-1024);
  thresholdFilterBelow->SetOutsideValue(-1024);
  thresholdFilterBelow->Update();

  spOutputShortImg = thresholdFilterBelow->GetOutput();
  std::cout << "Image file was loaded" << std::endl;

  return true;
}

void CbctRecon::ConvertShort2Ushort(
    ShortImageType::Pointer &spInputImgShort,
    UShortImageType::Pointer &spOutputImgUshort) {
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilterAbove =
      ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(spInputImgShort);
  thresholdFilterAbove->ThresholdAbove(3071);
  thresholdFilterAbove->SetOutsideValue(3071);

  ThresholdImageFilterType::Pointer thresholdFilterBelow =
      ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(-1024);
  thresholdFilterBelow->SetOutsideValue(-1024);
  thresholdFilterBelow->Update();

  using ImageCalculatorFilterType =
      itk::MinimumMaximumImageCalculator<ShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
      ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(thresholdFilterBelow->GetOutput());
  imageCalculatorFilter->Compute();
  auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  // Min value is always 3024 --> outside the FOV
  auto outputMinVal = static_cast<UShortImageType::PixelType>(minVal + 1024);
  auto outputMaxVal = static_cast<UShortImageType::PixelType>(maxVal + 1024);

  using RescaleFilterType =
      itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilterBelow->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();

  spOutputImgUshort = spRescaleFilter->GetOutput();
}

void CbctRecon::RotateImgBeforeFwd(UShortImageType::Pointer &spInputImgUS,
                                   UShortImageType::Pointer &spOutputImgUS) {
  if (spInputImgUS == nullptr) {
    std::cout << "ERROR! No 3D image file" << std::endl;
    return;
  }
  // 1) Transform
  UShortImageType::SizeType size_original =
      spInputImgUS->GetLargestPossibleRegion().GetSize();
  UShortImageType::SpacingType spacing_original = spInputImgUS->GetSpacing();

  // Same image type from original image -3D & float
  UShortImageType::IndexType start_trans{};
  start_trans[0] = 0;
  start_trans[1] = 0;
  start_trans[2] = 0;

  UShortImageType::SizeType size_trans{};
  size_trans[0] = size_original[1]; // X //512
  size_trans[1] = size_original[2]; // Y  //512
  size_trans[2] = size_original[0]; // Z //300

  UShortImageType::SpacingType spacing_trans;
  spacing_trans[0] = spacing_original[1];
  spacing_trans[1] = spacing_original[2];
  spacing_trans[2] = spacing_original[0];

  UShortImageType::PointType Origin_trans;
  Origin_trans[0] = -0.5 * size_trans[0] * spacing_trans[0];
  Origin_trans[1] = -0.5 * size_trans[1] * spacing_trans[1];
  Origin_trans[2] = -0.5 * size_trans[2] * spacing_trans[2];

  UShortImageType::RegionType region_trans;
  region_trans.SetSize(size_trans);
  region_trans.SetIndex(start_trans);

  using FilterType = itk::FlipImageFilter<UShortImageType>;
  FilterType::Pointer flipFilter = FilterType::New();
  using FlipAxesArrayType = FilterType::FlipAxesArrayType;

  FlipAxesArrayType arrFlipAxes;
  arrFlipAxes[0] = true;
  arrFlipAxes[1] = false;
  arrFlipAxes[2] = false;

  flipFilter->SetFlipAxes(arrFlipAxes);
  flipFilter->SetInput(spInputImgUS); // plan CT, USHORT image

  using TransformType = itk::Euler3DTransform<double>;
  TransformType::Pointer transform = TransformType::New();

  TransformType::ParametersType param;
  param.SetSize(6);
  param.put(0, itk::Math::pi / -2.0); // rot X // 0.5 = PI/2
  param.put(1, 0);                    // rot Y
  param.put(2, itk::Math::pi / 2.0);  // rot Z
  param.put(3, 0.0);                  // Trans X mm
  param.put(4, 0.0);                  // Trans Y mm
  param.put(5, 0.0);                  // Trans Z mm

  TransformType::ParametersType fixedParam(3); // rotation center
  fixedParam.put(0, 0);
  fixedParam.put(1, 0);
  fixedParam.put(2, 0);

  transform->SetParameters(param);
  transform->SetFixedParameters(fixedParam); // Center of the Transform

  /*std::cout << "Transform matrix:" << "	" << std::endl;
  std::cout << transform->GetMatrix() << std::endl;*/

  using ResampleFilterType =
      itk::ResampleImageFilter<UShortImageType, UShortImageType>;
  ResampleFilterType::Pointer resampler = ResampleFilterType::New();

  resampler->SetInput(flipFilter->GetOutput());
  resampler->SetSize(size_trans);
  resampler->SetOutputOrigin(Origin_trans);   // Lt Top Inf of Large Canvas
  resampler->SetOutputSpacing(spacing_trans); // 1 1 1
  resampler->SetOutputDirection(
      flipFilter->GetOutput()->GetDirection()); // image normal?
  resampler->SetTransform(transform);
  resampler->Update();

  spOutputImgUS = resampler->GetOutput();
}

void CbctRecon::ConvertUshort2AttFloat(UShortImageType::Pointer &spImgUshort,
                                       FloatImageType::Pointer &spAttImgFloat) {
  using CastFilterType =
      itk::CastImageFilter<UShortImageType,
                           FloatImageType>; // Maybe not inplace filter
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(spImgUshort);

  // Default value
  double calibF_A = 1.0;
  double calibF_B = 0.0;

  using MultiplyImageFilterType =
      itk::MultiplyImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  MultiplyImageFilterType::Pointer multiplyImageFilter =
      MultiplyImageFilterType::New();
  multiplyImageFilter->SetInput(castFilter->GetOutput());
  multiplyImageFilter->SetConstant(calibF_A / 65535.0);

  using AddImageFilterType =
      itk::AddImageFilter<FloatImageType, FloatImageType, FloatImageType>;
  AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
  addImageFilter->SetInput1(multiplyImageFilter->GetOutput());
  double addingVal = calibF_B / 65535.0;
  addImageFilter->SetConstant2(addingVal);
  addImageFilter->Update(); // will generate map of real_mu (att.coeff)

  // FloatImageType::Pointer spCTImg_mu;
  spAttImgFloat = multiplyImageFilter->GetOutput();
}

FloatImageType::Pointer
CbctRecon::ConvertUshort2WeplFloat(UShortImageType::Pointer &spImgUshort) {
  ShortImageType::Pointer hu_image_tmp;
  ConvertUshort2Short(spImgUshort, hu_image_tmp);

  using CastFilterType = itk::CastImageFilter<ShortImageType, FloatImageType>;
  CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(hu_image_tmp);
  castFilter->Update();

  using double_pair_list = std::array<std::pair<float, double>, 16>;
  const double_pair_list lookup = {
      {// Data from TRiP: 19990218.hlut
       std::pair<float, double>(NLMIN(float), 0),
       std::pair<float, double>(-1000.0f, 0.041),
       std::pair<float, double>(-798.0f, 0.244),
       std::pair<float, double>(-750.0f, 0.297),
       std::pair<float, double>(-108.0f, 0.943),
       std::pair<float, double>(-75.0f, 0.977),
       std::pair<float, double>(0.0f, 1.0),
       std::pair<float, double>(40.0f, 1.042),
       std::pair<float, double>(55.0f, 1.049),
       std::pair<float, double>(67.0f, 1.065),
       std::pair<float, double>(262.0f, 1.095),
       std::pair<float, double>(1033.0f, 1.468),
       std::pair<float, double>(1432.0f, 1.634),
       std::pair<float, double>(1974.0f, 1.778),
       std::pair<float, double>(3000.0f, 2.051),
       std::pair<float, double>(NLMAX(float), 2.051)}};
  /*{ // plastimatch data:
          std::pair<float, double>(NLMIN(float), 0),
                  std::pair<float, double>(-1000, 0.00106),
                  std::pair<float, double>(0, 1.0),
                  std::pair<float, double>(41.46, 1.048674),
                  std::pair<float, double>(NLMAX(float), 0.005011) // wtf?
  };*/

  // Linear interpolator (as first class function)
  std::function<float(float)> hu_to_dEdx = [&lookup](float val) {
    // Find first index in lookup that satisfies "val < lookup[i].first" :
    auto lookup_upper_ptr = std::find_if(
        lookup.begin(), lookup.end(), [val](std::pair<float, double> cur_pair) {
          return val < cur_pair.first;
        });
    const auto lookup_upper = *lookup_upper_ptr;

    // Get the previous index:
    const auto lookup_lower =
        lookup.at((lookup_upper_ptr - lookup.begin()) - 1);

    // Do linear interpolation between upper and lower data point:
    const double a =
        (lookup_upper.second - lookup_lower.second) /
        static_cast<double>(lookup_upper.first - lookup_lower.first);
    const double b =
        lookup_upper.second - a * static_cast<double>(lookup_upper.first);

    return a * val + b;
  };

  FloatImageType::Pointer wepl_image = castFilter->GetOutput();

  itk::ImageRegionIterator<FloatImageType> it(
      wepl_image, wepl_image->GetLargestPossibleRegion());

  // std::for_each(it.Begin(), it.End(), hu_to_dEdx); //ITK iterators doesn't
  // support <algorithm> (yet?)
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    it.Set(hu_to_dEdx(it.Get()));
  }

  using WriterType = itk::ImageFileWriter<FloatImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetInput(wepl_image);
  writer->SetFileName("wepl_image.mha");
  writer->Update();

  return wepl_image;
}

void CbctRecon::SLTM_CropMaskBatch() {
  // Specify mask file (USHORT)
  QString maskFilePath = QFileDialog::getOpenFileName(
      this, "Mask image (Ushort)", m_strPathDirDefault, "3D mask file (*.mha)",
      nullptr, nullptr);

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
      this, "Select one or more files to open", m_strPathDirDefault,
      "target files (*.mha)");

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
    m_pDlgRegistration->plm_mask_main(mask_option, input_fn, mask_fn, output_fn,
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

void CbctRecon::SLT_SaveCurrentSetting() {
  if (!SaveCurrentSetting(m_strPathDefaultConfigFile)) {
    std::cout << "Error! in SaveCurrentSetting" << std::endl;
    return;
  }
}

bool CbctRecon::SaveCurrentSetting(QString &strPathConfigFile) {
  QFileInfo fInfo(strPathConfigFile);
  if (!fInfo.exists()) {
    std::cout << "Config file not exist. will be created now" << std::endl;
  } else {
    std::cout << "Config file is found. it will be overwritten now"
              << std::endl;
  }

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

bool CbctRecon::LoadCurrentSetting(QString &strPathConfigFile) {
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

void CbctRecon::SLT_CropSupInf() {
  if (m_spCrntReconImg == nullptr) {
    return;
  }

  int bRaw = 0;

  if (m_spCrntReconImg == m_spRawReconImg) {
    bRaw = 1;
  }

  float dcmPosCutSup = ui.lineEdit_SupCutPos->text().toFloat(); // mm
  float dcmPosCutInf = ui.lineEdit_InfCutPos->text().toFloat(); // mm

  // CropFOV3D(m_spCrntReconImg, physPosX, physPosY, physRadius, physTablePosY);
  CropSupInf(m_spCrntReconImg, dcmPosCutInf, dcmPosCutSup);
  // QString strPath = m_strPathDirDefault + "/" + "TempSI_Cropped.mha";
  // QString strTmpFile = "C:/TmpSI_Cropped.mha";

  QString strPath =
      m_pDlgRegistration->m_strPathPlastimatch + "/" + "tmp_SI_cropped.mha";
  ExportReconSHORT_HU(m_spCrntReconImg, strPath);

  QString strName = "SI_Cropped";
  if (bRaw != 0) {
    if (!LoadShortImageToUshort(strPath, m_spRawReconImg)) {
      std::cout << "error! in LoadShortImageToUshort" << std::endl;
    }
    UpdateReconImage(m_spRawReconImg, strName);
  } else {
    if (!LoadShortImageToUshort(strPath, m_spRefCTImg)) {
      std::cout << "error! in LoadShortImageToUshort" << std::endl;
    }
    UpdateReconImage(m_spRefCTImg, strName);
  }

  ///*So buggy*/

  /*QString strName = "SI_Cropped";
  UpdateReconImage(m_spCrntReconImg, strName);*/

  // SLT_DrawReconImage();
}
