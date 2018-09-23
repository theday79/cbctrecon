/*All IO functions used with cbctrecon*/
#include "cbctrecon.h"

// std
#include <math.h>                                              // for sqrt
#include <algorithm>                                           // for copy, max
#include <stddef.h>                               // for size_t
#include <iostream>                               // for operator<<, basic_o...
#include <memory>                                 // for unique_ptr, allocat...
#include <string>                                 // for string
#include <vector>                                 // for vector

// Qt
#include <qbytearray.h>                           // for QByteArray
#include <qdir.h>                                 // for QDir
#include <qfile.h>                                // for QFile
#include <qfiledialog.h>                          // for QFileDialog, operator|
#include <qfileinfo.h>                            // for QFileInfo
#include <qglobal.h>                              // for qRound
#include <qlineedit.h>                            // for QLineEdit
#include <qmessagebox.h>
#include <qpushbutton.h>                          // for QPushButton
#include <qradiobutton.h>                         // for QRadioButton
#include <qspinbox.h>                             // for QSpinBox
#include <qstring.h>                              // for QString, operator+
#include <qstringlist.h>                          // for QStringList

// ITK
#include "itkAddImageFilter.h"
#include "itkArray.h"                                        // for Array:...
#include "itkBinaryFunctorImageFilter.h"                     // for Binary...
#include "itkCastImageFilter.h"
#include "itkConstNeighborhoodIterator.h"                    // for ConstN...
#include "itkConvertPixelBuffer.h"                           // for Conver...
#include "itkDefaultConvertPixelTraits.h"                      // for Defaul...
#include "itkExtractImageFilter.h"                           // for Extrac...
#include "itkFixedArray.h"                        // for FixedArray
#include "itkImage.h"                             // for Image<>::SizeType
#include "itkImageBase.h"                                    // for ImageB...
#include "itkImageConstIteratorWithIndex.h"                  // for ImageC...
#include "itkImageDuplicator.h"
#include "itkImageFileReader.h"                              // for ImageF...
#include "itkImageFileWriter.h"                              // for ImageF...
#include "itkImageFunction.h"                                // for ImageF...
#include "itkImageIORegion.h"                                  // for operat...
#include "itkImageRegion.h"                                  // for ImageR...
#include "itkImageRegionConstIteratorWithIndex.h"            // for ImageR...
#include "itkImageRegionIterator.h"                          // for ImageR...
#include "itkImageRegionIteratorWithIndex.h"                 // for ImageR...
#include "itkImageScanlineIterator.h"                        // for ImageS...
#include "itkImageSource.h"                       // for ImageSource<>::Outp...
#include "itkImageToImageFilter.h"              // for ImageToImageFilter:...
#include "itkImportImageContainer.h"                         // for Import...
#include "itkIndex.h"                                          // for operat...
#include "itkInPlaceImageFilter.h"              // for InPlaceImageFilter:...
#include "itkInterpolateImageFunction.h"                       // for Interp...
#include "itkLightObject.h"                       // for LightObject::Pointer
#include "itkLinearInterpolateImageFunction.h"                 // for Linear...
#include "itkLogImageFilter.h"
#include "itkMacro.h"                                          // for Except...
#include "itkMatrix.h"                                         // for operat...
#include "itkMultiplyImageFilter.h"
#include "itkMinimumMaximumImageCalculator.h"                // for Minimu...
#include "itkNeighborhood.h"                                 // for Neighb...
#include "itkOffset.h"                                         // for operat...
#include "itkPoint.h"                                        // for operat...
#include "itkProcessObject.h"                     // for ProcessObject::Data...
#include "itkRegion.h"                                         // for Region...
#include "itkRescaleIntensityImageFilter.h"
#include "itkSize.h"                              // for operator<<
#include "itkSmartPointer.h"                      // for SmartPointer, Smart...
#include "itkSimpleDataObjectDecorator.h"                    // for Simple...
#include "itkThresholdImageFilter.h"                         // for Thresh...
#include "itkUnaryFunctorImageFilter.h"                      // for UnaryF...
#include "itkVector.h"                                       // for operat...
#include "rtkLookupTableImageFilter.h"                       // for Lookup...
#include "vnl_matrix_fixed.h"                                  // for operat...

// RTK
#include "rtkMacro.h"                             // for TRY_AND_EXIT_ON_ITK...
#include "rtkElektaSynergyGeometryReader.h"
#include "rtkThreeDCircularProjectionGeometry.h"  // for ThreeDCircularProje...
#include "rtkThreeDCircularProjectionGeometryXMLFileReader.h"  // for ThreeDCircularProje...
#include "rtkThreeDCircularProjectionGeometryXMLFileWriter.h"  // for ThreeDCircularProje...
#include "rtkProjectionsReader.h"

// PLM
#include "nki_io.h"
#include "mha_io.h"

// local
#include "cbctrecon_io.h"
#include "YK16GrayImage.h"                        // for YK16GrayImage
#include "ui_cbctrecon.h"                         // for CbctReconClass

class Volume;


// Function for independent projection his images
void CbctRecon::LoadRawHisImages() {

  QStringList files = QFileDialog::getOpenFileNames(
      this, "Select one or more files to open", m_strPathDirDefault,
      "projection images (*.his,*.hnd,*.xim)");

  m_iImgCnt = files.size();
  std::vector<std::string> fileVector;

  for (auto &cur_file : files) {
    fileVector.push_back(cur_file.toStdString());
  }

  if (m_iImgCnt < 1) {
    return;
  }

  ReleaseMemory();

  m_arrYKImage.resize(m_iImgCnt);
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
  int sizeBuf =
      sizePix *
      sizeof(FloatImageType::PixelType);
  int bytesPerPix = qRound(sizeBuf / static_cast<double>(sizePix));

  size_t index = 0;
  for (auto &it : m_arrYKImage) {
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

  m_multiplyFactor = 1.0;

  ui.spinBoxImgIdx->setMinimum(0);
  ui.spinBoxImgIdx->setMaximum(m_iImgCnt - 1);
  ui.spinBoxImgIdx->setValue(0);

  SetMaxAndMinValueOfProjectionImage();
  SLT_InitializeGraphLim();

  SLT_DrawRawImages(); // Change FileName as well.. read spinbox value and draw
                       // image
}


void CbctRecon::SLT_LoadImageFloat3D() // Dose image for JPhillips
{
  using ReaderType = itk::ImageFileReader<FloatImageType>;
  ReaderType::Pointer reader = ReaderType::New();

  QString fileName = QFileDialog::getOpenFileName(
      this, "Open Image", m_strPathDirDefault, "3D dose float file (*.mha)",
      nullptr, nullptr);

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

  m_spRawReconImg = castFilter->GetOutput();
  m_spCrntReconImg = m_spRawReconImg;

  // Update UI
  UShortImageType::SizeType imgDim =
      m_spRawReconImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing = m_spRawReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  ui.lineEdit_Cur3DFileName->setText(fileName);

  m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}


void CbctRecon::SLT_Load3DImage() // mha reconstructed file, from external
                                  // source
{
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

  m_spRawReconImg = reader->GetOutput();
  m_spCrntReconImg = m_spRawReconImg;

  UShortImageType::SizeType imgDim =
      m_spRawReconImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing = m_spRawReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  ui.lineEdit_Cur3DFileName->setText(fileName);

  m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}


void CbctRecon::SLT_Load3DImageShort() {
  QString fileName =
      QFileDialog::getOpenFileName(this, "Open Image", m_strPathDirDefault,
                                   "short mha file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  if (!LoadShortImageToUshort(fileName, m_spRawReconImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  ImageCalculatorFilterType2::Pointer imageCalculatorFilter2 =
      ImageCalculatorFilterType2::New();
  // imageCalculatorFilter2->SetImage(m_spReconImg);
  imageCalculatorFilter2->SetImage(m_spRawReconImg);
  imageCalculatorFilter2->Compute();

  auto minVal2 = static_cast<double>(imageCalculatorFilter2->GetMinimum());
  auto maxVal2 = static_cast<double>(imageCalculatorFilter2->GetMaximum());
  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  m_spCrntReconImg = m_spRawReconImg;

  UShortImageType::SizeType imgDim =
      m_spCrntReconImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing = m_spCrntReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  ui.lineEdit_Cur3DFileName->setText(fileName);

  m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}


void CbctRecon::SLT_LoadNKIImage() {
  QString filePath =
      QFileDialog::getOpenFileName(this, "Open Image", m_strPathDirDefault,
                                   "NKI file (*.SCAN)", nullptr, nullptr);

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

  if (!LoadShortImageToUshort(newPath, m_spRawReconImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
      itk::MinimumMaximumImageCalculator<UShortImageType>;

  ImageCalculatorFilterType2::Pointer imageCalculatorFilter2 =
      ImageCalculatorFilterType2::New();
  // imageCalculatorFilter2->SetImage(m_spReconImg);
  imageCalculatorFilter2->SetImage(m_spRawReconImg);
  imageCalculatorFilter2->Compute();

  auto minVal2 = static_cast<double>(imageCalculatorFilter2->GetMinimum());
  auto maxVal2 = static_cast<double>(imageCalculatorFilter2->GetMaximum());
  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
            << std::endl;

  // Update UI
  m_spCrntReconImg = m_spRawReconImg;

  UShortImageType::SizeType imgDim =
      m_spCrntReconImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing = m_spCrntReconImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
            << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
            << "	" << spacing[2] << std::endl;

  ui.lineEdit_Cur3DFileName->setText(newFileName);

  m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(
      initVal); // DrawRecon Imge is called, but sometimes skipped
  ui.radioButton_graph_recon->setChecked(true);

  SLT_DrawReconImage();
}

QString CbctRecon::MakeElektaXML(const QString &filePath_ImageDBF,
                                 const QString &filePath_FrameDBF,
                                 const QString &DICOM_UID) {
  std::cout << "Elekta geometry XML file is being generated." << std::endl;
  // Define FRAME.DBF path
  rtk::ElektaSynergyGeometryReader::Pointer reader =
      rtk::ElektaSynergyGeometryReader::New();
  // string strDicomUID = DICOM_UID.toLocal8Bit().constData();
  // //DICOM_UID.toStdString()  string strDicomUID = DICOM_UID.toStdString();
  // string strDbfImg = filePath_ImageDBF.toStdString();
  // string strDbfFrame = filePath_FrameDBF.toStdString();

  QFileInfo info = QFileInfo(filePath_ImageDBF);
  QString dirPath = info.absolutePath();

  QString fileName = "ElektaGeom_" + DICOM_UID + ".xml";

  QString strOutput = dirPath + "/" + fileName;

  reader->SetDicomUID(DICOM_UID.toLocal8Bit().constData());
  reader->SetImageDbfFileName(filePath_ImageDBF.toLocal8Bit().constData());
  reader->SetFrameDbfFileName(filePath_FrameDBF.toLocal8Bit().constData());

  TRY_AND_EXIT_ON_ITK_EXCEPTION(reader->UpdateOutputData());

  // Write
  rtk::ThreeDCircularProjectionGeometryXMLFileWriter::Pointer xmlWriter =
      rtk::ThreeDCircularProjectionGeometryXMLFileWriter::New();
  xmlWriter->SetFilename(strOutput.toLocal8Bit().constData());
  xmlWriter->SetObject(reader->GetGeometry());
  TRY_AND_EXIT_ON_ITK_EXCEPTION(xmlWriter->WriteFile())

  std::cout << "Reading succeed" << std::endl;

  return strOutput;
}


void CbctRecon::SLT_ExportHis() {
  if (m_iImgCnt < 1) {
    std::cout << "Error: Load raw his images first" << std::endl;
    return;
  }

  // Get Folder Name!

  // For displaying Dir only..
  QString dir = QFileDialog::getExistingDirectory(
      this, "Open Directory", "/home",
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  // FileName should be same, only selected folder

  for (int i = 0; i < m_iImgCnt; i++) {
    QFileInfo tmpInfo = QFileInfo(m_arrYKImage[i].m_strFilePath);
    QString newPath = dir + "/" + tmpInfo.fileName();
    m_arrYKImage[i].SaveDataAsHis(newPath.toLocal8Bit().constData(), false);
  }

  std::cout << "File export was done successfully" << std::endl;
}


// Get the projection geometry
void CbctRecon::LoadRTKGeometryFile(const char *filePath) {
  rtk::ThreeDCircularProjectionGeometryXMLFileReader::Pointer geometryReader;
  geometryReader = rtk::ThreeDCircularProjectionGeometryXMLFileReader::New();
  geometryReader->SetFilename(filePath);
  TRY_AND_EXIT_ON_ITK_EXCEPTION(geometryReader->GenerateOutputInformation())
  std::cout << "Geometry reading succeed" << std::endl;

  m_spFullGeometry = geometryReader->GetOutputObject();

  // fullGeometry->GetGantryAngles();
  int geoDataSize =
      m_spFullGeometry->GetGantryAngles().size(); // This is MV gantry angle!!!
  std::cout << "Geometry data size(projection gantry angles): " << geoDataSize
            << std::endl;
  if (geoDataSize < 1) {
    return;
  }

  // CW: continuously ascending except 360 - 0 interface, no negative value
  // CCW: continuously descending except 0 - 360 interface, no negative value

  // m_bScanDirectionCW

  /*for (int i = 0 ; i<geoDataSize ; i++)
  {
  std::cout << "Projection gantry angle " << i << ": " <<
  m_spFullGeometry->GetGantryAngles().at(i) << std::endl;
  }
  std::cout << "Coordination is following IEC coordinate (e.g. 190 deg = RPO,
  170 LPO, 30 = LAO supposing HFS) no and kV source position (not MV gantry)" <<
  std::endl;*/

  std::vector<double> vTempConvAngles;

  auto itBegin = (m_spFullGeometry->GetGantryAngles()).begin();
  auto itEnd = (m_spFullGeometry->GetGantryAngles()).end();

  for (auto it = itBegin; it != itEnd; it++) {
    double tmpAngle = (*it);

    if (tmpAngle > 180.0) {
      tmpAngle = tmpAngle - 360.0;
    }

    vTempConvAngles.push_back(tmpAngle);
  }

  // compare 2 points in the middle of the angle list
  auto iLowerIdx = static_cast<int>(geoDataSize * 1.0 / 3.0);
  auto iUpperIdx = static_cast<int>(geoDataSize * 2.0 / 3.0);

  if (vTempConvAngles.at(iLowerIdx) <
      vTempConvAngles.at(iUpperIdx)) // ascending
  {
    m_bScanDirectionCW = true;
    std::cout << "The scan direction is CW" << std::endl;
  }

  else {
    m_bScanDirectionCW = false;
    std::cout << "The scan direction is CCW" << std::endl;
  }

  std::cout << "AngularGaps Size: "
            << m_spFullGeometry
                   ->GetAngularGaps(m_spFullGeometry->GetSourceAngles())
                   .size()
            << std::endl;

}

void CbctRecon::SLT_ReloadProjections() {
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
  m_spProjImg3DFloat = ImageReader->GetOutput();

  // Copied from SLT_LoadSelectedFiles:

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

  // m_spProjImgRaw3D is Ushort

  std::cout << "Projection reading succeeded." << m_vSelectedFileNames.size()
            << " files were read" << std::endl;

  // Because you can load projections from previous run:
  ui.pushButton_DoRecon->setEnabled(true);

  ui.spinBoxImgIdx->setMinimum(0);
  ui.spinBoxImgIdx->setMaximum(m_vSelectedFileNames.size() - 1);
  ui.spinBoxImgIdx->setValue(0);

  SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  SLT_InitializeGraphLim();

  this->SLT_DrawProjImages(); // Update Table is called
}

bool CbctRecon::LoadShortImageToUshort(QString &strPath,
  UShortImageType::Pointer &pUshortImage) {
  using ReaderType = itk::ImageFileReader<ShortImageType>;
  ReaderType::Pointer reader = ReaderType::New();

  // QString fileName = QFileDialog::getOpenFileName(this, "Open Image","",
  // "Plan CT file (*.mha)",0,0);

  if (strPath.length() < 1) {
    return false;
  }

  reader->SetFileName(strPath.toLocal8Bit().constData());
  reader->Update();

  // Figure out whether this is NKI
  using ImageCalculatorFilterType =
    itk::MinimumMaximumImageCalculator<ShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
    ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(reader->GetOutput());
  imageCalculatorFilter->Compute();

  auto minVal0 = static_cast<double>(imageCalculatorFilter->GetMinimum());
  auto maxVal0 = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Original Min and Max Values are	" << minVal0 << "	"
    << maxVal0 << std::endl;

  bool bNKI = false;
  if (minVal0 > -600) // impossible for normal Short image. IN NKI, always -512.
                      // don't know why
  {
    bNKI = true;
  }

  // Thresholding
  using ThresholdImageFilterType = itk::ThresholdImageFilter<ShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilter =
    ThresholdImageFilterType::New();

  if (!bNKI) {
    thresholdFilter->SetInput(reader->GetOutput());
    thresholdFilter->ThresholdOutside(-1024, 3072); //--> 0 ~ 4095
    thresholdFilter->SetOutsideValue(-1024);
    thresholdFilter->Update();
  }
  else {
    thresholdFilter->SetInput(reader->GetOutput());
    thresholdFilter->ThresholdOutside(0, 4095); //--> 0 ~ 4095
    thresholdFilter->SetOutsideValue(0);
    thresholdFilter->Update();
  }

  imageCalculatorFilter->SetImage(thresholdFilter->GetOutput());
  imageCalculatorFilter->Compute();

  auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  std::cout << "Current Min and Max Values are	" << minVal << "	"
    << maxVal << std::endl;

  USHORT_PixelType outputMinVal, outputMaxVal;
  if (!bNKI) {
    outputMinVal = static_cast<USHORT_PixelType>(minVal + 1024);
    outputMaxVal = static_cast<USHORT_PixelType>(maxVal + 1024);
  }
  else {
    outputMinVal = static_cast<USHORT_PixelType>(minVal);
    outputMaxVal = static_cast<USHORT_PixelType>(maxVal);
  }

  using RescaleFilterType =
    itk::RescaleIntensityImageFilter<ShortImageType, UShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilter->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);
  spRescaleFilter->Update();
  pUshortImage = spRescaleFilter->GetOutput();

  return true;
}

void CbctRecon::SLT_LoadPlanCT_mha() // m_spRecon -->m_spRefCT
{
  // typedef itk::ImageFileReader<ShortImageType> ReaderType;
  // ReaderType::Pointer reader = ReaderType::New();

  QString fileName =
    QFileDialog::getOpenFileName(this, "Open Image", m_strPathDirDefault,
      "Plan CT file (*.mha)", nullptr, nullptr);

  if (fileName.length() < 1) {
    return;
  }

  if (!LoadShortImageToUshort(fileName, m_spRefCTImg)) {
    std::cout << "error! in LoadShortImageToUshort" << std::endl;
  }

  using ImageCalculatorFilterType2 =
    itk::MinimumMaximumImageCalculator<UShortImageType>;

  ImageCalculatorFilterType2::Pointer imageCalculatorFilter2 =
    ImageCalculatorFilterType2::New();
  imageCalculatorFilter2->SetImage(m_spRefCTImg);
  imageCalculatorFilter2->Compute();

  auto minVal2 = static_cast<double>(imageCalculatorFilter2->GetMinimum());
  auto maxVal2 = static_cast<double>(imageCalculatorFilter2->GetMaximum());

  std::cout << "Min and Max Values are	" << minVal2 << "	" << maxVal2
    << std::endl;

  // Update UI
  UShortImageType::SizeType imgDim =
    m_spRefCTImg->GetBufferedRegion().GetSize();
  UShortImageType::SpacingType spacing = m_spRefCTImg->GetSpacing();

  std::cout << "Image Dimension:	" << imgDim[0] << "	" << imgDim[1]
    << "	" << imgDim[2] << std::endl;
  std::cout << "Image Spacing (mm):	" << spacing[0] << "	" << spacing[1]
    << "	" << spacing[2] << std::endl;

  RegisterImgDuplication(REGISTER_REF_CT, REGISTER_MANUAL_RIGID);
  m_spCrntReconImg = m_spRefCTImg;

  ui.lineEdit_Cur3DFileName->setText(fileName);
  m_dspYKReconImage->CreateImage(imgDim[0], imgDim[1], 0);

  ui.spinBoxReconImgSliceNo->setMinimum(0);
  ui.spinBoxReconImgSliceNo->setMaximum(imgDim[2] - 1);
  int initVal = qRound((imgDim[2] - 1) / 2.0);

  SLT_InitializeGraphLim();
  ui.spinBoxReconImgSliceNo->setValue(initVal); // DrawRecon Imge is called
  ui.radioButton_graph_recon->setChecked(true);

}


void CbctRecon::SLT_ExportReconUSHORT() {
  if (m_spCrntReconImg == nullptr) {
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
  writer->SetInput(m_spCrntReconImg);

  std::cout << "Writing is under progress...: "
    << strPath.toLocal8Bit().constData() << std::endl;
  writer->Update();
  std::cout << "Writing was successfully done" << std::endl;

  QString msgStr = QString("USHORT File Writing was successfully done");
  QMessageBox::information(this, "Procedure Done", msgStr);
}

void CbctRecon::ExportReconSHORT_HU(UShortImageType::Pointer &spUsImage,
  QString &outputFilePath) {
  if (spUsImage == nullptr) {
    std::cout << " no image to export" << std::endl;
    return;
  }

  using DuplicatorType = itk::ImageDuplicator<UShortImageType>;
  DuplicatorType::Pointer duplicator = DuplicatorType::New();
  duplicator->SetInputImage(spUsImage);
  duplicator->Update();
  UShortImageType::Pointer clonedReconImage = duplicator->GetOutput();
  ShortImageType::Pointer clonedReconImageSHORT;

  using ThresholdImageFilterType = itk::ThresholdImageFilter<UShortImageType>;
  ThresholdImageFilterType::Pointer thresholdFilterAbove =
    ThresholdImageFilterType::New();
  thresholdFilterAbove->SetInput(clonedReconImage);
  thresholdFilterAbove->ThresholdAbove(4095);
  thresholdFilterAbove->SetOutsideValue(4095);

  ThresholdImageFilterType::Pointer thresholdFilterBelow =
    ThresholdImageFilterType::New();
  thresholdFilterBelow->SetInput(thresholdFilterAbove->GetOutput());
  thresholdFilterBelow->ThresholdBelow(0);
  thresholdFilterBelow->SetOutsideValue(0);
  thresholdFilterBelow->Update();


  using ImageCalculatorFilterType =
    itk::MinimumMaximumImageCalculator<UShortImageType>;
  ImageCalculatorFilterType::Pointer imageCalculatorFilter =
    ImageCalculatorFilterType::New();
  imageCalculatorFilter->SetImage(thresholdFilterBelow->GetOutput());
  imageCalculatorFilter->Compute();
  auto minVal = static_cast<double>(imageCalculatorFilter->GetMinimum());
  auto maxVal = static_cast<double>(imageCalculatorFilter->GetMaximum());

  // Min value is always 3024 --> outside the FOV
  auto outputMinVal = static_cast<SHORT_PixelType>(minVal - 1024);
  auto outputMaxVal = static_cast<SHORT_PixelType>(maxVal - 1024);

  using RescaleFilterType =
    itk::RescaleIntensityImageFilter<UShortImageType, ShortImageType>;
  RescaleFilterType::Pointer spRescaleFilter = RescaleFilterType::New();
  spRescaleFilter->SetInput(thresholdFilterBelow->GetOutput());
  spRescaleFilter->SetOutputMinimum(outputMinVal);
  spRescaleFilter->SetOutputMaximum(outputMaxVal);

  spRescaleFilter->Update();
  // clonedReconImageSHORT = spRescaleFilter->GetOutput();

  // waterHU = 1024;

  using AddImageFilterType =
    itk::AddImageFilter<ShortImageType, ShortImageType, ShortImageType>;
  AddImageFilterType::Pointer addImageFilter = AddImageFilterType::New();
  addImageFilter->SetInput1(spRescaleFilter->GetOutput());

  int addingVal = 0; // 1024-680
  addImageFilter->SetConstant2(addingVal);
  addImageFilter->Update();

  using WriterType = itk::ImageFileWriter<ShortImageType>;
  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(outputFilePath.toLocal8Bit().constData());
  // writer->SetUseCompression(true);
  writer->SetUseCompression(false); // for plastimatch
  writer->SetInput(addImageFilter->GetOutput());

  std::cout << "Writing is under progress...: "
    << outputFilePath.toLocal8Bit().constData() << std::endl;
  writer->Update();
  std::cout << "Writing was successfully done" << std::endl;
}
