/*All IO functions used with cbctrecon*/


bool SaveDoseGrayImage(const char *filePath, int width, int height,
                       double spacingX, double spacingY, double originLeft_mm,
                       double originTop_mm, unsigned short *pData);


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

template <typename ImageType>
void saveImageAsMHA(typename ImageType::Pointer image) {
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(image);
  writer->SetFileName("Projections.mha");
  writer->Update();
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
