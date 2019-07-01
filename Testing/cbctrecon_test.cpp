// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

// For testing CbctRecon

#ifdef USE_TINYREFL
#include <tinyrefl/api.hpp> // MUST BE INCLUDED FIRST, FFS!

#include "cbctrecon.h"
#include "cbctrecon.h.tinyrefl"
#else
#include "cbctrecon.h"
#endif

#include <cassert>
#include <iostream>
#include <memory>
#include <thread>

#include <QDir>

#ifdef CUDA_FOUND
#undef CUDA_FOUND // Because both RTK and PLM defines this
#endif            // CUDA_FOUND
#include "cbctrecon_io.h"
#include "cbctrecon_test.hpp"
#include "cbctregistration.h"
#include "cbctregistration_test.hpp"

CbctReconTest::CbctReconTest() {
  m_dlgRegistration = std::make_unique<CbctRegistrationTest>(this);
  m_cbctregistration = std::make_unique<CbctRegistration>(m_cbctrecon.get());
  m_pTableModel = nullptr;

  QString tmp_folder("tmp");
  m_dlgRegistration->initCbctRegistrationTest(
      tmp_folder); // to Setup plastimatch folder. this is
                   // useful if registration will be only done
}

/* All of the following corresponds to a SLT_ function in mainwidget
 * all of which should have a test */
void CbctReconTest::test_LoadRawImages() {}
void CbctReconTest::test_Load3DImage() {}
void CbctReconTest::test_Load3DImageShort() {}
void CbctReconTest::test_LoadPlanCT_mha() {}
void CbctReconTest::test_LoadPlanCT_USHORT() {}
void CbctReconTest::test_LoadCBCTcorrMHA() {}
void CbctReconTest::test_LoadCTrigidMHA() {}
void CbctReconTest::test_LoadCTdeformMHA() {}
void CbctReconTest::test_LoadNKIImage() {}

QString getBowtiePath(const QDir &calDir) {
  return calDir.absolutePath() +
         "/AIR-Full-Bowtie-100KV/Current/FilterBowtie.xim";
}

FilterReaderType::Pointer CbctReconTest::ReadBowtieFileWhileProbing(
    const QString &proj_path, std::tuple<bool, bool> &answers) const {

  auto bowtiereader =
      FilterReaderType::New(); // we use is because we need the projections to
                               // be in the same unit (order of magnitude)

  // QDir guessDir(proj_path + QString("/../"));

  const auto calDir(proj_path + QString("/Calibrations/"));

  QString bowtiePath;
  answers = std::make_tuple(true, true);

  switch (this->m_cbctrecon->m_projFormat) {
  case XIM_FORMAT:
    bowtiePath = getBowtiePath(calDir);
    if (bowtiePath.length() > 1) {
      std::cerr << "loading bowtie-filter..."
                << "\n";
      std::vector<std::string> filepath;
      filepath.push_back(bowtiePath.toStdString());
      bowtiereader->SetFileNames(filepath);
      bowtiereader->Update();
    }
    break;
  default:
    break;
  }
  if (bowtiePath.length() > 1) {
    return bowtiereader;
  }
  return nullptr;
}

bool CbctReconTest::test_LoadSelectedProjFiles(const QString &proj_path,
                                               bool reconstruct) {
  // this->ui.pushButton_DoRecon->setDisabled(true);
  // 1) Get all projection file names
  auto dirPath = proj_path; // this->ui.lineEdit_HisDirPath->text();
  //.toLocal8Bit().constData();

  if (!QFile::exists(dirPath)) {
    std::cerr << "Projection file directory was not found. Retry."
              << "\n";
    return false;
  }

  auto names = this->m_cbctrecon->GetProjFileNames(dirPath);
  if (names.empty()) {
    std::cerr << "Regex didn't find any files in his, hns or xim format!"
              << "\n";
    return false;
  }

  if (this->m_cbctrecon->m_projFormat == HIS_FORMAT &&
      !this->m_cbctrecon->IsFileNameOrderCorrect(names)) {
    std::cerr << "Check the file name order"
              << "\n";
    return false;
  }

  std::cerr << "File name order was cross-checked and found to be OK!"
            << "\n";

  const auto fullCnt = names.size();
  if (fullCnt <= 0) {
    std::cerr << "No projection file was found. Retry."
              << "\n";
    return false;
  }

  std::cerr << fullCnt << "  projection files were found."
            << "\n";

  // 2) Elekta Geometry file
  QFileInfo geomFileInfo(this->m_cbctrecon->m_strPathGeomXML);

  if (!this->m_cbctrecon->LoadGeometry(geomFileInfo, names)) {
    if (!this->m_cbctrecon->m_strError.isEmpty()) {
      std::cerr << this->m_cbctrecon->m_strError.toStdString() << "\n";
      return false;
    }
  }

  const auto iFullGeoDataSize =
      this->m_cbctrecon->m_spFullGeometry->GetGantryAngles().size();
  if (iFullGeoDataSize < 1) {
    std::cerr << "Not enough projection image (should be > 0)"
              << "\n";
    return false;
  }

  if (iFullGeoDataSize != fullCnt) {
    if (this->m_cbctrecon->m_projFormat != XIM_FORMAT) {
      std::cerr << "Size of geometry data and file numbers are not same! Check "
                   "and retry"
                << "\n";
      return false;
    }
  }

  auto angle_gaps = this->m_cbctrecon->m_spFullGeometry->GetAngularGaps(
      this->m_cbctrecon->m_spFullGeometry->GetSourceAngles());

  auto sum_gap =
      std::accumulate(std::begin(angle_gaps), std::end(angle_gaps), 0.0);
  sum_gap /= itk::Math::pi * 180.0;
  const auto mean_gap = sum_gap / angle_gaps.size();

  std::cerr << "AngularGaps Sum (deg):" << sum_gap
            << ", Mean (deg): " << mean_gap << "\n";

  const auto gantryAngleInterval = 1.0;
  // this->ui.lineEdit_ManualProjAngleGap->text().toDouble();

  const auto exclude_ids = this->m_cbctrecon->GetExcludeProjFiles(
      false /*this->ui.Radio_ManualProjAngleGap->isChecked()*/,
      gantryAngleInterval);

  this->m_cbctrecon->LoadSelectedProj(exclude_ids, names);

  // Reads the cone beam projections
  using ReaderType = rtk::ProjectionsReader<FloatImageType>;
  auto reader = ReaderType::New();
  reader->SetFileNames(this->m_cbctrecon->m_vSelectedFileNames);
  std::thread calc_thread([&reader]() { reader->Update(); });

  std::cerr << "Reader detached from main thread"
            << "\n";

  // After reading the whole file,
  // HIS header should be saved
  this->m_cbctrecon->saveHisHeader();

  //  Insta Recon, Dcm read
  const auto geopath = geomFileInfo.absolutePath();
  std::tuple<bool, bool> answers;
  auto bowtie_reader = ReadBowtieFileWhileProbing(geopath, answers);

  calc_thread.join();
  std::cerr << "Reader re-attached to main thread"
            << "\n";

  this->m_cbctrecon->m_spProjImg3DFloat =
      reader->GetOutput(); // 1024 1024, line integ image

  auto &proj_ref = this->m_cbctrecon->m_spProjImg3DFloat;
  auto tmp_index = FloatImageType::IndexType();
  const auto tmp_size = proj_ref->GetLargestPossibleRegion().GetSize();
  tmp_index.SetElement(0, tmp_size[0] / 2);
  tmp_index.SetElement(1, tmp_size[1] / 2);
  tmp_index.SetElement(2, tmp_size[2] / 2);
  const auto test_value = proj_ref->GetPixel(tmp_index);
  std::cerr << "Test: " << test_value << " at index: " << tmp_index
            << "\n"; // this value seem to be different on different systems!
  // 5.26437... on Win, 5.10171... on Manjaro, 5.1791... on Ubuntu (travis-ci)

  if (bowtie_reader != nullptr) {
    const auto bowtie_ref = bowtie_reader->GetOutput();
    auto tmp_bt_index = FloatImage2DType::IndexType();
    const auto tmp_bt_size = bowtie_ref->GetLargestPossibleRegion().GetSize();
    tmp_bt_index.SetElement(0, tmp_bt_size[0] / 2);
    tmp_bt_index.SetElement(1, tmp_bt_size[1] / 2);
    const auto bowtie_test_value = bowtie_ref->GetPixel(tmp_bt_index);
    std::cerr << "Bowtie Test: " << bowtie_test_value << "\n";

    ApplyBowtie(proj_ref, bowtie_ref);

    const auto after_bowtie_test_value = proj_ref->GetPixel(tmp_index);
    std::cerr << "After Bowtie Test: " << after_bowtie_test_value << "\n";
    assert(fabs(after_bowtie_test_value - (test_value - bowtie_test_value)) <
           0.01f);
  }
  if (this->m_cbctrecon->m_projFormat == HND_FORMAT) {
    std::cerr << "Fitted bowtie-filter correction ongoing..."
              << "\n";
    test_DoBowtieCorrection();
  }

  auto res_factor = 0.5;
  if (!this->m_cbctrecon->ResampleProjections(res_factor)) { // 0.5
    // reset factor if image was not resampled
    std::cerr << "Could not resample projection size!\n";
  }

  this->m_cbctrecon->m_spProjImgRaw3D =
      this->m_cbctrecon->ConvertLineInt2Intensity(proj_ref);

  this->m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image

  test_InitializeGraphLim();

  this->test_DrawProjImages(); // Update Table is called

  if (!reconstruct) { //! std::get<0>(answers)) { // instaRecon
    std::cerr << "FINISHED!: Loading projection files."
              << "\n";
  } else {
    test_DoReconstruction();
  }

  if (std::get<1>(answers)) { // CT DCM dir was found
    test_ViewRegistration();
  }
  return true;
}

void CbctReconTest::test_ReloadProjections() {

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
  auto &p_projimg = this->m_cbctrecon->m_spProjImg3DFloat;
  // Copied from SLT_LoadSelectedFiles:

  if (this->m_cbctrecon->m_fResampleF != 1.0) {
    this->m_cbctrecon->ResampleItkImage(
        p_projimg, p_projimg,
        this->m_cbctrecon
            ->m_fResampleF); // was! BROKEN AF for .his where input size
                             // != 1024 (tested with 1016) -> outputs
                             // offset -inputoffset/refactor^2 and 4
                             // pixels too few in x and y
  }

  if (this->m_cbctrecon->m_projFormat == HND_FORMAT) { // -> hnd
    std::cout << "Fitted bowtie-filter correction ongoing..." << std::endl;
    test_DoBowtieCorrection();
  }

  this->m_cbctrecon->m_spProjImgRaw3D =
      this->m_cbctrecon->ConvertLineInt2Intensity(p_projimg);
  // if X not 1024 == input size: out_offset =
  // in_offset + (1024*res_f -
  // X*res_f)*out_spacing     <- will still
  // break down at fw_projection

  const auto p_projimgfloat = p_projimg;
  auto originPt = p_projimgfloat->GetOrigin();
  auto FloatImgSize = p_projimgfloat->GetBufferedRegion().GetSize();
  auto FloatImgSpacing = p_projimgfloat->GetSpacing();

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
  //this->ui.pushButton_DoRecon->setEnabled(true);

  //this->ui.spinBoxImgIdx->setMinimum(0);
  //this->ui.spinBoxImgIdx->setMaximum(
  //    this->m_cbctrecon->m_vSelectedFileNames.size() - 1);
  //this->ui.spinBoxImgIdx->setValue(0);

  this->m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  // SLT_InitializeGraphLim();

  // this->SLT_DrawProjImages(); // Update Table is called
}

void CbctReconTest::test_ExportHis() {}
void CbctReconTest::test_LoadImageFloat3D() {}

void CbctReconTest::test_LoadDICOMdir() const {
  auto dirPath = QString(this->m_cbctrecon->m_strPathDirDefault);

  if (dirPath.length() <= 1) {
    return;
  }

  if (this->m_cbctrecon->ReadDicomDir(dirPath)) {
    this->m_cbctrecon->RegisterImgDuplication(REGISTER_REF_CT,
                                              REGISTER_MANUAL_RIGID);
  }
}

void CbctReconTest::test_LoadRTKoutput() {}
void CbctReconTest::test_DrawRawImages() const {}
void CbctReconTest::test_DrawProjImages() {}
void CbctReconTest::test_DrawReconImage() {}
void CbctReconTest::test_FileNameHex2Dec() {}
void CbctReconTest::test_MakeElektaXML() {}
void CbctReconTest::test_OpenOffsetFile() {}
void CbctReconTest::test_OpenGainFile() {}
void CbctReconTest::test_OpenBadpixelFile() {}
void CbctReconTest::test_ApplyCalibration() const {}

void CbctReconTest::test_SetHisDir(QString &dirPath) {
  if (dirPath.length() <= 1) {
    return;
  }

  this->m_cbctrecon->SetProjDir(dirPath);
  // this->init_DlgRegistration(this->m_cbctrecon->m_strDCMUID);

  float kVp = 0.0;
  float mA = 0.0;
  float ms = 0.0;
  GetXrayParamFromINI(this->m_cbctrecon->m_strPathElektaINI, kVp, mA, ms);

  if (fabs(kVp * mA * ms) > 0.001) {
    // update GUI
    std::cout << "Updating current mAs setting from INI file: "
              << "kVp= " << kVp << ", mA= " << mA << ", ms= " << ms
              << std::endl;
  }

  VEC3D couch_trans = {-999, -999,
                       -999}; // mm. In the text file, these values are in cm.
  VEC3D couch_rot = {-999, -999,
                     -999}; // mm. In the text file, these values are in cm.

  /*
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
  }
  */

  this->m_cbctrecon->m_vSelectedFileNames.clear();
}

void CbctReconTest::test_OpenElektaGeomFile() {}
void CbctReconTest::test_SetOutputPath() {}

void CbctReconTest::test_DoReconstruction() {
  auto fdk_options = getFDKoptions();

  itk::TimeProbe reconTimeProbe;
  reconTimeProbe.Start();

#ifdef USE_CUDA
  const auto use_cuda = true;
#else
  const auto use_cuda = false;
#endif
  const auto use_opencl = true; //false; // CI OpenCL device doesn't support clImage

  if (use_cuda) {
    this->m_cbctrecon->DoReconstructionFDK<CUDA_DEVT>(REGISTER_RAW_CBCT,
                                                      fdk_options);
  } else if (use_opencl) {
    this->m_cbctrecon->DoReconstructionFDK<OPENCL_DEVT>(REGISTER_RAW_CBCT,
                                                        fdk_options);
  } else {
    this->m_cbctrecon->DoReconstructionFDK<CPU_DEVT>(REGISTER_RAW_CBCT,
                                                     fdk_options);
  }

  reconTimeProbe.Stop();
  std::cout << "It took " << reconTimeProbe.GetMean() << ' '
            << reconTimeProbe.GetUnit() << std::endl;

  test_DrawProjImages();
}

void CbctReconTest::test_CopyTableToClipBoard() const {}
void CbctReconTest::test_DataProbeProj() const {}
void CbctReconTest::test_DataProbeRecon() const {}
void CbctReconTest::test_DrawGraph() const {}
void CbctReconTest::test_InitializeGraphLim() const {}
void CbctReconTest::test_UpdateTable() {}
void CbctReconTest::test_CalculateROI_Recon() {}
void CbctReconTest::test_CalculateROI_Proj() {}
void CbctReconTest::test_GoForcedProbePos() {}
void CbctReconTest::test_PostApplyFOVDispParam() {}
void CbctReconTest::test_DoPostProcessing() {}
void CbctReconTest::test_PostProcCropInv() {}
void CbctReconTest::test_ExportReconUSHORT() {}
void CbctReconTest::test_ExportReconSHORT_HU() {}
void CbctReconTest::test_ExportALL_DCM_and_SHORT_HU_and_calc_WEPL() {}
void CbctReconTest::test_DoBHC() {}
void CbctReconTest::test_DoBowtieCorrection() {
  std::cerr << "HND BOWTIE CORRECTION TRIGGERED, BUT TESTS ONLY WRITTEN FOR XIM\n";
}
void CbctReconTest::test_Export2DDose_TIF() {}
void CbctReconTest::test_Export2DDoseMapAsMHA() {}
void CbctReconTest::test_ViewRegistration() const {}
void CbctReconTest::test_ViewHistogram() {}

void CbctReconTest::test_DoScatterCorrection_APRIORI() {

  const auto bExportProj_Fwd = false; // this->ui.checkBox_ExportFwd->isChecked();
  const auto bExportProj_Scat = false; // this->ui.checkBox_ExportScat->isChecked();
  const auto bExportProj_Cor = false; // this->ui.checkBox_ExportCor->isChecked();

  // ForwardProjection(m_spRefCTImg, m_spCustomGeometry, m_spProjImgCT3D,
  // false); //final moving image
  if (bExportProj_Fwd) {
    this->m_cbctrecon->m_strPathPatientDir = "";
  }
  auto &p_projimg = this->m_cbctrecon->m_spProjImgCT3D;

#ifdef USE_CUDA
  const auto use_cuda = true;
#else
  const auto use_cuda = false;
#endif

  if (m_dlgRegistration->m_spMoving != nullptr) {
    this->m_cbctrecon
        ->ForwardProjection_master<UShortImageType, UShortImageType>(
            m_dlgRegistration->m_spMoving,
            this->m_cbctrecon->m_spCustomGeometry,
            p_projimg, bExportProj_Fwd, use_cuda);
  } else if (this->m_cbctrecon->m_spRefCTImg != nullptr) {
    std::cerr << "No Moving image in Registration is found. Ref CT image will "
                 "be used instead"
              << "\n";
    this->m_cbctrecon
        ->ForwardProjection_master<UShortImageType, UShortImageType>(
            this->m_cbctrecon->m_spRefCTImg,
            this->m_cbctrecon->m_spCustomGeometry,
            p_projimg, bExportProj_Fwd, use_cuda);
  } else {
    std::cerr << "Error!: No ref image for forward projection is found."
              << "\n";
    return;
  }

  // double scaResam = this->ui.lineEdit_scaResam->text().toDouble();
  const auto scaMedian = 12.0; // this->ui.lineEdit_scaMedian->text().toDouble();
  const auto scaGaussian = 0.05; // this->ui.lineEdit_scaGaussian->text().toDouble();

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

  // this->ui.lineEdit_CurmAs->setText(this->m_cbctrecon->m_strCur_mAs);
  // this->ui.lineEdit_RefmAs->setText(this->m_cbctrecon->m_strRef_mAs);

  p_projimg->Initialize(); // memory saving

  std::cout << "Scatter correction is in progress..." << std::endl;

  const auto postScatMedianSize = 3.0; // this->ui.lineEdit_scaPostMedian->text().toInt();
  this->m_cbctrecon->ScatterCorr_PrioriCT(
      this->m_cbctrecon->m_spProjImgRaw3D, this->m_cbctrecon->m_spProjImgScat3D,
      this->m_cbctrecon->m_spProjImgCorr3D,
      this->m_cbctrecon->m_iFixedOffset_ScatterMap, postScatMedianSize,
      bExportProj_Cor);
  this->m_cbctrecon->m_spProjImgScat3D->Initialize(); // memory saving

  std::cout << "AfterCorrectionMacro is ongoing..." << std::endl;

  // Update UI
  // this->ui.pushButton_DoRecon->setEnabled(true);
  // this->ui.spinBoxImgIdx->setMinimum(0);
  const auto iSizeZ =
      this->m_cbctrecon->m_spProjImg3DFloat->GetRequestedRegion().GetSize()[2];
  // this->ui.spinBoxImgIdx->setMaximum(iSizeZ - 1);
  // this->ui.spinBoxImgIdx->setValue(0);
  this->m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image
  // SLT_InitializeGraphLim();
  // SLT_DrawProjImages(); // Update Table is called

  auto fdk_options = getFDKoptions();

  this->m_cbctrecon->AfterScatCorrectionMacro(
#ifdef USE_CUDA
            true,
#else
            false,
#endif
      true, // this->ui.radioButton_UseOpenCL->isChecked(),
      true, //this->ui.checkBox_ExportVolDICOM->isChecked(),
      fdk_options);

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

  // SLT_DrawProjImages();

  std::cout << "Updating ReconImage..";
  auto updated_text = QString("Scatter corrected CBCT");
  // UpdateReconImage(this->m_cbctrecon->m_spScatCorrReconImg,
  //                  updated_text); // main GUI update

  std::cout << "FINISHED!Scatter correction: CBCT DICOM files are saved"
            << std::endl;
}

void CbctReconTest::test_TempAudit() const {}
void CbctReconTest::test_CalcAndSaveAngularWEPL() {}
void CbctReconTest::test_DoScatterCorrectionUniform() {}
void CbctReconTest::test_FileExportShortDICOM_CurrentImg() {}
void CbctReconTest::test_AddConstHUToCurImg() {}
void CbctReconTest::test_SetCBCTSkinRSPath() {}
void CbctReconTest::test_CropSkinUsingRS() {}
void CbctReconTest::test_CropSkinUsingThreshold() {}
void CbctReconTest::test_ExportAngularWEPL_byFile() {}
void CbctReconTest::test_GeneratePOIData() const {}
void CbctReconTest::test_LoadPOIData() {}
void CbctReconTest::test_StartSyncFromSharedMem() {}
void CbctReconTest::test_StopSyncFromSharedMem() {}
void CbctReconTest::test_TimerEvent() {}
void CbctReconTest::test_ViewExternalCommand() const {}
void CbctReconTest::test_MedianFilterDoNow() {}
void CbctReconTest::test_ExportProjGeometryTXT() {}
void CbctReconTest::test_ForwardProjection() {}
void CbctReconTest::test_FineResolScatterCorrectrionMacro() {}
void CbctReconTest::test_FullScatterCorrectionMacroAP() {}
void CbctReconTest::test_BatchScatterCorrectionMacroAP() {}
void CbctReconTest::test_OpenPhaseData() {}
void CbctReconTest::test_Export4DCBCT() const {}
void CbctReconTest::test_DoCouchCorrection() {}
void CbctReconTest::test_WELPCalcMultipleFiles() {}
void CbctReconTest::test_ScatterCorPerProjRef() {}
void CbctReconTest::test_LoadPerProjRefList() {}
void CbctReconTest::test_CropMaskBatch() {}
void CbctReconTest::test_OutPathEdited() const {}
void CbctReconTest::test_SaveCurrentSetting() const {}
void CbctReconTest::test_CropSupInf() {}
