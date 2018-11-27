// For testing CbctRecon

#ifdef USE_TINYREFL
#include <tinyrefl/api.hpp> // MUST BE INCLUDED FIRST, FFS!

#include "cbctrecon.h"
#include "cbctrecon.h.tinyrefl"
#else
#include "cbctrecon.h"
#endif

#include <iostream>
#include <chrono>
#include <memory>
#include <thread>

#include <QDir>

#include "cbctrecon_test.hpp"
#include "cbctregistration.h"
#include "cbctrecon_io.h"

CbctReconTest::CbctReconTest() {
  m_cbctregistration = std::make_unique<CbctRegistration>(m_cbctrecon.get());
  m_pTableModel = nullptr;
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

FilterReaderType::Pointer
CbctReconTest::ReadBowtieFileWhileProbing(const QString &proj_path,
                                          std::tuple<bool, bool> &answers) const {

  auto bowtiereader =
      FilterReaderType::New(); // we use is because we need the projections to
                               // be in the same unit (order of magnitude)

  QDir guessDir(proj_path + QString("/../"));

  const auto calDir(proj_path + QString("/Calibrations/"));

  QString bowtiePath;
  answers = std::make_tuple(true, true);

  switch (this->m_cbctrecon->m_projFormat) {
  case XIM_FORMAT:
    bowtiePath = getBowtiePath(calDir);
    if (bowtiePath.length() > 1) {
      std::cerr << "loading bowtie-filter..." << "\n";
      std::vector<std::string> filepath;
      filepath.push_back(bowtiePath.toStdString());
      bowtiereader->SetFileNames(filepath);
      // std::thread calc_thread_bowtie(read_bowtie_projection, bowtiereader);
      std::thread calc_thread_bowtie(
          [&bowtiereader] { bowtiereader->Update(); });
      calc_thread_bowtie.join();
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

void CbctReconTest::test_LoadSelectedProjFiles(const QString &proj_path, const QString& geom_file)
{
  // this->ui.pushButton_DoRecon->setDisabled(true);
  // 1) Get all projection file names
  auto dirPath = proj_path; // this->ui.lineEdit_HisDirPath->text();
  //.toLocal8Bit().constData();

  if (!QFile::exists(dirPath)) {
    std::cerr << "Projection file directory was not found. Retry." << "\n";
    return;
  }

  auto names = this->m_cbctrecon->GetProjFileNames(dirPath);

  if (!this->m_cbctrecon->IsFileNameOrderCorrect(names) &&
      this->m_cbctrecon->m_projFormat != XIM_FORMAT) {
    std::cerr << "Check the file name order" << "\n";
    return;
  }

  std::cerr << "File name order was cross-checked and found to be OK!"
            << "\n";

  const auto fullCnt = names.size();
  if (fullCnt <= 0) {
    std::cerr << "No projection file was found. Retry." << "\n";
    return;
  }

  std::cerr << fullCnt << "  projection files were found." << "\n";

  // 2) Elekta Geometry file
  QFileInfo geomFileInfo(geom_file);

  if (!this->m_cbctrecon->LoadGeometry(geomFileInfo, names)) {
    if (!this->m_cbctrecon->m_strError.isEmpty()) {
      std::cerr << this->m_cbctrecon->m_strError.toStdString() << "\n";
    }
  }

  const auto iFullGeoDataSize =
      this->m_cbctrecon->m_spFullGeometry->GetGantryAngles().size();
  if (iFullGeoDataSize < 1) {
    std::cerr << "Not enough projection image (should be > 0)" << "\n";
    return;
  }

  if (iFullGeoDataSize != fullCnt) {
    if (this->m_cbctrecon->m_projFormat != XIM_FORMAT) {
      std::cerr << "Size of geometry data and file numbers are not same! Check "
                   "and retry"
                << "\n";
      return;
    }

    const auto reply = true; /*QMessageBox::question(
        this, "Mismatch in number of files and Geometry information!",
        "Mismatch in number of files and Geometry information!\nHowever, Xim "
        "detected, so it may be safe to continue anyway?",
        QMessageBox::Yes | QMessageBox::No);*/
        if (reply == true) { // QMessageBox::Yes) {
      std::cerr << "continuing despite warning..." << "\n";
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

  std::cerr << "AngularGaps Sum (deg):" << sum_gap
            << ", Mean (deg): " << mean_gap << "\n";

  const auto gantryAngleInterval = 1.0;
      // this->ui.lineEdit_ManualProjAngleGap->text().toDouble();

  if (false){ // this->ui.Radio_ManualProjAngleGap->isChecked()) {

    if (gantryAngleInterval < mean_gap) {
      std::cerr << "Angle gap size is too small. Terminating the app"
                << "\n";
      return;
    }
  }

  const auto exclude_ids = this->m_cbctrecon->GetExcludeProjFiles(
      false /*this->ui.Radio_ManualProjAngleGap->isChecked()*/, gantryAngleInterval);

  this->m_cbctrecon->LoadSelectedProj(exclude_ids, names);

  // Reads the cone beam projections
  using ReaderType = rtk::ProjectionsReader<FloatImageType>;
  auto reader = ReaderType::New();
  reader->SetFileNames(this->m_cbctrecon->m_vSelectedFileNames);
  std::thread calc_thread([&reader]() { reader->Update(); });

  std::cerr << "Reader detached from main thread" << "\n";

  // After reading the whole file,
  // HIS header should be saved
  this->m_cbctrecon->saveHisHeader();

  //  Insta Recon, Dcm read
  const auto geopath = geomFileInfo.absolutePath();
  std::tuple<bool, bool> answers;
  auto bowtie_reader = ReadBowtieFileWhileProbing(geopath, answers);

  calc_thread.join();
  std::cerr << "Reader re-attached to main thread" << "\n";

  this->m_cbctrecon->m_spProjImg3DFloat =
      reader->GetOutput(); // 1024 1024, line integ image

  if (bowtie_reader != nullptr) {
    ApplyBowtie(reader, bowtie_reader);
  }
  if (this->m_cbctrecon->m_projFormat == HND_FORMAT) {
    std::cerr << "Fitted bowtie-filter correction ongoing..." << "\n";
    test_DoBowtieCorrection();
  }

  saveImageAsMHA<FloatImageType>(this->m_cbctrecon->m_spProjImg3DFloat);
  auto res_factor = 0.5;
   // this->ui.lineEdit_DownResolFactor->text().toDouble();
  if (!this->m_cbctrecon->ResampleProjections(res_factor)) { // 0.5
    // reset factor if image was not resampled
    std::cerr << "Could not resample projection size!\n";
    // this->ui.lineEdit_DownResolFactor->setText("1.0");
  }

  this->m_cbctrecon->ConvertLineInt2Intensity(
      this->m_cbctrecon->m_spProjImg3DFloat,
      this->m_cbctrecon->m_spProjImgRaw3D,
      65535); // if X not 1024 == input size: out_offset =
              // in_offset + (1024*res_f -
              // X*res_f)*out_spacing     <- will still
              // break down at fw_projection

  this->m_cbctrecon
      ->SetMaxAndMinValueOfProjectionImage(); // update min max projection image

  test_InitializeGraphLim();

  this->test_DrawProjImages(); // Update Table is called

  if (!std::get<0>(answers)) { // instaRecon
    std::cerr
        << "FINISHED!: Loading projection files. Proceed to reconstruction"
        << "\n";
  } else {
    test_DoReconstruction();
  }

  if (std::get<1>(answers)) { // CT DCM dir was found
    test_ViewRegistration();
  }
}

void CbctReconTest::test_ReloadProjections() {}
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
void CbctReconTest::test_SetHisDir() {}
void CbctReconTest::test_OpenElektaGeomFile() {}
void CbctReconTest::test_SetOutputPath() {}
void CbctReconTest::test_DoReconstruction() {}
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
void CbctReconTest::test_DoBowtieCorrection() {}
void CbctReconTest::test_Export2DDose_TIF() {}
void CbctReconTest::test_Export2DDoseMapAsMHA() {}
void CbctReconTest::test_ViewRegistration() const {}
void CbctReconTest::test_ViewHistogram() {}
void CbctReconTest::test_DoScatterCorrection_APRIORI() {}
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

int main(const int argc, char *argv[]) {

  if (argc < 3) {
    std::cerr << "Usage:\n" << argv[0] << " ./dicom/directory ./CB_proj/directory\n";
    return -1;
  }

  std::cerr << "Running cbctrecon_test!\n";
  auto cbctrecon_test = std::make_unique<CbctReconTest>();

  auto dcm_dir = QDir(argv[1]);
  auto dcm_path = dcm_dir.absolutePath();
  if (!dcm_dir.exists()) {
    std::cerr << "Directory didn't exist: " << dcm_path.toStdString() << "\n";
    return -2;
  }
  if (dcm_dir.isEmpty(QDir::AllEntries | QDir::NoDotAndDotDot)) {
    std::cerr << "Directory was empty: " << dcm_path.toStdString() << "\n";
    return -3;
  }

  cbctrecon_test->m_cbctrecon->m_strPathDirDefault = dcm_path;
  cbctrecon_test->test_LoadDICOMdir();
  if (cbctrecon_test->m_cbctrecon->m_spManualRigidCT.IsNull()) {
    std::cerr << "Manual Rigid CT was NULL -> Dicom dir was not read!\n";
    return -4;
  }
  auto ss = cbctrecon_test->m_cbctrecon->m_structures->get_ss(PLAN_CT);
  for (auto &structure : ss->slist) {
    std::cerr << structure.name << "\n";
  }

  const auto voi = std::string("CTV1");
  auto start_time = std::chrono::steady_clock::now();
  cbctrecon_test->m_cbctregistration->CalculateWEPLtoVOI(
      voi, 45, 45, cbctrecon_test->m_cbctrecon->m_spManualRigidCT);
  auto end_time = std::chrono::steady_clock::now();
  std::cerr << "WEPL was calculated in: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                     start_time)
                   .count()
            << " ms" << "\n";

  /* Some verification of the WEPL results should go here */

  /* Load projections (Needs to be uploaded to girder first) */
  auto cbct_dir = QDir(argv[2]);
  auto cbct_path = cbct_dir.absolutePath();
  if (!cbct_dir.exists()) {
    std::cerr << "Directory didn't exist: " << cbct_path.toStdString() << "\n";
    return -2;
  }
  if (cbct_dir.isEmpty(QDir::AllEntries | QDir::NoDotAndDotDot)) {
    std::cerr << "Directory was empty: " << cbct_path.toStdString() << "\n";
    return -3;
  }
  const auto geom_file = cbct_path + "/Scan.xml";
  start_time = std::chrono::steady_clock::now();
  cbctrecon_test->test_LoadSelectedProjFiles(cbct_path, geom_file);
  end_time = std::chrono::steady_clock::now();
  std::cerr << "Proj. was loaded and reconstructed in: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
                                                                     start_time)
                   .count()
            << " ms" << "\n";


  /* Scatter correction algorithm "Batch" style */

  return 0;
}
