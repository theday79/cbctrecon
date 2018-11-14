// For testing CbctRecon

#ifdef USE_TINYREFL
#include <tinyrefl/api.hpp> // MUST BE INCLUDED FIRST, FFS!

#include "cbctrecon.h"
#include "cbctrecon.h.tinyrefl"
#else
#include "cbctrecon.h"
#endif

#include <memory>

#include <QDir>

#include "cbctrecon_test.hpp"
#include "cbctregistration.h"

CbctReconTest::CbctReconTest() {
  m_cbctregistration = std::make_unique<CbctRegistration>(m_cbctrecon.get());
  m_pTableModel = nullptr;
}

void CbctReconTest::test_LoadRawImages() {}
void CbctReconTest::test_Load3DImage() {}
void CbctReconTest::test_Load3DImageShort() {}
void CbctReconTest::test_LoadPlanCT_mha() {}
void CbctReconTest::test_LoadPlanCT_USHORT() {}
void CbctReconTest::test_LoadCBCTcorrMHA() {}
void CbctReconTest::test_LoadCTrigidMHA() {}
void CbctReconTest::test_LoadCTdeformMHA() {}
void CbctReconTest::test_LoadNKIImage() {}
void CbctReconTest::test_LoadSelectedProjFiles() {}
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

  if (argc < 2) {
    std::cerr << "Usage:\n" << argv[0] << " ./dicom/directory\n";
    return -1;
  }

  std::cout << "Running cbctrecon_test!\n";
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
  /*try { // This will have to wait, unfortunately
    cbctrecon->ReadDicomDir(dcm_path);
  }
  catch (std::exception& e) {
    std::cerr << "Couldn't read DICOM: " << e.what() << std::endl;
    return 1;
  }

  std::cout << "image was read" << std::endl;
  */
  /*cbctrecon->m_pDlgRegistration->UpdateVOICombobox(PLAN_CT);
  cbctrecon->m_pDlgRegistration->UpdateListOfComboBox(1);
  cbctrecon->m_pDlgRegistration->LoadImgFromComboBox(1, QString("REF_CT"));
  cbctrecon->m_pDlgRegistration->ui.comboBox_VOI->setCurrentIndex(1);

  std::cout << "DlgRegi. is ready for calculating WEPL for "
    <<
  cbctrecon->m_pDlgRegistration->ui.comboBox_VOI->currentText().toStdString()
    << std::endl;
  auto start_time = std::chrono::steady_clock::now(); //clock();
  cbctrecon->m_pDlgRegistration->SLT_WEPLcalc();
  auto end_time = std::chrono::steady_clock::now();

  std::cout << "WEPL was calculated in: "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end_time -
  start_time).count()
    << " ms"
    << std::endl;*/

  return 0;
}
