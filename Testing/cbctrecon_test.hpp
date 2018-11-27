#ifndef CBCTRECON_TEST_HPP
#define CBCTRECON_TEST_HPP

/* Defines necessary for using TinyRefl */

// #define TINYREFL_API_CODEGEN_VERSION_MAJOR 0
// #define TINYREFL_API_CODEGEN_VERSION_MINOR 1
// #define TINYREFL_API_CODEGEN_VERSION_FIX   1

#include <memory>
#include <tuple>

// Qt
#include <QStandardItemModel>

#include "cbctrecon.h"
#include "cbctregistration.h"

class CbctReconTest {
public:
  CbctReconTest();
  //~CbctReconTest() = default;

  static FDK_options getFDKoptions() {
    return FDK_options{}; // Uses default options
  }

private:
  bool FullScatterCorrectionMacroSingle(QString &outputDirPath,
                                        enREGI_IMAGES enFwdRefImg,
                                        bool bFullResolRecon,
                                        bool bExportImages,
                                        bool bCBCT_IntensityShift);

public:
  std::unique_ptr<CbctRecon> m_cbctrecon = std::make_unique<CbctRecon>();
  std::unique_ptr<CbctRegistration> m_cbctregistration; // just for convienience
  std::unique_ptr<QStandardItemModel> m_pTableModel;

  // still public:
  void test_LoadRawImages(); // independent 2d projection files //not used in
                             // clinical case
  void test_Load3DImage();   // indenepndent 3D mha file. UshortFormat. Do
                             // reconstruction is an antoher way to make
                             // m_spReconImg
  void test_Load3DImageShort();
  void test_LoadPlanCT_mha();
  void test_LoadPlanCT_USHORT();
  void test_LoadCBCTcorrMHA();
  void test_LoadCTrigidMHA();
  void test_LoadCTdeformMHA();
  void test_LoadNKIImage();
  void test_LoadSelectedProjFiles(); // based on presetting values on GUI,
                                     // including geometry files
  void test_ReloadProjections();
  void test_ExportHis();

  void test_LoadImageFloat3D(); // Dose file
  void test_LoadDICOMdir() const;
  void test_LoadRTKoutput();

  void test_DrawRawImages() const; // external *.his images
  void test_DrawProjImages(); // draw images from HIS FILE READER or filtered
                              // image before going into recon.
  void test_DrawReconImage();

  // tools
  void test_FileNameHex2Dec();
  void test_MakeElektaXML();

  // Gain/ Offset correction
  void test_OpenOffsetFile();
  void test_OpenGainFile();
  void test_OpenBadpixelFile();
  void test_ApplyCalibration() const;

  // Gain/ Offset correction
  void test_SetHisDir();
  void test_OpenElektaGeomFile();
  void test_SetOutputPath();
  void test_DoReconstruction();
  // Profile table
  // void test_GetProjectionProfile();
  // void test_GetReconImgProfile();
  void test_CopyTableToClipBoard() const;
  void test_DataProbeProj() const;
  void test_DataProbeRecon() const;
  void test_DrawGraph() const;
  void test_InitializeGraphLim() const;
  void test_UpdateTable();
  void test_CalculateROI_Recon();
  void test_CalculateROI_Proj();
  void test_GoForcedProbePos();
  void test_PostApplyFOVDispParam();
  void test_DoPostProcessing(); // cropping Circle
  void test_PostProcCropInv();
  void test_ExportReconUSHORT();
  void test_ExportReconSHORT_HU();
  void test_ExportALL_DCM_and_SHORT_HU_and_calc_WEPL();
  void test_DoBHC();
  void test_DoBowtieCorrection();
  void test_Export2DDose_TIF();
  void test_Export2DDoseMapAsMHA();
  void test_ViewRegistration() const;
  void test_ViewHistogram();
  void test_DoScatterCorrection_APRIORI();
  void test_TempAudit() const;
  void test_CalcAndSaveAngularWEPL();
  void test_DoScatterCorrectionUniform();
  void test_FileExportShortDICOM_CurrentImg();
  void test_AddConstHUToCurImg();
  void test_SetCBCTSkinRSPath();
  void test_CropSkinUsingRS();
  void test_CropSkinUsingThreshold();
  void test_ExportAngularWEPL_byFile();
  void test_GeneratePOIData() const;
  void test_LoadPOIData();
  void test_StartSyncFromSharedMem();
  static void test_StopSyncFromSharedMem();
  void test_TimerEvent();
  void test_ViewExternalCommand() const;
  void test_MedianFilterDoNow();
  void test_ExportProjGeometryTXT();
  void test_ForwardProjection();
  void test_FineResolScatterCorrectrionMacro(); // projection: full, scatter
                                                // map:512x512

  void test_FullScatterCorrectionMacroAP();
  void test_BatchScatterCorrectionMacroAP();
  void test_OpenPhaseData();      // fill lineEdit_PhaseTxtPath
  void test_Export4DCBCT() const; // phase resorting
  void test_DoCouchCorrection();
  void test_WELPCalcMultipleFiles();
  void test_ScatterCorPerProjRef();
  void test_LoadPerProjRefList();
  void test_CropMaskBatch();
  void test_OutPathEdited() const;
  void test_SaveCurrentSetting() const;
  void test_CropSupInf();
};

#endif // CBCTRECON_TEST_HPP
