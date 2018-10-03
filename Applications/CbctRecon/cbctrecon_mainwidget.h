#ifndef CBCTRECON_MAINWIDGET_H
#define CBCTRECON_MAINWIDGET_H

#include "cbctrecon_config.h"

#include <memory>

// Qt
#include <QtWidgets/QMainWindow>
#include <qtimer.h>

#include "cbctrecon.h"


class CBCTRECON_API CbctReconWidget : public QMainWindow {
  Q_OBJECT

  CbctReconWidget(QWidget *parent = nullptr, Qt::WindowFlags flags = nullptr);
  ~CbctReconWidget() = default;

private:
  std::tuple<bool, bool> probeUser(const QString &guessDir);
  FilterReaderType::Pointer ReadBowtieFileWhileProbing(const QString &proj_path, std::tuple<bool, bool> &answers);
  void UpdateReconImage(UShortImageType::Pointer &spNewImg, QString &fileName);
  bool FullScatterCorrectionMacroSingle(QString &outputDirPath,
    enREGI_IMAGES enFwdRefImg,
    bool bFullResolRecon,
    bool bExportImages,
    bool bCBCT_IntensityShift);
  void ForwardProjection(UShortImageType::Pointer &spVolImg3D,
    GeometryType::Pointer &spGeometry,
    UShortImageType::Pointer &spProjCT3D,
    bool bSave, bool use_cuda);

  template<enREGI_IMAGES imagetype>
  void LoadMHAfileAs();
  void FileExportByGUI();
  bool SaveCurrentSetting(QString &strPathConfigFile);
  bool LoadCurrentSetting(QString &strPathConfigFile);

public: 
  std::unique_ptr<CbctRecon> m_cbctrecon;
  std::unique_ptr<QTimer> m_Timer;
  bool m_busyTimer;

public slots:
  void SLT_LoadRawImages(); // independent 2d projection files //not used in
                            // clinical case
  void SLT_Load3DImage();   // indenepndent 3D mha file. UshortFormat. Do
                            // reconstruction is an antoher way to make
                            // m_spReconImg
  void SLT_Load3DImageShort();
  void SLT_LoadPlanCT_mha();
  void SLT_LoadPlanCT_USHORT();
  void SLT_LoadCBCTcorrMHA();
  void SLT_LoadCTrigidMHA();
  void SLT_LoadCTdeformMHA();
  void SLT_LoadNKIImage();
  void SLT_LoadSelectedProjFiles(); // based on presetting values on GUI,
                                    // including geometry files
  void SLT_ReloadProjections();
  void SLT_ExportHis();

  void SLT_LoadImageFloat3D(); // Dose file
  void SLTM_LoadDICOMdir();
  void SLTM_LoadRTKoutput();

  void SLT_DrawRawImages();  // external *.his images
  void SLT_DrawProjImages(); // draw images from HIS FILE READER or filtered
                             // image before going into recon.
  void SLT_DrawReconImage();

  // tools
  void SLT_FileNameHex2Dec();
  void SLT_MakeElektaXML();

  // Gain/ Offset correction
  void SLT_OpenOffsetFile();
  void SLT_OpenGainFile();
  void SLT_OpenBadpixelFile();
  void SLT_ApplyCalibration();

  // Gain/ Offset correction
  void SLT_SetHisDir();
  void SLT_OpenElektaGeomFile();
  void SLT_SetOutputPath();
  void SLT_DoReconstruction();
  // Profile table
  // void SLT_GetProjectionProfile();
  // void SLT_GetReconImgProfile();
  void SLT_CopyTableToClipBoard();
  void SLT_DataProbeProj();
  void SLT_DataProbeRecon();
  void SLT_DrawGraph();
  void SLT_InitializeGraphLim();
  void SLT_UpdateTable();
  void SLT_CalculateROI_Recon();
  void SLT_CalculateROI_Proj();
  void SLT_GoForcedProbePos();
  void SLT_PostApplyFOVDispParam();
  void SLT_DoPostProcessing(); // cropping Circle
  void SLT_PostProcCropInv();
  void SLT_ExportReconUSHORT();
  void SLT_ExportReconSHORT_HU();
  void SLT_ExportALL_DCM_and_SHORT_HU_and_calc_WEPL();
  void SLT_DoBHC();
  void SLT_DoBowtieCorrection();
  void SLT_Export2DDose_TIF();
  void SLTM_Export2DDoseMapAsMHA();
  void SLT_ViewRegistration();
  void SLT_ViewHistogram();
  void SLT_DoScatterCorrection_APRIORI();
  void SLT_TempAudit();
  void SLT_CalcAndSaveAngularWEPL();
  void SLT_DoScatterCorrectionUniform();
  void SLT_FileExportShortDICOM_CurrentImg();
  void SLT_AddConstHUToCurImg();
  void SLT_SetCBCTSkinRSPath();
  void SLT_CropSkinUsingRS();
  void SLT_CropSkinUsingThreshold();
  void SLT_ExportAngularWEPL_byFile();
  void SLT_GeneratePOIData();
  void SLT_LoadPOIData();
  void SLT_StartSyncFromSharedMem();
  void SLT_StopSyncFromSharedMem();
  void SLT_TimerEvent();
  void SLTM_ViewExternalCommand();
  void SLT_MedianFilterDoNow();
  void SLTM_ExportProjGeometryTXT();
  void SLTM_ForwardProjection();
  void SLTM_FineResolScatterCorrectrionMacro(); // projection: full, scatter
                                                // map:512x512

  void SLTM_FullScatterCorrectionMacroAP();
  void SLTM_BatchScatterCorrectionMacroAP();
  void SLT_OpenPhaseData(); // fill lineEdit_PhaseTxtPath
  void SLT_Export4DCBCT();  // phase resorting
  void SLT_DoCouchCorrection();
  void SLTM_WELPCalcMultipleFiles();
  void SLTM_ScatterCorPerProjRef();
  void SLTM_LoadPerProjRefList();
  void SLTM_CropMaskBatch();
  void SLT_OutPathEdited();
  void SLT_SaveCurrentSetting();
  void SLT_CropSupInf();

public:
  Ui::CbctReconClass ui{};

};

#endif // CBCTRECON_MAINWIDGET_H