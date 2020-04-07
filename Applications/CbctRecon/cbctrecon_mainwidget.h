#ifndef CBCTRECON_MAINWIDGET_H
#define CBCTRECON_MAINWIDGET_H

#include <filesystem>
#include <memory>
#include <tuple>

// Qt
#include <QStandardItemModel>
#include <qstring.h>
#include <qtimer.h>

#include "cbctrecon.h"
#include "cbctregistration.h"

#include "DlgExternalCommand.h"
#include "DlgHistogram.h"
#include "DlgRegistration.h"

#include "ui_cbctrecon.h"

namespace fs = std::filesystem;

class CbctReconWidget : public QMainWindow {
  Q_OBJECT

public:
  explicit CbctReconWidget(QWidget *parent = nullptr,
                           Qt::WindowFlags flags = nullptr);
  //~CbctReconWidget() = default;
  void UpdateReconImage(UShortImageType::Pointer &spNewImg, const QString& fileName);
  FDK_options getFDKoptions() const;

private:
  std::tuple<bool, bool> probeUser(const fs::path &guessDir);
  FilterReaderType::Pointer
  ReadBowtieFileWhileProbing(const fs::path &proj_path,
                             std::tuple<bool, bool> &answers);
  bool FullScatterCorrectionMacroSingle(const fs::path &outputDirPath,
                                        enREGI_IMAGES enFwdRefImg,
                                        bool bFullResolRecon,
                                        bool bExportImages,
                                        bool bCBCT_IntensityShift);

  template <enREGI_IMAGES imagetype> void LoadMHAfileAs();
  bool SaveCurrentSetting(const fs::path &strPathConfigFile) const;
  bool LoadCurrentSetting(const fs::path &strPathConfigFile) const;
  void init_DlgRegistration(std::string &str_dcm_uid) const;

public:
  std::unique_ptr<CbctRecon> m_cbctrecon;
  std::unique_ptr<DlgRegistration> m_dlgRegistration;
  CbctRegistration *m_cbctregistration; // just for convienience
  std::unique_ptr<DlgExternalCommand> m_dlgExternalCommand;
  std::unique_ptr<DlgHistogram> m_dlgHistogram;
  std::unique_ptr<QTimer> m_Timer;
  std::unique_ptr<QStandardItemModel> m_pTableModel;
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

  void SLT_DrawRawImages() const; // external *.his images
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
  void SLT_ApplyCalibration() const;

  // Gain/ Offset correction
  void SLT_SetHisDir();
  void SLT_OpenElektaGeomFile();
  void SLT_SetOutputPath();
  void SLT_DoReconstruction();
  // Profile table
  // void SLT_GetProjectionProfile();
  // void SLT_GetReconImgProfile();
  void SLT_CopyTableToClipBoard() const;
  void SLT_DataProbeProj() const;
  void SLT_DataProbeRecon() const;
  void SLT_DrawGraph() const;
  void SLT_InitializeGraphLim() const;
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
  void SLT_ViewRegistration() const;
  void SLT_ViewHistogram() const;
  void SLT_DoScatterCorrection_APRIORI();
  void SLT_CalcAndSaveAngularWEPL();
  void SLT_DoScatterCorrectionUniform();
  void SLT_FileExportShortDICOM_CurrentImg();
  void SLT_AddConstHUToCurImg();
  void SLT_CropSkinUsingRS();
  void SLT_CropSkinUsingThreshold();
  void SLT_ExportAngularWEPL_byFile();
  void SLT_GeneratePOIData() const;
  void SLT_LoadPOIData();
  static void SLT_StartSyncFromSharedMem();
  static void SLT_StopSyncFromSharedMem();
  void SLT_TimerEvent();
  void SLTM_ViewExternalCommand() const;
  void SLT_MedianFilterDoNow();
  void SLTM_ExportProjGeometryTXT();
  void SLTM_ForwardProjection();
  void SLTM_FineResolScatterCorrectrionMacro(); // projection: full, scatter
                                                // map:512x512

  void SLTM_FullScatterCorrectionMacroAP();
  void SLTM_BatchScatterCorrectionMacroAP();
  void SLT_OpenPhaseData();      // fill lineEdit_PhaseTxtPath
  void SLT_Export4DCBCT() const; // phase resorting
  void SLT_DoCouchCorrection();
  void SLTM_WELPCalcMultipleFiles();
  void SLTM_ScatterCorPerProjRef();
  void SLTM_LoadPerProjRefList();
  void SLTM_CropMaskBatch();
  void SLT_OutPathEdited() const;
  void SLT_SaveCurrentSetting() const;
  void SLT_CropSupInf();

public:
  Ui::CbctReconClass ui{};
};

#endif // CBCTRECON_MAINWIDGET_H
