#ifndef CBCTRECON_H
#define CBCTRECON_H
#include "cbctrecon_config.h"

// std
#include <memory> // unique_, shared_ and weak_ptr

// Qt
#include <QDir>

// Local
#include "AG17RGBAImage.h"
#include "WEPL.h"
#include "YK16GrayImage.h"
#include "cbctrecon_types.h"

class QFileInfo;
class QString;
class QStringList;
class QXmlStreamReader;

class StructureSet;

class CBCTRECON_API CbctRecon {

public:
  CbctRecon();
  ~CbctRecon() = default;

  ITK_DISALLOW_COPY_AND_ASSIGN(CbctRecon);

  // void DoRecon();
  void ReleaseMemory();

  bool FillProjForDisplay(int slice_number);
  void LoadCalibData(std::string &filepath, enCalibType calib_type);
  void RenameFromHexToDecimal(QStringList &filenameList) const;
  QString HexStr2IntStr(QString &str_hex) const;

  std::unique_ptr<YK16GrayImage>
  ApplyCalibrationMaps(YK16GrayImage *const &rawImg, bool DarkCorr,
                       bool GainCorr, bool DefectCorr);
  QString CorrectSingleFile(const char *filePath, bool DarkCorr, bool GainCorr,
                            bool DefectCorr);
  void CorrectSingleFile(YK16GrayImage *pYKRawImg, bool DarkCorr, bool GainCorr,
                         bool DefectCorr);

  void LoadBadPixelMap(const char *filePath);
  // void BadPixReplacement(YK16GrayImage *targetImg);
  std::unique_ptr<YK16GrayImage>
  BadPixReplacement(std::unique_ptr<YK16GrayImage> targetImg);

  void LoadRTKGeometryFile(const char *filePath);

  std::vector<std::string> GetProjFileNames(QString &dirPath);
  bool LoadGeometry(QFileInfo &geomFileInfo, std::vector<std::string> &names);
  std::vector<size_t> GetExcludeProjFiles(bool bManAngleGap,
                                          double gantryAngleInterval);
  void LoadSelectedProj(const std::vector<size_t> &exclude_ids,
                        const std::vector<std::string> &names);
  void saveHisHeader();
  void NormalizeProjections(const FloatImageType::Pointer &reader_output);
  bool ResampleProjections(double &resample_factor);
  void BowtieByFit(bool fullfan, const QStringList &params) const;
  int CropSkinUsingThreshold(int threshold, int erode_radius,
                             int dilate_radius);
  void GeneratePOIData(bool AnteriorToPosterior, double table_posY);
  void Export2DDoseMapAsMHA(QString &strPath) const;
  void ExportProjGeometryTXT(QString &strPath) const;
  void ScatterCorPerProjRef(double scaMedian, double scaGaussian,
                            int postScatMedianSize, bool use_cuda,
                            bool use_opencl, bool save_dicom,
                            FDK_options &fdk_options);

  // void GetSelectedIndices(const std::vector<double>& vFullAngles,
  // std::vector<double>& vNormAngles, std::vector<int>& vTargetIdx, bool bCW);
  void GetExcludeIndexByNames(const QString &outlierListPath,
                              std::vector<std::string> &vProjFileFullPath,
                              std::vector<int> &vExcludeIdx) const;
  void GetSelectedIndices(const std::vector<double> &vFullAngles,
                          std::vector<double> &vNormAngles,
                          std::vector<size_t> &vTargetIdx, bool bCW,
                          std::vector<size_t> &vExcludingIdx) const;

  void SetMaxAndMinValueOfProjectionImage(); // scan m_spProjImg3D and update
                                             // m_fProjImgValueMin, max

  bool IsFileNameOrderCorrect(std::vector<std::string> &vFileNames) const;

  void PostApplyFOVDispParam(float physPosX, float physPosY, float physRadius,
                             float physTablePosY) const;

  // void ExportDICOM_SHORT(SHORT_ImageType::Pointer& sp3DshortImage);//NOT
  // COMPLETED YET!! Export DICOM without Source DICOM is not possible
  void
  CopyDictionary(itk::MetaDataDictionary &fromDict,
                 itk::MetaDataDictionary &toDict) const; // NOT COMPLETED YET!!
                                                         // Export DICOM without
                                                         // Source DICOM is not
                                                         // possible

  void DoBeamHardeningCorrection() const;

  void Draw2DFrom3D(UShortImageType::Pointer &pImg, enPLANE direction,
                    double pos, YK16GrayImage &Output2D) const;
  void Draw2DFrom3DDouble(UShortImageType::Pointer &spFixedImg,
                          UShortImageType::Pointer &spMovingImg,
                          enPLANE enPlane, double pos, YK16GrayImage &YKFixed,
                          YK16GrayImage &YKMoving) const;
  void Draw2DFrom3DDouble(UShortImageType::Pointer &spFixedImg,
                          UShortImageType::Pointer &spMovingImg,
                          enPLANE enPlane, double pos, AG17RGBAImage &YKFixed,
                          AG17RGBAImage &YKMoving) const;

  void RegisterImgDuplication(enREGI_IMAGES src, enREGI_IMAGES target);

  // plastimatch skin / bubble-remover
  // QString getPathCTDir(enMachineType enType);//DICOM Dir
  // QString getPathRS(enMachineType enType);//RS path
  // bool loadPlanCTFromDCM(QString& strCTDirPath, QString& strRSPath);//using
  // plastimatch, prepare m_spRefCTImg. Remove air, RS is needed  Skin will be
  // removed, bubble will be filled

  void FindAllRelevantPaths(const QString &pathProjHisDir);

  template <typename CTImageType, typename ProjImageType>
  void ForwardProjection_master(typename CTImageType::Pointer &spVolImg3D,
                                GeometryType::Pointer &spGeometry,
                                typename ProjImageType::Pointer &spProjCT3D,
                                bool bSave, bool use_cuda);
  // using RTK forward projection algorithm, generate 2D projection image files
  // (as line integral, mu_t)
  template <typename DevFloatImageType>
  void ForwardProjection(UShortImageType::Pointer &spVolImg3D,
                         GeometryType::Pointer &spGeometry,
                         UShortImageType::Pointer &spProjCT3D) const;

  // to be implemented: Save projection3D to *.his files
  void GenScatterMap_PriorCT(UShortImageType::Pointer &spProjRaw3D,
                             UShortImageType::Pointer &spProjCT3D,
                             UShortImageType::Pointer &spProjScat3D,
                             double medianRadius, double gaussianSigma,
                             int nonNegativeScatOffset, bool bSave);
  void ScatterCorr_PrioriCT(UShortImageType::Pointer &spProjRaw3D,
                            UShortImageType::Pointer &spProjScat3D,
                            UShortImageType::Pointer &m_spProjCorr3D,
                            int nonNegativeScatOffset, int postMedian,
                            bool bSave);
  void AfterScatCorrectionMacro(bool use_cuda, bool use_opencl, bool save_dicom,
                                FDK_options &fdk_options);

  // His file export from 3D proj file
  void SaveProjImageAsHIS(UShortImageType::Pointer &spProj3D,
                          std::vector<YK16GrayImage> arrYKImage,
                          QString &strSavingFolder,
                          double resampleF)
      const; // arrYKImage include HIS header and original file name

  static UShortImageType::Pointer
  ConvertLineInt2Intensity(FloatImageType::Pointer &spProjLineInt3D);

  static FloatImageType::Pointer
  ConvertIntensity2LineInt(UShortImageType::Pointer &spProjIntensity3D);

  void Set2DTo3D(FloatImage2DType::Pointer &spSrcImg2D,
                 UShortImageType::Pointer &spTargetImg3D, int idx,
                 enPLANE iDirection) const;

  // void ResampleItkImage(OutputImageType::Pointer& spImgFloat, double
  // resampleF);  Resample proj images
  void ResampleItkImage(FloatImageType::Pointer &spSrcImg,
                        FloatImageType::Pointer &spTarImg,
                        double resFactor) const;
  void ResampleItkImage(UShortImageType::Pointer &spSrcImg,
                        UShortImageType::Pointer &spTarImg,
                        double resFactor) const;
  void ResampleItkImage2D(FloatImage2DType::Pointer &spSrcImg2D,
                          FloatImage2DType::Pointer &spTarImg2D,
                          double resFactor) const; // using slice iterator

  template <enDeviceType Tdev>
  void DoReconstructionFDK(enREGI_IMAGES target,
                           const FDK_options &fdk_options);

  template <enDeviceType Tdev, typename ImageType, typename DDFType,
            typename PSSFType, typename FDKType>
  void DoReconstructionFDK(enREGI_IMAGES target,
                           const FDK_options &fdk_options);

  // void CudaDoReconstructionFDK(enREGI_IMAGES target);
  // void OpenCLDoReconstructionFDK(enREGI_IMAGES target);

  // void SaveUSHORTAsSHORT_DICOM (USHORT_ImageType::Pointer& spImg, QString&
  // strPatientID, QString& strPatientName);//ushort image --> short image -->

  // void GetAngularWEPL_SinglePoint(USHORT_ImageType::Pointer& spImage, int
  // angleGap, VEC3D calcPt, int curPtIdx, std::vector<WEPLData>&
  // vOutputWEPLData, bool bAppend);//output std::vector: append
  void GetAngularWEPL_SinglePoint(UShortImageType::Pointer &spUshortImage,
                                  float fAngleGap, float fAngleStart,
                                  float fAngleEnd, const VEC3D &calcPt,
                                  size_t curPtIdx,
                                  std::vector<WEPLData> &vOutputWEPLData,
                                  bool bAppend) const;
  void GetAngularWEPL_window(UShortImageType::Pointer &spUshortImage,
                             float fAngleGap, float fAngleStart,
                             float fAngleEnd,
                             std::vector<WEPLData> &vOutputWEPLData,
                             bool bAppend);

  //	void UpdateUIAfterLoading(QString& imgName);

  void LoadExternalFloatImage(QString &strPath, bool bConversion);

  void MedianFilterByGUI(const UShortImageType::SizeType
                             &indexRadius); // params are given at the UI

  /*Temporary implementation for XVI5 xml*/
  bool
  LoadXVIGeometryFile(const char *filePath); // temporary implenetation
                                             // using QT XML. This is for
                                             // XVI v >5.0.2. _Frames.xml is
                                             // in every projection folder

  void SetProjDir(QString &strProjPath);

  void ExportAngularWEPL_byFile(QString &strPathOutput, double fAngleStart,
                                double fAngleEnd, double fAngleGap);

  /*Temporary implementation for XVI5 xml*/

  // using itkimage and cuda do median filter faster. Output will be overwritten
  //	void cudaMedianFilter2DITK(OutputImageType2D::Pointer& spUshortImage,
  // int wndSizeX, int wndSizeY);

  // double CropSkinUsingRS(USHORT_ImageType::Pointer& spImgUshort, QString&
  // strPathRS, double cropMargin);

  // void AuditMemory();

  void CropSupInf(UShortImageType::Pointer &sp_Img, float physPosInfCut,
                  float physPosSupCut);
  void CropFOV3D(UShortImageType::Pointer &sp_Img, float physPosX,
                 float physPosY, float physRadius, float physTablePosY) const;

  void GenerateCylinderMask(UShortImageType::Pointer &spImgCanvas,
                            float fDcmPosX, float fDcmPosY,
                            float fRadius) const;

  float GetMeanIntensity(UShortImageType::Pointer &spImg, float sphereR,
                         float *sdIntensity = nullptr) const;

  bool ResortCBCTProjection(std::vector<int> &vIntPhaseBinSelected,
                            QString &strPathForXML, QString &strPathProjRoot,
                            QString &strUID,
                            std::vector<float> &vFloatPhaseFull,
                            GeometryType::Pointer &spGeomFull,
                            std::vector<std::string> &vProjPathsFull) const;

  void AppendInPhaseIndex(int iPhase, std::vector<float> &vFloatPhaseFull,
                          std::vector<size_t> &vOutputIndex,
                          int margin = 5) const;

  void LoadShort3DImage(QString &filePath, enREGI_IMAGES enTarget);

  void GetWEPLDataFromSingleFile(const QString &filePath,
                                 std::vector<VEC3D> &vPOI,
                                 std::vector<WEPLData> &vOutputWEPL,
                                 double fAngleStart, double fAngleEnd) const;

  template <typename DevImageType>
  void SingleForwardProjection(FloatImageType::Pointer &spVolImgFloat,
                               float fMVGanAngle, float panelOffsetX,
                               float panelOffsetY,
                               UShortImageType::Pointer &spProjImg3D,
                               int iSliceIdx) const;

  bool ReadDicomDir(QString &dirPath);

  // using RTK forward projection algorithm, generate 2D projection image files
  // (as line integral, mu_t)

  // still public:
  std::unique_ptr<StructureSet> m_structures;
  // independent raw images
  std::vector<YK16GrayImage> m_arrYKImage;
  int m_iImgCnt{}; // for independent raw images --> no relation to Directroy
                   // based projections

  std::vector<YK16GrayImage> m_arrYKBufProj;
  size_t m_iCntSelectedProj;

  std::unique_ptr<YK16GrayImage> m_pImgOffset;
  std::unique_ptr<YK16GrayImage> m_pImgGain;
  // Badpixmap;
  std::vector<BADPIXELMAP> m_vPixelReplMap{};

  // RTK recon
  GeometryType::Pointer m_spFullGeometry; // sp = smart pointer
  GeometryType::Pointer m_spCustomGeometry;

  enProjFormat m_projFormat = HIS_FORMAT;
  bool m_bScanDirectionCW;

  FloatImageType::Pointer
      m_spProjImg3DFloat; // This is float image loaded by
                          // RTK. line integral (mu_t value).
                          // To convert this to Intensity,
                          // use mu t = ln(65535/I)

  UShortImageType::Pointer m_spProjImgRaw3D; // raw intensity value converted
                                             // from line integral (mu_t).
                                             // 0-65535
  UShortImageType::Pointer
      m_spProjImgCT3D; // 1.5G // release this after Scatter Generation
  UShortImageType::Pointer m_spProjImgScat3D; // scatter map proj file using any
                                              // scatter-estimation method//1.5G
                                              // //release this after Cor gen
  UShortImageType::Pointer m_spProjImgCorr3D; // proj file-scatter corrected one
                                              // using either priori CT or any
                                              // other method//1.5G

  UShortImageType::Pointer m_spCrntReconImg; // fixed image // ID: RawCBCT
  UShortImageType::Pointer m_spRawReconImg;  // just added --> when file is
                                             // loaded
  UShortImageType::Pointer
      m_spScatCorrReconImg; // just added --> after scatter correction

  UShortImageType::Pointer
      m_spRefCTImg; // filled by SLT_LoadPlanCT_mha(); ID: RefCT_Original
  UShortImageType::Pointer m_spManualRigidCT; // copied from RefCTImg; ID: RefCT
                                              // --> Moving Img, cloned
  UShortImageType::Pointer m_spAutoRigidCT;   // ID: AutoRigidCT
  UShortImageType::Pointer m_spDeformedCT1; // Deformmation will be carried out
                                            // based on Moving IMage of GUI
                                            // //AutoDeformCT1
  UShortImageType::Pointer m_spDeformedCT2; // AutoDeformCT2
  UShortImageType::Pointer m_spDeformedCT3; // AutoDeformCT3
  UShortImageType::Pointer m_spDeformedCT_Final; // AutoDeformCT3

  std::unique_ptr<YK16GrayImage> m_dspYKReconImage{};
  std::unique_ptr<YK16GrayImage> m_dspYKImgProj{};
  int m_iTmpIdx;

  double m_fProjImgValueMax; // value of float image
  double m_fProjImgValueMin;

  double m_multiplyFactor{};

  // std::unique_ptr<DlgRegistration> m_pDlgRegistration;
  // DlgHistogram* m_pDlgHistogram;
  // std::unique_ptr<DlgExternalCommand> m_pDlgExternalCommand;

  // Automatically detected relavant file/Dir path when load the projection
  // files (SLT_SetHisDir)

  // Belows are ELEKTA specific
  // 1) Find DICOM UID by subtracting "img_" from whole proj path name

  // Below paths will be decided after the Find... Func.
  QString m_strDCMUID;
  QString m_strPathPatientDir; // full path of patient Directory
  QString m_strPatientDirName; // just the name --> later I can extract the
                               // patient ID from here
  QString m_strPathFRAME_DBF;
  QString m_strPathIMAGE_DBF;
  QString m_strPathGeomXML; // after Generation of the XML from DBF files
  QString m_strPathPlanCTDir;
  QString m_strPathRS;            // for body and lung contours
  QString m_strPathPlan;          // for isocenter position
  QString m_strPathDirDefault;    // QFileDialog default starting point
  QString m_strPathRS_CBCT;       // QFileDialog default starting point
  QString m_strPathElektaINI;     // for mAs values
  QString m_strPathIMAGES;        // upper folder of projection files (His)
  QString m_strPathElektaINIXVI2; // this includes couch shift values. longer
                                  // INI.XVI file
  QString m_strCur_mAs;           // QString("20,20")
  QString m_strRef_mAs;           // QString("64,40")
  QString m_strError;

  QDir m_dcm_dir;

  int m_iFixedOffset_ScatterMap; // fixed! allows negative value of scatter
  double m_fResampleF; // typically 0.5. this is updated during LoadSelectedProj
                       // image and ui.lineEdit_DownResolFactor.//also affects
                       // all other scatter correction method
  // double m_fProjSpacingX; // updated from SelectedProjLoad
  // double m_fProjSpacingY; // updated from SelectedProjLoad

  std::vector<VEC3D> m_vPOI_DCM{}; // initialized by file Load

  std::vector<std::string> m_vSelectedFileNames{};

  bool m_bMacroContinue;

  std::vector<float> m_vPhaseFloat{};

  QStringList m_strListPerProjRefVol;

  QString m_strPathDefaultConfigFile;

  std::vector<int> m_vExcludeProjIdx{}; // if kVON (exposed_ tag is false

  // private:
};

#ifndef ITK_MANUAL_INSTANTIATION
#include "cbctrecon_fdk.hxx"
#endif

#endif // CBCTRECON_H
