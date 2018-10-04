#ifndef CBCTRECON_H
#define CBCTRECON_H
#include "cbctrecon_config.h"

// std
#include <array>
#include <memory> // unique_, shared_ and weak_ptr


// configs
#include <plm_config.h> // first due to lack of ifdef guards
#include <rtkConfiguration.h> // also includes itkConfiguration.h

// ITK
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#ifdef USE_CUDA
#include <itkCudaImage.h>
#endif

// RTK
#include <rtkThreeDCircularProjectionGeometry.h>

// PLM
#include <itk_mask.h>

// Local
#include "ui_cbctrecon.h"
#include "YK16GrayImage.h"
#include "WEPL.h"
#include "cbctrecon_io.h"
#include "cbctrecon_compute.h"

#if USE_CUDA
using CUDAFloatImageType = itk::CudaImage<float, 3U>;
#endif // USE_CUDA

using FloatPixelType = float;
using USHORT_PixelType = unsigned short;
using SHORT_PixelType = short;

using FloatImageType = itk::Image<FloatPixelType, 3U>;
using FloatImage2DType = itk::Image<FloatPixelType, 2U>;
using UShortImageType = itk::Image<USHORT_PixelType, 3U>;
using ShortImageType = itk::Image<SHORT_PixelType, 3U>;

using FloatReaderType = itk::ImageFileReader<FloatImageType>;
using FloatWriterType = itk::ImageFileWriter<FloatImageType>;

using GeometryType = rtk::ThreeDCircularProjectionGeometry;

using FilterReaderType =
rtk::ProjectionsReader<FloatImage2DType>; // This one does a bit more than
                                          // we need, but mainly
                                          // statistical filters
using ProjReaderType = rtk::ProjectionsReader<FloatImageType>;

#define DEFAULT_ELEKTA_PROJ_WIDTH 1024
#define DEFAULT_ELEKTA_PROJ_HEIGHT 1024
#define DEFAULT_VARIAN_PROJ_WIDTH 1024
#define DEFAULT_VARIAN_PROJ_HEIGHT 768
#define MAX_LINE_LENGTH 1024

// lineEdit_scaMedian, when downsampling =1.0
#define DEFAULT_SCA_MEDIAN 25.0
// lineEdit_scaGaussian, when downsampling =1.0
#define DEFAULT_SCA_GAUSSIAN 3.0
// lineEdit_sca, when downsampling =1.0
#define DEFAULT_SCA_POST_PROJ_MEDIAN 6.0

enum enPLANE {
  PLANE_AXIAL = 0,
  PLANE_FRONTAL,
  PLANE_SAGITTAL,
};

enum enREGI_IMAGES {
  REGISTER_RAW_CBCT = 0,
  REGISTER_REF_CT,       // manual moving image
  REGISTER_MANUAL_RIGID, // manual moving image
  REGISTER_AUTO_RIGID,
  REGISTER_DEFORM1,
  REGISTER_DEFORM2,
  REGISTER_DEFORM3,
  REGISTER_DEFORM_FINAL,
  REGISTER_COR_CBCT,
  REGISTER_DEFORM_SKIP_AUTORIGID,
};

enum enMachineType {
  MACHINE_ELEKTA = 0,
  MACHINE_VARIAN,
};

enum enProjFormat {
  HIS_FORMAT,
  HND_FORMAT,
  XIM_FORMAT,
};

enum enCalibType {
  GAIN_CALIB,
  OFFSET_CALIB,
  BADPIXEL_CALIB,
};

enum FWD_METHOD {
  en_Joseph = 0,
  en_CudaRayCast,
  // en_RayCastInterpolator, Deprecated in rtk 1.4
};

struct WEPLData {
  double fWEPL;
  int ptIndex;
  double fGanAngle;
};
struct VEC3D {
  double x;
  double y;
  double z;
};

//
struct FLEXDATA {
  float fGanAngle;     // MV beam gantry angle
  float fPanelOffsetX; // MV beam gantry angle
  float fPanelOffsetY; // MV beam gantry angle
  bool bKV_On;
  bool bMV_On;
};

class QStandardItemModel;

//class DlgRegistration;
// class DlgHistogram;
//class DlgExternalCommand;
class StructureSet;
class QTimer;
class QXmlStreamReader;

// using namespace std;

class CBCTRECON_API CbctRecon {

  //friend class TestCbctRecon;

public:
  CbctRecon();
  ~CbctRecon();
  // void DoRecon();
  void ReleaseMemory();

  bool FillProjForDisplay(const int slice_number);
  void LoadCalibData(std::string filepath, enCalibType calib_type);
  void RenameFromHexToDecimal(QStringList &filenameList);
  QString HexStr2IntStr(QString &strHex);

  std::unique_ptr<YK16GrayImage> ApplyCalibrationMaps(YK16GrayImage* const &rawImg, bool DarkCorr, bool GainCorr, bool DefectCorr);
  QString CorrectSingleFile(const char *filePath, bool DarkCorr, bool GainCorr, bool DefectCorr);
  void CorrectSingleFile(YK16GrayImage *pYKRawImg, bool DarkCorr, bool GainCorr, bool DefectCorr);

  void LoadBadPixelMap(const char *filePath);
  //void BadPixReplacement(YK16GrayImage *targetImg);
  std::unique_ptr<YK16GrayImage> BadPixReplacement(std::unique_ptr<YK16GrayImage> targetImg);

  void LoadRTKGeometryFile(const char *filePath);

  std::vector<std::string> GetProjFileNames(QString dirPath);
  bool LoadGeometry(QFileInfo geomFileInfo, std::vector<std::string> names);
  std::vector<int> GetExcludeProjFiles(bool bManAngleGap, double gantryAngleInterval);
  void LoadSelectedProj(const std::vector<int> &exclude_ids, const std::vector<std::string> &names);
  void saveHisHeader();
  void ApplyBowtie(ProjReaderType::Pointer &reader,
    FilterReaderType::Pointer &bowtie_reader);
  void NormalizeProjections(ProjReaderType::Pointer &reader);
  bool ResampleProjections(double resample_factor);
  void BowtieByFit(bool fullfan, QStringList params);
  int CropSkinUsingThreshold(int threshold, int erode_radius, int dilate_radius);
  void GeneratePOIData(bool AnteriorToPosterior, double table_posY);
  void Export2DDoseMapAsMHA(QString strPath);
  void ExportProjGeometryTXT(QString strPath);
  void ScatterCorPerProjRef(double scaMedian, double scaGaussian, int postScatMedianSize, bool use_cuda, bool save_dicom);

  // void GetSelectedIndices(const std::vector<double>& vFullAngles,
  // std::vector<double>& vNormAngles, std::vector<int>& vTargetIdx, bool bCW);
  void GetExcludeIndexByNames(const QString &outlierListPath,
                              std::vector<std::string> &vProjFileFullPath,
                              std::vector<int> &vExcludeIdx);
  void GetSelectedIndices(const std::vector<double> &vFullAngles,
                          std::vector<double> &vNormAngles,
                          std::vector<int> &vTargetIdx, bool bCW,
                          std::vector<int> &vExcludingIdx);

  void SetMaxAndMinValueOfProjectionImage(); // scan m_spProjImg3D and update
                                             // m_fProjImgValueMin, max
  double GetMaxAndMinValueOfProjectionImage(
      double &fProjImgValueMax, double &fProjImgValueMin,
      const FloatImageType::Pointer &projImage); // , double theoreticalMin);

  double GetValueFrom3DImageFloat(int reqX, int reqY, int reqZ,
                                  FloatImageType::Pointer &sp3DFloatImage);
  double GetValueFrom3DImageUshort(int reqX, int reqY, int reqZ,
                                   UShortImageType::Pointer &sp3DUshortImage);

  bool IsFileNameOrderCorrect(std::vector<std::string> &vFileNames);

  void PostApplyFOVDispParam(float physPosX, float physPosY, float physRadius, float physTablePosY);

  // void ExportDICOM_SHORT(SHORT_ImageType::Pointer& sp3DshortImage);//NOT
  // COMPLETED YET!! Export DICOM without Source DICOM is not possible
  void CopyDictionary(itk::MetaDataDictionary &fromDict,
                      itk::MetaDataDictionary &toDict); // NOT COMPLETED YET!!
                                                        // Export DICOM without
                                                        // Source DICOM is not
                                                        // possible

  void DoBeamHardeningCorrection();

  void Draw2DFrom3D(UShortImageType::Pointer &pImg, enPLANE direction,
                    double pos, YK16GrayImage &Output2D);
  void Draw2DFrom3DDouble(UShortImageType::Pointer &spFixedImg,
                          UShortImageType::Pointer &spMovingImg,
                          enPLANE enPlane, double pos, YK16GrayImage &YKFixed,
                          YK16GrayImage &YKMoving);
  void Draw2DFrom3DDouble(UShortImageType::Pointer &spFixedImg,
                          UShortImageType::Pointer &spMovingImg,
                          enPLANE enPlane, double pos, AG17RGBAImage &YKFixed,
                          AG17RGBAImage &YKMoving);

  void RegisterImgDuplication(enREGI_IMAGES src, enREGI_IMAGES target);

  // plastimatch skin / bubble-remover
  // QString getPathCTDir(enMachineType enType);//DICOM Dir
  // QString getPathRS(enMachineType enType);//RS path
  // bool loadPlanCTFromDCM(QString& strCTDirPath, QString& strRSPath);//using
  // plastimatch, prepare m_spRefCTImg. Remove air, RS is needed  Skin will be
  // removed, bubble will be filled

  void FindAllRelevantPaths(const QString &pathProjHisDir);
  QString MakeElektaXML(const QString &filePath_ImageDBF,
                        const QString &filePath_FrameDBF,
                        const QString &DICOM_UID);

  bool LoadShortImageToUshort(QString &strPath,
                              UShortImageType::Pointer &pUshortImage);
  void init_DlgRegistration(QString &strDCM_UID);

  // using RTK forward projection algorithm, generate 2D projection image files
  // (as line integral, mu_t)
  void CPU_ForwardProjection(UShortImageType::Pointer &spVolImg3D,
                             GeometryType::Pointer &spGeometry,
                             UShortImageType::Pointer &spProjCT3D);
#ifdef USE_CUDA
  void CUDA_ForwardProjection(UShortImageType::Pointer &spVolImg3D,
                              GeometryType::Pointer &spGeometry,
                              UShortImageType::Pointer &spProjCT3D);
#endif
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
  void AfterScatCorrectionMacro(bool use_cuda, bool save_dicom);

  // His file export from 3D proj file
  void SaveProjImageAsHIS(
      UShortImageType::Pointer &spProj3D, std::vector<YK16GrayImage> arrYKImage,
      QString &strSavingFolder,
      double resampleF); // arrYKImage include HIS header and original file name

  void ConvertLineInt2Intensity(FloatImageType::Pointer &spProjLineInt3D,
                                UShortImageType::Pointer &spProjIntensity3D,
                                int bkIntensity);
  void ConvertIntensity2LineInt(UShortImageType::Pointer &spProjIntensity3D,
                                FloatImageType::Pointer &spProjLineInt3D,
                                int bkIntensity);

  void Get2DFrom3D(UShortImageType::Pointer &spSrcImg3D,
                   FloatImage2DType::Pointer &spTargetImg2D, int idx,
                   enPLANE iDirection);
  void Set2DTo3D(FloatImage2DType::Pointer &spSrcImg2D,
                 UShortImageType::Pointer &spTargetImg3D, int idx,
                 enPLANE iDirection);

  void AllocateByRef(UShortImageType::Pointer &spRefImg3D,
                     UShortImageType::Pointer &spTarImg3D);
  void AllocateByRef(FloatImage2DType::Pointer &spRefImg2D,
                     FloatImage2DType::Pointer &spTarImg2D);
  void AllocateByRef(FloatImageType::Pointer &spRefImg3D,
                     FloatImageType::Pointer &spTarImg3D);
  void AllocateByRef(UShortImageType::Pointer &spRefImg3D,
                     FloatImageType::Pointer &spTarImg3D);
  void AllocateByRef(FloatImageType::Pointer &spRefImg3D,
                     UShortImageType::Pointer &spTarImg3D);

  // void ResampleItkImage(OutputImageType::Pointer& spImgFloat, double
  // resampleF);  Resample proj images
  void ResampleItkImage(FloatImageType::Pointer &spSrcImg,
                        FloatImageType::Pointer &spTarImg, double resFactor);
  void ResampleItkImage(UShortImageType::Pointer &spSrcImg,
                        UShortImageType::Pointer &spTarImg, double resFactor);
  void ResampleItkImage2D(FloatImage2DType::Pointer &spSrcImg2D,
                          FloatImage2DType::Pointer &spTarImg2D,
                          double resFactor); // using slice iterator

  void DoReconstructionFDK(enREGI_IMAGES target);
  void CudaDoReconstructionFDK(enREGI_IMAGES target);
  void OpenCLDoReconstructionFDK(enREGI_IMAGES target);

  // void SaveUSHORTAsSHORT_DICOM (USHORT_ImageType::Pointer& spImg, QString&
  // strPatientID, QString& strPatientName);//ushort image --> short image -->
  void SaveUSHORTAsSHORT_DICOM(UShortImageType::Pointer &spImg,
                               QString &strPatientID, QString &strPatientName,
                               QString &strPathTargetDir);

  static void ConvertUshort2Short(UShortImageType::Pointer &spImgUshort,
                           ShortImageType::Pointer &spImgShort);

  double
  CalculateIntensityScaleFactorFromMeans(UShortImageType::Pointer &spProjRaw3D,
                                         UShortImageType::Pointer &spProjCT3D);
  double GetRawIntensityScaleFactor(QString strRef_mAs, QString strCur_mAs);

  // void GetAngularWEPL_SinglePoint(USHORT_ImageType::Pointer& spImage, int
  // angleGap, VEC3D calcPt, int curPtIdx, std::vector<WEPLData>&
  // vOutputWEPLData, bool bAppend);//output std::vector: append
  void GetAngularWEPL_SinglePoint(UShortImageType::Pointer &spUshortImage,
                                  float fAngleGap, float fAngleStart,
                                  float fAngleEnd, VEC3D calcPt, int curPtIdx,
                                  std::vector<WEPLData> &vOutputWEPLData,
                                  bool bAppend);
  void GetAngularWEPL_MultiPoint(UShortImageType::Pointer &spUshortImage,
                                 float fAngleGap, float fAngleStart,
                                 float fAngleEnd,
                                 std::vector<WEPLData> &vOutputWEPLData,
                                 bool bAppend);
  void GetAngularWEPL_window(UShortImageType::Pointer &spUshortImage,
                             float fAngleGap, float fAngleStart,
                             float fAngleEnd,
                             std::vector<WEPLData> &vOutputWEPLData,
                             bool bAppend);

  //	void UpdateUIAfterLoading(QString& imgName);

  void LoadExternalFloatImage(QString &strPath, bool bConversion);
  void TransformationRTK2IEC(FloatImageType::Pointer &spSrcTarg);

  void MedianFilterByGUI(UShortImageType::SizeType indexRadius); // params are given at the UI

  /*Temporary implementation for XVI5 xml*/
  bool
  LoadXVIGeometryFile(const char *filePath); // temporary implenetation
                                             // using QT XML. This is for
                                             // XVI v >5.0.2. _Frames.xml is
                                             // in every projection folder

  FLEXDATA XML_parseFrameForXVI5(QXmlStreamReader &xml);
  QString XML_GetSingleItemString(QXmlStreamReader &xml);

  bool GetXrayParamFromINI(QString &strPathINI, float &kVp, float &mA,
                           float &ms);

  void SetProjDir(QString &strProjPath);

  void ExportAngularWEPL_byFile(QString &strPathOutput, double fAngleStart, double fAngleEnd, double fAngleGap);

  void ExportReconSHORT_HU(UShortImageType::Pointer &spUsImage,
                           QString &outputFilePath);
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
                 float physPosY, float physRadius, float physTablePosY);

  void GenerateCylinderMask(UShortImageType::Pointer &spImgCanvas,
                            float fDcmPosX, float fDcmPosY, float fRadius);

  float GetMeanIntensity(UShortImageType::Pointer &spImg, float sphereR,
                         float *sdIntensity = nullptr);

  void AddConstHU(UShortImageType::Pointer &spImg, int HUval);

  bool ResortCBCTProjection(std::vector<int> &vIntPhaseBinSelected,
                            QString &strPathForXML, QString &strPathProjRoot,
                            QString &strUID,
                            std::vector<float> &vFloatPhaseFull,
                            GeometryType::Pointer &spGeomFull,
                            std::vector<std::string> &vProjPathsFull);

  void AppendInPhaseIndex(int iPhase, std::vector<float> &vFloatPhaseFull,
                          std::vector<int> &vOutputIndex, int margin = 5);

  void LoadShort3DImage(QString &filePath, enREGI_IMAGES enTarget);
  // Read long INIXVI text file and read couch shift values. apply cm -> mm
  // conversion (multiply 10). NO sign changes.
  bool GetCouchShiftFromINIXVI(QString &strPathINIXVI, VEC3D *pTrans,
                               VEC3D *pRot);

  // This function came from the tracking project. trans values are all in mm,
  // DICOM x, y, z
  void
  ImageTransformUsingCouchCorrection(UShortImageType::Pointer &spUshortInput,
                                     UShortImageType::Pointer &spUshortOutput,
                                     VEC3D couch_trans, VEC3D couch_rot);

  void GetWEPLDataFromSingleFile(const QString &filePath,
                                 std::vector<VEC3D> &vPOI,
                                 std::vector<WEPLData> &vOutputWEPL, double fAngleStart, double fAngleEnd);

  void SingleForwardProjection(FloatImageType::Pointer &spVolImgFloat,
                               float fMVGanAngle, float panelOffsetX,
                               float panelOffsetY,
                               UShortImageType::Pointer &spProjImg3D,
                               int iSliceIdx);

  bool LoadShortImageDirOrFile(QString &strPathDir,
                               ShortImageType::Pointer &spOutputShortImg);
  static void ConvertShort2Ushort(ShortImageType::Pointer &spInputImgShort,
                           UShortImageType::Pointer &spOutputImgUshort);

  void RotateImgBeforeFwd(UShortImageType::Pointer &spInputImgUS,
                          UShortImageType::Pointer &spOutputImgUS);

  static void ConvertUshort2AttFloat(UShortImageType::Pointer &spImgUshort,
                              FloatImageType::Pointer &spAttImgFloat);

  bool ReadDicomDir(QString &dirPath);

  // using RTK forward projection algorithm, generate 2D projection image files
  // (as line integral, mu_t)


public:
  std::unique_ptr<StructureSet> m_structures;
  //independent raw images
  std::vector<YK16GrayImage> m_arrYKImage;
  int m_iImgCnt{}; // for independent raw images --> no relation to Directroy
                   // based projections

  std::vector<YK16GrayImage> m_arrYKBufProj;
  int m_iCntSelectedProj;

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
  std::unique_ptr<QStandardItemModel> m_pTableModel;
  
  //std::unique_ptr<DlgRegistration> m_pDlgRegistration;
  // DlgHistogram* m_pDlgHistogram;
  //std::unique_ptr<DlgExternalCommand> m_pDlgExternalCommand;

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
  QString m_strCur_mAs; // QString("20,20")
  QString m_strRef_mAs; // QString("64,40")
  QString m_strError;

  int m_iFixedOffset_ScatterMap; // fixed! allows negative value of scatter
  double m_fResampleF; // typically 0.5. this is updated during LoadSelectedProj
                       // image and ui.lineEdit_DownResolFactor.//also affects
                       // all other scatter correction method
  double m_fProjSpacingX; // updated from SelectedProjLoad
  double m_fProjSpacingY; // updated from SelectedProjLoad

  std::vector<VEC3D> m_vPOI_DCM{}; // initialized by file Load

  std::vector<std::string> m_vSelectedFileNames{};

  bool m_bMacroContinue;

  std::vector<float> m_vPhaseFloat{};

  QStringList m_strListPerProjRefVol;

  QString m_strPathDefaultConfigFile;

  std::vector<int> m_vExcludeProjIdx{}; // if kVON (exposed_ tag is false

  // private:
};

#endif // CBCTRECON_H
