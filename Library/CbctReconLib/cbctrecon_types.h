//
// Created by andreas on 10/20/18.
//

#ifndef CBCTRECON_TYPES_H
#define CBCTRECON_TYPES_H

#include <QString>

#include "itkImage.h"
#ifdef USE_CUDA
#include <itkCudaImage.h>
#endif
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>

#include "rtkProjectionsReader.h"
#include "rtkThreeDCircularProjectionGeometry.h"

using FloatPixelType = float;
using USHORT_PixelType = unsigned short;
using SHORT_PixelType = short;

using FloatImageType = itk::Image<FloatPixelType, 3U>;
using FloatImage2DType = itk::Image<FloatPixelType, 2U>;
using UShortImageType = itk::Image<USHORT_PixelType, 3U>;
using UShortImage2DType = itk::Image<USHORT_PixelType, 2U>;
using ShortImageType = itk::Image<SHORT_PixelType, 3U>;
#ifdef USE_CUDA
using CUDAFloatImageType = itk::CudaImage<FloatPixelType, 3U>;
#endif // USE_CUDA

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

#define DEFAULT_WINLEVEL_MID 10000
#define DEFAULT_WINLEVEL_WIDTH 20000

#define DEFAULT_ELEKTA_HIS_HEADER_SIZE 100

struct BADPIXELMAP {
  int BadPixX;
  int BadPixY;
  int ReplPixX;
  int ReplPixY;
};

enum enProfileDirection {
  DIRECTION_HOR = 0,
  DIRECTION_VER,
};

enum enSplitOption {
  PRI_LEFT_TOP = 0, // Primary Left Top
  PRI_RIGHT_TOP,    // Primary Left Top
  PRI_LEFT,
  PRI_RIGHT,
  PRI_TOP,
  PRI_BOTTOM,
};

enum enPLANE {
  PLANE_AXIAL = 0,
  PLANE_FRONTAL,
  PLANE_SAGITTAL,
};

enum ctType {
  PLAN_CT = 0,
  RIGID_CT = 1,
  DEFORM_CT = 2,
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

enum enDeviceType {
  CUDA_DEVT,
  CPU_DEVT,
  OPENCL_DEVT,
};

enum DCM_MODALITY {
    RTIMAGE,
    RTDOSE,
    RTSTRUCT,
    RTPLAN,
    RTRECORD,
    RTUNKNOWN,
};

//
struct FLEXDATA {
  float fGanAngle;     // MV beam gantry angle
  float fPanelOffsetX; // MV beam gantry angle
  float fPanelOffsetY; // MV beam gantry angle
  bool bKV_On;
  bool bMV_On;
};

struct TIFIFD {
  unsigned short TagID;
  unsigned short DataType;
  int DataCnt;
  int DataOrOffset;
};

struct RATIONAL {
  long a;
  long b;
};

struct FDK_options {
  bool displacedDetectorFilter = true;
  double HannCutX = 0.0;
  double HannCutY = 0.0;
  double CosCut = 0.0;
  double HammCut = 0.0;
  double TruncCorFactor = 0.0;
  bool updateAfterDDF = false; // ui.checkBox_UpdateAfterFiltering->isChecked()
  bool ParkerShortScan = true; // ui.checkBox_UsePSSF->isChecked()
  double ct_spacing[3]; // ui.lineEdit_outImgSp_[AP, SI, LR]->text().toDouble();
  int ct_size[3];       // ui.lineEdit_outImgDim_[AP, SI, LR]->text().toInt()
  int medianRadius[3];
  bool medianFilter = true;
  // indexRadius[0] = ui.lineEdit_PostMedSizeX->text().toInt(); // radius along
  // x
  QString outputFilePath; // ui.lineEdit_OutputFilePath->text();
};

using ItkVectorType = itk::Vector<float, 3U>;
using VectorFieldType = itk::Image<ItkVectorType, 3U>;
using PointType = itk::Point<double, 3U>;
// Sorry, I can't control myself, I just love std::function
using TransformType = std::function<PointType(PointType)>;

struct FloatVector {
  float x;
  float y;
  float z;
};

struct WEPLVector {
  double WEPL;
  FloatVector point;
};

struct DoubleVector {
  double x;
  double y;
  double z;
};

struct IntVector {
  int x;
  int y;
  int z;
};

#endif // CBCTRECON_TYPES_H
