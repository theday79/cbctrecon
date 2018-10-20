//
// Created by andreas on 10/20/18.
//

#ifndef CBCTRECON_TYPES_H
#define CBCTRECON_TYPES_H

#include <QString>

#include "itkImage.h"
#include "rtkThreeDCircularProjectionGeometry.h"
#include "rtkProjectionsReader.h"

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

enum enDeviceType {
  CUDA_DEVT,
  CPU_DEVT,
  OPENCL_DEVT,
};

//
struct FLEXDATA {
  float fGanAngle;     // MV beam gantry angle
  float fPanelOffsetX; // MV beam gantry angle
  float fPanelOffsetY; // MV beam gantry angle
  bool bKV_On;
  bool bMV_On;
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

#endif // CBCTRECON_TYPES_H
