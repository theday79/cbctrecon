#ifndef YK16GRAYIMAGE_H
#define YK16GRAYIMAGE_H

// v20130830 : his header buffer, itk compatible
#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

class QPixmap;
class QLabel;
class QPainter;

#include <QImage>

class CBCTRECON_API YK16GrayImage {
public:
  YK16GrayImage();
  YK16GrayImage(int width, int height);
  ~YK16GrayImage();
  // YK16GrayImage(const YK16GrayImage &) = delete;
  // void operator=(const YK16GrayImage &) = delete;
  // YK16GrayImage(YK16GrayImage &&) = delete;
  // void operator=(YK16GrayImage &&) = delete;

  unsigned short *m_pData; // 0 - 65535

  int m_iWidth;
  int m_iHeight;
  // added: 20140206
  double m_fSpacingX; //[mm/px]
  double m_fSpacingY;

  QPixmap *m_pPixmap; // Actually, no need!
  QImage m_QImage;
  // QPainter* m_pPainter;

  bool LoadRawImage(const char *filePath, int width, int height);
  bool CopyFromBuffer(const unsigned short *p_image_buf, int width,
                      int height) const;
  bool CloneImage(YK16GrayImage &other);

  bool CreateImage(int width, int height, unsigned short usVal);

  bool FillPixMap(int winMid, int winWidth);
  bool FillPixMapMinMax(int winMin, int winMax); // 0-65535 window level

  bool FillPixMapDual(int winMid1, int winMid2, int winWidth1, int winWidth2);
  bool FillPixMapMinMaxDual(int winMin1, int winMin2, int winMax1,
                            int winMax2); // 0-65535 window level

  bool SaveDataAsRaw(const char *filePath) const;
  // bool DrawToLabel(QLabel* lbDisplay);

  bool IsEmpty() const;
  bool ReleaseBuffer();

  // bool CalcImageInfo (double& meanVal, double& STDV, double& minVal, double&
  // maxVal);
  bool CalcImageInfo();
  double CalcAveragePixelDiff(YK16GrayImage &other) const;

  // bool DoPixelReplacement(std::vector<BADPIXELMAP>& vPixelMapping); //based
  // on pixel mapping information, some bad pixels will be replaced with median
  // pixel value near by

  static void CopyYKImage2ItkImage(YK16GrayImage *pYKImage,
                                   UShortImage2DType::Pointer &spTarImage);
  static void CopyItkImage2YKImage(UShortImage2DType::Pointer &spSrcImage,
                                   YK16GrayImage *pYKImage);
  static std::unique_ptr<YK16GrayImage>
  CopyItkImage2YKImage(UShortImage2DType::Pointer &spSrcImage,
                       std::unique_ptr<YK16GrayImage> pYKImage);

  QString m_strFilePath;

  double m_fPixelMean;
  double m_fPixelSD;
  double m_fPixelMin;
  double m_fPixelMax;

  static void Swap(YK16GrayImage *pImgA, YK16GrayImage *pImgB);

  QRect m_rtROI;
  bool setROI(int left, int top, int right,
              int bottom); // if there is error, go to default: entire image
  bool CalcImageInfo_ROI();
  double m_fPixelMean_ROI;
  double m_fPixelSD_ROI;
  double m_fPixelMin_ROI;
  double m_fPixelMax_ROI;
  bool m_bDrawROI;

  void DrawROIOn(bool bROI_Draw); // only rectangle

  // Elekta CBCT recon
  char *m_pElektaHisHeader;
  void CopyHisHeader(const char *hisFilePath);
  // bool SaveDataAsHis (const char *filePath);
  bool SaveDataAsHis(const char *filePath, bool bInverse) const;
  bool m_bShowInvert;

  void MultiplyConstant(double multiplyFactor) const;

  void SetSpacing(const double spacingX, const double spacingY) {
    m_fSpacingX = spacingX;
    m_fSpacingY = spacingY;
  }

  QPoint m_ptProfileProbe; // Mouse Clicked Position --> Data
  bool m_bDrawProfileX;
  bool m_bDrawProfileY;

  QPoint m_ptFOVCenter; // data pos
  int m_iFOVRadius;     // data pos (pixel)
  bool m_bDrawFOVCircle;

  int m_iTableTopPos; // data pos
  bool m_bDrawTableLine;

  QPoint m_ptCrosshair; // data position
  bool m_bDrawCrosshair;

  ////ZOOM and PAN function. Using these information below, prepare the m_QImage
  /// for displaying
  // in qlabel in FillPixMap function
  int m_iOffsetX; // for Pan function.. this is data based offset
  int m_iOffsetY;
  void SetOffset(const int offsetX, const int offsetY) {
    m_iOffsetX = offsetX;
    m_iOffsetY = offsetY;
  }
  double m_fZoom;
  void SetZoom(double fZoom);

  // SPLIT VIEW
  QPoint m_ptSplitCenter; // Fixed image with Moving image. center is based on
                          // dataPt.//Fixed Image: Left Top + Right Bottom,
                          // Moving: Right Top + Left Bottom
  enSplitOption m_enSplitOption;
  // This cetner is moved while Left Dragging //All split and crosshair are data
  // point based!
  void SetSplitOption(const enSplitOption option) { m_enSplitOption = option; }
  void SetSplitCenter(QPoint &ptSplitCenter); // From mouse event, data point
  // void SetSplitCenter(int centerX, int centerY)
  // {m_ptSplitCenter.setX(centerX); m_ptSplitCenter.setY(centerY);}//From mouse
  // event, data point
  bool ConstituteFromTwo(YK16GrayImage &YKImg1,
                         YK16GrayImage &YKImg2); // YKImg1 and two should be in
                                                 // exactly same dimension and
                                                 // spacing
  bool isPtInFirstImage(int dataX, int dataY) const;

  void SetProfileProbePos(int dataX, int dataY);
  unsigned short GetProfileProbePixelVal() const;
  void GetProfileData(int dataX, int dataY, QVector<double> &vTarget,
                      enProfileDirection direction) const;
  void GetProfileData(QVector<double> &vTarget,
                      enProfileDirection direction) const;

  void EditImage_Flip() const;
  void EditImage_Mirror() const;

  void MedianFilter(int iMedianSizeX, int iMedianSizeY);

  double
      m_fResampleFactor; // if it is not the 1.0, the data is already resampled.

  UShortImage2DType::Pointer CloneItkImage() const;
  void ResampleImage(double fResampleFactor);

  void UpdateFromItkImage(UShortImage2DType::Pointer &spRefItkImg);
  void UpdateFromItkImageFloat(FloatImage2DType::Pointer &spRefItkImg);

  void InvertImage() const;

  // will be added later
  /*void EditImage_CW90();
  void EditImage_CCW90();
  void EditImage_Rotation(double angle);*/

  // should be implemented later

  // Flip
  // mirror
};

#endif
