#ifndef AG17RGBAIMAGE_H
#define AG17RGBAIMAGE_H

// Qt (gets these from YK16 )
// #include <QImage>
// #include <QVector>
// #include <QString>

// Local
#include "YK16GrayImage.h" // for enumerators and types without overloading
#include "cbctrecon_config.h"

class CBCTRECON_API AG17RGBAImage {
public:
  AG17RGBAImage();
  AG17RGBAImage(int width, int height);
  ~AG17RGBAImage();

  int m_iWidth;
  int m_iHeight;
  // added: 20140206
  double m_fSpacingX; //[mm/px]
  double m_fSpacingY;

  std::valarray<unsigned short> m_pData; // 0 - 65535

  QPixmap *m_pPixmap; // Actually, no need!
  QImage m_QImage;
  // QPainter* m_pPainter;

  bool CopyFromBuffer(const std::valarray<unsigned short> &pImageBuf, int width,
                      int height);
  bool CloneImage(AG17RGBAImage &other);

  bool CreateImage(int width, int height, unsigned short usVal);

  bool FillPixMap(int winMid, int winWidth);
  bool FillPixMapMinMax(int winMin, int winMax); // 0-65535 Сп window level

  bool FillPixMapDual(int winMid1, int winMid2, int winWidth1, int winWidth2);
  bool FillPixMapMinMaxDual(int winMin1, int winMin2, int winMax1,
                            int winMax2); // 0-65535 Сп window level

  bool IsEmpty() const;
  bool ReleaseBuffer();

  // bool CalcImageInfo (double& meanVal, double& STDV, double& minVal, double&
  // maxVal);
  bool CalcImageInfo();
  double CalcAveragePixelDiff(AG17RGBAImage &other);

  static void CopyYKImage2ItkImage(AG17RGBAImage *pYKImage,
                                   UnsignedShortImageType::Pointer &spTarImage);
  static void CopyItkImage2YKImage(UnsignedShortImageType::Pointer &spSrcImage,
                                   AG17RGBAImage *pYKImage);

  QString m_strFilePath;

  double m_fPixelMean;
  double m_fPixelSD;
  double m_fPixelMin;
  double m_fPixelMax;

  static void Swap(AG17RGBAImage *pImgA, AG17RGBAImage *pImgB);

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

  bool m_bShowInvert;

  void MultiplyConstant(double multiplyFactor);

  void SetSpacing(double spacingX, double spacingY) {
    m_fSpacingX = spacingX;
    m_fSpacingY = spacingY;
  };

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
  void SetOffset(int offsetX, int offsetY) {
    m_iOffsetX = offsetX;
    m_iOffsetY = offsetY;
  }
  double m_fZoom;
  void SetZoom(double fZoom);

  // SPLIT VIEW
  QPoint m_ptSplitCenter; // Fixed image with Moving image. center is based on
                          // dataPt.//Fixed Image: Left Top + Right Bottom,
                          // Moving: Right Top + Left Bottom
  int m_enSplitOption{};
  // This cetner is moved while Left Dragging //All split and crosshair are data
  // point based!
  void SetSplitOption(enSplitOption option) { m_enSplitOption = option; }
  void SetSplitCenter(QPoint &ptSplitCenter); // From mouse event, data point
  // void SetSplitCenter(int centerX, int centerY)
  // {m_ptSplitCenter.setX(centerX); m_ptSplitCenter.setY(centerY);}//From mouse
  // event, data point
  bool ConstituteFromTwo(AG17RGBAImage &YKImg1,
                         AG17RGBAImage &YKImg2); // YKImg1 and two should be in
                                                 // exactly same dimension and
                                                 // spacing
  bool isPtInFirstImage(int dataX, int dataY) const;

  void SetProfileProbePos(int dataX, int dataY);
  unsigned short GetProfileProbePixelVal();
  void GetProfileData(int dataX, int dataY, QVector<double> &vTarget,
                      enProfileDirection direction);
  void GetProfileData(QVector<double> &vTarget, enProfileDirection direction);

  void EditImage_Flip();
  void EditImage_Mirror();

  void MedianFilter(int iMedianSizeX, int iMedianSizeY);

  double m_fResampleFactor{}; // if it is not the 1.0, the data is already
                              // resampled.

  UnsignedShortImageType::Pointer CloneItkImage() const;
  void ResampleImage(double fResampleFactor);

  void UpdateFromItkImage(UnsignedShortImageType::Pointer &spRefItkImg);
  void UpdateFromItkImageFloat(FloatImageType2D::Pointer &spRefItkImg);

  void InvertImage();
};

#endif
