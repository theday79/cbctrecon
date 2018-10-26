#include <cstdio>
#include <fstream>

#include <QPixmap>

#include "itkAffineTransform.h"
#include "itkMedianImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

#include "YK16GrayImage.h"

YK16GrayImage::YK16GrayImage() {
  m_iWidth = 0;
  m_iHeight = 0;

  m_pData = nullptr;
  m_pPixmap = nullptr; // for display
  m_pElektaHisHeader = nullptr;

  // m_pQImage = nullptr;

  // m_pPainter = nullptr;

  m_fPixelMean = 0.0;
  m_fPixelSD = 0.0;
  m_fPixelMin = 0.0;
  m_fPixelMax = 0.0;

  m_fPixelMean_ROI = 0.0;
  m_fPixelSD_ROI = 0.0;
  m_fPixelMin_ROI = 0.0;
  m_fPixelMax_ROI = 0.0;
  m_bDrawROI = false;
  m_bShowInvert = false;

  m_fSpacingX = 0.0;
  m_fSpacingY = 0.0;

  m_ptProfileProbe.setX(0); // Mouse Clicked Position --> Data
  m_ptProfileProbe.setY(0);

  m_bDrawProfileX = false;
  m_bDrawProfileY = false;

  m_ptFOVCenter.setX(0);
  m_ptFOVCenter.setY(0);
  m_iFOVRadius = 0; // data pixel
  m_bDrawFOVCircle = false;

  m_iTableTopPos = 0;
  m_bDrawTableLine = false;

  m_bDrawCrosshair = false;

  m_fZoom = 1.0;
  m_iOffsetX = 0;
  m_iOffsetY = 0;

  m_enSplitOption = PRI_LEFT_TOP;
  m_fResampleFactor = 1.0;
}

YK16GrayImage::YK16GrayImage(const int width, const int height) {
  m_iWidth = width;
  m_iHeight = height;

  m_pData = nullptr;
  m_pPixmap = nullptr; // for display
  m_pElektaHisHeader = nullptr;

  m_fPixelMean = 0.0;
  m_fPixelSD = 0.0;
  m_fPixelMin = 0.0;
  m_fPixelMax = 0.0;

  m_fPixelMean_ROI = 0.0;
  m_fPixelSD_ROI = 0.0;
  m_fPixelMin_ROI = 0.0;
  m_fPixelMax_ROI = 0.0;
  m_bDrawROI = false;
  m_bShowInvert = false;

  m_fSpacingX = 0.0;
  m_fSpacingY = 0.0;

  m_ptProfileProbe.setX(0); // Mouse Clicked Position --> Data
  m_ptProfileProbe.setY(0);

  m_bDrawProfileX = false;
  m_bDrawProfileY = false;

  m_ptFOVCenter.setX(0);
  m_ptFOVCenter.setY(0);
  m_iFOVRadius = 0; // data pixel
  m_bDrawFOVCircle = false;

  m_iTableTopPos = 0;
  m_bDrawTableLine = false;

  m_bDrawCrosshair = false;

  m_fZoom = 1.0;
  m_iOffsetX = 0;
  m_iOffsetY = 0;

  // m_pQImage = nullptr;
  m_enSplitOption = PRI_LEFT_TOP;
  m_fResampleFactor = 1.0;

  CreateImage(width, height, 0);
}

YK16GrayImage::~YK16GrayImage() { ReleaseBuffer(); }

bool YK16GrayImage::ReleaseBuffer() {
  if (m_pData != nullptr) {
    delete[] m_pData;
    m_pData = nullptr;
  }
  if (m_pPixmap != nullptr) {
    delete m_pPixmap; // delete  m_pPixmap??
    m_pPixmap = nullptr;
  }

  m_iWidth = 0;
  m_iHeight = 0;

  if (m_pElektaHisHeader != nullptr) {
    delete[] m_pElektaHisHeader;
    m_pElektaHisHeader = nullptr;
  }

  return true;
}

bool YK16GrayImage::IsEmpty() const { return m_pData == nullptr; }

bool YK16GrayImage::CreateImage(const int width, const int height,
                                unsigned short usVal) {
  if (width < 1 || height < 1)
    return false;

  if (usVal < 0 || usVal > 65535)
    usVal = 0;

  // if (m_pData != nullptr)
  // delete [] m_pData;

  ReleaseBuffer();

  m_iWidth = width;
  m_iHeight = height;

  const auto imgSize = width * height;
  m_pData = new unsigned short[imgSize];

  for (auto i = 0; i < imgSize; i++) {
    m_pData[i] = usVal;
  }

  // m_QImage = QImage(m_iWidth, m_iHeight,QImage::Format_RGB888);

  return true;
}

bool YK16GrayImage::LoadRawImage(const char *filePath, const int width,
                                 const int height) {
  if (width < 1 || height < 1)
    return false;

  // if (m_pData != nullptr)
  //	delete [] m_pData;
  ReleaseBuffer();

  m_strFilePath = filePath;

  m_iWidth = width;
  m_iHeight = height;

  const auto img_size = static_cast<size_t>(width * height);
  m_pData = new unsigned short[img_size];

  // aqprintf("ImageInfo in LoadRawImage, w: %d  h: %d   %d  %d \n",width,
  // height, m_iWidth, m_iHeight);
  FILE *fd = nullptr;
  if (fopen_s(&fd, filePath, "rb") == 0) {
    std::cerr << "Could not open file: " << filePath << " for writing!"
              << std::endl;
    return false;
  }

  if (fd == nullptr)
    return false;

  auto buf = std::valarray<unsigned short>(img_size);
  if (fread(&buf[0], 2, img_size, fd) != img_size) {
    std::cerr << "Could not read Raw Image" << std::endl;
    return false;
  }

  fclose(fd);

  std::copy_n(std::begin(buf), img_size, &m_pData[0]);

  return true;
}

bool YK16GrayImage::CopyFromBuffer(const unsigned short *p_image_buf,
                                   const int width, const int height) const {
  if (m_pData == nullptr)
    return false;
  if (p_image_buf == nullptr)
    return false;
  if (width != m_iWidth || height != m_iHeight)
    return false;

  const auto imgSize = static_cast<size_t>(m_iWidth * m_iHeight);
  std::copy_n(&p_image_buf[0], imgSize, &m_pData[0]);

  // CalcImageInfo();

  return true;
}

template <typename T>
quint32 fill_pixel(T data, size_t lowVal, size_t uppVal, double d_win_width) {
  if (data >= uppVal) {
    return 0xffffffff;
  }
  if (data <= lowVal) {
    return 0xff000000;
  }
  return qRgba(static_cast<uchar>((data - lowVal) / d_win_width * 255.0),
               static_cast<uchar>((data - lowVal) / d_win_width * 255.0),
               static_cast<uchar>((data - lowVal) / d_win_width * 255.0), 255);
}

template <typename T>
quint32 fill_pixel_invert(T data, size_t low_val, size_t upp_val,
                          double d_win_width) {
  if (data >= upp_val) {
    return 0xff000000;
  }
  if (data <= low_val) {
    return 0xffffffff;
  }
  return qRgba(255 - static_cast<uchar>((data - low_val) / d_win_width * 255.0),
               255 - static_cast<uchar>((data - low_val) / d_win_width * 255.0),
               255 - static_cast<uchar>((data - low_val) / d_win_width * 255.0),
               255);
}

void fill_array(std::valarray<quint32> &tmpData, unsigned short *m_pData,
                size_t lowVal, size_t uppVal, size_t d_win_width,
                const bool m_bShowInvert) {
  if (!m_bShowInvert) {
    std::transform(&m_pData[0], &m_pData[tmpData.size() - 1],
                   std::begin(tmpData),
                   [&lowVal, &uppVal, &d_win_width](auto data) {
                     return fill_pixel(data, lowVal, uppVal, d_win_width);
                   });
  } else {
    std::transform(
        &m_pData[0], &m_pData[tmpData.size() - 1], std::begin(tmpData),
        [&lowVal, &uppVal, &d_win_width](auto data) {
          return fill_pixel_invert(data, lowVal, uppVal, d_win_width);
        });
  }
}

bool YK16GrayImage::FillPixMap(const int winMid,
                               const int winWidth) // 0-65535 �� window level
{
  if (m_pData == nullptr)
    return false;

  if (m_pPixmap != nullptr) {
    delete m_pPixmap;
    m_pPixmap = nullptr;
  }
  m_pPixmap = new QPixmap(QSize(
      m_iWidth, m_iHeight)); // something happened here!!!: w: 4289140  h: 0

  // 8 bit gray buffer preparing

  const auto size = static_cast<size_t>(m_iWidth * m_iHeight);

  // uchar* tmpData = new uchar [size*3];//RGB
  auto tmpData = std::valarray<quint32>(size); // rgb represented as 0xffRRGGBB
  const auto uppVal = static_cast<unsigned short>(winMid + winWidth / 2.0);
  const auto lowVal = static_cast<unsigned short>(winMid - winWidth / 2.0);
  const auto d_win_width = static_cast<double>(winWidth);
  // It takes 0.4 s in Release mode
  fill_array(tmpData, m_pData, lowVal, uppVal, d_win_width, m_bShowInvert);
  /*
  for (int i = 0; i < m_iHeight; i++) // So long time.... - YKP (not really -
  AGA)
  {
    for (int j = 0; j < m_iWidth; j++) {
      const int tmpIdx = (i * m_iWidth + j); // *3

      if (!m_bShowInvert) {

        if (m_pData[tmpIdx] >= uppVal) {
          tmpData[tmpIdx] = 0xffffffff;

        } else if (m_pData[tmpIdx] <= lowVal) {
          tmpData[tmpIdx] = 0xff000000;
        } else {
          tmpData[tmpIdx] =
              qRgba(static_cast<uchar>((m_pData[tmpIdx] - lowVal) /
                                       d_win_width * 255.0),
                    static_cast<uchar>((m_pData[tmpIdx] - lowVal) /
                                       d_win_width * 255.0),
                    static_cast<uchar>((m_pData[tmpIdx] - lowVal) /
                                       d_win_width * 255.0),
                    255);
        }
      } else {
        if (m_pData[tmpIdx] >= uppVal) {
          tmpData[tmpIdx] = 0xff000000;

        } else if (m_pData[tmpIdx] <= lowVal) {
          tmpData[tmpIdx] = 0xffffffff;
        } else {
          tmpData[tmpIdx] =
              qRgba(255 - static_cast<uchar>((m_pData[tmpIdx] - lowVal) /
                                             d_win_width * 255.0),
                    255 - static_cast<uchar>((m_pData[tmpIdx] - lowVal) /
                                             d_win_width * 255.0),
                    255 - static_cast<uchar>((m_pData[tmpIdx] - lowVal) /
                                             d_win_width * 255.0),
                    255);
        }
      }
    }
  }*/
  /// m_QImage.save("C:\\FromFillPixmap.png");//it works well

  // int iBytesPerLine = m_iWidth*3;

  // QImage tmpQImage = QImage((unsigned char*)tmpData,m_iWidth,
  // m_iHeight,iBytesPerLine, QImage::Format_RGB32); //not deep copy!
  auto tmpQImage =
      QImage(reinterpret_cast<unsigned char *>(&tmpData[0]), m_iWidth,
             m_iHeight, 4 * m_iWidth, QImage::Format_ARGB32); // not deep copy!

  // image ratio (width / height) should be kept constant.
  // PAN: center of image is moving

  // m_QImage = tmpQImage.copy(0,0,m_iWidth, m_iHeight); //memory allocated
  // here!!!

  const auto newWidth = qRound(m_iWidth / m_fZoom);
  const auto newHeight = qRound(m_iHeight / m_fZoom);

  const auto centerX = m_iOffsetX + qRound(m_iWidth / 2.0);
  const auto centerY = m_iOffsetY + qRound(m_iHeight / 2.0);

  const auto newLeftTopX = centerX - qRound(newWidth / 2.0);
  const auto newLeftTopY = centerY - qRound(newHeight / 2.0);
  m_QImage = tmpQImage.copy(newLeftTopX, newLeftTopY, newWidth,
                            newHeight); // memory allocated here!!!

  // copy (0,0, width, height)

  //*m_pPixmap = QPixmap::fromImage(m_QImage); //copy data to pre-allcated
  // pixmap buffer

  // delete[] tmpData;
  return true;
}

bool YK16GrayImage::FillPixMapMinMax(int winMin,
                                     int winMax) // 0-65535 �� window level
{
  if (winMin < 0 || winMax > 65535 || winMin > winMax) {
    winMin = 0;
    winMax = 65535;
  }

  const auto midVal = static_cast<int>((winMin + winMax) / 2.0);
  const auto widthVal = winMax - winMin;

  return FillPixMap(midVal, widthVal);
}

bool YK16GrayImage::SaveDataAsRaw(const char *filePath) const
// save 16 bit gray raw file
{
  if (m_pData == nullptr)
    return false;

  const auto imgSize = m_iWidth * m_iHeight;

  FILE *fd = nullptr;
  if (fopen_s(&fd, filePath, "wb") == 0) {
    std::cerr << "Could not open file: " << filePath << " for writing!"
              << std::endl;
    return false;
  }

  for (auto i = 0; i < imgSize; i++) {
    fwrite(&m_pData[i], 2, 1, fd);
  }

  fclose(fd);
  return true;
}

bool YK16GrayImage::CalcImageInfo() {
  if (m_pData == nullptr)
    return false;

  int i;

  const auto npixels = m_iWidth * m_iWidth;

  if (npixels <= 0)
    return false;

  auto nTotal = 0;
  // minPixel = 4095;
  long minPixel = 65535;
  long maxPixel = 0;
  auto sumPixel = 0.0;

  for (i = 0; i < npixels; i++) {
    const auto pixel = static_cast<double>(m_pData[i]);
    sumPixel += pixel;
    if (m_pData[i] > maxPixel)
      maxPixel = m_pData[i];
    if (m_pData[i] < minPixel)
      minPixel = m_pData[i];
    nTotal++;
  }

  auto meanPixelval = 0.0;

  if (nTotal > 0)
    meanPixelval = sumPixel / static_cast<double>(nTotal);

  auto sqrSum = 0.0;
  for (i = 0; i < npixels; i++) {
    sqrSum = sqrSum + pow(static_cast<double>(m_pData[i]) - meanPixelval, 2.0);
  }
  const auto SD = sqrt(sqrSum / static_cast<double>(nTotal));

  m_fPixelMean = meanPixelval;
  m_fPixelSD = SD;
  m_fPixelMin = minPixel;
  m_fPixelMax = maxPixel;

  return true;
}

double YK16GrayImage::CalcAveragePixelDiff(YK16GrayImage &other) const {
  if (m_pData == nullptr || other.m_pData == nullptr)
    return 0.0;

  const auto totalPixCnt = m_iWidth * m_iHeight;
  auto tmpSum = 0.0;
  for (auto i = 0; i < totalPixCnt; i++) {
    tmpSum = tmpSum + fabs(static_cast<double>(m_pData[i]) -
                           static_cast<double>(other.m_pData[i]));
  }

  return tmpSum / static_cast<double>(totalPixCnt);
}

void YK16GrayImage::Swap(YK16GrayImage *pImgA, YK16GrayImage *pImgB) {
  if (pImgA == nullptr || pImgB == nullptr)
    return;

  if (pImgA->IsEmpty() || pImgB->IsEmpty())
    return;

  YK16GrayImage tmpImg(pImgA->m_iWidth, pImgB->m_iHeight);
  tmpImg.CopyFromBuffer(pImgA->m_pData, pImgA->m_iWidth, pImgA->m_iHeight);
  tmpImg.m_strFilePath = pImgA->m_strFilePath;

  pImgA->CopyFromBuffer(pImgB->m_pData, pImgB->m_iWidth, pImgB->m_iHeight);
  pImgA->m_strFilePath = pImgB->m_strFilePath;

  if (!pImgB->CopyFromBuffer(tmpImg.m_pData, tmpImg.m_iWidth,
                             tmpImg.m_iHeight)) {
    std::cerr << "Copy from buffer failed!" << std::endl;
  }
  pImgB->m_strFilePath = tmpImg.m_strFilePath;
}

bool YK16GrayImage::SaveDataAsHis(const char *filePath,
                                  const bool bInverse) const {
  if (m_pData == nullptr)
    return false;

  if (m_pElektaHisHeader == nullptr)
    return false;

  FILE *fd = nullptr;
  if (fopen_s(&fd, filePath, "wb") == 0) {
    std::cerr << "Could not open file: " << filePath << " for writing!"
              << std::endl;
    return false;
  }

  fwrite(m_pElektaHisHeader, 100, 1, fd);

  const auto imgSize = m_iWidth * m_iHeight;

  for (auto i = 0; i < imgSize; i++) {
    unsigned short tmpVal = 0;

    if (bInverse)
      tmpVal = 65535 - m_pData[i];
    else
      tmpVal = m_pData[i];

    fwrite(&tmpVal, 2, 1, fd);
  }

  fclose(fd);

  return true;
}

void YK16GrayImage::CopyHisHeader(const char *hisFilePath) {
  // open file
  std::ifstream file(hisFilePath, std::ios::in | std::ios::binary);

  if (file.fail())
    std::cout << "Fail to open"
              << "	" << hisFilePath << std::endl;

  // read header

  delete[] m_pElektaHisHeader;

  m_pElektaHisHeader = new char
      [DEFAULT_ELEKTA_HIS_HEADER_SIZE]; // DEFAULT_ELEKTA_HIS_HEADER_SIZE
                                        // = 100
  file.read(m_pElektaHisHeader, DEFAULT_ELEKTA_HIS_HEADER_SIZE);
}
//
// void YK16GrayImage::Swap(YK16GrayImage* pImgA, YK16GrayImage* pImgB)
//{
//	if (pImgA == NULL || pImgB == NULL )
//		return;
//
//	if (pImgA->IsEmpty() || pImgB->IsEmpty() )
//		return;
//
//	YK16GrayImage tmpImg(pImgA->m_iWidth, pImgB->m_iHeight);
//	tmpImg.CopyFromBuffer(pImgA->m_pData,pImgA->m_iWidth, pImgA->m_iHeight);
//	tmpImg.m_strFilePath = pImgA->m_strFilePath;
//
//	pImgA->CopyFromBuffer(pImgB->m_pData,pImgB->m_iWidth, pImgB->m_iHeight);
//	pImgA->m_strFilePath = pImgB->m_strFilePath;
//
//
//	pImgB->CopyFromBuffer(tmpImg.m_pData,tmpImg.m_iWidth, tmpImg.m_iHeight);
//	pImgB->m_strFilePath = tmpImg.m_strFilePath;
//}

void YK16GrayImage::CopyYKImage2ItkImage(
    YK16GrayImage *pYKImage, UShortImage2DType::Pointer &spTarImage) {
  if (pYKImage == nullptr)
    return;
  // Raw File open
  // UShortImage2DType::SizeType tmpSize =
  auto region = spTarImage->GetRequestedRegion();
  auto tmpSize = region.GetSize();

  const int sizeX = tmpSize[0];
  const int sizeY = tmpSize[1];

  if (sizeX < 1 || sizeY < 1)
    return;

  itk::ImageRegionIterator<UShortImage2DType> it(spTarImage, region);

  auto i = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    it.Set(pYKImage->m_pData[i]);
    i++;
  }
  // int totCnt = i;
  // writerType::Pointer writer = writerType::New();
  // writer->SetInput(spTarImage);
  // writer->SetFileName("C:\\ThisImageIs_spSrcImage.png");	//It works!
  // writer->Update();
}
std::unique_ptr<YK16GrayImage>
YK16GrayImage::CopyItkImage2YKImage(UShortImage2DType::Pointer &spSrcImage,
                                    std::unique_ptr<YK16GrayImage> pYKImage) {
  if (pYKImage == nullptr)
    return pYKImage;
  auto region = spSrcImage->GetRequestedRegion();
  auto tmpSize = region.GetSize();

  const int sizeX = tmpSize[0];
  const int sizeY = tmpSize[1];

  if (sizeX < 1 || sizeY < 1)
    return pYKImage;

  itk::ImageRegionIterator<UShortImage2DType> it(spSrcImage, region);

  auto i = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    pYKImage->m_pData[i] = it.Get();
    i++;
  }

  return pYKImage;
}

void YK16GrayImage::CopyItkImage2YKImage(UShortImage2DType::Pointer &spSrcImage,
                                         YK16GrayImage *pYKImage) {
  if (pYKImage == nullptr)
    return;
  auto region = spSrcImage->GetRequestedRegion();
  auto tmpSize = region.GetSize();

  const int sizeX = tmpSize[0];
  const int sizeY = tmpSize[1];

  if (sizeX < 1 || sizeY < 1)
    return;

  itk::ImageRegionIterator<UShortImage2DType> it(spSrcImage, region);

  auto i = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    pYKImage->m_pData[i] = it.Get();
    i++;
  }
}

bool YK16GrayImage::CalcImageInfo_ROI() {
  if (m_pData == nullptr) {
    m_fPixelMean_ROI = -1.0;
    m_fPixelSD_ROI = -1.0;
    m_fPixelMin_ROI = -1.0;
    m_fPixelMax_ROI = -1.0;
    return false;
  }

  if (m_rtROI.width() < 1 || m_rtROI.height() < 1) {
    m_fPixelMean_ROI = -1.0;
    m_fPixelSD_ROI = -1.0;
    m_fPixelMin_ROI = -1.0;
    m_fPixelMax_ROI = -1.0;
    return false;
  }

  // int npixels = m_iWidth * m_iWidth;
  auto nTotal = 0;
  // minPixel = 4095;
  long minPixel = 65535;
  long maxPixel = 0;
  auto sumPixel = 0.0;

  int i, j;

  for (i = m_rtROI.top(); i < m_rtROI.bottom(); i++) {
    for (j = m_rtROI.left(); j < m_rtROI.right(); j++) {
      const auto idx = m_iWidth * i + j;
      const auto pixel = static_cast<double>(m_pData[idx]);
      sumPixel += pixel;
      if (m_pData[idx] > maxPixel)
        maxPixel = m_pData[idx];
      if (m_pData[idx] < minPixel)
        minPixel = m_pData[idx];
      nTotal++;
    }
  }

  const auto meanPixelval = sumPixel / static_cast<double>(nTotal);

  auto sqrSum = 0.0;
  /*for (i = 0; i < nTotal; i++)
  {
          sqrSum = sqrSum + pow(((double)m_pData[i] - meanPixelval),2.0);
  }*/

  for (i = m_rtROI.top(); i < m_rtROI.bottom(); i++) {
    for (j = m_rtROI.left(); j < m_rtROI.right(); j++) {
      const auto idx = m_iWidth * i + j;
      sqrSum =
          sqrSum + pow(static_cast<double>(m_pData[idx]) - meanPixelval, 2.0);
    }
  }

  const auto SD = sqrt(sqrSum / static_cast<double>(nTotal));

  m_fPixelMean_ROI = meanPixelval;
  m_fPixelSD_ROI = SD;
  m_fPixelMin_ROI = minPixel;
  m_fPixelMax_ROI = maxPixel;

  return true;
}

bool YK16GrayImage::setROI(const int left, const int top, const int right,
                           const int bottom) {

  if (left >= right || top >= bottom || left < 0 || right > m_iWidth - 1 ||
      top < 0 || bottom > m_iHeight - 1) {
    m_rtROI.setLeft(0);
    m_rtROI.setTop(0);
    m_rtROI.setRight(m_iWidth - 1);
    m_rtROI.setBottom(m_iHeight - 1);
    return false;
  }
  m_rtROI.setLeft(left);
  m_rtROI.setTop(top);
  m_rtROI.setRight(right);
  m_rtROI.setBottom(bottom);
  return true;
}

void YK16GrayImage::DrawROIOn(const bool bROI_Draw) {
  if (bROI_Draw)
    m_bDrawROI = true;
  else
    m_bDrawROI = false;
}

bool YK16GrayImage::CloneImage(YK16GrayImage &other) {
  if (other.m_pData == nullptr)
    return false;

  // first, delete existing things
  ReleaseBuffer(); // Redundancy. CreateImage will call this too. also Create is
                   // calling this func.

  const auto width = other.m_iWidth;
  const auto height = other.m_iHeight;

  CreateImage(width, height, 0);
  if (!CopyFromBuffer(other.m_pData, width, height)) {
    std::cerr << "Copy from buffer failed!" << std::endl;
  }

  if (other.m_pElektaHisHeader != nullptr) {
    m_pElektaHisHeader = new char[DEFAULT_ELEKTA_HIS_HEADER_SIZE];
    for (auto i = 0; i < DEFAULT_ELEKTA_HIS_HEADER_SIZE; i++) {
      m_pElektaHisHeader[i] = other.m_pElektaHisHeader[i];
    }
  }

  m_strFilePath = other.m_strFilePath;

  m_fPixelMean = other.m_fPixelMean;
  m_fPixelSD = other.m_fPixelSD;
  m_fPixelMin = other.m_fPixelMin;
  m_fPixelMax = other.m_fPixelMax;

  m_rtROI = other.m_rtROI;
  m_fPixelMean_ROI = other.m_fPixelMean_ROI;
  m_fPixelSD_ROI = other.m_fPixelSD_ROI;
  m_fPixelMin_ROI = other.m_fPixelMin_ROI;
  m_fPixelMax_ROI = other.m_fPixelMax_ROI;
  m_bDrawROI = other.m_bDrawROI;
  m_QImage = other.m_QImage;
  m_bShowInvert = other.m_bShowInvert;

  SetSpacing(other.m_fSpacingX, other.m_fSpacingY);

  m_ptProfileProbe = other.m_ptProfileProbe;
  m_bDrawProfileX = other.m_bDrawProfileX;
  m_bDrawProfileY = other.m_bDrawProfileY;

  m_ptFOVCenter = other.m_ptFOVCenter; // data pos
  m_iFOVRadius = other.m_iFOVRadius;   // data pos (pixel)
  m_bDrawFOVCircle = other.m_bDrawFOVCircle;

  m_iTableTopPos = other.m_iTableTopPos; // data pos
  m_bDrawTableLine = other.m_bDrawTableLine;

  m_ptCrosshair = other.m_ptCrosshair; // data position
  m_bDrawCrosshair = other.m_bDrawCrosshair;

  m_ptSplitCenter = other.m_ptSplitCenter; // Fixed image with Moving image.
                                           // center is based on dataPt.//Fixed
                                           // Image: Left Top + Right Bottom,
                                           // Moving: Right Top + Left Bottom
  m_enSplitOption = other.m_enSplitOption;

  m_iOffsetX = other.m_iOffsetX;
  m_iOffsetY = other.m_iOffsetY;

  m_fZoom = other.m_fZoom;

  return true;
}

void YK16GrayImage::MultiplyConstant(const double multiplyFactor) const {
  if (m_pData == nullptr)
    return;

  for (auto i = 0; i < m_iHeight; i++) {
    for (auto j = 0; j < m_iWidth; j++) {
      m_pData[m_iWidth * i + j] = static_cast<unsigned short>(
          static_cast<double>(m_pData[m_iWidth * i + j]) * multiplyFactor);
    }
  }
}

void YK16GrayImage::SetProfileProbePos(const int dataX, const int dataY) {
  if (m_ptProfileProbe.y() >= 0 && m_ptProfileProbe.y() < m_iHeight &&
      m_ptProfileProbe.x() >= 0 && m_ptProfileProbe.x() < m_iWidth) {
    m_ptProfileProbe.setX(dataX);
    m_ptProfileProbe.setY(dataY);
  } else {
    m_ptProfileProbe.setX(0);
    m_ptProfileProbe.setY(0);
  }
}

unsigned short YK16GrayImage::GetProfileProbePixelVal() const {
  unsigned short resultVal = 0;
  if (m_pData == nullptr)
    return 0;

  if (m_ptProfileProbe.y() >= 0 && m_ptProfileProbe.y() < m_iHeight &&
      m_ptProfileProbe.x() >= 0 && m_ptProfileProbe.x() < m_iWidth) {
    resultVal = m_pData[m_iWidth * m_ptProfileProbe.y() + m_ptProfileProbe.x()];
  }

  return resultVal;
}

void YK16GrayImage::GetProfileData(const int dataX, const int dataY,
                                   QVector<double> &vTarget,
                                   const enProfileDirection direction) const {
  if (m_pData == nullptr)
    return;

  if (dataY < 0 || dataY >= m_iHeight || dataX < 0 || dataX >= m_iWidth)
    return;

  vTarget.clear();

  if (direction == DIRECTION_HOR) {
    const auto fixedY = dataY;
    for (auto j = 0; j < m_iWidth; j++) {
      vTarget.push_back(m_pData[m_iWidth * fixedY + j]);
    }
  } else if (direction == DIRECTION_VER) {
    // Upper to Lower profile

    const auto fixedX = dataX;
    for (auto i = 0; i < m_iHeight; i++) {
      vTarget.push_back(m_pData[m_iWidth * i + fixedX]);
    }
  }
}
void YK16GrayImage::GetProfileData(QVector<double> &vTarget,
                                   const enProfileDirection direction) const {
  if (m_pData == nullptr)
    return;

  const auto dataX = m_ptProfileProbe.x();
  const auto dataY = m_ptProfileProbe.y();

  if (dataY < 0 || dataY >= m_iHeight || dataX < 0 || dataX >= m_iWidth)
    return;

  vTarget.clear();

  if (direction == DIRECTION_HOR) {
    const auto fixedY = dataY;
    for (auto j = 0; j < m_iWidth; j++) {
      vTarget.push_back(static_cast<double>(m_pData[m_iWidth * fixedY + j]));
    }
  } else if (direction == DIRECTION_VER) {
    // Upper to Lower profile

    const auto fixedX = dataX;
    for (auto i = 0; i < m_iHeight; i++) {
      vTarget.push_back(static_cast<double>(m_pData[m_iWidth * i + fixedX]));
    }
  }
}

bool YK16GrayImage::ConstituteFromTwo(YK16GrayImage &YKImg1,
                                      YK16GrayImage &YKImg2) {
  // Filtering
  if (YKImg1.IsEmpty() || YKImg2.IsEmpty() ||
      YKImg1.m_iWidth != YKImg2.m_iWidth ||
      YKImg1.m_iHeight != YKImg2.m_iHeight ||
      YKImg1.m_iWidth * YKImg1.m_iHeight == 0)
    return false;

  const auto width = YKImg1.m_iWidth;
  const auto height = YKImg1.m_iHeight;

  CreateImage(width, height, 0);

  const auto centerX = m_ptSplitCenter.x(); // data index
  const auto centerY = m_ptSplitCenter.y();

  switch (m_enSplitOption) {
  case PRI_LEFT_TOP:
    for (auto i = 0; i < centerY; i++) {
      for (auto j = 0; j < centerX; j++) {
        m_pData[width * i + j] = YKImg1.m_pData[width * i + j];
      }
    }

    for (auto i = centerY; i < height; i++) {
      for (auto j = centerX; j < width; j++) {
        m_pData[width * i + j] = YKImg1.m_pData[width * i + j];
      }
    }

    for (auto i = 0; i < centerY; i++) {
      for (auto j = centerX; j < width; j++) {
        m_pData[width * i + j] = YKImg2.m_pData[width * i + j];
      }
    }

    for (auto i = centerY; i < height; i++) {
      for (auto j = 0; j < centerX; j++) {
        m_pData[width * i + j] = YKImg2.m_pData[width * i + j];
      }
    }

    break;
  default:
    break;
  }

  return true;
}

void YK16GrayImage::EditImage_Flip() const {
  if (m_pData == nullptr)
    return;

  const auto imgSize = static_cast<size_t>(m_iWidth * m_iHeight);

  if (imgSize <= 0)
    return;

  //������ �ӽ� buffer ����
  auto pPrevImg = std::valarray<unsigned short>(imgSize);

  int i, j;

  for (i = 0; i < m_iHeight; i++) {
    for (j = 0; j < m_iWidth; j++) {
      pPrevImg[i * m_iWidth + j] = m_pData[i * m_iWidth + j];
    }
  }

  for (i = 0; i < m_iHeight; i++) {
    for (j = 0; j < m_iWidth; j++) {
      const auto tmpX = j;
      const auto tmpY = m_iHeight - i - 1;

      m_pData[i * m_iWidth + j] = pPrevImg[tmpY * m_iWidth + tmpX];
    }
  }
}

void YK16GrayImage::EditImage_Mirror() const {
  if (m_pData == nullptr)
    return;

  const auto imgSize = m_iWidth * m_iHeight;

  if (imgSize <= 0)
    return;

  auto pPrevImg = std::valarray<unsigned short>(imgSize);

  //��ȯ �� �̹����� copy

  for (auto i = 0; i < m_iHeight; i++) {
    for (auto j = 0; j < m_iWidth; j++) {
      pPrevImg[i * m_iWidth + j] = m_pData[i * m_iWidth + j];
    }
  }

  for (auto i = 0; i < m_iHeight; i++) {
    for (auto j = 0; j < m_iWidth; j++) {
      const auto tmpX = m_iWidth - j - 1;
      const auto tmpY = i;

      m_pData[i * m_iWidth + j] = pPrevImg[tmpY * m_iWidth + tmpX];
    }
  }

  // delete[] pPrevImg;
}

inline void fill_index(const size_t i, const size_t j, const size_t iWidth,
                       const size_t uppVal, const size_t lowVal,
                       const bool bShowInvert, const unsigned short *data,
                       std::valarray<quint32> &tmpData, const int winWidth) {

  const auto tmpIdx = static_cast<size_t>(i * iWidth + j); // *3

  const auto d_win_width = static_cast<double>(winWidth);
  if (!bShowInvert) {
    tmpData[tmpIdx] = fill_pixel(data[tmpIdx], lowVal, uppVal, d_win_width);
  } else {
    tmpData[tmpIdx] =
        fill_pixel_invert(data[tmpIdx], lowVal, uppVal, d_win_width);
  }
}

bool YK16GrayImage::FillPixMapDual(const int winMid1, const int winMid2,
                                   const int winWidth1, const int winWidth2) {
  if (m_pData == nullptr)
    return false;

  if (m_pPixmap != nullptr) {
    delete m_pPixmap;
    m_pPixmap = nullptr;
  }
  m_pPixmap = new QPixmap(QSize(
      m_iWidth, m_iHeight)); // something happened here!!!: w: 4289140  h: 0

  // 8 bit gray buffer preparing
  const auto size = m_iWidth * m_iHeight;

  // uchar* tmpData = new uchar [size*3];//RGB
  auto tmpData = std::valarray<quint32>(size); // RGB

  auto uppVal1 = static_cast<int>(winMid1 + winWidth1 / 2.0);
  auto lowVal1 = static_cast<int>(winMid1 - winWidth1 / 2.0);

  auto uppVal2 = static_cast<int>(winMid2 + winWidth2 / 2.0);
  auto lowVal2 = static_cast<int>(winMid2 - winWidth2 / 2.0);

  if (uppVal1 > 65535)
    uppVal1 = 65535;
  if (uppVal2 > 65535)
    uppVal2 = 65535;

  if (lowVal1 <= 0)
    lowVal1 = 0;
  if (lowVal2 <= 0)
    lowVal2 = 0;

  // It takes 0.4 s in Release mode

  if (m_enSplitOption != PRI_LEFT_TOP)
    return false;

  const auto splitX = m_ptSplitCenter.x();
  const auto splitY = m_ptSplitCenter.y();

  // 1/4 sector
  for (auto i = 0; i < splitY; i++) // So long time....
  {
    for (auto j = 0; j < splitX; j++) {
      fill_index(i, j, m_iWidth, uppVal1, lowVal1, m_bShowInvert, m_pData,
                 tmpData, winWidth1);
    }
  }

  // 2/4 sector
  for (auto i = 0; i < splitY; i++) // So long time....
  {
    for (auto j = splitX; j < m_iWidth; j++) {
      fill_index(i, j, m_iWidth, uppVal2, lowVal2, m_bShowInvert, m_pData,
                 tmpData, winWidth2);
    }
  }

  // 3/4 sector
  for (auto i = splitY; i < m_iHeight; i++) // So long time....
  {
    for (auto j = 0; j < splitX; j++) {
      fill_index(i, j, m_iWidth, uppVal2, lowVal2, m_bShowInvert, m_pData,
                 tmpData, winWidth2);
    }
  }

  // 4/4 sector
  for (auto i = splitY; i < m_iHeight; i++) // So long time....
  {
    for (auto j = splitX; j < m_iWidth; j++) {
      fill_index(i, j, m_iWidth, uppVal1, lowVal1, m_bShowInvert, m_pData,
                 tmpData, winWidth1);
    }
  }

  // int iBytesPerLine = m_iWidth*3;
  auto tmpQImage =
      QImage(reinterpret_cast<unsigned char *>(&tmpData[0]), m_iWidth,
             m_iHeight, 4 * m_iWidth, QImage::Format_ARGB32); // not deep copy!

  // Copy only a ROI region. All below are data-point, rather than display
  // points  Outside region wiill be filled with Black by QImage inherent
  // function

  //
  const auto newWidth = qRound(m_iWidth / m_fZoom);
  const auto newHeight = qRound(m_iHeight / m_fZoom);

  const auto centerX = m_iOffsetX + qRound(m_iWidth / 2.0);
  const auto centerY = m_iOffsetY + qRound(m_iHeight / 2.0);

  const auto newLeftTopX = centerX - qRound(newWidth / 2.0);  // data position
  const auto newLeftTopY = centerY - qRound(newHeight / 2.0); // data position
  m_QImage = tmpQImage.copy(newLeftTopX, newLeftTopY, newWidth,
                            newHeight); // memory allocated here!!!
  //                        ^~~~~~~~ and ^~~~~~~~~~~ is already initialized as
  //                        int, no need to qRound

  // m_QImage = tmpQImage.copy(0,0,m_iWidth, m_iHeight); //memory allocated
  // here!!!  YKTEMP: is it needed? no it worked without below: *m_pPixmap =
  // QPixmap::fromImage(m_QImage); //copy data to pre-allcated pixmap buffer

  // delete[] tmpData;
  return true;
}

bool YK16GrayImage::FillPixMapMinMaxDual(int winMin1, int winMin2, int winMax1,
                                         int winMax2) {
  if (winMin1 < 0 || winMax1 > 65535 || winMin1 > winMax1) {
    winMin1 = 0;
    winMax1 = 65535;
  }

  if (winMin2 < 0 || winMax2 > 65535 || winMin2 > winMax2) {
    winMin2 = 0;
    winMax2 = 65535;
  }

  const auto midVal1 = static_cast<int>((winMin1 + winMax1) / 2.0);
  const auto midVal2 = static_cast<int>((winMin2 + winMax2) / 2.0);

  const auto widthVal1 = static_cast<int>(winMax1 - winMin1);
  const auto widthVal2 = static_cast<int>(winMax2 - winMin2);

  return FillPixMapDual(midVal1, midVal2, widthVal1, widthVal2);
}

bool YK16GrayImage::isPtInFirstImage(const int dataX, const int dataY) const {
  if (dataX < 0 || dataX >= m_iWidth || dataY < 0 || dataY >= m_iHeight) {
    std::cout
        << "Fatal error in isPtInFirstImage! Given point is out of image point"
        << std::endl;
    return false;
  }

  if (m_enSplitOption == PRI_LEFT_TOP && !IsEmpty()) {
    return (dataX >= 0 && dataX < m_ptSplitCenter.x() && dataY >= 0 &&
            dataY < m_ptSplitCenter.y()) ||
           (dataX >= m_ptSplitCenter.x() && dataX < m_iWidth &&
            dataY >= m_ptSplitCenter.y() && dataY < m_iHeight);
  }
  return false;
}

void YK16GrayImage::SetSplitCenter(QPoint &ptSplitCenter) {
  if (IsEmpty())
    return;

  if (ptSplitCenter.x() < 0 || ptSplitCenter.x() >= m_iWidth ||
      ptSplitCenter.y() < 0 || ptSplitCenter.y() >= m_iHeight) {
    m_ptSplitCenter.setX(0);
    m_ptSplitCenter.setY(0);
  } else
    m_ptSplitCenter = ptSplitCenter;
}

void YK16GrayImage::SetZoom(const double fZoom) {
  if (fZoom <= 1)
    m_fZoom = 1.0;
  else
    m_fZoom = fZoom;
}

void YK16GrayImage::MedianFilter(const int iMedianSizeX,
                                 const int iMedianSizeY) {
  if (m_pData == nullptr)
    return;

  auto spTmpItkImg = UShortImage2DType::New();

  UShortImage2DType::SizeType size;
  size[0] = m_iWidth;
  size[1] = m_iHeight;
  UShortImage2DType::IndexType idxStart;
  idxStart[0] = 0;
  idxStart[1] = 0;
  UShortImage2DType::SpacingType spacing;
  if (m_fSpacingX * m_fSpacingY == 0) {
    spacing[0] = 1.0;
    spacing[1] = 1.0;
  } else {
    spacing[0] = m_fSpacingX;
    spacing[1] = m_fSpacingY;
  }

  UShortImage2DType::PointType origin;
  origin[0] = size[0] * spacing[0] / -2.0;
  origin[1] = size[1] * spacing[1] / -2.0;

  UShortImage2DType::RegionType region;
  region.SetSize(size);
  region.SetIndex(idxStart);

  spTmpItkImg->SetRegions(region);
  spTmpItkImg->SetSpacing(spacing);
  spTmpItkImg->SetOrigin(origin);
  spTmpItkImg->Allocate();

  CopyYKImage2ItkImage(this, spTmpItkImg);

  typedef itk::MedianImageFilter<UShortImage2DType, UShortImage2DType>
      MedianFilterType;
  auto medianFilter = MedianFilterType::New();
  // medianFilter->SetInput(spTmpItkImg);

  MedianFilterType::InputSizeType radius;
  radius[0] = qRound(iMedianSizeX / 2.0);
  radius[1] = qRound(iMedianSizeY / 2.0);

  medianFilter->SetRadius(radius);
  medianFilter->SetInput(spTmpItkImg);
  medianFilter->Update();

  spTmpItkImg = medianFilter->GetOutput();

  CopyItkImage2YKImage(spTmpItkImg, this);
}

UShortImage2DType::Pointer YK16GrayImage::CloneItkImage() const {
  if (m_pData == nullptr) {
    std::cerr << "Could not CloneItkImage" << std::endl;
    return nullptr;
  }

  auto spTmpItkImg = UShortImage2DType::New();

  UShortImage2DType::SizeType size;
  size[0] = m_iWidth;
  size[1] = m_iHeight;
  UShortImage2DType::IndexType idxStart;
  idxStart[0] = 0;
  idxStart[1] = 0;
  UShortImage2DType::SpacingType spacing;
  if (m_fSpacingX * m_fSpacingY == 0) {
    spacing[0] = 1.0;
    spacing[1] = 1.0;
  } else {
    spacing[0] = m_fSpacingX;
    spacing[1] = m_fSpacingY;
  }
  UShortImage2DType::PointType origin;
  origin[0] = size[0] * spacing[0] / -2.0;
  origin[1] = size[1] * spacing[1] / -2.0;

  UShortImage2DType::RegionType region;
  region.SetSize(size);
  region.SetIndex(idxStart);

  spTmpItkImg->SetRegions(region);
  spTmpItkImg->SetSpacing(spacing);
  spTmpItkImg->SetOrigin(origin);
  spTmpItkImg->Allocate();

  // Raw File open
  // UShortImage2DType::SizeType tmpSize =
  // UShortImage2DType::RegionType region =
  // spTmpItkImg->GetRequestedRegion();  UShortImage2DType::SizeType
  // tmpSize = region.GetSize();

  // int sizeX = tmpSize[0];
  // int sizeY = tmpSize[1];

  // if (sizeX < 1 || sizeY <1)
  //	return;

  itk::ImageRegionIterator<UShortImage2DType> it(
      spTmpItkImg, spTmpItkImg->GetRequestedRegion());

  auto i = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    it.Set(m_pData[i]);
    i++;
  }

  return spTmpItkImg;
}

void YK16GrayImage::ResampleImage(const double fResampleFactor) {
  if (m_pData == nullptr)
    return;

  if (fResampleFactor <= 0)
    return;

  m_fResampleFactor = fResampleFactor;

  UShortImage2DType::SizeType inputSize;
  inputSize[0] = m_iWidth;
  inputSize[1] = m_iHeight;

  UShortImage2DType::SizeType outputSize;
  outputSize[0] = qRound(m_iWidth * fResampleFactor);
  outputSize[1] = qRound(m_iHeight * fResampleFactor);
  // m_iWidth = outputSize[0];
  // m_iHeight = outputSize[1];

  UShortImage2DType::SpacingType outputSpacing;

  if (m_fSpacingX <= 0 || m_fSpacingY <= 0) {
    m_fSpacingX = 1.0;
    m_fSpacingY = 1.0;
  }

  outputSpacing[0] = m_fSpacingX * (static_cast<double>(inputSize[0]) /
                                    static_cast<double>(outputSize[0]));
  outputSpacing[1] = m_fSpacingY * (static_cast<double>(inputSize[1]) /
                                    static_cast<double>(outputSize[1]));

  const auto input = CloneItkImage(); // returns Ushort itk image from data buf
  const auto outputOrigin = input->GetOrigin(); //-204.6 - 204.6  0

  //// Resample the image
  // typedef itk::IdentityTransform<float, 2> TransformType;
  typedef itk::ResampleImageFilter<UShortImage2DType, UShortImage2DType, float>
      ResampleImageFilterType;
  auto resample = ResampleImageFilterType::New();

  typedef itk::AffineTransform<float, 2> TransformType;
  auto transform = TransformType::New();
  typedef itk::NearestNeighborInterpolateImageFunction<UShortImage2DType, float>
      InterpolatorType;
  const auto interpolator = InterpolatorType::New();
  transform->SetIdentity();

  resample->SetInput(input);
  resample->SetOutputDirection(input->GetDirection());
  resample->SetInterpolator(interpolator);
  resample->SetDefaultPixelValue(50);
  resample->SetSize(outputSize);
  resample->SetOutputSpacing(outputSpacing);
  resample->SetOutputOrigin(outputOrigin);
  resample->SetTransform(transform);
  resample->Update();

  UShortImage2DType::Pointer outputImg = resample->GetOutput();
  UpdateFromItkImage(outputImg); // itk --> YKImage
}

void YK16GrayImage::UpdateFromItkImage(
    UShortImage2DType::Pointer &spRefItkImg) {
  if (!spRefItkImg)
    return;

  if (m_pData != nullptr) {
    delete[] m_pData;
    m_pData = nullptr;
  }
  if (m_pPixmap != nullptr) {
    delete m_pPixmap;
    m_pPixmap = nullptr;
  }

  auto size = spRefItkImg->GetRequestedRegion().GetSize();
  // UShortImage2DType::SpacingType spacing = spRefItkImg->GetSpacing();

  m_iWidth = static_cast<int>(size[0]);
  m_iHeight = static_cast<int>(size[1]);

  m_pData = new unsigned short[m_iWidth * m_iHeight];

  itk::ImageRegionIterator<UShortImage2DType> it(
      spRefItkImg, spRefItkImg->GetRequestedRegion());

  auto i = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    m_pData[i] = it.Get();
    i++;
  }
}

void YK16GrayImage::UpdateFromItkImageFloat(
    FloatImage2DType::Pointer &spRefItkImg) {
  if (!spRefItkImg)
    return;

  if (m_pData != nullptr) {
    delete[] m_pData;
    m_pData = nullptr;
  }
  if (m_pPixmap != nullptr) {
    delete m_pPixmap;
    m_pPixmap = nullptr;
  }

  auto size = spRefItkImg->GetRequestedRegion().GetSize();
  // FloatImageType2D::SpacingType spacing = spRefItkImg->GetSpacing();

  m_iWidth = static_cast<int>(size[0]);
  m_iHeight = static_cast<int>(size[1]);

  m_pData = new unsigned short[m_iWidth * m_iHeight];

  itk::ImageRegionIterator<FloatImage2DType> it(
      spRefItkImg, spRefItkImg->GetRequestedRegion());

  auto i = 0;
  for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
    const auto curVal = it.Get();
    unsigned short outVal;

    if (curVal < 0.0)
      outVal = 0;
    else if (curVal > 65535.0)
      outVal = 65535;
    else
      outVal = static_cast<unsigned short>(qRound(curVal));

    m_pData[i] = outVal;
    i++;
  }
}

void YK16GrayImage::InvertImage() const {
  if (m_pData == nullptr)
    return;

  const auto imgSize = m_iWidth * m_iHeight;
  // Data inversion: Default for Elekta XVI system
  for (auto i = 0; i < imgSize; i++) {
    m_pData[i] = static_cast<unsigned short>(65535U - m_pData[i]);
  }
}
