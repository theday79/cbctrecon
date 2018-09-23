#ifndef IO_H
#define IO_H
/*All IO functions used with cbctrecon*/
#include "cbctrecon.h"


bool SaveDoseGrayImage(const char *filePath, int width, int height,
                       double spacingX, double spacingY, double originLeft_mm,
                       double originTop_mm, unsigned short *pData);


template <typename ImageType> // image by value, because we call this from std::thread
void saveImageAsMHA(typename ImageType::Pointer const image) {
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(image);
  writer->SetFileName("Projections.mha");
  writer->Update();
}

#endif // IO_H
