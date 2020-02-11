#ifndef CBCTRECON_IO_H
#define CBCTRECON_IO_H
/*All IO functions used with cbctrecon*/
#ifdef OF
#undef OF
#endif // OF
#include "gdcmAttribute.h"
#include <itkDOMNode.h>

#include "cbctrecon_config.h"
#include "cbctrecon_types.h"

#include <filesystem>
#include <string>

struct Rtss_modern;

namespace fs = std::filesystem;

namespace crl {

CBCTRECON_API fs::path MakeElektaXML(const fs::path &filePath_ImageDBF,
                                        const fs::path &filePath_FrameDBF,
                                        const std::string &DICOM_UID);

CBCTRECON_API FLEXDATA XML_parseFrameForXVI5(itk::DOMNode::Pointer xml);

CBCTRECON_API [[nodiscard]] rtk::ThreeDCircularProjectionGeometry::Pointer
LoadRTKGeometryFile(const fs::path &filePath);

CBCTRECON_API bool
LoadShortImageToUshort(fs::path &strPath,
                       UShortImageType::Pointer &pUshortImage);

CBCTRECON_API void ExportReconSHORT_HU(UShortImageType::Pointer &spUsImage,
                                       std::string &outputFilePath);
CBCTRECON_API bool
LoadShortImageDirOrFile(fs::path &strPathDir,
                        ShortImageType::Pointer &spOutputShortImg);
bool CBCTRECON_API SaveDoseGrayImage(const fs::path& filePath, int width,
                                     int height, double spacingX,
                                     double spacingY, double originLeft_mm,
                                     double originTop_mm,
                                     unsigned short *pData);

template <typename ImageType>
void saveImageAsMHA(typename ImageType::Pointer const &image,
                    std::string filename = "Projections.mha") {
  using ImageWriterType = itk::ImageFileWriter<ImageType>;
  typename ImageWriterType::Pointer writer = ImageWriterType::New();
  writer->SetInput(image);
  writer->SetFileName(filename);
  // writer->SetUseCompression(true); // quite expensive in CPU time
  // ITK's zlib and native-compiled zlib are almost equally fast:
  // Being about 75 seconds slower than without compression (in test with
  // alderson head scan)
  // (maybe an option to change compression rate will appear some day)
  writer->Update();
}

template <typename ImageType> auto loadMHAImageAs(const std::string &filename) {
  auto reader = itk::ImageFileReader<ImageType>::New();
  reader->SetFileName(filename);
  reader->Update();
  return reader->GetOutput();
}

std::vector<std::string> CBCTRECON_API get_dcm_image_files(fs::path &dir);

template <int group, int element, typename T>
auto gdcm_attribute_from(T &parent) {
  // auto attribute = gdcm::Attribute<group, element>();
  const auto tag = gdcm::Tag(group, element);
  const auto &data_element = parent->GetDataElement(tag);
  auto at = gdcm::Attribute<group, element>();
  at.SetFromDataElement(data_element);
  return at;
}

std::unique_ptr<Rtss_modern>
    CBCTRECON_API load_rtstruct(const std::string &filename);

bool CBCTRECON_API AlterData_RTStructureSetStorage(
    const fs::path &input_file, const Rtss_modern *input_rt_struct,
    const fs::path &output_file);

void CBCTRECON_API
ConvertUshort2Short(const UShortImageType::Pointer &spImgUshort,
                    ShortImageType::Pointer &spImgShort);

void CBCTRECON_API
ConvertShort2Ushort(ShortImageType::Pointer &spInputImgShort,
                    UShortImageType::Pointer &spOutputImgUshort);

fs::path CBCTRECON_API SaveUSHORTAsSHORT_DICOM(
    UShortImageType::Pointer &spImg, std::string &strPatientID,
    std::string &strPatientName, fs::path &strPathTargetDir);
fs::path CBCTRECON_API SaveUSHORTAsSHORT_DICOM_gdcmITK(
    UShortImageType::Pointer &spImg, std::string &strPatientID,
    std::string &strPatientName, fs::path &strPathTargetDir);

std::string CBCTRECON_API
get_output_options(const UShortImageType::Pointer &m_spFixed);

} // namespace crl

#endif // CBCTRECON_IO_H
