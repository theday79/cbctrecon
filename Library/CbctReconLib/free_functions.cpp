// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#include <charconv>
#include <memory>
#include <string>
#include <string_view>

#include "free_functions.h"

#include "YK16GrayImage.h"
#include "cbctrecon.h"

// CbctReconLib
namespace crl {

// code inspired by
// https://marcoarena.wordpress.com/2017/01/03/string_view-odi-et-amo and B.
// Filipek's book C++17 in Detail
std::vector<std::string_view> split_string(std::string_view strv,
                                           std::string_view delims) {
  std::vector<std::string_view> output;
  auto first = strv.begin();

  while (first != strv.end()) {
    const auto second = std::find_first_of(
        first, std::cend(strv), std::cbegin(delims), std::cend(delims));

    if (first != second) {
      output.emplace_back(strv.substr(std::distance(strv.begin(), first),
                                      std::distance(first, second)));
    }

    if (second == strv.end())
      break;

    first = std::next(second);
  }

  return output;
}

std::string HexStr2IntStr(std::string_view str_hex) {
  int val = 0;
  if (auto [ptr, ec] = std::from_chars(
          str_hex.data(), str_hex.data() + str_hex.size(), val, 16);
      ec == std::errc()) {
    return std::to_string(val);
  }
  return {};
}

std::unique_ptr<YK16GrayImage>
ApplyCalibrationMaps(CbctRecon *p_cr, YK16GrayImage *const &rawImg,
                     const bool DarkCorr, const bool GainCorr,
                     const bool DefectCorr) {
  auto corrImg = std::make_unique<YK16GrayImage>(DEFAULT_ELEKTA_PROJ_WIDTH,
                                                 DEFAULT_ELEKTA_PROJ_HEIGHT);

  // m_pParent->m_pCurrImageRaw->m_pData[i]

  if (!DarkCorr && !GainCorr) {
    for (auto i = 0; i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT;
         i++) {
      corrImg->m_pData[i] = rawImg->m_pData[i]; // raw image
    }
  } else if (DarkCorr && !GainCorr) {
    if (p_cr->m_pImgOffset->IsEmpty()) {
      for (auto i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg->m_pData[i] = rawImg->m_pData[i]; // raw image
      }
    } else {
      for (auto i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        if (rawImg->m_pData[i] > p_cr->m_pImgOffset->m_pData[i]) {
          corrImg->m_pData[i] =
              rawImg->m_pData[i] - p_cr->m_pImgOffset->m_pData[i];
        } else {
          corrImg->m_pData[i] = 0;
        }
      }
    }
  } else if (!DarkCorr && GainCorr) {
    if (p_cr->m_pImgGain->IsEmpty()) {
      for (auto i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg->m_pData[i] = rawImg->m_pData[i]; // raw image
      }
    } else {
      // get a mean value for m_pGainImage
      auto sum = 0.0;
      for (auto i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        sum = sum + p_cr->m_pImgGain->m_pData[i];
      }
      const auto MeanVal =
          sum / static_cast<double>(DEFAULT_ELEKTA_PROJ_WIDTH *
                                    DEFAULT_ELEKTA_PROJ_HEIGHT);

      for (auto i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        if (p_cr->m_pImgGain->m_pData[i] == 0) {
          corrImg->m_pData[i] = rawImg->m_pData[i];
        } else {
          corrImg->m_pData[i] = static_cast<unsigned short>(
              static_cast<double>(rawImg->m_pData[i]) /
              static_cast<double>(p_cr->m_pImgGain->m_pData[i]) * MeanVal);
        }
      }
    }
  }

  else if (DarkCorr && GainCorr) {
    auto bRawImage = false;
    if (p_cr->m_pImgOffset->IsEmpty()) {
      bRawImage = true;
    }
    if (p_cr->m_pImgGain->IsEmpty()) {
      bRawImage = true;
    }

    if (bRawImage) {
      for (auto i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        corrImg->m_pData[i] = rawImg->m_pData[i]; // raw image
      }
    } else // if not raw image
    {
      // get a mean value for m_pGainImage
      auto sum = 0.0;
      for (auto i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        sum = sum +
              (p_cr->m_pImgGain->m_pData[i] - p_cr->m_pImgOffset->m_pData[i]);
      }
      const auto MeanVal =
          sum / static_cast<double>(DEFAULT_ELEKTA_PROJ_WIDTH *
                                    DEFAULT_ELEKTA_PROJ_HEIGHT);

      auto iDenomLessZero = 0;
      auto iDenomLessZero_RawIsGreaterThanDark = 0;
      auto iDenomLessZero_RawIsSmallerThanDark = 0;
      auto iDenomOK_RawValueMinus = 0;
      auto iValOutOfRange = 0;

      for (auto i = 0;
           i < DEFAULT_ELEKTA_PROJ_WIDTH * DEFAULT_ELEKTA_PROJ_HEIGHT; i++) {
        const auto denom = static_cast<double>(p_cr->m_pImgGain->m_pData[i] -
                                               p_cr->m_pImgOffset->m_pData[i]);

        if (denom <= 0) {
          iDenomLessZero++;

          if (rawImg->m_pData[i] > p_cr->m_pImgOffset->m_pData[i]) {
            corrImg->m_pData[i] =
                rawImg->m_pData[i] - p_cr->m_pImgOffset->m_pData[i];
            iDenomLessZero_RawIsGreaterThanDark++;
          } else {
            corrImg->m_pData[i] = 0;
            iDenomLessZero_RawIsSmallerThanDark++;
          }
        } else {
          const auto tmpVal =
              (rawImg->m_pData[i] - p_cr->m_pImgOffset->m_pData[i]) / denom *
              MeanVal;

          if (tmpVal < 0) {
            corrImg->m_pData[i] = 0;
            iDenomOK_RawValueMinus++;
          } else {
            if (tmpVal > 65535) { // 16bit max value
              iValOutOfRange++;
            }

            corrImg->m_pData[i] = static_cast<unsigned short>(tmpVal);
          }
        }
      } // end of for
    }   // end if not bRawImage

  } // else if (m_bDarkCorrApply && m_bGainCorrApply)

  if (DefectCorr && !p_cr->m_vPixelReplMap.empty()) // pixel replacement
  {
    corrImg = BadPixReplacement(p_cr->m_vPixelReplMap, std::move(corrImg));
  }

  return corrImg;
}

std::filesystem::path CorrectSingleFile(CbctRecon *p_cr,
                                        const std::filesystem::path &filePath,
                                        const bool DarkCorr,
                                        const bool GainCorr,
                                        const bool DefectCorr) {
  // Load raw file
  auto rawImg = std::make_unique<YK16GrayImage>(DEFAULT_ELEKTA_PROJ_WIDTH,
                                                DEFAULT_ELEKTA_PROJ_HEIGHT);
  rawImg->LoadRawImage(filePath, DEFAULT_ELEKTA_PROJ_WIDTH,
                       DEFAULT_ELEKTA_PROJ_HEIGHT);

  const auto corrImg =
      ApplyCalibrationMaps(p_cr, rawImg.get(), DarkCorr, GainCorr, DefectCorr);
  // filePath
  // std::string exportName = filePath;
  // corrImg.SaveDataAsRaw();
  const auto endFix = std::string();

  auto dir = std::filesystem::absolute(filePath);
  auto baseName = filePath.stem();
  const auto extName = filePath.extension();

  const auto newFileName = baseName.string() + "_CORR." + extName.string();
  auto newPath = std::filesystem::absolute(dir) / newFileName;

  if (!corrImg->SaveDataAsRaw(newPath)) {
    std::cerr << "Could not save as Raw in: " << newPath << std::endl;
  }

  return newPath;
  // corrImg.ReleaseBuffer();
}

void CorrectSingleFile(CbctRecon *p_cr, YK16GrayImage *pYKRawImg,
                       const bool DarkCorr, const bool GainCorr,
                       const bool DefectCorr) {
  if (pYKRawImg == nullptr) {
    return;
  }

  const auto corrImg =
      ApplyCalibrationMaps(p_cr, pYKRawImg, DarkCorr, GainCorr, DefectCorr);

  // Replace old buffer with new one.
  pYKRawImg->CopyFromBuffer(corrImg->m_pData, corrImg->m_iWidth,
                            corrImg->m_iHeight);
}

void LoadBadPixelMap(std::vector<BADPIXELMAP> &vPixelReplMap,
                     const std::filesystem::path &filePath) {
  vPixelReplMap.clear();

  std::ifstream fin;
  fin.open(filePath.c_str());

  if (fin.fail()) {
    return;
  }

  char str[MAX_LINE_LENGTH];
  // memset(str, 0, MAX_LINE_LENGTH);

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    auto tmpStr = std::string(&str[0]);

    if (tmpStr.find("#ORIGINAL_X") != std::string::npos) {
      break;
    }
  }

  while (!fin.eof()) {
    memset(&str[0], 0, MAX_LINE_LENGTH);
    fin.getline(&str[0], MAX_LINE_LENGTH);
    std::string_view tmpStr{&str[0], MAX_LINE_LENGTH};

    auto strList = split_string(tmpStr, "	");

    if (strList.size() == 4) {
      BADPIXELMAP tmpData{};

      if (from_sv(strList.at(0), tmpData.BadPixX) &&
          from_sv(strList.at(1), tmpData.BadPixY) &&
          from_sv(strList.at(2), tmpData.ReplPixX) &&
          from_sv(strList.at(3), tmpData.ReplPixY)) {
        vPixelReplMap.push_back(tmpData);
      }
    }
  }
  fin.close();
}

std::unique_ptr<YK16GrayImage>
BadPixReplacement(std::vector<BADPIXELMAP> &vPixelReplMap,
                  std::unique_ptr<YK16GrayImage> targetImg) {
  if (vPixelReplMap.empty()) {
    return targetImg;
  }

  for (auto &it : vPixelReplMap) {
    const auto tmpData = it;
    const auto oriIdx =
        static_cast<size_t>(tmpData.BadPixY) * DEFAULT_ELEKTA_PROJ_WIDTH +
        tmpData.BadPixX;
    const auto replIdx =
        static_cast<size_t>(tmpData.ReplPixY) * DEFAULT_ELEKTA_PROJ_WIDTH +
        tmpData.ReplPixX;
    targetImg->m_pData[oriIdx] = targetImg->m_pData[replIdx];
  }

  return targetImg;
}

std::vector<std::string>
GetProjFileNames(CbctRecon *p_cr,
                 std::filesystem::path &dirPath) // main loading fuction for
                                                 // projection images
{

  p_cr->m_iImgCnt = 0; // should be reset
  p_cr->m_iCntSelectedProj = 0;
  p_cr->ReleaseMemory(); // only reset mem for indepent projection images

  std::string regexp;
  switch (p_cr->m_projFormat) {
  case enProjFormat::HIS_FORMAT:
    regexp = "(.[0-9a-fA-F]).his";
    break;
  case enProjFormat::HND_FORMAT:
    regexp = "Proj_(.*).hnd";
    break;
  case enProjFormat::XIM_FORMAT:
    regexp = "Proj_(.*).xim";
    break;
  }
  auto regexpnames = itk::RegularExpressionSeriesFileNames::New();
  regexpnames->SetDirectory(dirPath.string());
  // regexpnames->SetNumericSort(false);
  regexpnames->SetNumericSort(true); // doesn't work with hexadecimal. and
                                     // [true/false] doesn't mean ascending or
                                     // descending
  regexpnames->SetRegularExpression(regexp);
  const auto submatch = 1;
  regexpnames->SetSubMatch(
      submatch); // SetSubMatch(0) led to sorting from last digit of the name

  auto names = regexpnames->GetFileNames();

  rtk::RegisterIOFactories();
  std::vector<size_t> idxtopop;
  for (auto &fn : names) {
    auto imageio = itk::ImageIOFactory::CreateImageIO(
        fn.c_str(), itk::ImageIOFactory::ReadMode);

    if (imageio.IsNull()) {
      idxtopop.push_back(&fn - &names[0]);
    }
  }
  std::reverse(idxtopop.begin(), idxtopop.end());
  for (auto &id : idxtopop) {
    names.erase(names.begin() + id);
  }

  return names;
}

std::vector<size_t> GetExcludeProjFiles(
    const rtk::ThreeDCircularProjectionGeometry::Pointer spFullGeometry,
    const bool bScanDirectionCW, const bool bManAngleGap,
    const double gantryAngleInterval) {
  ///////////////////////////////////Exclude outlier projection files
  auto angle_gaps =
      spFullGeometry->GetAngularGaps(spFullGeometry->GetSourceAngles());

  auto sum_gap =
      std::accumulate(std::begin(angle_gaps), std::end(angle_gaps), 0.0);
  sum_gap /= itk::Math::pi * 180.0;

  auto &gantry_angles = spFullGeometry->GetGantryAngles();
  std::vector<size_t> vSelectedIdx;
  std::vector<size_t> vExcludeIdx;

  if (bManAngleGap) {
    // Select indices for recon
    // Generate norminal gantry values from the first angle
    const auto firstAngle = gantry_angles.at(0);
    const auto lastAngle = gantry_angles.at(gantry_angles.size() - 1);

    std::vector<double> vNormAngles;

    const auto multiSize = std::lround(sum_gap / gantryAngleInterval) + 2;

    // CW only (179.xx -> 181.xx -> 359.xx --> 1.xx --> 179.xx), CCW should be
    // checked later
    for (auto i = 0; i < multiSize; i++) {
      auto curAngle = 0.0;

      if (bScanDirectionCW) {
        curAngle = firstAngle + i * gantryAngleInterval;
        if (curAngle >= 360.0) {
          curAngle = curAngle - 360.0;
        }
      } else {
        curAngle = firstAngle - i * gantryAngleInterval;
        if (curAngle < 0.0) {
          curAngle = curAngle + 360.0;
        }
      }
      // Don't add last gantry angle if their intervals are too small.

      // Last data will be added at the last part
      if (i > multiSize - 5) // last parts of the data
      {
        if (bScanDirectionCW) {
          if (curAngle <=
              lastAngle - gantryAngleInterval / 2.0) // from 5 latest indices,
          {
            vNormAngles.push_back(curAngle);
          }
        } else {
          if (curAngle >=
              lastAngle - gantryAngleInterval / 2.0) // from 5 latest indices,
          {
            vNormAngles.push_back(curAngle);
          }
        }
        // gantryAngleInterval/2.0 is given to remove "very near" value to the
        // last value
      } else {
        vNormAngles.push_back(curAngle);
      }
    }
    vNormAngles.push_back(lastAngle);

    for (auto vNormAngle : vNormAngles) {
      std::cout << "Nominal proj. angle: ";
      std::cout << vNormAngle << std::endl;
    }

    // Collect appropriate indices
    GetSelectedIndices(gantry_angles, vNormAngles, vSelectedIdx,
                       bScanDirectionCW, vExcludeIdx);

    for (auto &it_idx : vSelectedIdx) {
      std::cout << "Index: " << it_idx << "     "
                << "GantryAngle: " << gantry_angles.at(it_idx) << std::endl;
    }
  } else // not manual
  {
    for (size_t i = 0; i < gantry_angles.size(); i++) {
      if (std::find(vExcludeIdx.begin(), vExcludeIdx.end(), i) ==
          vExcludeIdx.end()) { // if i is not included in vExcludeIdx
        vSelectedIdx.push_back(i);
      }
    }
  }
  // Another exlusion for kV off images

  // std::vector<int>::iterator itExclude;
  // for (itExclude = m_vExcludeProjIdx.begin(); itExclude !=
  // m_vExcludeProjIdx.end(); ++itExclude)
  //{
  //    int idx = (*itExclude);
  //    //if (std::find(m_vExcludeProjIdx.begin(), m_vExcludeProjIdx.end(),
  //    curIdx) == m_vExcludeProjIdx.end()) // if i is not included in
  //    vExcludeIdx
  //    //    vSelectedIdx_final.push_back(curIdx);

  //    std::cout << "Exclude " << idx << std::endl;
  //}

  // m_vExcludeProjIdx.clear();

  // for (auto &it_final : vSelectedIdx) {
  //   if (std::find(m_vExcludeProjIdx.begin(), m_vExcludeProjIdx.end(),
  //                 it_final) ==
  //       m_vExcludeProjIdx.end()) { // if i is not included in vExcludeIdx
  //     vSelectedIdx_final.push_back(it_final);
  //   }
  // }

  std::cout << "Total proj count: " << vSelectedIdx.size() << std::endl;

  return vSelectedIdx;
}

void GetSelectedIndices(const std::vector<double> &vFullAngles,
                        std::vector<double> &vNormAngles,
                        std::vector<size_t> &vTargetIdx, const bool bCW,
                        std::vector<size_t> &vExcludingIdx) {
  // projection time. Begins with 179.xxx (CW)
  size_t latest_Idx = 0;

  const auto sizeNom = vNormAngles.size();
  const auto sizeFull = vFullAngles.size();

  for (size_t i = 0; i < sizeNom; ++i) {
    const auto tmpNominalValue = vNormAngles.at(i);

    for (auto j = latest_Idx + 1; j < sizeFull - 1; j++) {
      auto enExcludingMode = 0; // 0: safe,1: right is outlier, 2: left is
                                // outlier, 3: both are outlier

      // 1) Left point is outlier
      if (std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j) !=
              vExcludingIdx.end() &&
          std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j + 1) ==
              vExcludingIdx.end()) {
        enExcludingMode = 2;
      }
      // 2) Right point is outlier
      else if (std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j) ==
                   vExcludingIdx.end() &&
               std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j + 1) !=
                   vExcludingIdx.end()) {
        enExcludingMode = 1;
      }
      // 2) No outlier
      else if (std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j) ==
                   vExcludingIdx.end() &&
               std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j + 1) ==
                   vExcludingIdx.end()) {
        enExcludingMode = 0;
      }
      // 3) Both are outliers
      else if (std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j) !=
                   vExcludingIdx.end() &&
               std::find(vExcludingIdx.begin(), vExcludingIdx.end(), j + 1) !=
                   vExcludingIdx.end()) {
        enExcludingMode = 3;
      }

      auto cur_val = vFullAngles.at(j);
      auto next_val = vFullAngles.at(j + 1);

      if (bCW) {
        // for full gantry angle value of 359.0 - 1.0 interface in CW
        if (cur_val >
            next_val + 0.2) // e.g.359)  - 0.5,, 0.2-->minimum angle diff.
        {
          if (tmpNominalValue < 100) {
            cur_val = cur_val - 360.0;
          } else if (tmpNominalValue > 260) {
            next_val = next_val + 360.0;
          }
        }
        if (tmpNominalValue >= cur_val && tmpNominalValue <= next_val) {

          // Add filtering
          // if j is among the excluding index list (e.g. outlier), just pass
          // it.

          const auto diff_cur = fabs(tmpNominalValue - cur_val);
          const auto diff_next = fabs(tmpNominalValue - next_val);

          if (diff_cur <= diff_next || enExcludingMode == 0 ||
              enExcludingMode == 1) {
            latest_Idx = j;
            vTargetIdx.push_back(latest_Idx);
          } else if (j != sizeFull - 2 && enExcludingMode == 2) {
            latest_Idx = j + 1;
            vTargetIdx.push_back(latest_Idx);
          } else {
            latest_Idx = j;
            // Skip to pushback
          }

          break;
        }
      } else // in CCW case
      {
        // for full gantry angle value of 1.0 - 359.0 interface in CCW
        if (cur_val <
            next_val + 0.01) // e.g.359)  - 0.5,, 0.2-->minimum angle diff.
        {
          if (tmpNominalValue < 100) { // for redundant check
            next_val = next_val - 360.0;
          } else if (tmpNominalValue > 260) { // for redundant check
            cur_val = cur_val + 360.0;
          }
        }

        // in CCW, next value should be smaller than curVal
        if (tmpNominalValue >= next_val && tmpNominalValue <= cur_val) {
          const auto diffCur = fabs(tmpNominalValue - cur_val);
          const auto diffNext = fabs(tmpNominalValue - next_val);

          if (diffCur <= diffNext || enExcludingMode == 0 ||
              enExcludingMode == 1) {
            latest_Idx = j;
            vTargetIdx.push_back(latest_Idx);
          } else if (j != sizeFull - 2 && enExcludingMode == 2) {
            latest_Idx = j + 1;
            vTargetIdx.push_back(latest_Idx);
          } else {
            latest_Idx = j;
            // Skip to pushback
          }
          break;
        }
      }
    }
  }
  // vTargetIdx.push_back(sizeFull-1); //omit the last image --> should be same
  // as first...
}

bool IsFileNameOrderCorrect(std::vector<std::string> &vFileNames) {
  // regardless of whether number or hexa codes,
  // we can convert it from number to hexa number and compare the order

  const auto size = vFileNames.size();

  if (size < 2) {
    return false;
  }

  std::vector<int> arrNum(size);

  size_t index = 0;
  while (index < size) {
    std::filesystem::path crntFilePath{vFileNames.at(index++)};
    auto dir = std::filesystem::absolute(crntFilePath);
    auto file_basename = crntFilePath.stem();
    auto newBaseName = crl::HexStr2IntStr(file_basename.string());
    if (auto o_iName = crl::from_string<int>(newBaseName);
        o_iName.has_value()) {
      arrNum.push_back(o_iName.value());
    }
  }

  /*bool bOrderOK = true;
  for (int i = 0; i < size - 1; i++)
  {
      if (arrNum[i] >= arrNum[i + 1]) {
          bOrderOK = false;
              }
  }*/

  const auto index_of_nonascending =
      std::adjacent_find(arrNum.begin(), arrNum.end(), std::greater<>());

  return index_of_nonascending == arrNum.end();
}

void CopyDictionary(itk::MetaDataDictionary &fromDict,
                    itk::MetaDataDictionary &toDict) {
  using DictionaryType = itk::MetaDataDictionary;

  DictionaryType::ConstIterator itr = fromDict.Begin();
  const DictionaryType::ConstIterator end = fromDict.End();
  using MetaDataStringType = itk::MetaDataObject<std::string>;

  while (itr != end) {
    auto entry = itr->second;

    MetaDataStringType::Pointer entryvalue =
        dynamic_cast<MetaDataStringType *>(entry.GetPointer());
    if (entryvalue != nullptr) {
      auto tagkey = itr->first;
      auto tagvalue = entryvalue->GetMetaDataObjectValue();
      itk::EncapsulateMetaData<std::string>(toDict, tagkey, tagvalue);
    }
    ++itr;
  }
}

} // namespace crl