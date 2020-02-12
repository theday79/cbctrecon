// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

// Qt
#include <QFileDialog>

// Local
#include "DlgExternalCommand.h"
#include "DlgRegistration.h"
#include "cbctrecon_io.h"
#include "cbctrecon_mainwidget.h"
#include "cbctregistration.h"
#include "free_functions.h"
#include "qtwrap.h"

using namespace std::literals;

DlgExternalCommand::DlgExternalCommand() {
  /* Sets up the GUI */
  ui.setupUi(this);
}

DlgExternalCommand::DlgExternalCommand(QWidget *parent) : QDialog(parent) {
  /* Sets up the GUI */
  ui.setupUi(this);
  m_pParent = dynamic_cast<CbctReconWidget *>(parent);

  // int len = BuildRTKCommandFilter();
}

DlgExternalCommand::~DlgExternalCommand() = default;

void DlgExternalCommand::SLT_SetRTKPath() {
  auto dirPath = to_path(QFileDialog::getExistingDirectory(
      this, tr("Open RTK bin Directory"), "",
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks));

  if (!dirPath.empty()) {
    SetRTKBinPath(dirPath);
  }
}

void DlgExternalCommand::SetRTKBinPath(const fs::path &strDirPath) {
  // If DIr is exist

  if (!fs::exists(strDirPath)) {
    std::cout << "ERROR! " << strDirPath << " doesn't exist." << std::endl;
    return;
  }

  ui.comboBoxRTKOption->clear();

  m_strDirRTKBin = strDirPath;
  ui.plainTextRTKPath->setPlainText(to_qstr(m_strDirRTKBin));

  // Search for available commands
  if (m_strDirRTKBin.empty()) {
    return;
  }

  const auto cnt = m_listRTKCommandFilter.length();

  for (auto i = 0; i < cnt; i++) {
    const auto strCommandFilter = m_listRTKCommandFilter.at(i);
    auto tmpStrPath = m_strDirRTKBin;
    tmpStrPath = tmpStrPath / (strCommandFilter.toStdString() + ".exe");

    std::cout << m_strDirRTKBin << std::endl;
    std::cout << tmpStrPath << std::endl;

    if (fs::exists(tmpStrPath)) // add combo
    {
      ui.comboBoxRTKOption->addItem(strCommandFilter);
    }
  }
}

void DlgExternalCommand::SLT_GenRTKCommand() {
  auto crntPath = ui.plainTextRTKPath->toPlainText();
  auto crntCommand = ui.comboBoxRTKOption->currentText();

  fs::path strFinalCommandText;

  if (crntCommand.length() > 0 && crntPath.length() > 0) {
    strFinalCommandText = to_path(crntPath) / crntCommand.toStdString();
  } else {
    return;
  }

  // QString strHardware;
  // QString strTruncation;
  // QString strHann;
  // QString strHannY;

  // QString strProjPath;
  // QString strRegExp;
  // QString strOutOrigin; //default: centered
  // QString strOutDimension;
  // QString strOutSpacing;
  // QString strOutDirection; //no clue about what it is

  // geometry
  auto str_mainGeometry =
      m_pParent->ui.lineEdit_ElektaGeomPath->text().toStdString();

  if (str_mainGeometry.size() < 1) {
    std::cout << "Command will not be valid. set geometry file path in the "
                 "main UI first."
              << std::endl;
  }

  std::string str_mainHardware;
  if (m_pParent->ui.radioButton_UseCPU->isChecked()) {
    str_mainHardware = "cpu";
  } else if (m_pParent->ui.radioButton_UseCUDA->isChecked()) {
    str_mainHardware = "cuda";
  } else if (m_pParent->ui.radioButton_UseOpenCL->isChecked()) {
    str_mainHardware = "opencl";
  }

  // QString str_mainTruncation;
  const auto f_mainTrunc =
      m_pParent->ui.lineEdit_Ramp_TruncationCorrection->text().toDouble();
  const auto f_mainHann =
      m_pParent->ui.lineEdit_Ramp_HannCut->text().toDouble();
  const auto f_mainHannY =
      m_pParent->ui.lineEdit_Ramp_HannCutY->text().toDouble();

  const auto str_mainProjPath =
      m_pParent->ui.lineEdit_HisDirPath->text().toStdString();

  const auto str_mainProjRegExp = ".*.his"s;

  const auto str_mainDimension = crl::make_sep_str<','>(
      m_pParent->ui.lineEdit_outImgDim_AP->text().toStdString(),
      m_pParent->ui.lineEdit_outImgDim_SI->text().toStdString(),
      m_pParent->ui.lineEdit_outImgDim_LR->text().toStdString());

  const auto str_mainSpacing = crl::make_sep_str<','>(
      m_pParent->ui.lineEdit_outImgSp_AP->text().toStdString(),
      m_pParent->ui.lineEdit_outImgSp_SI->text().toStdString(),
      m_pParent->ui.lineEdit_outImgSp_LR->text().toStdString());

  auto curTime = QTime::currentTime();
  const auto strTimeStamp = curTime.toString("hhmmss");
  const auto &tmpPlmDir =
      m_pParent->m_dlgRegistration->m_cbctregistration->m_strPathPlastimatch;

  if (!fs::exists(tmpPlmDir)) {
    std::cout << "Error! No tmp plm path is available."
              << fs::absolute(tmpPlmDir) << std::endl;
    return;
  }

  const auto strOutput =
      fs::absolute(tmpPlmDir) / ("ExternalRtk_" + crntCommand.toStdString() +
                                 "_" + strTimeStamp.toStdString() + ".mha");

  if (crntCommand == "rtkfdk") {
    // For FDK, this is sure!
    strFinalCommandText = crl::make_sep_str<' '>(
        strFinalCommandText, "--geometry", str_mainGeometry, "--path",
        str_mainProjPath, "--verbose", "--regexp", str_mainProjRegExp,
        "--output", strOutput, "--spacing", str_mainSpacing, "--dimension",
        str_mainDimension, "--hardware", str_mainHardware, "--pad", f_mainTrunc,
        "--hann", f_mainHann, "--hannY", f_mainHannY);
  }
  // ui
  /*lineEditIteration
          lineEditSARTlamda
          lineEditSARTpositivity
          lineEditSARTsubsetproj*/

  else if (crntCommand == "rtksart") {
    const auto strIteration = ui.lineEditIteration->text()
                                  .trimmed()
                                  .toStdString(); // niterations default 5
    const auto strLamda =
        ui.lineEditSARTlamda->text()
            .trimmed()
            .toStdString(); // Convergence factor : default 0.3
    const auto strPositivity = ui.lineEditSARTpositivity->text()
                                   .trimmed()
                                   .toStdString(); // Enforces positivity
                                                   // during the reconstruction
                                                   // (default=off)",
    const auto strNprojpersubset = ui.lineEditSARTsubsetproj->text()
                                       .trimmed()
                                       .toStdString(); // Number of projections
                                                       // processed between each
                                                       // update of the
                                                       // reconstructed volume
                                                       // (1 for SART, several
                                                       // for OSSART, all for
                                                       // SIRT) (default=`1')",
    auto strFwdMethod = "Joseph"s;
    if (str_mainHardware == "cuda") {
      strFwdMethod = "CudaRayCast"s;
    }
    //			--fp = ENUM Forward projection method(possible values =
    //\"Joseph\",
    //\"RayCastInterpolator\", "CudaRayCast\" default=`Joseph')",

    auto strBackMethod = "VoxelBasedBackProjection"s;
    if (str_mainHardware == "cuda") {
      strBackMethod = "CudaVoxelBased"s;
    }
    strFinalCommandText = crl::make_sep_str<' '>(
        strFinalCommandText, "--geometry", str_mainGeometry, "--path",
        str_mainProjPath, "--verbose", "--regexp", str_mainProjRegExp,
        "--spacing", str_mainSpacing, "--dimension", str_mainDimension,
        "--niterations", strIteration, "--lambda", strLamda, "--positivity",
        strPositivity, "--nprojpersubset", strNprojpersubset, "--fp",
        strFwdMethod, "--bp", strBackMethod, "--time", "on", "--output",
        strOutput);
  } else if (crntCommand == "rtkadmmtotalvariation") {
    // lineEditIteration
    // lineEditTValpha
    //			lineEditTVbeta
    //		lineEditTVCGiter* /

    const auto strIteration = ui.lineEditIteration->text()
                                  .trimmed()
                                  .toStdString(); // niterations default 5
    const auto strTValpha =
        ui.lineEditTValpha->text()
            .trimmed()
            .toStdString(); // Convergence factor : default 0.3
    const auto strTVbeta = ui.lineEditTVbeta->text()
                               .trimmed()
                               .toStdString(); // Enforces positivity during the
                                               // reconstruction (default=off)",
    const auto strTVCGiter =
        ui.lineEditTVCGiter->text()
            .trimmed()
            .toStdString(); // Enforces positivity during the
                            // reconstruction (default=off)",
    auto strFwdMethod = "Joseph"s;
    if (str_mainHardware == "cuda") {
      strFwdMethod = "CudaRayCast"s;
    }
    //			--fp = ENUM Forward projection method(possible values =
    //\"Joseph\",
    //\"RayCastInterpolator\", "CudaRayCast\" default=`Joseph')",

    auto strBackMethod = "VoxelBasedBackProjection"s;
    if (str_mainHardware == "cuda") {
      strBackMethod = "CudaVoxelBased"s;
    }
    strFinalCommandText = crl::make_sep_str<' '>(
        strFinalCommandText, "--geometry", str_mainGeometry, "--path",
        str_mainProjPath, "--verbose", "--regexp", str_mainProjRegExp,
        "--spacing", str_mainSpacing, "--dimension", str_mainDimension,
        "--niterations", strIteration, "--alpha", strTValpha, "--beta",
        strTVbeta, "--CGiter", strTVCGiter, "--fp", strFwdMethod, "--bp",
        strBackMethod, "--time", "on", "--output", strOutput);
  } else if (crntCommand == "rtkadmmwavelets") {
    const auto strIteration = ui.lineEditIteration->text()
                                  .trimmed()
                                  .toStdString(); // niterations default 5
    const auto strTValpha =
        ui.lineEditTValpha->text()
            .trimmed()
            .toStdString(); // Convergence factor : default 0.3
    const auto strTVbeta = ui.lineEditTVbeta->text()
                               .trimmed()
                               .toStdString(); // Enforces positivity during the
                                               // reconstruction (default=off)",
    const auto strTVCGiter =
        ui.lineEditTVCGiter->text()
            .trimmed()
            .toStdString(); // Enforces positivity during the
                            // reconstruction (default=off)",

    const auto strWVorder = ui.lineEditWVorder->text().trimmed().toStdString();
    const auto strWVlevel = ui.lineEditWVlevel->text().trimmed().toStdString();
    //"      --order=INT         The order of the Daubechies wavelets
    //(default=`3')", "      --levels=INT        The number of decomposition
    // levels in the wavelets \n                            transform
    //(default=`5')",

    auto strFwdMethod = "Joseph"s;
    if (str_mainHardware == "cuda") {
      strFwdMethod = "CudaRayCast"s;
    }
    //			--fp = ENUM Forward projection method(possible values =
    //\"Joseph\",
    //\"RayCastInterpolator\", "CudaRayCast\" default=`Joseph')",

    auto strBackMethod = "VoxelBasedBackProjection"s;
    if (str_mainHardware == "cuda") {
      strBackMethod = "CudaVoxelBased"s;
    }
    strFinalCommandText = crl::make_sep_str<' '>(
        strFinalCommandText, "--geometry", str_mainGeometry, "--path",
        str_mainProjPath, "--verbose", "--regexp", str_mainProjRegExp,
        "--spacing", str_mainSpacing, "--dimension", str_mainDimension,
        "--niterations", strIteration, "--alpha", strTValpha, "--beta",
        strTVbeta, "--CGiter", strTVCGiter, "--fp", strFwdMethod, "--bp",
        strBackMethod, "--order", strWVorder, "--levels", strWVlevel, "--time",
        "on", "--output", strOutput);
  }

  m_strRecentOutputPath = strOutput;
  ui.plainTextRTKCommandLine->setPlainText(to_qstr(strFinalCommandText));
}

void DlgExternalCommand::SLT_RunRTKCommand() {
  const auto strFinalExternalCommand =
      ui.plainTextRTKCommandLine->toPlainText().toStdString();
  if (std::system(strFinalExternalCommand.c_str()) < 0) {
    std::cerr << "Failed to run\n";
  }

  std::cout << "External RTK reconstruction is done" << std::endl;
  std::cout << "File is being loaded" << std::endl;

  m_pParent->m_cbctrecon->LoadExternalFloatImage(
      m_strRecentOutputPath,
      true); // true: conversion (float, direction

  if (m_pParent->ui.checkBox_PostMedianOn->isChecked()) {
    UShortImageType::SizeType indexRadius{};
    indexRadius[0] =
        m_pParent->ui.lineEdit_PostMedSizeX->text().toInt(); // radius along x
    indexRadius[1] =
        m_pParent->ui.lineEdit_PostMedSizeY->text().toInt(); // radius along y
    indexRadius[2] =
        m_pParent->ui.lineEdit_PostMedSizeZ->text().toInt(); // radius along y
    m_pParent->m_cbctrecon->MedianFilterByGUI(
        indexRadius); // applied to raw image
  }

  auto outputFilePath =
      to_path(this->m_pParent->ui.lineEdit_OutputFilePath->text());
  auto outFileDir = fs::absolute(outputFilePath);

  // bool b = outFileDir.exists();
  // QString tmpPath = outFileDir.absolutePath();

  if (outputFilePath.empty() || !fs::exists(outFileDir)) {
    std::cout << "No available output path. Should be exported later"
              << std::endl;
  } else {
    crl::saveImageAsMHA<UShortImageType>(
        this->m_pParent->m_cbctrecon->m_spRawReconImg, outputFilePath.string());

    std::cout << "Wrote the image to: " << outputFilePath << std::endl;
  }
}

int DlgExternalCommand::BuildRTKCommandFilter() // called when it is created
{
  // Search for some preset option file, later.
  // Temporarily it is hard coded.
  m_listRTKCommandFilter.clear();

  m_listRTKCommandFilter.push_back("rtksart");
  m_listRTKCommandFilter.push_back("rtkfdk");
  m_listRTKCommandFilter.push_back("rtkadmmtotalvariation");
  m_listRTKCommandFilter.push_back("rtkadmmwavelets");
  // m_listRTKCommandFilter.push_back("rtktotalvariationdenoising");
  // m_listRTKCommandFilter.push_back("rtkwaveletsdenoising");

  return m_listRTKCommandFilter.length();
}

void DlgExternalCommand::SLT_SetRTKPathManual() // apply button
{
  auto tmpPlainText = ui.plainTextRTKPath->toPlainText();
  SetRTKBinPath(to_path(tmpPlainText));
}
