// Qt
#include <QFileDialog>
#include <QProcess>

// Local
#include "DlgExternalCommand.h"
#include "DlgRegistration.h"
#include "cbctrecon_mainwidget.h"
#include "cbctregistration.h"

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
  QString dirPath = QFileDialog::getExistingDirectory(
      this, tr("Open RTK bin Directory"), "",
      QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

  if (dirPath.length() > 1) {
    SetRTKBinPath(dirPath);
  }
}

void DlgExternalCommand::SetRTKBinPath(QString &strDirPath) {
  // If DIr is exist
  QDir dirInfo = QDir(strDirPath);

  if (!dirInfo.exists()) {
    std::cout << "ERROR! " << strDirPath.toLocal8Bit().constData()
              << " doesn't exist." << std::endl;
    return;
  }

  ui.comboBoxRTKOption->clear();

  m_strDirRTKBin = strDirPath;
  ui.plainTextRTKPath->setPlainText(m_strDirRTKBin);

  // Search for available commands
  if (m_strDirRTKBin.isEmpty()) {
    return;
  }

  const auto cnt = m_listRTKCommandFilter.length();

  for (int i = 0; i < cnt; i++) {
    const auto strCommandFilter = m_listRTKCommandFilter.at(i);
    auto tmpStrPath = m_strDirRTKBin;
    tmpStrPath = tmpStrPath.append("/").append(strCommandFilter).append(".exe");

    std::cout << m_strDirRTKBin.toLocal8Bit().constData() << std::endl;
    std::cout << tmpStrPath.toLocal8Bit().constData() << std::endl;

    QFileInfo fInfo = QFileInfo(tmpStrPath);
    if (fInfo.exists()) // add combo
    {
      ui.comboBoxRTKOption->addItem(strCommandFilter);
    }
  }
}

void DlgExternalCommand::SLT_GenRTKCommand() {
  QString crntPath = ui.plainTextRTKPath->toPlainText();
  QString crntCommand = ui.comboBoxRTKOption->currentText();

  QString strFinalCommandText;

  if (crntCommand.length() > 0 && crntPath.length() > 0) {
    strFinalCommandText = crntPath.append("/").append(crntCommand);
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
  QString str_mainGeometry = m_pParent->ui.lineEdit_ElektaGeomPath->text();

  if (str_mainGeometry.length() < 1) {
    std::cout << "Command will not be valid. set geometry file path in the "
                 "main UI first."
              << std::endl;
  }

  QString str_mainHardware;
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

  const auto str_mainProjPath = m_pParent->ui.lineEdit_HisDirPath->text();

  const auto str_mainProjRegExp = ".*.his";

  const auto str_mainDimension =
      QString("%1,%2,%3")
          .arg(m_pParent->ui.lineEdit_outImgDim_AP->text())
          .arg(m_pParent->ui.lineEdit_outImgDim_SI->text())
          .arg(m_pParent->ui.lineEdit_outImgDim_LR->text());

  const auto str_mainSpacing =
      QString("%1,%2,%3")
          .arg(m_pParent->ui.lineEdit_outImgSp_AP->text())
          .arg(m_pParent->ui.lineEdit_outImgSp_SI->text())
          .arg(m_pParent->ui.lineEdit_outImgSp_LR->text());

  QTime curTime = QTime::currentTime();
  const auto strTimeStamp = curTime.toString("hhmmss");
  QDir tmpPlmDir = QDir(
      m_pParent->m_dlgRegistration->m_cbctregistration->m_strPathPlastimatch);

  if (!tmpPlmDir.exists()) {
    std::cout << "Error! No tmp plm path is available."
              << tmpPlmDir.absolutePath().toLocal8Bit().constData()
              << std::endl;
    return;
  }

  const auto strOutput = tmpPlmDir.absolutePath() + "/" + "ExternalRtk_" +
                         crntCommand + "_" + strTimeStamp + ".mha";

  if (crntCommand == "rtkfdk") {
    // clang-format off
		//For FDK, this is sure!	
	strFinalCommandText = strFinalCommandText +
		" --geometry " + str_mainGeometry +
		" --path " + str_mainProjPath +
		" --verbos" +
		" --regexp " + str_mainProjRegExp +
		" --output " + strOutput +
		" --spacing " + str_mainSpacing +
		" --dimension " + str_mainDimension +

		" --hardware " + str_mainHardware +
		" --pad " + QString("%1").arg(f_mainTrunc) +
		" --hann " + QString("%1").arg(f_mainHann) +
		" --hannY " + QString("%1").arg(f_mainHannY);
    // clang-format on
  }
  // ui
  /*lineEditIteration
          lineEditSARTlamda
          lineEditSARTpositivity
          lineEditSARTsubsetproj*/

  else if (crntCommand == "rtksart") {
    const auto strIteration =
        ui.lineEditIteration->text().trimmed(); // niterations default 5
    const auto strLamda = ui.lineEditSARTlamda->text()
                           .trimmed(); // Convergence factor : default 0.3
    const auto strPositivity = ui.lineEditSARTpositivity->text()
                                .trimmed(); // Enforces positivity
                                            // during the reconstruction
                                            // (default=off)",
    const auto strNprojpersubset =
        ui.lineEditSARTsubsetproj->text().trimmed(); // Number of projections
                                                     // processed between each
                                                     // update of the
                                                     // reconstructed volume (1
                                                     // for SART, several for
                                                     // OSSART, all for SIRT)
                                                     // (default=`1')",
    QString strFwdMethod = "Joseph";
    if (str_mainHardware == "cuda") {
      strFwdMethod = "CudaRayCast";
    }
    //			--fp = ENUM Forward projection method(possible values =
    //\"Joseph\",
    //\"RayCastInterpolator\", "CudaRayCast\" default=`Joseph')",

    QString strBackMethod = "VoxelBasedBackProjection";
    if (str_mainHardware == "cuda") {
      strBackMethod = "CudaVoxelBased";
    }
    // clang-format off
	strFinalCommandText = strFinalCommandText +
		" --geometry " + str_mainGeometry +
		" --path " + str_mainProjPath +
		" --verbos" +
		" --regexp " + str_mainProjRegExp +			
		" --spacing " + str_mainSpacing +
		" --dimension " + str_mainDimension +

		" --niterations " + strIteration +
		" --lambda " + strLamda +
		" --positivity " + strPositivity +
		" --nprojpersubset " + strNprojpersubset +
		" --fp " + strFwdMethod +
		" --bp " + strBackMethod +
		" --time " + "on" +
		" --output " + strOutput;
    // clang-format on
  } else if (crntCommand == "rtkadmmtotalvariation") {
    // lineEditIteration
    // lineEditTValpha
    //			lineEditTVbeta
    //		lineEditTVCGiter* /

    const auto strIteration =
        ui.lineEditIteration->text().trimmed(); // niterations default 5
    const auto strTValpha = ui.lineEditTValpha->text()
                             .trimmed(); // Convergence factor : default 0.3
    const auto strTVbeta =
        ui.lineEditTVbeta->text().trimmed(); // Enforces positivity during the
                                             // reconstruction (default=off)",
    const auto strTVCGiter =
        ui.lineEditTVCGiter->text().trimmed(); // Enforces positivity during the
                                               // reconstruction (default=off)",
    QString strFwdMethod = "Joseph";
    if (str_mainHardware == "cuda") {
      strFwdMethod = "CudaRayCast";
    }
    //			--fp = ENUM Forward projection method(possible values =
    //\"Joseph\",
    //\"RayCastInterpolator\", "CudaRayCast\" default=`Joseph')",

    QString strBackMethod = "VoxelBasedBackProjection";
    if (str_mainHardware == "cuda") {
      strBackMethod = "CudaVoxelBased";
    }
    // clang-format off
	strFinalCommandText = strFinalCommandText +
		" --geometry " + str_mainGeometry +
		" --path " + str_mainProjPath +
		" --verbos" +
		" --regexp " + str_mainProjRegExp +			
		" --spacing " + str_mainSpacing +
		" --dimension " + str_mainDimension +
		" --niterations " + strIteration +
		" --alpha " + strTValpha +
		" --beta " + strTVbeta +
		" --CGiter " + strTVCGiter +
		" --fp " + strFwdMethod +
		" --bp " + strBackMethod +
		" --time " + "on" +
		" --output " + strOutput;
    // clang-format on
  } else if (crntCommand == "rtkadmmwavelets") {
    const auto strIteration =
        ui.lineEditIteration->text().trimmed(); // niterations default 5
    const auto strTValpha = ui.lineEditTValpha->text()
                             .trimmed(); // Convergence factor : default 0.3
    const auto strTVbeta =
        ui.lineEditTVbeta->text().trimmed(); // Enforces positivity during the
                                             // reconstruction (default=off)",
    const auto strTVCGiter =
        ui.lineEditTVCGiter->text().trimmed(); // Enforces positivity during the
                                               // reconstruction (default=off)",

    const auto strWVorder = ui.lineEditWVorder->text().trimmed();
    const auto strWVlevel = ui.lineEditWVlevel->text().trimmed();
    //"      --order=INT         The order of the Daubechies wavelets
    //(default=`3')", "      --levels=INT        The number of decomposition
    // levels in the wavelets \n                            transform
    //(default=`5')",

    QString strFwdMethod = "Joseph";
    if (str_mainHardware == "cuda") {
      strFwdMethod = "CudaRayCast";
    }
    //			--fp = ENUM Forward projection method(possible values =
    //\"Joseph\",
    //\"RayCastInterpolator\", "CudaRayCast\" default=`Joseph')",

    QString strBackMethod = "VoxelBasedBackProjection";
    if (str_mainHardware == "cuda") {
      strBackMethod = "CudaVoxelBased";
    }
    // clang-format off
	strFinalCommandText = strFinalCommandText +
		" --geometry " + str_mainGeometry +
		" --path " + str_mainProjPath +
		" --verbos" +
		" --regexp " + str_mainProjRegExp +		
		" --spacing " + str_mainSpacing +
		" --dimension " + str_mainDimension +
		" --niterations " + strIteration +
		" --alpha " + strTValpha +
		" --beta " + strTVbeta +
		" --CGiter " + strTVCGiter +
		" --fp " + strFwdMethod +
		" --bp " + strBackMethod +
		" --order " + strWVorder +
		" --levels " + strWVlevel +
		" --time " + "on" +
		" --output " + strOutput;
    // clang-format on
  }

  m_strRecentOutputPath = strOutput;
  ui.plainTextRTKCommandLine->setPlainText(strFinalCommandText);
}

void DlgExternalCommand::SLT_RunRTKCommand() {
  const auto strFinalExternalCommand =
      ui.plainTextRTKCommandLine->toPlainText();
  if (QProcess::execute(strFinalExternalCommand) < 0) {
    qDebug() << "Failed to run";
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

  m_pParent->FileExportByGUI(); // applied to raw image
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

  return (m_listRTKCommandFilter.length());
}

void DlgExternalCommand::SLT_SetRTKPathManual() // apply button
{
  QString tmpPlainText = ui.plainTextRTKPath->toPlainText();
  SetRTKBinPath(tmpPlainText);
}
