#ifndef DLGEXTERNALCOMMAND_H
#define DLGEXTERNALCOMMAND_H

#include <filesystem>

// Qt
#include <QDialog>
#include <QStringList>

// Local
#include "ui_DlgExternalCommand.h"

namespace fs = std::filesystem;

class CbctReconWidget;

class DlgExternalCommand : public QDialog, public Ui::DlgExternalCommandClass {
  Q_OBJECT

public slots:
  void SLT_SetRTKPath();
  void SLT_GenRTKCommand();
  void SLT_RunRTKCommand();
  void SLT_SetRTKPathManual();

public:
  DlgExternalCommand();
  explicit DlgExternalCommand(QWidget *parent);
  ~DlgExternalCommand() override;
  DlgExternalCommand(const DlgExternalCommand &) = delete;
  void operator=(const DlgExternalCommand &) = delete;
  DlgExternalCommand(DlgExternalCommand &&) = delete;
  void operator=(DlgExternalCommand &&) = delete;

  int BuildRTKCommandFilter(); // pull predefined command names that you want.
                               // From some preset option file

  void SetRTKBinPath(const fs::path &strDirPath);

  CbctReconWidget *m_pParent{}; // to pull 3D images
  fs::path m_strDirRTKBin;

  QStringList m_listRTKCommandFilter;

  fs::path m_strRecentOutputPath;

private:
  Ui::DlgExternalCommandClass ui{};
};

#endif
