#ifndef DLGEXTERNALCOMMAND_H
#define DLGEXTERNALCOMMAND_H

// Qt
#include <QDialog>
#include <QString>
#include <QStringList>

// Local
#include "ui_DlgExternalCommand.h"

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
  DlgExternalCommand(QWidget *parent);
  ~DlgExternalCommand() override;

  int BuildRTKCommandFilter(); // pull predefined command names that you want.
                               // From some preset option file

  void SetRTKBinPath(QString &strDirPath);

public:
  CbctReconWidget *m_pParent{}; // to pull 3D images
  QString m_strDirRTKBin;

  QStringList m_listRTKCommandFilter;

  QString m_strRecentOutputPath;

private:
  Ui::DlgExternalCommandClass ui{};
};

#endif
