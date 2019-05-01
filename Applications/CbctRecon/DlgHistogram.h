#ifndef DLGHISTOGRAM_H
#define DLGHISTOGRAM_H

#include "ui_DlgHistogram.h"
#include <QDialog>

class CbctReconWidget;

class DlgHistogram : public QDialog, public Ui::DlgHistogramClass {
  Q_OBJECT

public slots:
  void SLT_DrawGraph() const;
  void SLT_ReDrawGraph_dial() const;
  void SLT_ReDrawGraph_limits() const;
  void SLT_ReturnCF() const;

public:
  DlgHistogram();
  explicit DlgHistogram(CbctReconWidget *parent);
  ~DlgHistogram() = default;

  CbctReconWidget *m_pParent{}; // to pull 3D images

private:
  Ui::DlgHistogramClass ui{};
};

#endif
