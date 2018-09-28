// For testing CbctRecon
#if USE_TINYREFL
#include <tinyrefl/api.hpp>
#include "cbctrecon.h"
#include "cbctrecon.h.tinyrefl"
#else
#include "cbctrecon.h"
#endif

#include <QtTest/QtTest>

#include "DlgRegistration.h"

class TestCbctRecon : public QObject
{
  Q_OBJECT

private slots:
  void TestWEPLContours();

private:
  CbctRecon app;
};

void TestCbctRecon::TestWEPLContours() {
  std::cout << "Running cbctrecon_test!" << std::endl;
  auto dcm_dir = QString("lol");
  if (!app.ReadDicomDir(dcm_dir)) {
    std::cerr << "Couldn't read DICOM" << std::endl;
    return;
  }

  std::cout << "image was read" << std::endl;

  app.m_pDlgRegistration->UpdateVOICombobox(PLAN_CT);
  app.m_pDlgRegistration->UpdateListOfComboBox(1);
  auto img_type = QString("REF_CT");
  app.m_pDlgRegistration->LoadImgFromComboBox(1, img_type);
  app.m_pDlgRegistration->ui.comboBox_VOI->setCurrentIndex(1);

  std::cout << "DlgRegi. is ready for calculating WEPL for "
    << app.m_pDlgRegistration->ui.comboBox_VOI->currentText().toStdString()
    << std::endl;
  auto start_time = std::chrono::steady_clock::now(); //clock();
  app.m_pDlgRegistration->SLT_WEPLcalc();
  auto end_time = std::chrono::steady_clock::now();

  std::cout << "WEPL was calculated in: "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
    << " ms"
    << std::endl;
}

QTEST_MAIN(TestCbctRecon)

#include "cbctrecon_test.moc"
