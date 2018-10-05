// For testing CbctRecon
#if USE_TINYREFL
#include <tinyrefl/api.hpp>
#include "cbctrecon.h"
#include "cbctrecon.h.tinyrefl"
#else
#include "cbctrecon.h"
#endif

//#include <QtTest/QtTest>

#include "DlgRegistration.h"

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage:\n"
      << argv[0]
      << " ./dicom/directory"
      << std::endl;
    return -1;
  }

  //QApplication a(argc, argv);

  std::cout << "Running cbctrecon_test!" << std::endl;
  auto cbctrecon = std::make_unique<CbctRecon>();
  auto dcm_dir = QString(argv[1]);
  if (!cbctrecon->ReadDicomDir(dcm_dir)) {
    std::cerr << "Couldn't read DICOM" << std::endl;
    return 1;
  }

  std::cout << "image was read" << std::endl;

  /*cbctrecon->m_pDlgRegistration->UpdateVOICombobox(PLAN_CT);
  cbctrecon->m_pDlgRegistration->UpdateListOfComboBox(1);
  cbctrecon->m_pDlgRegistration->LoadImgFromComboBox(1, QString("REF_CT"));
  cbctrecon->m_pDlgRegistration->ui.comboBox_VOI->setCurrentIndex(1);

  std::cout << "DlgRegi. is ready for calculating WEPL for "
    << cbctrecon->m_pDlgRegistration->ui.comboBox_VOI->currentText().toStdString()
    << std::endl;
  auto start_time = std::chrono::steady_clock::now(); //clock();
  cbctrecon->m_pDlgRegistration->SLT_WEPLcalc();
  auto end_time = std::chrono::steady_clock::now();

  std::cout << "WEPL was calculated in: "
    << std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()
    << " ms"
    << std::endl;*/

  return 0;
}


/*
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
*/
