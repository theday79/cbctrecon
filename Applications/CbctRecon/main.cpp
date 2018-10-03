#include "cbctrecon_mainwidget.h"

int main(int argc, char *argv[]) {
  QApplication a(argc, argv);
  // Load color theme if it exists:
  QString app_path = QCoreApplication::applicationDirPath();
  QFile f(app_path + "/ui/style.qss");
  if (!f.exists())
  {
      std::cout << "Unable to set stylesheet, file: "
                << f.fileName().toStdString()
                << " not found\n"
                << std::endl;
  }
  else
  {
      f.open(QFile::ReadOnly | QFile::Text);
      QTextStream ts(&f);
      qApp->setStyleSheet(ts.readAll());
  }

  CbctReconWidget w;
  w.show();
  return a.exec();
}
