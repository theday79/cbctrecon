#include "cbctrecon.h"
#include <QtWidgets/QApplication>

int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	CbctRecon w;
	w.show();
	return a.exec();
}
