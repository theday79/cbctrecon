#pragma once
#include <QDialog>
#include <QString>
#include "ui_DlgExternalCommand.h"
#include "cbctrecon.h"
#include <QStringList>


using namespace std;

class DlgExternalCommand : public QDialog,
	public Ui::DlgExternalCommandClass
{
    Q_OBJECT
    ;

public slots:  
void SLT_SetRTKPath();
void SLT_GenRTKCommand();
void SLT_RunRTKCommand();
void SLT_SetRTKPathManual();

public:
	DlgExternalCommand();
	DlgExternalCommand(QWidget *parent);
	~DlgExternalCommand();

	int BuildRTKCommandFilter(); //pull predefined command names that you want. From some preset option file
	
	void SetRTKBinPath(QString& strDirPath);
public: 
    CbctRecon* m_pParent; //to pull 3D images 
	QString m_strDirRTKBin;

	QStringList m_listRTKCommandFilter;

	QString m_strRecentOutputPath;

private:
    Ui::DlgExternalCommandClass ui;
	
};
