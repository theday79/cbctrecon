#pragma once
#include <QDialog>
#include <QString>
#include "ui_DlgHistogram.h"
#include "cbctrecon.h"



class DlgHistogram : public QDialog,
	public Ui::DlgHistogramClass
{
    Q_OBJECT
    ;

public slots:
	void SLT_DrawGraph();
	void SLT_ReDrawGraph_dial();
	void SLT_ReDrawGraph_limits();
	void SLT_ReturnCF();

public:
	DlgHistogram();
	DlgHistogram(QWidget *parent);
	~DlgHistogram();


public: 
    CbctRecon* m_pParent; //to pull 3D images 


private:
    Ui::DlgHistogramClass ui;
	
};
