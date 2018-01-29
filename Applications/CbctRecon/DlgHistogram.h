#pragma once
#include "cbctrecon.h"
#include "ui_DlgHistogram.h"
#include <QDialog>
#include <QString>



class DlgHistogram : public QDialog,
	public Ui::DlgHistogramClass
{
    Q_OBJECT

public slots:
	void SLT_DrawGraph();
	void SLT_ReDrawGraph_dial();
	void SLT_ReDrawGraph_limits();
	void SLT_ReturnCF();

public:
	DlgHistogram();
	DlgHistogram(QWidget *parent);
	~DlgHistogram() override;


public: 
    CbctRecon* m_pParent; //to pull 3D images 


private:
    Ui::DlgHistogramClass ui;
	
};
