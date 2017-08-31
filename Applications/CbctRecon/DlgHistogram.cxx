#include "DlgHistogram.h"
#include "cbctrecon.h"

#define MAX_USHORT 65535U

DlgHistogram::DlgHistogram() : QDialog()
{
    /* Sets up the GUI */
    ui.setupUi (this);
}

DlgHistogram::DlgHistogram(QWidget *parent) : QDialog(parent)
{
    /* Sets up the GUI */
    ui.setupUi (this);
    m_pParent = (CbctRecon*)(parent);
}

DlgHistogram::~DlgHistogram()
{
}

void threaded_calculate_histogram(unsigned int* histogram_out, const int n_bins, const size_t bin_size,
	const unsigned short* raw_data, const size_t size_raw_data) {
	assert(n_bins > 0);

	size_t n_threads = 16U;
	while ((size_raw_data % n_threads != 0U) && (n_threads >= 1U))
		n_threads /= 2U;

#pragma omp parallel for shared(histogram_out)
	for (int i = 0; i < (int)n_threads; i++) {

		unsigned int* histogram_local = new unsigned int[(size_t)n_bins];

#pragma omp parallel for // initialize to 0
		for (int j = 0; j < n_bins; j++)
			histogram_local[j] = 0;

		for (size_t j = (i*size_raw_data / n_threads); j < ((i + 1)*size_raw_data / n_threads); j++)
			histogram_local[(size_t)(raw_data[j]) / bin_size]++;

#pragma omp critical
		{
#pragma omp parallel for // initialize to 0
			for (int j = 0; j < n_bins; j++)
				histogram_out[j] += histogram_local[j];
		}

		delete[] histogram_local;
	}

}

void DlgHistogram::SLT_DrawGraph() //based on profile
{
	if (!m_pParent->m_spProjImgRaw3D || !m_pParent->m_spProjImgCT3D)
		return;

	const unsigned short* raw_data = m_pParent->m_spProjImgRaw3D->GetBufferPointer();
	const size_t size_raw_data = 
		m_pParent->m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize()[0] *
		m_pParent->m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize()[1] *
		m_pParent->m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize()[2];

	const unsigned short* ct_data = m_pParent->m_spProjImgCT3D->GetBufferPointer();
	const size_t size_ct_data =
		m_pParent->m_spProjImgCT3D->GetLargestPossibleRegion().GetSize()[0] *
		m_pParent->m_spProjImgCT3D->GetLargestPossibleRegion().GetSize()[1] *
		m_pParent->m_spProjImgCT3D->GetLargestPossibleRegion().GetSize()[2];

	const size_t n_bins = (size_t) ui.spinBox_nBins->value();
	const size_t bin_size = MAX_USHORT / n_bins;

	unsigned int* histogram_raw = new unsigned int[n_bins];
	unsigned int* histogram_ct = new unsigned int[n_bins]; // same size should be enough, because raw projections are normalised at read time.

	// initialize to 0
#pragma omp parallel sections
	{
#pragma omp section
		{
			for (size_t i = 0; i < n_bins; i++)
				histogram_raw[i] = 0;
		}
#pragma omp section
		{
			for (size_t i = 0; i < n_bins; i++)
				histogram_ct[i] = 0;
		}
	}

#pragma omp parallel sections
	{
#pragma omp section
		{
			threaded_calculate_histogram(histogram_raw, (int)n_bins, bin_size, raw_data, size_raw_data);
		}
#pragma omp section
		{
			threaded_calculate_histogram(histogram_ct, (int)n_bins, bin_size, ct_data, size_raw_data);
		}
	}

	QVector<double> vAxisX; //can be rows or columns
	QVector<double> vAxisY_raw;
	QVector<double> vAxisY_ct;

	ui.customPlot->clearGraphs();

	const double minX = 0.0;
	double maxX = 0.0;
	const double minY = 0.0;
	double maxY = 0.0;

	for (int i = 0; i< (int)n_bins; i++)
	{
		if (maxY < histogram_raw[i])
			maxY = histogram_raw[i];
		if (histogram_raw[i] > 0)
			maxX = i * bin_size;
		vAxisX.push_back(bin_size * i);
		vAxisY_raw.push_back(histogram_raw[i]);
		vAxisY_ct.push_back(histogram_ct[i]);
	}

	ui.customPlot->addGraph();
	ui.customPlot->graph(0)->setData(vAxisX, vAxisY_raw);
	ui.customPlot->graph(0)->setPen(QPen(Qt::red));
	ui.customPlot->graph(0)->setName("Raw proj. histogram");

	ui.customPlot->addGraph();
	ui.customPlot->graph(1)->setData(vAxisX, vAxisY_ct);
	ui.customPlot->graph(1)->setPen(QPen(Qt::blue));
	ui.customPlot->graph(1)->setName("CT fwd. proj. histogram");

	ui.lineEditXMin->setText(QString("%1").arg(minX));
	ui.lineEditXMax->setText(QString("%1").arg(maxX));
	ui.lineEditYMin->setText(QString("%1").arg(minY));
	ui.lineEditYMax->setText(QString("%1").arg(maxY));

	ui.customPlot->xAxis->setRange(minX, maxX);
	ui.customPlot->yAxis->setRange(minY, maxY);

	ui.customPlot->xAxis->setLabel("value");
	ui.customPlot->yAxis->setLabel("N");
	ui.customPlot->setTitle("Intensity histograms");
	QFont titleFont = font();
	titleFont.setPointSize(10);
	ui.customPlot->setTitleFont(titleFont);

	ui.customPlot->legend->setVisible(true);
	QFont legendFont = font();  // start out with MainWindow's font..
	legendFont.setPointSize(9); // and make a bit smaller for legend
	ui.customPlot->legend->setFont(legendFont);
	ui.customPlot->legend->setPositionStyle(QCPLegend::psTopRight);
	ui.customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 200)));

	ui.customPlot->replot();
	
	delete[] histogram_raw;
	delete[] histogram_ct;
}


void DlgHistogram::SLT_ReDrawGraph_limits() //based on profile
{

	double tmpXMin = ui.lineEditXMin->text().toDouble();
	double tmpXMax = ui.lineEditXMax->text().toDouble();
	double tmpYMin = ui.lineEditYMin->text().toDouble();
	double tmpYMax = ui.lineEditYMax->text().toDouble();

	ui.customPlot->xAxis->setRange(tmpXMin, tmpXMax);
	ui.customPlot->yAxis->setRange(tmpYMin, tmpYMax);

	ui.customPlot->replot();
}

void threaded_calculate_histogram_with_scaling(unsigned int* histogram_out, const int n_bins, const size_t bin_size, const float scaling, 
	const unsigned short* raw_data, const size_t size_raw_data) {

	size_t n_threads = 16U;
	while ((size_raw_data % n_threads != 0U) && (n_threads >= 1U))
		n_threads /= 2U;

	// 8 threads is probably ok for most systems
#pragma omp parallel for
	for (int i = 0; i < (int)n_threads; i++) {

		unsigned int* histogram_local = new unsigned int[(size_t)n_bins];
#pragma omp parallel for // initialize to 0
		for (int j = 0; j < n_bins; j++)
			histogram_local[j] = 0;

		for (size_t j = (i*size_raw_data / n_threads); j < ((i + 1)*size_raw_data / n_threads); j++)
			histogram_local[(size_t)(raw_data[j] * scaling) / bin_size] += 1;

#pragma omp critical
		{
#pragma omp parallel for // initialize to 0
			for (int j = 0; j < n_bins; j++)
				histogram_out[j] += histogram_local[j];
		}

		delete[] histogram_local;
	}

}


void DlgHistogram::SLT_ReDrawGraph_dial() //based on profile
{
	if (!m_pParent->m_spProjImgRaw3D || !m_pParent->m_spProjImgCT3D)
		return;

	float scaling = (ui.dial->value() * 20.0f) / 1023.0f; // 0 to 20
	// psudo code for histogram
	const unsigned short* raw_data = m_pParent->m_spProjImgRaw3D->GetBufferPointer();
	const size_t size_raw_data =
		m_pParent->m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize()[0] *
		m_pParent->m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize()[1] *
		m_pParent->m_spProjImgRaw3D->GetLargestPossibleRegion().GetSize()[2];

	const int n_bins = ui.spinBox_nBins->value();
	const size_t bin_size = MAX_USHORT / n_bins;

	unsigned int* histogram_raw = new unsigned int[(size_t)n_bins];
#pragma omp parallel for // initialize to 0
	for (int i = 0; i < n_bins; i++)
		histogram_raw[i] = 0;

	threaded_calculate_histogram_with_scaling(histogram_raw, n_bins, bin_size, scaling, raw_data, size_raw_data);

	QVector<double> vAxisX; //can be rows or columns
	QVector<double> vAxisY_raw;

	for (int i = 0; i< n_bins; i++)
	{
		vAxisX.push_back(bin_size * i);
		vAxisY_raw.push_back(histogram_raw[i]);
	}

	ui.customPlot->graph(0)->setData(vAxisX, vAxisY_raw);

	ui.customPlot->replot();

	ui.lcdNumber->display(scaling);

	delete[] histogram_raw;
}

// CF = mAs_ref / mAs (ref = CT)
void DlgHistogram::SLT_ReturnCF()
{
	float scaling = (ui.dial->value() * 20.0f) / 1023.0f; // 0 to 20

	m_pParent->ui.lineEdit_CurmAs->setText(QString("%1,20").arg((64 * 40 / 20) / scaling));

	m_pParent->ui.lineEdit_RefmAs->setText(QString("64,40"));
}