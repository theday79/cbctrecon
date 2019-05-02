// This is an open source non-commercial project. Dear PVS-Studio, please check
// it. PVS-Studio Static Code Analyzer for C, C++, C#, and Java:
// http://www.viva64.com

#if USE_OPENMP
#include <omp.h>
#endif

#include "DlgHistogram.h"
#include "cbctrecon_mainwidget.h"
#include "qcustomplot.h"

DlgHistogram::DlgHistogram() {
  /* Sets up the GUI */
  ui.setupUi(this);
}

DlgHistogram::DlgHistogram(CbctReconWidget *parent) : QDialog(parent) {
  /* Sets up the GUI */
  ui.setupUi(this);
  m_pParent = parent;
}

struct HistogramType {
  HistogramType(const size_t n_bins) {
    hist_data = std::valarray<unsigned int>(n_bins);
    bin_size =
        static_cast<float>(std::numeric_limits<unsigned short>::max()) / n_bins;
  }
  std::valarray<unsigned int> hist_data;
  float bin_size = 0.0f;
};

void threaded_calculate_histogram(
    HistogramType &histogram_out,
    const std::valarray<unsigned short> &raw_data) {
  const auto bin_scaling = 1.f / histogram_out.bin_size;

  for (auto &val : raw_data) {
    ++histogram_out.hist_data[std::lroundf(val * bin_scaling)];
  }
}

void DlgHistogram::SLT_DrawGraph() const {
  auto &raw_img = m_pParent->m_cbctrecon->m_spProjImgRaw3D;
  auto &ct_img = m_pParent->m_cbctrecon->m_spProjImgCT3D;
  if (!raw_img || !ct_img) {
    return;
  }

  const auto &raw_size = raw_img->GetLargestPossibleRegion().GetSize();
  const auto size_raw_data = raw_size[0] * raw_size[1] * raw_size[2];
  const std::valarray<unsigned short> raw_data(raw_img->GetBufferPointer(),
                                               size_raw_data);

  const auto &ct_size = raw_img->GetLargestPossibleRegion().GetSize();
  const auto size_ct_data = ct_size[0] * ct_size[1] * ct_size[2];
  const std::valarray<unsigned short> ct_data(ct_img->GetBufferPointer(),
                                              size_ct_data);

  const auto n_bins = static_cast<size_t>(ui.spinBox_nBins->value());

  HistogramType histogram_raw(n_bins);
  HistogramType histogram_ct(n_bins);
  // same size should be enough, because raw
  // projections are normalised at read time.

#pragma omp parallel sections
  {
#pragma omp section
    { threaded_calculate_histogram(histogram_raw, raw_data); }
#pragma omp section
    { threaded_calculate_histogram(histogram_ct, ct_data); }
  }

  QVector<double> vAxisX(n_bins, histogram_raw.bin_size);
  std::transform(vAxisX.begin(), vAxisX.end(), vAxisX.begin(),
                 [i = 0](const double val) mutable { return val * i++; });

  QVector<double> vAxisY_raw(n_bins);
  std::copy(std::begin(histogram_raw.hist_data),
            std::end(histogram_raw.hist_data), std::begin(vAxisY_raw));

  ui.customPlot->clearGraphs();

  ui.customPlot->addGraph();
  ui.customPlot->graph(0)->setData(vAxisX, vAxisY_raw);
  ui.customPlot->graph(0)->setPen(QPen(Qt::red));
  ui.customPlot->graph(0)->setName("Raw proj. histogram");

  QVector<double> vAxisY_ct;
  std::copy(std::begin(histogram_ct.hist_data),
            std::end(histogram_ct.hist_data), std::back_inserter(vAxisY_ct));
  ui.customPlot->addGraph();
  ui.customPlot->graph(1)->setData(vAxisX, vAxisY_ct);
  ui.customPlot->graph(1)->setPen(QPen(Qt::blue));
  ui.customPlot->graph(1)->setName("CT fwd. proj. histogram");

  // Set limits of graph
  const auto minX = 0.0;
  ui.lineEditXMin->setText(QString("%1").arg(minX));

  auto maxX = 0.0;
  size_t index = 0;
  for (auto &val : histogram_raw.hist_data) {
    if (val > 0) {
      maxX = index * histogram_raw.bin_size;
    }
    ++index;
  }
  ui.lineEditXMax->setText(QString("%1").arg(maxX));

  const auto minY = 0.0;
  ui.lineEditYMin->setText(QString("%1").arg(minY));

  const auto maxY = histogram_raw.hist_data.max();
  ui.lineEditYMax->setText(QString("%1").arg(maxY));

  ui.customPlot->xAxis->setRange(minX, maxX);
  ui.customPlot->yAxis->setRange(minY, maxY);

  // Set labels
  ui.customPlot->xAxis->setLabel("value");
  ui.customPlot->yAxis->setLabel("N");
  ui.customPlot->setWindowTitle("Intensity histograms");

  ui.customPlot->legend->setVisible(true);
  auto legendFont = font();   // start out with MainWindow's font..
  legendFont.setPointSize(9); // and make a bit smaller for legend
  ui.customPlot->legend->setFont(legendFont);
  ui.customPlot->legend->setBrush(QBrush(QColor(255, 255, 255, 200)));

  ui.customPlot->replot();
}

void DlgHistogram::SLT_ReDrawGraph_limits() const {

  const auto tmpXMin = ui.lineEditXMin->text().toDouble();
  const auto tmpXMax = ui.lineEditXMax->text().toDouble();
  const auto tmpYMin = ui.lineEditYMin->text().toDouble();
  const auto tmpYMax = ui.lineEditYMax->text().toDouble();

  ui.customPlot->xAxis->setRange(tmpXMin, tmpXMax);
  ui.customPlot->yAxis->setRange(tmpYMin, tmpYMax);

  ui.customPlot->replot();
}

void threaded_calculate_histogram_with_scaling(
    HistogramType &histogram_out, const float scaling,
    const std::valarray<unsigned short> &raw_data) {

  const auto bin_scaling = scaling / histogram_out.bin_size;
  for (auto &val : raw_data) {
    ++histogram_out.hist_data[std::lroundf(val * bin_scaling)];
  }
}

void DlgHistogram::SLT_ReDrawGraph_dial() const {
  auto &raw_img = m_pParent->m_cbctrecon->m_spProjImgRaw3D;
  if (!raw_img) {
    return;
  }

  const auto scaling =
      (ui.dial->value() * 20.0f) / ui.dial->maximum(); // 0 to 20

  const auto &raw_size = raw_img->GetLargestPossibleRegion().GetSize();
  const auto size_raw_data = raw_size[0] * raw_size[1] * raw_size[2];
  const std::valarray<unsigned short> raw_data(raw_img->GetBufferPointer(),
                                               size_raw_data);

  const auto n_bins = ui.spinBox_nBins->value();

  HistogramType histogram_raw(n_bins);

  threaded_calculate_histogram_with_scaling(histogram_raw, scaling, raw_data);

  QVector<double> vAxisX(n_bins, histogram_raw.bin_size);
  std::transform(vAxisX.begin(), vAxisX.end(), vAxisX.begin(),
                 [i = 0](const double val) mutable { return val * i++; });

  QVector<double> vAxisY_raw(histogram_raw.hist_data.size());
  std::copy(std::begin(histogram_raw.hist_data),
            std::end(histogram_raw.hist_data), vAxisY_raw.begin());

  ui.customPlot->graph(0)->setData(vAxisX, vAxisY_raw);

  ui.customPlot->replot();

  ui.lcdNumber->display(scaling);
}

// CF = mAs_ref / mAs (ref = CT)
void DlgHistogram::SLT_ReturnCF() const {
  const auto scaling = ui.dial->value() * 20.0f / 1023.0f; // 0 to 20

  m_pParent->ui.lineEdit_CurmAs->setText(
      QString("%1,20").arg((64.0 * 40.0 / 20.0) / scaling));

  m_pParent->ui.lineEdit_RefmAs->setText(QString("64,40"));
}
