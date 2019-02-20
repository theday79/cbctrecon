// This is an open source non-commercial project. Dear PVS-Studio, please check it.
// PVS-Studio Static Code Analyzer for C, C++, C#, and Java: http://www.viva64.com

#include "qyklabel.h"
#include "AG17RGBAImage.h"
#include "YK16GrayImage.h"
#include <QPainter>

qyklabel::qyklabel(QWidget *parent) : QLabel(parent) {
  m_pYK16Image = nullptr;
  m_pRGBAImage = nullptr;
  // this->width();
  // m_Rt = this->rect();

  // m_Rt.setRect()
  m_bDrawPoints = true;

  m_iMouseWheelDelta = 0;

  m_bFocusIn = false;
}

qyklabel::~qyklabel() = default;

void qyklabel::mouseMoveEvent(QMouseEvent *ev) {
  if (m_pYK16Image == nullptr) {
    return;
  }

  this->x = ev->x();
  this->y = ev->y();
  emit Mouse_Move();
}

void qyklabel::mousePressEvent(QMouseEvent *ev) {
  if (m_pYK16Image == nullptr) {
    return;
  }
  this->x = ev->x();
  this->y = ev->y();

  if (ev->button() == Qt::LeftButton) {
    emit Mouse_Pressed_Left();
  }
  if (ev->button() == Qt::RightButton) {
    emit Mouse_Pressed_Right();
  }
}

void qyklabel::mouseReleaseEvent(QMouseEvent *ev) {
  if (m_pYK16Image == nullptr) {
    return;
  }

  // std::cout << "Mouse Released from qlabel" << std::endl;//it worked

  if (ev->button() == Qt::LeftButton) {
    emit Mouse_Released_Left();
  }
  if (ev->button() == Qt::RightButton) {
    emit Mouse_Released_Right();
  }
}

void qyklabel::wheelEvent(QWheelEvent *event) {
  if (m_pYK16Image == nullptr) {
    return;
  }

  m_iMouseWheelDelta = static_cast<int>(event->delta() / 100.0);
  emit Mouse_Wheel();
}

void qyklabel::enterEvent(QEvent * /*event*/) {
  if (m_pYK16Image == nullptr) {
    return;
  }

  m_bFocusIn = true;

  emit FocusIn();
}

void qyklabel::leaveEvent(QEvent * /*event*/) {
  if (m_pYK16Image == nullptr) {
    return;
  }

  m_bFocusIn = false;

  emit FocusOut();
}

// void focusInEvent ( QFocusEvent * ev );
// void focusOutEvent ( QFocusEvent * ev );
//
// void qyklabel::focusInEvent( QFocusEvent * )
//{
//  std::cout <<"focus in" << std::endl;
//  emit FocusIn();
//}
//
// void qyklabel::focusOutEvent( QFocusEvent * )
//{
//  std::cout <<"focus out" << std::endl;
//  emit FocusOut();
//}

QImage *compose_image_with_overlay(QRect *rect, QImage *src, QImage *overlay) {
  QPainter painter(src);
  painter.drawImage(*rect, *overlay, *rect, Qt::AutoColor);
  return src;
}

void qyklabel::paintEvent(QPaintEvent * /*unused*/) {
  QPainter painter(this);
  painter.setPen(QPen(Qt::black, 2));
  const auto target_rt = this->rect();
  painter.drawRect(target_rt);

  if (m_pYK16Image == nullptr) {
    return;
  }

  if (m_pYK16Image->m_iWidth < 1 || m_pYK16Image->m_iHeight < 1)
    return;

  // bool bHorLonger = false;

  auto phys_hor = static_cast<double>(m_pYK16Image->m_iWidth);
  auto phys_ver = static_cast<double>(m_pYK16Image->m_iHeight);

  // int labelNewFixedWidth = 0;
  // int labelNewFixedHeight = 0;

  if (fabs(m_pYK16Image->m_fSpacingX * m_pYK16Image->m_fSpacingY) > 0.0001) {
    phys_hor *= m_pYK16Image->m_fSpacingX;
    phys_ver *= m_pYK16Image->m_fSpacingY;
  }

  const auto VH_ratio =
      phys_ver / phys_hor; // if horizontal is longer than vertical

  if (phys_hor > phys_ver) {
    // bHorLonger = true;
    const auto newFixedHeight = qRound(this->width() * VH_ratio);
    this->setFixedHeight(newFixedHeight);
  } else {
    // bHorLonger = false;
    const auto newFixedWidth = qRound(this->height() / VH_ratio);
    this->setFixedWidth(newFixedWidth);
  }

  // Calculate ver / hor ratio.
  /*VH_ratio = m_pYK16Image->m_iHeight / (double) m_pYK16Image->m_iWidth;
  int newFixedHeight = qRound(this->width() * VH_ratio);
  this->setFixedHeight(newFixedHeight);*/
  if (m_pRGBAImage != nullptr && m_pYK16Image != nullptr) {
    painter.setCompositionMode(QPainter::CompositionMode_Source);
    painter.fillRect(target_rt, Qt::transparent);
    painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
    painter.drawImage(target_rt, m_pYK16Image->m_QImage.convertToFormat(
                                     QImage::Format_ARGB32));
    painter.setCompositionMode(QPainter::CompositionMode_SourceOver);
    painter.drawImage(target_rt, m_pRGBAImage->m_QImage.convertToFormat(
                                     QImage::Format_ARGB32));
    if (m_pRGBAImage->m_QImage.allGray())
      std::cerr << "RGBA was grayscale!?\n";
  } else if (m_pYK16Image != nullptr) {
    painter.drawImage(target_rt, m_pYK16Image->m_QImage); // it Works!YKTEMP
  }

  if (m_bDrawPoints) {
    painter.setPen(QPen(Qt::red, 2));
    for (auto &it : m_vPt) {
      painter.drawPoint(it.x(), it.y());
    }
    painter.setPen(QPen(Qt::green, 2));
    for (auto &it : m_vPt_green) {
      painter.drawPoint(it.x(), it.y());
    }
  }
  if (m_pYK16Image == nullptr) {
    return;
  }

  if (m_pYK16Image->m_bDrawROI) {
    painter.setPen(QPen(Qt::red, 2));
    QRect rtDraw;
    rtDraw.setTopLeft(Data2View(
        QPoint(m_pYK16Image->m_rtROI.left(), m_pYK16Image->m_rtROI.top()),
        this->width(), this->height(), m_pYK16Image->m_iWidth,
        m_pYK16Image->m_iHeight));
    rtDraw.setBottomRight(Data2View(
        QPoint(m_pYK16Image->m_rtROI.right(), m_pYK16Image->m_rtROI.bottom()),
        this->width(), this->height(), m_pYK16Image->m_iWidth,
        m_pYK16Image->m_iHeight));

    painter.drawRect(rtDraw);
  }
  // if x, y >0
  if (m_pYK16Image->m_bDrawProfileX) {
    auto crntViewPt =
        Data2View(m_pYK16Image->m_ptProfileProbe, this->width(), this->height(),
                  m_pYK16Image->m_iWidth, m_pYK16Image->m_iHeight);
    painter.setPen(QPen(Qt::red, 1, Qt::DotLine));
    painter.drawLine(0, crntViewPt.y(), this->width() - 1, crntViewPt.y());
  }
  if (m_pYK16Image->m_bDrawProfileY) {
    auto crntViewPt =
        Data2View(m_pYK16Image->m_ptProfileProbe, this->width(), this->height(),
                  m_pYK16Image->m_iWidth, m_pYK16Image->m_iHeight);
    painter.setPen(QPen(Qt::red, 1, Qt::DotLine));
    painter.drawLine(crntViewPt.x(), 0, crntViewPt.x(), this->height() - 1);
  }

  if (m_pYK16Image->m_bDrawFOVCircle) {
    painter.setPen(QPen(Qt::yellow, 1, Qt::SolidLine));
    const auto crnt_center_pt =
        Data2View(m_pYK16Image->m_ptFOVCenter, this->width(), this->height(),
                  m_pYK16Image->m_iWidth, m_pYK16Image->m_iHeight);
    const auto crnt_radius = static_cast<int>(
        m_pYK16Image->m_iFOVRadius /
        static_cast<double>(m_pYK16Image->m_iWidth) * this->width());
    painter.drawEllipse(crnt_center_pt, crnt_radius, crnt_radius);
  }

  if (m_pYK16Image->m_bDrawTableLine) {
    painter.setPen(QPen(Qt::yellow, 1, Qt::SolidLine));
    const auto crntTablePosY = static_cast<int>(
        m_pYK16Image->m_iTableTopPos /
        static_cast<double>(m_pYK16Image->m_iHeight) * this->height());
    // int crntRadius = (int)(m_pYK16Image->m_iFOVRadius /
    // (double)m_pYK16Image->m_iWidth * this->width());
    painter.drawLine(0, crntTablePosY, this->width() - 1, crntTablePosY);
    // painter.drawEllipse(crntCenterPt,crntRadius, crntRadius);
  }

  if (m_pYK16Image->m_bDrawCrosshair) // objects should be addressed one by one
  {
    painter.setPen(QPen(Qt::yellow, 1, Qt::SolidLine));
    // QPoint crosshair;

    // int dispCrossX = (int)(m_pYK16Image->m_ptCrosshair.x() /
    //                       (double)m_pYK16Image->m_iWidth * this->width());
    // int dispCrossY = (int)(m_pYK16Image->m_ptCrosshair.y() /
    //                       (double)m_pYK16Image->m_iHeight * this->height());

    auto ptDispCrosshair = GetViewPtFromDataPt(m_pYK16Image->m_ptCrosshair.x(),
                                               m_pYK16Image->m_ptCrosshair.y());
    // crosshair.setX(crossX);
    // crosshair.setY(crossY);

    // int crntRadius = (int)(m_pYK16Image->m_iFOVRadius /
    // (double)m_pYK16Image->m_iWidth * this->width());
    painter.drawLine(0, ptDispCrosshair.y(), this->width() - 1,
                     ptDispCrosshair.y());
    painter.drawLine(ptDispCrosshair.x(), 0, ptDispCrosshair.x(),
                     this->height() - 1);
  }
}

void qyklabel::SetBaseImage(YK16GrayImage *pYKImage) {
  if (pYKImage->m_pData != nullptr && !pYKImage->m_QImage.isNull()) // YKTEMP
    m_pYK16Image = pYKImage;
}

void qyklabel::SetOverlayImage(AG17RGBAImage *pRGBAImage) {
  if (pRGBAImage->m_pData.size() != 0 && !pRGBAImage->m_QImage.isNull()) {
    m_pRGBAImage = pRGBAImage;
  }
}

void qyklabel::ConvertAndCopyPoints(std::vector<QPoint> &vSrcPoint,
                                    const int iDataWidth,
                                    const int iDataHeight) {
  m_vPt.clear();

  const auto dspWidth = this->width();
  const auto dspHeight = this->height();

  for (auto &it : vSrcPoint) {
    auto tmp_dsp_pt =
        Data2View(it, dspWidth, dspHeight, iDataWidth, iDataHeight);
    m_vPt.push_back(tmp_dsp_pt);
  }
}

QPoint qyklabel::View2Data(QPoint view_pt, const int view_width,
                           const int view_height, const int data_width,
                           const int data_height) const {
  const auto f_zoom_x = view_width / static_cast<double>(data_width);
  const auto f_zoom_y = view_height / static_cast<double>(data_height);

  QPoint data_pt;
  data_pt.setX(qRound(view_pt.x() / f_zoom_x));
  data_pt.setY(qRound(view_pt.y() / f_zoom_y));

  return data_pt;
}

QPoint qyklabel::Data2View(QPoint data_pt, const int view_width,
                           const int view_height, const int data_width,
                           const int data_height) const {
  const auto f_zoom_x = view_width / static_cast<double>(data_width);
  const auto f_zoom_y = view_height / static_cast<double>(data_height);

  QPoint view_pt;
  view_pt.setX(qRound(data_pt.x() * f_zoom_x));
  view_pt.setY(qRound(data_pt.y() * f_zoom_y));

  return view_pt;
}

void qyklabel::SetDrawPointToggle(const bool bToggle) {
  if (bToggle) {
    m_bDrawPoints = true;
  } else {
    m_bDrawPoints = false;
  }

  update();
}

QPoint qyklabel::GetDataPtFromMousePos() const {
  if (m_pYK16Image == nullptr) {
    return {0, 0};
  }

  if (fabs(m_pYK16Image->m_fZoom - 1.0) < 0.001 && m_pYK16Image->m_iOffsetX == 0 &&
      m_pYK16Image->m_iOffsetY == 0)
    return View2Data(QPoint(this->x, this->y), width(), height(),
                     m_pYK16Image->m_iWidth, m_pYK16Image->m_iHeight);
  return View2DataExt(
      QPoint(this->x, this->y), width(), height(), m_pYK16Image->m_iWidth,
      m_pYK16Image->m_iHeight,
      QPoint(m_pYK16Image->m_iOffsetX, m_pYK16Image->m_iOffsetY),
      m_pYK16Image->m_fZoom);
}

QPoint qyklabel::GetDataPtFromViewPt(const int viewPtX, const int viewPtY) const {
  if (m_pYK16Image == nullptr) {
    return {0, 0};
  }

  if (fabs(m_pYK16Image->m_fZoom - 1.0) < 0.001 &&
      m_pYK16Image->m_iOffsetX == 0 &&
      m_pYK16Image->m_iOffsetY == 0)
    return View2Data(QPoint(viewPtX, viewPtY), width(), height(),
                     m_pYK16Image->m_iWidth, m_pYK16Image->m_iHeight);
  return View2DataExt(
      QPoint(viewPtX, viewPtY), width(), height(), m_pYK16Image->m_iWidth,
      m_pYK16Image->m_iHeight,
      QPoint(m_pYK16Image->m_iOffsetX, m_pYK16Image->m_iOffsetY),
      m_pYK16Image->m_fZoom);
}

QPoint qyklabel::GetViewPtFromDataPt(const int dataPtX, const int dataPtY) const {
  if (m_pYK16Image == nullptr) {
    return {0, 0};
  }

  if (fabs(m_pYK16Image->m_fZoom - 1.0) < 0.001 &&
      m_pYK16Image->m_iOffsetX == 0 &&
      m_pYK16Image->m_iOffsetY == 0)
    return Data2View(QPoint(dataPtX, dataPtY), width(), height(),
                     m_pYK16Image->m_iWidth, m_pYK16Image->m_iHeight);
  return Data2ViewExt(
      QPoint(dataPtX, dataPtY), width(), height(), m_pYK16Image->m_iWidth,
      m_pYK16Image->m_iHeight,
      QPoint(m_pYK16Image->m_iOffsetX, m_pYK16Image->m_iOffsetY),
      m_pYK16Image->m_fZoom);
}

QPoint qyklabel::View2DataExt(QPoint viewPt, const int viewWidth,
                              const int viewHeight, const int dataWidth,
                              const int dataHeight, QPoint ptDataOffset,
                              const double fUserZoom) const {
  // double fZoomX = viewWidth / (double)dataWidth * fUserZoom; // < 1
  // double fZoomY = viewHeight / (double)dataHeight * fUserZoom;

  //  dataPt.setX(qRound(viewPt.x() * fZoomX  + ptDataOffset.x()));
  //  dataPt.setY(qRound(viewPt.y() * fZoomY + ptDataOffset.y()));
  const auto newWidth = qRound(dataWidth / fUserZoom);
  const auto newHeight = qRound(dataHeight / fUserZoom);

  const auto dataCenterX = ptDataOffset.x() + qRound(dataWidth / 2.0);
  const auto dataCenterY = ptDataOffset.y() + qRound(dataHeight / 2.0);

  const auto dataLeftTopX =
      dataCenterX - qRound(newWidth / 2.0); // data position
  const auto dataLeftTopY =
      dataCenterY - qRound(newHeight / 2.0); // data position

  QPoint dataPt;
  dataPt.setX(qRound(viewPt.x() * newWidth / static_cast<double>(viewWidth) +
                     dataLeftTopX));
  dataPt.setY(qRound(viewPt.y() * newHeight / static_cast<double>(viewHeight) +
                     dataLeftTopY));

  return dataPt;
}

QPoint qyklabel::Data2ViewExt(QPoint dataPt, const int viewWidth,
                              const int viewHeight, const int dataWidth,
                              const int dataHeight, QPoint ptDataOffset,
                              const double fUserZoom) const {
  const auto newWidth = qRound(dataWidth / fUserZoom);
  const auto newHeight = qRound(dataHeight / fUserZoom);

  const auto dataCenterX = ptDataOffset.x() + qRound(dataWidth / 2.0);
  const auto dataCenterY = ptDataOffset.y() + qRound(dataHeight / 2.0);

  const auto dataLeftTopX =
      dataCenterX - qRound(newWidth / 2.0); // data position
  const auto dataLeftTopY =
      dataCenterY - qRound(newHeight / 2.0); // data position

  QPoint viewPt;
  viewPt.setX(qRound((dataPt.x() - dataLeftTopX) * viewWidth /
                     static_cast<double>(newWidth)));
  viewPt.setY(qRound((dataPt.y() - dataLeftTopY) * viewHeight /
                     static_cast<double>(newHeight)));

  return viewPt;
}

// void qyklabel::keyPressEvent( QKeyEvent* ev )
//{
//  if (!isFocusIn())
//	return;
//
//  int enMovingKey = -1;
//
//  if (ev->key() == Qt::Key_Left)
//	enMovingKey = 0;
//  else if (ev->key() == Qt::Key_Right)
//	enMovingKey = 1;
//  else if (ev->key() == Qt::Key_Up)
//	enMovingKey = 2;
//  else if (ev->key() == Qt::Key_Down)
//	enMovingKey = 3;
//  else if (ev->key() == Qt::Key_PageUp)
//	enMovingKey = 4;
//  else if (ev->key() == Qt::Key_PageUp)
//	enMovingKey = 5;
//
//  std::cout << enMovingKey << std::endl;
//
//  emit ArrowPressed(enMovingKey);
//}
