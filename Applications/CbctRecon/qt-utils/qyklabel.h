#ifndef QYKLABEL_H
#define QYKLABEL_H

#include <QLabel>
#include <QMouseEvent>
#include <QRect>
#include <vector>

class YK16GrayImage;
class AG17RGBAImage;

class qyklabel : public QLabel {
  Q_OBJECT

public:
  YK16GrayImage *m_pYK16Image;
  AG17RGBAImage *m_pRGBAImage;
  QRect m_Rt;
  std::vector<QPoint> m_vPt;
  std::vector<QPoint> m_vPt_green;
  bool m_bDrawPoints;
  int m_iMouseWheelDelta;

  bool m_bFocusIn;

  explicit qyklabel(QWidget *parent);
  ~qyklabel() override;
  qyklabel(const qyklabel &) = delete;
  void operator=(const qyklabel &) = delete;
  qyklabel(qyklabel &&) = delete;
  void operator=(qyklabel &&) = delete;

  bool isFocusIn() const { return m_bFocusIn; }

  // virtual function reimplementation
  void
  mousePressEvent(QMouseEvent *ev) override; // ev->buttons() == Qt::LeftButton
  void mouseMoveEvent(QMouseEvent *ev) override;
  void mouseReleaseEvent(QMouseEvent *ev) override;
  void wheelEvent(
      QWheelEvent *event) override; // this->setText("Delta Value:
                                    // "+QString::number(event->delta()));

  // void keyPressEvent ( QKeyEvent *ev);
  // void focusInEvent ( QFocusEvent * ev );
  // void focusOutEvent ( QFocusEvent * ev );

  void enterEvent(QEvent * /*event*/) override;
  void leaveEvent(QEvent * /*event*/) override;
  int x{}, y{};

  void SetBaseImage(YK16GrayImage *pYKImage);
  void SetOverlayImage(AG17RGBAImage *pRGBAImage);
  // void ConvertAndCopyPoints(std::vector<QPoint>& vSrcPoint);
  void ConvertAndCopyPoints(std::vector<QPoint> &vSrcPoint, int iDataWidth,
                            int iDataHeight);

  QPoint View2Data(QPoint view_pt, int view_width, int view_height,
                   int data_width, int data_height) const;
  QPoint View2DataExt(QPoint viewPt, int viewWidth, int viewHeight,
                      int dataWidth, int dataHeight, QPoint ptDataOffset,
                      double fUserZoom) const;
  QPoint Data2View(QPoint data_pt, int view_width, int view_height,
                   int data_width, int data_height) const;
  QPoint Data2ViewExt(QPoint dataPt, int viewWidth, int viewHeight,
                      int dataWidth, int dataHeight, QPoint ptDataOffset,
                      double fUserZoom) const;

  QPoint GetDataPtFromMousePos() const; // Return data position of the mouse
                                        // position.m_pYK16 image is mandatory
  QPoint GetDataPtFromViewPt(int viewPtX, int viewPtY) const;

  QPoint GetViewPtFromDataPt(int dataPtX, int dataPtY) const;

protected:
  void paintEvent(QPaintEvent * /*unused*/) override;

signals:
  void Mouse_Pressed_Left();
  void Mouse_Pressed_Right();
  void Mouse_Move();

  void Mouse_Released_Left();
  void Mouse_Released_Right();
  void Mouse_Wheel();

  void FocusIn();
  void FocusOut();

  // void ArrowPressed(int arrowDirection);// 0: Left, 1: Upward, 2: Right, 3:
  // Downward  void OutFromWindow();

public slots:
  void SetDrawPointToggle(bool bToggle);
};

#endif // QYKLABEL_H
