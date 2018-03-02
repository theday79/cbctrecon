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
  bool m_bDrawPoints;
  int m_iMouseWheelDelta;

  bool m_bFocusIn;

public:
  qyklabel(QWidget *parent);
  ~qyklabel() override;

  bool isFocusIn() { return m_bFocusIn; }

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

  QPoint View2Data(QPoint viewPt, int viewWidth, int viewHeight, int dataWidth,
                   int dataHeight);
  QPoint View2DataExt(QPoint viewPt, int viewWidth, int viewHeight,
                      int dataWidth, int dataHeight, QPoint ptDataOffset,
                      double fUserZoom);
  QPoint Data2View(QPoint dataPt, int viewWidth, int viewHeight, int dataWidth,
                   int dataHeight);
  QPoint Data2ViewExt(QPoint dataPt, int viewWidth, int viewHeight,
                      int dataWidth, int dataHeight, QPoint ptDataOffset,
                      double fUserZoom);

  QPoint GetDataPtFromMousePos(); // Return data position of the mouse
                                  // position.m_pYK16 image is mandatory
  QPoint GetDataPtFromViewPt(int viewPtX, int viewPtY);

  QPoint GetViewPtFromDataPt(int dataPtX, int dataPtY);

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

private:
};

#endif // QYKLABEL_H
