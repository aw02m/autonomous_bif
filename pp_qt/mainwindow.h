#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include "dynamical_system.hpp"
#include "qcustomplot.h" // the header file of QCustomPlot. Don't forget to add it to your project, if you use an IDE, so it gets compiled.
#include <QMainWindow>
#include <QTimer>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow {
  Q_OBJECT

public:
  explicit MainWindow(QWidget *parent = 0);
  ~MainWindow();

  void ppsetup();
  void ppcalc(QCustomPlot *customPlot);

  void setupPlayground(QCustomPlot *customPlot);

private slots:
  void ppSlot();
  void mousePress(QMouseEvent *event);

private:
  Ui::MainWindow *ui;
  QTimer dataTimer;
  QCPCurve *trajectory;
  QCPGraph *poincare;
  dynamical_system ds;
  bool show_trajectory;
  bool show_poincare;

protected:
  void keyPressEvent(QKeyEvent *event);
  // void mousePressEvent(QMouseEvent *event);
};

#endif // MAINWINDOW_H
