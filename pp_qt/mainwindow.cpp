#include "mainwindow.h"
#include "ds_derivatives.hpp"
#include "runge_kutta.cpp"
#include "ui_mainwindow.h"
#include <QDebug>
#if QT_VERSION < QT_VERSION_CHECK(6, 0, 0)
#include <QDesktopWidget>
#endif
#include <QMessageBox>
#include <QMetaEnum>
#include <QScreen>
#include <iostream>

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent), ui(new Ui::MainWindow),
      ds(QApplication::arguments().at(1).toStdString()) {

  // ds = dynamical_system(QApplication::arguments().at(1).toStdString());

  ui->setupUi(this);
  setGeometry(400, 250, 768, 768);

  ppsetup();
}

void MainWindow::ppsetup() {
  ppcalc(ui->customPlot);
  setWindowTitle("Phase Portrait : pp");
  statusBar()->clearMessage();
  ui->customPlot->replot();
}

void MainWindow::ppcalc(QCustomPlot *customPlot) {
  trajectory = new QCPCurve(customPlot->xAxis, customPlot->yAxis);
  poincare = new QCPGraph(customPlot->xAxis, customPlot->yAxis);

  // QVector<QCPCurveData> sol = integrate(f, 0, ds.x0, ds.tick, ds);
  // trajectory->data()->set(sol, true);

  trajectory->setPen(QPen(Qt::blue));
  poincare->setLineStyle(QCPGraph::lsNone);
  poincare->setScatterStyle(QCPScatterStyle::ssDisc);
  poincare->setPen(QPen(QBrush(Qt::red), 0.1));

  customPlot->xAxis->setLabel("x");
  customPlot->yAxis->setLabel("y");
  customPlot->xAxis->setRange(ds.xrange[0], ds.xrange[1]);
  customPlot->yAxis->setRange(ds.yrange[0], ds.yrange[1]);

  customPlot->setMouseTracking(true);

  connect(&dataTimer, SIGNAL(timeout()), this, SLOT(ppSlot()));
  connect(ui->customPlot, SIGNAL(mousePress(QMouseEvent *)), this,
          SLOT(mousePress(QMouseEvent *)));
  dataTimer.start(0); // Interval 0 means to refresh as fast as possible
}

void MainWindow::ppSlot() {
  double secs =
      QCPAxisTickerDateTime::dateTimeToKey(QDateTime::currentDateTime());

  ds.integrate(0, ds.last_state, ds.tick);
  trajectory->data()->set(ds.QCPCsol, true);
  // ds.tau += ds.QCPCsol.last().t;
  unsigned int sol_size = ds.QCPCsol.size();
  if (sol_size >= ds.max_plot) {
    ds.QCPCsol.remove(0, sol_size - ds.max_plot);
  }

  if (ds.is_hit_section()) {
    ds.QCPGpoincare << QCPGraphData(ds.last_state(ds.axis[0]),
                                    ds.last_state(ds.axis[1]));
    poincare->data()->set(ds.QCPGpoincare, true);
    std::cout << ds.tau << std::endl;
    ds.tau = 0;
  }
  unsigned int poincare_size = ds.QCPGpoincare.size();
  if (poincare_size >= ds.max_poincare_plot) {
    ds.QCPGpoincare.remove(0, poincare_size - ds.max_poincare_plot);
  }

  ui->customPlot->replot();

  // calculate frames per second:
  double key = secs;
  static double lastFpsKey;
  static int frameCount;
  ++frameCount;
  if (key - lastFpsKey > 2) // average fps over 2 seconds
  {
    ui->statusBar->showMessage(
        QString("%1 FPS, Total Data points: %2")
            .arg(frameCount / (key - lastFpsKey), 0, 'f', 0)
            .arg(trajectory->data()->size()),
        0);
    lastFpsKey = key;
    frameCount = 0;
  }
}

void MainWindow::mousePress(QMouseEvent *event) {
  if (event->button() == Qt::LeftButton) {
    ds.last_state(ds.axis[0]) =
        ui->customPlot->xAxis->pixelToCoord(event->pos().x());
    ds.last_state(ds.axis[1]) =
        ui->customPlot->yAxis->pixelToCoord(event->pos().y());
    ds.QCPCsol.clear();
    ds.QCPGpoincare.clear();
  }
}

void MainWindow::keyPressEvent(QKeyEvent *event) {
  static unsigned int param_index = 0;
  static Eigen::IOFormat Comma(8, 0, ", ", "\n", "[", "]");
  switch (event->key()) {
  case Qt::Key_Left:
    if (param_index == 0) {
      param_index = ds.params.size() - 1;
    } else {
      param_index--;
    }
    std::cout << "parameter index changed : " << param_index << std::endl;
    break;
  case Qt::Key_Right:
    if (param_index == ds.params.size() - 1) {
      param_index = 0;
    } else {
      param_index++;
    }
    std::cout << "parameter index changed : " << param_index << std::endl;
    break;
  case Qt::Key_Up:
    ds.params(param_index) += ds.dparams[param_index];
    std::cout << param_index << ":" << ds.params.transpose().format(Comma)
              << std::endl;
    break;
  case Qt::Key_Down:
    ds.params(param_index) -= ds.dparams[param_index];
    std::cout << param_index << ":" << ds.params.transpose().format(Comma)
              << std::endl;
    break;
  case Qt::Key_Space:
    std::cout << "state : " << ds.x0.transpose().format(Comma)
              << std::endl;
    std::cout << "param : " << ds.params.transpose().format(Comma) << std::endl;
    break;
  case Qt::Key_Q:
    QApplication::exit();
    break;
  }
}

void MainWindow::setupPlayground(QCustomPlot *customPlot) {
  Q_UNUSED(customPlot)
}

MainWindow::~MainWindow() { delete ui; }