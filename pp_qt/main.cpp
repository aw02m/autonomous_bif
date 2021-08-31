#include "mainwindow.h"
#include <QApplication>
#include <iostream>

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cerr << "usage: ./main [json_file]" << std::endl;
    std::exit(1);
  }

#if QT_VERSION < QT_VERSION_CHECK(5, 0, 0)
  QApplication::setGraphicsSystem("raster");
#endif

  QApplication a(argc, argv);
  MainWindow w;
  w.show();

  return a.exec();
}