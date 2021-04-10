#include "mainwindow.h"
#include "glvisual.h"
#include <QtCore/qglobal.h>
#if QT_VERSION >= 0x050000
#include <QtWidgets/QApplication>
#else
#include <QtGui/QApplication>
#endif

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
    w.show();

    //[2]
    GLVisual v;
    v.rotate_x+=5;
    v.rotate_y-=5;
    v.show();
    v.setWindowTitle("3D Data Visualization");
    v.updateGL();

    return a.exec();
}
