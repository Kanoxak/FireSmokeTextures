QT       += core gui opengl

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

CONFIG += c++11

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

mac {
  PKG_CONFIG = /usr/local/bin/pkg-config
}

INCLUDEPATH += -I/opt/local/include
INCLUDEPATH += -I/opt/local/include/opencv
INCLUDEPATH += -I/opt/local/include/opencv2

INCLUDEPATH += /usr/local/Cellar/opencv/4.1.2/include
INCLUDEPATH += /usr/local/Cellar/opencv/4.1.2/include/opencv4
INCLUDEPATH += /usr/local/Cellar/opencv/4.1.2/include/opencv4/opencv2

LIBS += -L/usr/local/lib -lopencv_core -lopencv_imgcodecs -lopencv_highgui
LIBS += -L/opt/local/lib
        -lopencv_shape
        -lopencv_stitching
        -lopencv_objdetect
        -lopencv_superres
        -lopencv_videostab
        -lopencv_calib3d
        -lopencv_features2d
        -lopencv_highgui
        -lopencv_videoio
        -lopencv_videoio.3.2.0
        -lopencv_videoio.3.2
        -lopencv_imgcodecs
        -lopencv_video
        -lopencv_video.3.2.0
        -lopencv_video.3.2
        -lopencv_videostab
        -lopencv_videostab.3.2.0
        -lopencv_videostab.3.2
        -lopencv_photo
        -lopencv_ml
        -lopencv_imgproc
        -lopencv_flann
        -lopencv_core
LIBS += -L/usr/lib
        -lopencv_shape
        -lopencv_stitching
        -lopencv_objdetect
        -lopencv_superres
        -lopencv_videostab
        -lopencv_calib3d
        -lopencv_features2d
        -lopencv_highgui
        -lopencv_videoio
        -lopencv_imgcodecs
        -lopencv_video
        -lopencv_photo
        -lopencv_ml
        -lopencv_imgproc
        -lopencv_flann
        -lopencv_core

SOURCES += \
    glvisual.cpp \
    main.cpp \
    mainwindow.cpp \
    slic.cpp

HEADERS += \
    glvisual.h \
    mainwindow.h \
    slic.h

FORMS += \
    mainwindow.ui

# User-defined Configuration
CONFIG += link_pkgconfig
CONFIG += sdk_no_version_check
PKGCONFIG += opencv

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target
