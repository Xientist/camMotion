QT -= gui

CONFIG += c++11 console
CONFIG -= app_bundle

# The following define makes your compiler emit warnings if you use
# any Qt feature that has been marked deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if it uses deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0

SOURCES += \
    getCorrespondences/GetCorrespondences.cpp \
    getCorrespondences/integralimage.cpp \
    getCorrespondences/isprocessor2.cpp \
        main.cpp

INCLUDEPATH += /usr/local/include/opencv4/
INCLUDEPATH += /usr/local/include/eigen3/

LIBS += /usr/local/lib/libopencv_core.so
LIBS += /usr/local/lib/libopencv_highgui.so
LIBS += /usr/local/lib/libopencv_imgcodecs.so
LIBS += /usr/local/lib/libopencv_imgproc.so
LIBS += /usr/local/lib/libopencv_features2d.so
LIBS += /usr/local/lib/libopencv_calib3d.so

# Default rules for deployment.
qnx: target.path = /tmp/$${TARGET}/bin
else: unix:!android: target.path = /opt/$${TARGET}/bin
!isEmpty(target.path): INSTALLS += target

HEADERS += \
    basicGeometry.h \
    epipolarGeometry.h \
    getCorrespondences/cdcorner.h \
    getCorrespondences/gsmean.h \
    getCorrespondences/integralimage.h \
    getCorrespondences/isprocessor2.h \
    diff.h

QMAKE_CXXFLAGS += -mpopcnt
QMAKE_CXXFLAGS_RELEASE += -Ofast
