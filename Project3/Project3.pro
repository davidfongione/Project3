TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    planet.cpp \
    protoplanet.cpp \
    System.cpp

HEADERS += \
    planet.h \
    protoplanet.h \
    System.h
