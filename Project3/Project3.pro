TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    planet.cpp \
    protoplanet.cpp \
    system.cpp

HEADERS += \
    planet.h \
    protoplanet.h \
    system.h
