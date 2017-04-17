QT += core
QT -= gui

CONFIG += c++11

TARGET = GraphGenerator
CONFIG += console
CONFIG -= app_bundle

TEMPLATE = app

SOURCES += main.cpp \
    dcel/anchor.cpp \
    dcel/barcode_template.cpp \
    dcel/dcel.cpp \
    dcel/mesh.cpp \
    point.cpp

HEADERS += \
    dcel/anchor.h \
    dcel/barcode_template.h \
    dcel/dcel.h \
    dcel/mesh.h \
    point.h
