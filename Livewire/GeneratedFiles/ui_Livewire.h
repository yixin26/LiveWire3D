/********************************************************************************
** Form generated from reading UI file 'Livewire.ui'
**
** Created by: Qt User Interface Compiler version 5.12.4
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_LIVEWIRE_H
#define UI_LIVEWIRE_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QWidget>
#include "MyGLCanvas.h"

QT_BEGIN_NAMESPACE

class Ui_LivewireClass
{
public:
    QAction *OpenMesh;
    QAction *SaveCurves;
    QWidget *centralWidget;
    MyGLCanvas *openGLWidget;
    QWidget *verticalLayoutWidget;
    QVBoxLayout *verticalLayout;
    QCheckBox *Geodesy_checkBox;
    QCheckBox *AnisGeoMin_checkBox;
    QCheckBox *AnisGeoMax_checkBox;
    QMenuBar *menuBar;
    QMenu *menuOpen;
    QMenu *menuSave;

    void setupUi(QMainWindow *LivewireClass)
    {
        if (LivewireClass->objectName().isEmpty())
            LivewireClass->setObjectName(QString::fromUtf8("LivewireClass"));
        LivewireClass->resize(1600, 1000);
        OpenMesh = new QAction(LivewireClass);
        OpenMesh->setObjectName(QString::fromUtf8("OpenMesh"));
        SaveCurves = new QAction(LivewireClass);
        SaveCurves->setObjectName(QString::fromUtf8("SaveCurves"));
        centralWidget = new QWidget(LivewireClass);
        centralWidget->setObjectName(QString::fromUtf8("centralWidget"));
        openGLWidget = new MyGLCanvas(centralWidget);
        openGLWidget->setObjectName(QString::fromUtf8("openGLWidget"));
        openGLWidget->setGeometry(QRect(0, 0, 1600, 930));
        verticalLayoutWidget = new QWidget(centralWidget);
        verticalLayoutWidget->setObjectName(QString::fromUtf8("verticalLayoutWidget"));
        verticalLayoutWidget->setGeometry(QRect(10, 10, 251, 80));
        verticalLayout = new QVBoxLayout(verticalLayoutWidget);
        verticalLayout->setSpacing(6);
        verticalLayout->setContentsMargins(11, 11, 11, 11);
        verticalLayout->setObjectName(QString::fromUtf8("verticalLayout"));
        verticalLayout->setContentsMargins(0, 0, 0, 0);
        Geodesy_checkBox = new QCheckBox(verticalLayoutWidget);
        Geodesy_checkBox->setObjectName(QString::fromUtf8("Geodesy_checkBox"));
        QSizePolicy sizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
        sizePolicy.setHorizontalStretch(0);
        sizePolicy.setVerticalStretch(0);
        sizePolicy.setHeightForWidth(Geodesy_checkBox->sizePolicy().hasHeightForWidth());
        Geodesy_checkBox->setSizePolicy(sizePolicy);

        verticalLayout->addWidget(Geodesy_checkBox);

        AnisGeoMin_checkBox = new QCheckBox(verticalLayoutWidget);
        AnisGeoMin_checkBox->setObjectName(QString::fromUtf8("AnisGeoMin_checkBox"));
        sizePolicy.setHeightForWidth(AnisGeoMin_checkBox->sizePolicy().hasHeightForWidth());
        AnisGeoMin_checkBox->setSizePolicy(sizePolicy);
        AnisGeoMin_checkBox->setChecked(true);

        verticalLayout->addWidget(AnisGeoMin_checkBox);

        AnisGeoMax_checkBox = new QCheckBox(verticalLayoutWidget);
        AnisGeoMax_checkBox->setObjectName(QString::fromUtf8("AnisGeoMax_checkBox"));
        sizePolicy.setHeightForWidth(AnisGeoMax_checkBox->sizePolicy().hasHeightForWidth());
        AnisGeoMax_checkBox->setSizePolicy(sizePolicy);
        AnisGeoMax_checkBox->setChecked(true);

        verticalLayout->addWidget(AnisGeoMax_checkBox);

        LivewireClass->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(LivewireClass);
        menuBar->setObjectName(QString::fromUtf8("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 1600, 18));
        menuOpen = new QMenu(menuBar);
        menuOpen->setObjectName(QString::fromUtf8("menuOpen"));
        menuSave = new QMenu(menuBar);
        menuSave->setObjectName(QString::fromUtf8("menuSave"));
        LivewireClass->setMenuBar(menuBar);

        menuBar->addAction(menuOpen->menuAction());
        menuBar->addAction(menuSave->menuAction());
        menuOpen->addAction(OpenMesh);
        menuSave->addAction(SaveCurves);

        retranslateUi(LivewireClass);
        QObject::connect(Geodesy_checkBox, SIGNAL(clicked()), LivewireClass, SLOT(Geodesy()));
        QObject::connect(AnisGeoMin_checkBox, SIGNAL(clicked()), LivewireClass, SLOT(AnisGeodesyMin()));
        QObject::connect(AnisGeoMax_checkBox, SIGNAL(clicked()), LivewireClass, SLOT(AnisGeodesyMax()));
        QObject::connect(OpenMesh, SIGNAL(triggered()), LivewireClass, SLOT(OpenMeshFile()));
        QObject::connect(SaveCurves, SIGNAL(triggered()), LivewireClass, SLOT(SaveCurveFile()));

        QMetaObject::connectSlotsByName(LivewireClass);
    } // setupUi

    void retranslateUi(QMainWindow *LivewireClass)
    {
        LivewireClass->setWindowTitle(QApplication::translate("LivewireClass", "Livewire", nullptr));
        OpenMesh->setText(QApplication::translate("LivewireClass", "Mesh", nullptr));
        SaveCurves->setText(QApplication::translate("LivewireClass", "Curves", nullptr));
        Geodesy_checkBox->setText(QApplication::translate("LivewireClass", "Geodesy", nullptr));
        AnisGeoMin_checkBox->setText(QApplication::translate("LivewireClass", "Anisotropic Geodesy Min", nullptr));
        AnisGeoMax_checkBox->setText(QApplication::translate("LivewireClass", "Anisotropic Geodesy Max", nullptr));
        menuOpen->setTitle(QApplication::translate("LivewireClass", "Open", nullptr));
        menuSave->setTitle(QApplication::translate("LivewireClass", "Save", nullptr));
    } // retranslateUi

};

namespace Ui {
    class LivewireClass: public Ui_LivewireClass {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_LIVEWIRE_H
