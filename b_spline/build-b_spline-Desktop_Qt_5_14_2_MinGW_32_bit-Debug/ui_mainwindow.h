/********************************************************************************
** Form generated from reading UI file 'mainwindow.ui'
**
** Created by: Qt User Interface Compiler version 5.14.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtWidgets/QApplication>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QTextBrowser>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QWidget *layoutWidget;
    QGridLayout *gridLayout;
    QHBoxLayout *horizontalLayout_3;
    QComboBox *comboBox_lineType;
    QHBoxLayout *horizontalLayout;
    QLabel *label;
    QTextBrowser *textBrowser_X;
    QHBoxLayout *horizontalLayout_4;
    QPushButton *pushButton_clear;
    QHBoxLayout *horizontalLayout_2;
    QLabel *label_2;
    QTextBrowser *textBrowser_Y;
    QTextBrowser *textBrowser_msg;
    QGraphicsView *canvas;
    QMenuBar *menubar;
    QStatusBar *statusbar;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1112, 566);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        layoutWidget = new QWidget(centralwidget);
        layoutWidget->setObjectName(QString::fromUtf8("layoutWidget"));
        layoutWidget->setGeometry(QRect(850, 10, 254, 471));
        gridLayout = new QGridLayout(layoutWidget);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        gridLayout->setContentsMargins(0, 0, 0, 0);
        horizontalLayout_3 = new QHBoxLayout();
        horizontalLayout_3->setObjectName(QString::fromUtf8("horizontalLayout_3"));
        comboBox_lineType = new QComboBox(layoutWidget);
        comboBox_lineType->setObjectName(QString::fromUtf8("comboBox_lineType"));
        comboBox_lineType->setMinimumSize(QSize(111, 41));
        comboBox_lineType->setMaximumSize(QSize(111, 41));

        horizontalLayout_3->addWidget(comboBox_lineType);

        horizontalLayout = new QHBoxLayout();
        horizontalLayout->setObjectName(QString::fromUtf8("horizontalLayout"));
        label = new QLabel(layoutWidget);
        label->setObjectName(QString::fromUtf8("label"));
        label->setMinimumSize(QSize(54, 12));
        label->setMaximumSize(QSize(54, 12));

        horizontalLayout->addWidget(label);

        textBrowser_X = new QTextBrowser(layoutWidget);
        textBrowser_X->setObjectName(QString::fromUtf8("textBrowser_X"));
        textBrowser_X->setMinimumSize(QSize(71, 31));
        textBrowser_X->setMaximumSize(QSize(71, 31));

        horizontalLayout->addWidget(textBrowser_X);


        horizontalLayout_3->addLayout(horizontalLayout);


        gridLayout->addLayout(horizontalLayout_3, 0, 0, 1, 1);

        horizontalLayout_4 = new QHBoxLayout();
        horizontalLayout_4->setObjectName(QString::fromUtf8("horizontalLayout_4"));
        pushButton_clear = new QPushButton(layoutWidget);
        pushButton_clear->setObjectName(QString::fromUtf8("pushButton_clear"));
        pushButton_clear->setMinimumSize(QSize(111, 41));
        pushButton_clear->setMaximumSize(QSize(111, 41));

        horizontalLayout_4->addWidget(pushButton_clear);

        horizontalLayout_2 = new QHBoxLayout();
        horizontalLayout_2->setObjectName(QString::fromUtf8("horizontalLayout_2"));
        label_2 = new QLabel(layoutWidget);
        label_2->setObjectName(QString::fromUtf8("label_2"));
        label_2->setMinimumSize(QSize(54, 12));
        label_2->setMaximumSize(QSize(54, 12));

        horizontalLayout_2->addWidget(label_2);

        textBrowser_Y = new QTextBrowser(layoutWidget);
        textBrowser_Y->setObjectName(QString::fromUtf8("textBrowser_Y"));
        textBrowser_Y->setMinimumSize(QSize(71, 31));
        textBrowser_Y->setMaximumSize(QSize(71, 31));

        horizontalLayout_2->addWidget(textBrowser_Y);


        horizontalLayout_4->addLayout(horizontalLayout_2);


        gridLayout->addLayout(horizontalLayout_4, 1, 0, 1, 1);

        textBrowser_msg = new QTextBrowser(layoutWidget);
        textBrowser_msg->setObjectName(QString::fromUtf8("textBrowser_msg"));
        textBrowser_msg->setMinimumSize(QSize(231, 371));
        textBrowser_msg->setMaximumSize(QSize(231, 381));

        gridLayout->addWidget(textBrowser_msg, 2, 0, 1, 1);

        canvas = new QGraphicsView(centralwidget);
        canvas->setObjectName(QString::fromUtf8("canvas"));
        canvas->setGeometry(QRect(10, 10, 831, 511));
        MainWindow->setCentralWidget(centralwidget);
        menubar = new QMenuBar(MainWindow);
        menubar->setObjectName(QString::fromUtf8("menubar"));
        menubar->setGeometry(QRect(0, 0, 1112, 23));
        MainWindow->setMenuBar(menubar);
        statusbar = new QStatusBar(MainWindow);
        statusbar->setObjectName(QString::fromUtf8("statusbar"));
        MainWindow->setStatusBar(statusbar);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QCoreApplication::translate("MainWindow", "MainWindow", nullptr));
        label->setText(QCoreApplication::translate("MainWindow", "X\345\235\220\346\240\207", nullptr));
        pushButton_clear->setText(QCoreApplication::translate("MainWindow", "\346\270\205\351\231\244", nullptr));
        label_2->setText(QCoreApplication::translate("MainWindow", "Y\345\235\220\346\240\207", nullptr));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
