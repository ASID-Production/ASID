# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'search_dialogLMfxvM.ui'
##
## Created by: Qt User Interface Compiler version 6.7.2
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QBrush, QColor, QConicalGradient, QCursor,
    QFont, QFontDatabase, QGradient, QIcon,
    QImage, QKeySequence, QLinearGradient, QPainter,
    QPalette, QPixmap, QRadialGradient, QTransform)
from PySide6.QtWidgets import (QApplication, QCheckBox, QComboBox, QDialog,
    QFormLayout, QFrame, QHBoxLayout, QLineEdit,
    QPushButton, QSizePolicy, QVBoxLayout, QWidget, QLabel)

class Ui_Dialog(object):
    def setupUi(self, Dialog):
        if not Dialog.objectName():
            Dialog.setObjectName(u"Dialog")
        Dialog.resize(400, 260)
        self.horizontalLayout = QHBoxLayout(Dialog)
        self.horizontalLayout.setObjectName(u"horizontalLayout")
        self.frame_2 = QFrame(Dialog)
        self.frame_2.setObjectName(u"frame_2")
        self.frame_2.setFrameShape(QFrame.Shape.StyledPanel)
        self.frame_2.setFrameShadow(QFrame.Shadow.Raised)
        self.verticalLayout = QVBoxLayout(self.frame_2)
        self.verticalLayout.setSpacing(0)
        self.verticalLayout.setObjectName(u"verticalLayout")
        self.verticalLayout.setContentsMargins(0, 0, 0, 0)
        self.comboBox = QComboBox(self.frame_2)
        self.comboBox.addItem("")
        self.comboBox.addItem("")
        self.comboBox.setObjectName(u"comboBox")

        self.verticalLayout.addWidget(self.comboBox)

        self.frame = QFrame(self.frame_2)
        self.frame.setObjectName(u"frame")
        self.frame.setFrameShape(QFrame.Shape.StyledPanel)
        self.frame.setFrameShadow(QFrame.Shadow.Raised)
        self.formLayout = QFormLayout(self.frame)
        self.formLayout.setObjectName(u"formLayout")
        self.formLayout.setContentsMargins(0, 0, 0, 0)
        self.checkBox = QCheckBox(self.frame)
        self.checkBox.setObjectName(u"checkBox")
        self.checkBox.setChecked(True)
        self.checkBox.setAutoExclusive(True)

        self.formLayout.setWidget(0, QFormLayout.LabelRole, self.checkBox)

        self.lineEdit = QLineEdit(self.frame)
        self.lineEdit.setObjectName(u"lineEdit")

        self.formLayout.setWidget(0, QFormLayout.FieldRole, self.lineEdit)

        self.checkBox_8 = QCheckBox(self.frame)
        self.checkBox_8.setObjectName(u"checkBox_8")
        self.checkBox_8.setAutoExclusive(True)

        self.formLayout.setWidget(1, QFormLayout.LabelRole, self.checkBox_8)

        self.lineEdit_8 = QLineEdit(self.frame)
        self.lineEdit_8.setObjectName(u"lineEdit_8")

        self.formLayout.setWidget(1, QFormLayout.FieldRole, self.lineEdit_8)

        self.checkBox_9 = QCheckBox(self.frame)
        self.checkBox_9.setObjectName(u"checkBox_9")
        self.checkBox_9.setAutoExclusive(True)

        self.formLayout.setWidget(2, QFormLayout.LabelRole, self.checkBox_9)

        self.lineEdit_9 = QLineEdit(self.frame)
        self.lineEdit_9.setObjectName(u"lineEdit_9")

        self.formLayout.setWidget(2, QFormLayout.FieldRole, self.lineEdit_9)

        self.checkBox_10 = QCheckBox(self.frame)
        self.checkBox_10.setObjectName(u"checkBox_10")
        self.checkBox_10.setAutoExclusive(True)

        self.formLayout.setWidget(3, QFormLayout.LabelRole, self.checkBox_10)

        self.lineEdit_10 = QLineEdit(self.frame)
        self.lineEdit_10.setObjectName(u"lineEdit_10")

        self.formLayout.setWidget(3, QFormLayout.FieldRole, self.lineEdit_10)

        self.checkBox_11 = QCheckBox(self.frame)
        self.checkBox_11.setObjectName(u"checkBox_11")
        self.checkBox_11.setAutoExclusive(True)

        self.formLayout.setWidget(4, QFormLayout.LabelRole, self.checkBox_11)

        self.lineEdit_11 = QLineEdit(self.frame)
        self.lineEdit_11.setObjectName(u"lineEdit_11")

        self.formLayout.setWidget(4, QFormLayout.FieldRole, self.lineEdit_11)

        self.checkBox_12 = QCheckBox(self.frame)
        self.checkBox_12.setObjectName(u"checkBox_12")
        self.checkBox_12.setAutoExclusive(True)

        self.formLayout.setWidget(5, QFormLayout.LabelRole, self.checkBox_12)

        self.lineEdit_12 = QLineEdit(self.frame)
        self.lineEdit_12.setObjectName(u"lineEdit_12")

        self.formLayout.setWidget(5, QFormLayout.FieldRole, self.lineEdit_12)

        self.pushButton = QPushButton(self.frame)
        self.pushButton.setObjectName(u"pushButton")
        self.pushButton.setAutoDefault(False)

        self.formLayout.setWidget(6, QFormLayout.FieldRole, self.pushButton)

        self.checkBox_13 = QCheckBox(self.frame)
        self.checkBox_13.setObjectName(u"checkBox_13")
        self.checkBox_13.setAutoExclusive(True)

        self.lineEdit_13 = QLineEdit(self.frame)
        self.formLayout.insertRow(3, self.checkBox_13, self.lineEdit_13)
        self.lineEdit_13.setObjectName(u"lineEdit_13")

        self.verticalLayout.addWidget(self.frame)


        self.horizontalLayout.addWidget(self.frame_2)


        self.retranslateUi(Dialog)

        QMetaObject.connectSlotsByName(Dialog)
    # setupUi

    def retranslateUi(self, Dialog):
        Dialog.setWindowTitle(QCoreApplication.translate("Dialog", u"Dialog", None))
        self.comboBox.setItemText(0, QCoreApplication.translate("Dialog", u"Crystal structures db", None))
        self.comboBox.setItemText(1, QCoreApplication.translate("Dialog", u"QM structures db", None))

        self.checkBox.setText(QCoreApplication.translate("Dialog", u"Refcode", None))
        self.checkBox_8.setText(QCoreApplication.translate("Dialog", u"Name", None))
        self.checkBox_9.setText(QCoreApplication.translate("Dialog", u"Elements", None))
        self.checkBox_13.setText(QCoreApplication.translate("Dialog", u"Formula", None))
        self.checkBox_10.setText(QCoreApplication.translate("Dialog", u"Doi", None))
        self.checkBox_11.setText(QCoreApplication.translate("Dialog", u"Authors", None))
        self.checkBox_12.setText(QCoreApplication.translate("Dialog", u"Cell", None))
        self.pushButton.setText(QCoreApplication.translate("Dialog", u"Search", None))
    # retranslateUi

