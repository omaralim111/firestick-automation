# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '_settingsdialog.ui'
#
# Created: Sun Jun 20 01:18:29 2010
#      by: PyQt4 UI code generator 4.7.2
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

class Ui_SettingsDialog(object):
    def setupUi(self, SettingsDialog):
        SettingsDialog.setObjectName("SettingsDialog")
        SettingsDialog.setWindowModality(QtCore.Qt.ApplicationModal)
        SettingsDialog.resize(432, 256)
        self.gridLayout = QtGui.QGridLayout(SettingsDialog)
        self.gridLayout.setObjectName("gridLayout")
        self.groupBox_2 = QtGui.QGroupBox(SettingsDialog)
        self.groupBox_2.setObjectName("groupBox_2")
        self.horizontalLayout_2 = QtGui.QHBoxLayout(self.groupBox_2)
        self.horizontalLayout_2.setObjectName("horizontalLayout_2")
        self.nameFormatLineEdit = QtGui.QLineEdit(self.groupBox_2)
        self.nameFormatLineEdit.setReadOnly(True)
        self.nameFormatLineEdit.setObjectName("nameFormatLineEdit")
        self.horizontalLayout_2.addWidget(self.nameFormatLineEdit)
        self.nameFormatButton = QtGui.QPushButton(self.groupBox_2)
        self.nameFormatButton.setObjectName("nameFormatButton")
        self.horizontalLayout_2.addWidget(self.nameFormatButton)
        self.gridLayout.addWidget(self.groupBox_2, 1, 0, 1, 2)
        self.groupBox_3 = QtGui.QGroupBox(SettingsDialog)
        self.groupBox_3.setObjectName("groupBox_3")
        self.horizontalLayout_3 = QtGui.QHBoxLayout(self.groupBox_3)
        self.horizontalLayout_3.setObjectName("horizontalLayout_3")
        self.numThreadsSpinBox = QtGui.QSpinBox(self.groupBox_3)
        self.numThreadsSpinBox.setMinimum(1)
        self.numThreadsSpinBox.setMaximum(5)
        self.numThreadsSpinBox.setObjectName("numThreadsSpinBox")
        self.horizontalLayout_3.addWidget(self.numThreadsSpinBox)
        self.gridLayout.addWidget(self.groupBox_3, 2, 0, 1, 1)
        self.groupBox_4 = QtGui.QGroupBox(SettingsDialog)
        self.groupBox_4.setObjectName("groupBox_4")
        self.horizontalLayout_4 = QtGui.QHBoxLayout(self.groupBox_4)
        self.horizontalLayout_4.setObjectName("horizontalLayout_4")
        self.toolkitComboBox = QtGui.QComboBox(self.groupBox_4)
        self.toolkitComboBox.setObjectName("toolkitComboBox")
        self.toolkitComboBox.addItem("")
        self.toolkitComboBox.addItem("")
        self.horizontalLayout_4.addWidget(self.toolkitComboBox)
        self.gridLayout.addWidget(self.groupBox_4, 2, 1, 1, 1)
        self.buttonBox = QtGui.QDialogButtonBox(SettingsDialog)
        self.buttonBox.setOrientation(QtCore.Qt.Horizontal)
        self.buttonBox.setStandardButtons(QtGui.QDialogButtonBox.Cancel|QtGui.QDialogButtonBox.Ok)
        self.buttonBox.setObjectName("buttonBox")
        self.gridLayout.addWidget(self.buttonBox, 3, 1, 1, 1)
        self.groupBox = QtGui.QGroupBox(SettingsDialog)
        self.groupBox.setObjectName("groupBox")
        self.horizontalLayout = QtGui.QHBoxLayout(self.groupBox)
        self.horizontalLayout.setObjectName("horizontalLayout")
        self.saveDirLineEdit = QtGui.QLineEdit(self.groupBox)
        self.saveDirLineEdit.setReadOnly(True)
        self.saveDirLineEdit.setObjectName("saveDirLineEdit")
        self.horizontalLayout.addWidget(self.saveDirLineEdit)
        self.saveDirButton = QtGui.QPushButton(self.groupBox)
        self.saveDirButton.setObjectName("saveDirButton")
        self.horizontalLayout.addWidget(self.saveDirButton)
        self.gridLayout.addWidget(self.groupBox, 0, 0, 1, 2)

        self.retranslateUi(SettingsDialog)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("accepted()"), SettingsDialog.accept)
        QtCore.QObject.connect(self.buttonBox, QtCore.SIGNAL("rejected()"), SettingsDialog.reject)
        QtCore.QMetaObject.connectSlotsByName(SettingsDialog)

    def retranslateUi(self, SettingsDialog):
        SettingsDialog.setWindowTitle(QtGui.QApplication.translate("SettingsDialog", "Preferences", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_2.setTitle(QtGui.QApplication.translate("SettingsDialog", "Name Format", None, QtGui.QApplication.UnicodeUTF8))
        self.nameFormatButton.setText(QtGui.QApplication.translate("SettingsDialog", "Change", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_3.setTitle(QtGui.QApplication.translate("SettingsDialog", "Downloader Threads", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_4.setTitle(QtGui.QApplication.translate("SettingsDialog", "Preferred GUI Backend", None, QtGui.QApplication.UnicodeUTF8))
        self.toolkitComboBox.setItemText(0, QtGui.QApplication.translate("SettingsDialog", "qt4", None, QtGui.QApplication.UnicodeUTF8))
        self.toolkitComboBox.setItemText(1, QtGui.QApplication.translate("SettingsDialog", "gtk", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setTitle(QtGui.QApplication.translate("SettingsDialog", "Save Directory", None, QtGui.QApplication.UnicodeUTF8))
        self.saveDirButton.setText(QtGui.QApplication.translate("SettingsDialog", "Change", None, QtGui.QApplication.UnicodeUTF8))

