# -*- coding: utf-8 -*-

################################################################################
## Form generated from reading UI file 'bbarolowindow.ui'
##
## Created by: Qt User Interface Compiler version 6.6.1
##
## WARNING! All changes made in this file will be lost when recompiling UI file!
################################################################################

from PySide6.QtCore import (QCoreApplication, QDate, QDateTime, QLocale,
    QMetaObject, QObject, QPoint, QRect,
    QSize, QTime, QUrl, Qt)
from PySide6.QtGui import (QAction, QBrush, QColor, QConicalGradient,
    QCursor, QFont, QFontDatabase, QGradient,
    QIcon, QImage, QKeySequence, QLinearGradient,
    QPainter, QPalette, QPixmap, QRadialGradient,
    QTransform)
from PySide6.QtWidgets import (QAbstractSpinBox, QApplication, QCheckBox, QComboBox,
    QDoubleSpinBox, QFrame, QGraphicsView, QGridLayout,
    QGroupBox, QHBoxLayout, QLabel, QLayout,
    QLineEdit, QListView, QListWidget, QListWidgetItem,
    QMainWindow, QMenu, QMenuBar, QPushButton,
    QSizePolicy, QSpacerItem, QSpinBox, QStackedWidget,
    QStatusBar, QTextEdit, QToolButton, QVBoxLayout,
    QWidget)
import resources_rc

class Ui_BBaroloWindow(object):
    def setupUi(self, BBaroloWindow):
        if not BBaroloWindow.objectName():
            BBaroloWindow.setObjectName(u"BBaroloWindow")
        BBaroloWindow.setEnabled(True)
        BBaroloWindow.resize(874, 860)
        sizePolicy = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Minimum)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(BBaroloWindow.sizePolicy().hasHeightForWidth())
        BBaroloWindow.setSizePolicy(sizePolicy)
        BBaroloWindow.setLayoutDirection(Qt.LeftToRight)
        BBaroloWindow.setUnifiedTitleAndToolBarOnMac(False)
        self.actionOpen_FITS_file = QAction(BBaroloWindow)
        self.actionOpen_FITS_file.setObjectName(u"actionOpen_FITS_file")
        self.actionOpen_parameter_file = QAction(BBaroloWindow)
        self.actionOpen_parameter_file.setObjectName(u"actionOpen_parameter_file")
        self.actionExport_parameter_file = QAction(BBaroloWindow)
        self.actionExport_parameter_file.setObjectName(u"actionExport_parameter_file")
        self.actionOpen_list_of_FITS_files = QAction(BBaroloWindow)
        self.actionOpen_list_of_FITS_files.setObjectName(u"actionOpen_list_of_FITS_files")
        self.actionReset_GUI = QAction(BBaroloWindow)
        self.actionReset_GUI.setObjectName(u"actionReset_GUI")
        self.centralWidget = QWidget(BBaroloWindow)
        self.centralWidget.setObjectName(u"centralWidget")
        self.gridLayout = QGridLayout(self.centralWidget)
        self.gridLayout.setSpacing(6)
        self.gridLayout.setContentsMargins(11, 11, 11, 11)
        self.gridLayout.setObjectName(u"gridLayout")
        self.gridLayout_4 = QGridLayout()
        self.gridLayout_4.setSpacing(6)
        self.gridLayout_4.setObjectName(u"gridLayout_4")
        self.gridLayout_4.setSizeConstraint(QLayout.SetFixedSize)
        self.gridLayout_4.setVerticalSpacing(5)
        self.ParampushButton = QPushButton(self.centralWidget)
        self.ParampushButton.setObjectName(u"ParampushButton")

        self.gridLayout_4.addWidget(self.ParampushButton, 0, 2, 1, 1)

        self.BoxcheckBox = QCheckBox(self.centralWidget)
        self.BoxcheckBox.setObjectName(u"BoxcheckBox")

        self.gridLayout_4.addWidget(self.BoxcheckBox, 2, 0, 1, 1)

        self.Inputlabel = QLabel(self.centralWidget)
        self.Inputlabel.setObjectName(u"Inputlabel")

        self.gridLayout_4.addWidget(self.Inputlabel, 0, 0, 1, 1)

        self.FitspushButton = QPushButton(self.centralWidget)
        self.FitspushButton.setObjectName(u"FitspushButton")

        self.gridLayout_4.addWidget(self.FitspushButton, 1, 2, 1, 1)

        self.FitslineEdit = QLineEdit(self.centralWidget)
        self.FitslineEdit.setObjectName(u"FitslineEdit")

        self.gridLayout_4.addWidget(self.FitslineEdit, 1, 1, 1, 1)

        self.horizontalSpacer_25 = QSpacerItem(40, 20, QSizePolicy.Ignored, QSizePolicy.Minimum)

        self.gridLayout_4.addItem(self.horizontalSpacer_25, 3, 2, 1, 1)

        self.RunpushButton = QPushButton(self.centralWidget)
        self.RunpushButton.setObjectName(u"RunpushButton")
        sizePolicy1 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Maximum)
        sizePolicy1.setHorizontalStretch(0)
        sizePolicy1.setVerticalStretch(0)
        sizePolicy1.setHeightForWidth(self.RunpushButton.sizePolicy().hasHeightForWidth())
        self.RunpushButton.setSizePolicy(sizePolicy1)
        self.RunpushButton.setMinimumSize(QSize(80, 0))

        self.gridLayout_4.addWidget(self.RunpushButton, 2, 4, 2, 1)

        self.verticalSpacer = QSpacerItem(20, 40, QSizePolicy.Minimum, QSizePolicy.Expanding)

        self.gridLayout_4.addItem(self.verticalSpacer, 0, 3, 3, 1)

        self.ParamlineEdit = QLineEdit(self.centralWidget)
        self.ParamlineEdit.setObjectName(u"ParamlineEdit")
        sizePolicy1.setHeightForWidth(self.ParamlineEdit.sizePolicy().hasHeightForWidth())
        self.ParamlineEdit.setSizePolicy(sizePolicy1)

        self.gridLayout_4.addWidget(self.ParamlineEdit, 0, 1, 1, 1)

        self.AutocheckBox = QCheckBox(self.centralWidget)
        self.AutocheckBox.setObjectName(u"AutocheckBox")
        sizePolicy2 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Minimum)
        sizePolicy2.setHorizontalStretch(0)
        sizePolicy2.setVerticalStretch(0)
        sizePolicy2.setHeightForWidth(self.AutocheckBox.sizePolicy().hasHeightForWidth())
        self.AutocheckBox.setSizePolicy(sizePolicy2)
        self.AutocheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_4.addWidget(self.AutocheckBox, 3, 0, 1, 2)

        self.Fitslabel = QLabel(self.centralWidget)
        self.Fitslabel.setObjectName(u"Fitslabel")

        self.gridLayout_4.addWidget(self.Fitslabel, 1, 0, 1, 1)

        self.horizontalLayout_3 = QHBoxLayout()
        self.horizontalLayout_3.setSpacing(5)
        self.horizontalLayout_3.setObjectName(u"horizontalLayout_3")
        self.horizontalLayout_3.setContentsMargins(-1, -1, -1, 0)
        self.Xboxlabel = QLabel(self.centralWidget)
        self.Xboxlabel.setObjectName(u"Xboxlabel")
        self.Xboxlabel.setEnabled(False)
        sizePolicy3 = QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Maximum)
        sizePolicy3.setHorizontalStretch(0)
        sizePolicy3.setVerticalStretch(0)
        sizePolicy3.setHeightForWidth(self.Xboxlabel.sizePolicy().hasHeightForWidth())
        self.Xboxlabel.setSizePolicy(sizePolicy3)

        self.horizontalLayout_3.addWidget(self.Xboxlabel)

        self.XminspinBox = QSpinBox(self.centralWidget)
        self.XminspinBox.setObjectName(u"XminspinBox")
        self.XminspinBox.setEnabled(False)
        font = QFont()
        font.setPointSize(12)
        self.XminspinBox.setFont(font)
        self.XminspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.XminspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.XminspinBox.setMaximum(99999)

        self.horizontalLayout_3.addWidget(self.XminspinBox)

        self.label_19 = QLabel(self.centralWidget)
        self.label_19.setObjectName(u"label_19")
        self.label_19.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.label_19.sizePolicy().hasHeightForWidth())
        self.label_19.setSizePolicy(sizePolicy3)

        self.horizontalLayout_3.addWidget(self.label_19)

        self.XmaxspinBox = QSpinBox(self.centralWidget)
        self.XmaxspinBox.setObjectName(u"XmaxspinBox")
        self.XmaxspinBox.setEnabled(False)
        self.XmaxspinBox.setFont(font)
        self.XmaxspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.XmaxspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.XmaxspinBox.setMaximum(99999)

        self.horizontalLayout_3.addWidget(self.XmaxspinBox)

        self.Yboxlabel = QLabel(self.centralWidget)
        self.Yboxlabel.setObjectName(u"Yboxlabel")
        self.Yboxlabel.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.Yboxlabel.sizePolicy().hasHeightForWidth())
        self.Yboxlabel.setSizePolicy(sizePolicy3)

        self.horizontalLayout_3.addWidget(self.Yboxlabel)

        self.YminspinBox = QSpinBox(self.centralWidget)
        self.YminspinBox.setObjectName(u"YminspinBox")
        self.YminspinBox.setEnabled(False)
        self.YminspinBox.setFont(font)
        self.YminspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.YminspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.YminspinBox.setMaximum(99999)
        self.YminspinBox.setValue(0)

        self.horizontalLayout_3.addWidget(self.YminspinBox)

        self.label_23 = QLabel(self.centralWidget)
        self.label_23.setObjectName(u"label_23")
        self.label_23.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.label_23.sizePolicy().hasHeightForWidth())
        self.label_23.setSizePolicy(sizePolicy3)

        self.horizontalLayout_3.addWidget(self.label_23)

        self.YmaxspinBox = QSpinBox(self.centralWidget)
        self.YmaxspinBox.setObjectName(u"YmaxspinBox")
        self.YmaxspinBox.setEnabled(False)
        self.YmaxspinBox.setFont(font)
        self.YmaxspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.YmaxspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.YmaxspinBox.setMaximum(99999)
        self.YmaxspinBox.setValue(0)

        self.horizontalLayout_3.addWidget(self.YmaxspinBox)

        self.Zboxlabel = QLabel(self.centralWidget)
        self.Zboxlabel.setObjectName(u"Zboxlabel")
        self.Zboxlabel.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.Zboxlabel.sizePolicy().hasHeightForWidth())
        self.Zboxlabel.setSizePolicy(sizePolicy3)

        self.horizontalLayout_3.addWidget(self.Zboxlabel)

        self.ZminspinBox = QSpinBox(self.centralWidget)
        self.ZminspinBox.setObjectName(u"ZminspinBox")
        self.ZminspinBox.setEnabled(False)
        self.ZminspinBox.setFont(font)
        self.ZminspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ZminspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.ZminspinBox.setMaximum(99999)
        self.ZminspinBox.setValue(0)

        self.horizontalLayout_3.addWidget(self.ZminspinBox)

        self.label_29 = QLabel(self.centralWidget)
        self.label_29.setObjectName(u"label_29")
        self.label_29.setEnabled(False)
        sizePolicy3.setHeightForWidth(self.label_29.sizePolicy().hasHeightForWidth())
        self.label_29.setSizePolicy(sizePolicy3)

        self.horizontalLayout_3.addWidget(self.label_29)

        self.ZmaxspinBox = QSpinBox(self.centralWidget)
        self.ZmaxspinBox.setObjectName(u"ZmaxspinBox")
        self.ZmaxspinBox.setEnabled(False)
        self.ZmaxspinBox.setFont(font)
        self.ZmaxspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ZmaxspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.ZmaxspinBox.setMaximum(99999)
        self.ZmaxspinBox.setValue(0)

        self.horizontalLayout_3.addWidget(self.ZmaxspinBox)


        self.gridLayout_4.addLayout(self.horizontalLayout_3, 2, 1, 1, 1)

        self.ResetpushButton = QPushButton(self.centralWidget)
        self.ResetpushButton.setObjectName(u"ResetpushButton")

        self.gridLayout_4.addWidget(self.ResetpushButton, 2, 2, 1, 1)

        self.gridLayout_25 = QGridLayout()
        self.gridLayout_25.setSpacing(6)
        self.gridLayout_25.setObjectName(u"gridLayout_25")
        self.ExtenspinBox = QSpinBox(self.centralWidget)
        self.ExtenspinBox.setObjectName(u"ExtenspinBox")
        self.ExtenspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ExtenspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.gridLayout_25.addWidget(self.ExtenspinBox, 1, 1, 1, 1)

        self.Extenlabel = QLabel(self.centralWidget)
        self.Extenlabel.setObjectName(u"Extenlabel")

        self.gridLayout_25.addWidget(self.Extenlabel, 1, 0, 1, 1)

        self.Threadslabel = QLabel(self.centralWidget)
        self.Threadslabel.setObjectName(u"Threadslabel")
        sizePolicy4 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Preferred)
        sizePolicy4.setHorizontalStretch(0)
        sizePolicy4.setVerticalStretch(0)
        sizePolicy4.setHeightForWidth(self.Threadslabel.sizePolicy().hasHeightForWidth())
        self.Threadslabel.setSizePolicy(sizePolicy4)

        self.gridLayout_25.addWidget(self.Threadslabel, 0, 0, 1, 1)

        self.ThreadspinBox = QSpinBox(self.centralWidget)
        self.ThreadspinBox.setObjectName(u"ThreadspinBox")
        sizePolicy5 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Fixed)
        sizePolicy5.setHorizontalStretch(0)
        sizePolicy5.setVerticalStretch(0)
        sizePolicy5.setHeightForWidth(self.ThreadspinBox.sizePolicy().hasHeightForWidth())
        self.ThreadspinBox.setSizePolicy(sizePolicy5)
        self.ThreadspinBox.setMinimumSize(QSize(60, 20))
        self.ThreadspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ThreadspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.ThreadspinBox.setMinimum(1)
        self.ThreadspinBox.setMaximum(4096)

        self.gridLayout_25.addWidget(self.ThreadspinBox, 0, 1, 1, 1)


        self.gridLayout_4.addLayout(self.gridLayout_25, 0, 4, 2, 1)


        self.gridLayout.addLayout(self.gridLayout_4, 0, 1, 1, 1)

        self.gridLayout_26 = QGridLayout()
        self.gridLayout_26.setSpacing(6)
        self.gridLayout_26.setObjectName(u"gridLayout_26")
        self.AdvancedOptcheckBox = QCheckBox(self.centralWidget)
        self.AdvancedOptcheckBox.setObjectName(u"AdvancedOptcheckBox")
        self.AdvancedOptcheckBox.setChecked(False)

        self.gridLayout_26.addWidget(self.AdvancedOptcheckBox, 1, 0, 1, 1)

        self.horizontalSpacer_2 = QSpacerItem(100, 20, QSizePolicy.Fixed, QSizePolicy.Minimum)

        self.gridLayout_26.addItem(self.horizontalSpacer_2, 0, 3, 1, 1)

        self.KillpushButton = QPushButton(self.centralWidget)
        self.KillpushButton.setObjectName(u"KillpushButton")
        sizePolicy6 = QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Fixed)
        sizePolicy6.setHorizontalStretch(0)
        sizePolicy6.setVerticalStretch(0)
        sizePolicy6.setHeightForWidth(self.KillpushButton.sizePolicy().hasHeightForWidth())
        self.KillpushButton.setSizePolicy(sizePolicy6)

        self.gridLayout_26.addWidget(self.KillpushButton, 0, 4, 1, 1)

        self.OutfolderlineEdit = QLineEdit(self.centralWidget)
        self.OutfolderlineEdit.setObjectName(u"OutfolderlineEdit")
        self.OutfolderlineEdit.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.gridLayout_26.addWidget(self.OutfolderlineEdit, 0, 1, 1, 1)

        self.OutfolderpushButton = QPushButton(self.centralWidget)
        self.OutfolderpushButton.setObjectName(u"OutfolderpushButton")

        self.gridLayout_26.addWidget(self.OutfolderpushButton, 0, 2, 1, 1)

        self.Outfolderlabel = QLabel(self.centralWidget)
        self.Outfolderlabel.setObjectName(u"Outfolderlabel")

        self.gridLayout_26.addWidget(self.Outfolderlabel, 0, 0, 1, 1)

        self.AdvancedOptlineEdit = QLineEdit(self.centralWidget)
        self.AdvancedOptlineEdit.setObjectName(u"AdvancedOptlineEdit")
        self.AdvancedOptlineEdit.setEnabled(False)

        self.gridLayout_26.addWidget(self.AdvancedOptlineEdit, 1, 1, 1, 2)


        self.gridLayout.addLayout(self.gridLayout_26, 3, 1, 1, 1)

        self.MaingridLayout = QGridLayout()
        self.MaingridLayout.setSpacing(6)
        self.MaingridLayout.setObjectName(u"MaingridLayout")
        self.MaingridLayout.setHorizontalSpacing(6)
        self.MaingridLayout.setVerticalSpacing(0)
        self.stackedWidget = QStackedWidget(self.centralWidget)
        self.stackedWidget.setObjectName(u"stackedWidget")
        self.stackedWidget.setEnabled(True)
        self.stackedWidget.setFrameShape(QFrame.Panel)
        self.stackedWidget.setFrameShadow(QFrame.Plain)
        self.page = QWidget()
        self.page.setObjectName(u"page")
        self.verticalLayout_7 = QVBoxLayout(self.page)
        self.verticalLayout_7.setSpacing(6)
        self.verticalLayout_7.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_7.setObjectName(u"verticalLayout_7")
        self.textEdit = QTextEdit(self.page)
        self.textEdit.setObjectName(u"textEdit")
        self.textEdit.setEnabled(True)
        palette = QPalette()
        brush = QBrush(QColor(255, 255, 255, 0))
        brush.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Active, QPalette.Base, brush)
        palette.setBrush(QPalette.Inactive, QPalette.Base, brush)
        brush1 = QBrush(QColor(237, 237, 237, 255))
        brush1.setStyle(Qt.SolidPattern)
        palette.setBrush(QPalette.Disabled, QPalette.Base, brush1)
        self.textEdit.setPalette(palette)
        self.textEdit.setFrameShape(QFrame.NoFrame)
        self.textEdit.setReadOnly(True)

        self.verticalLayout_7.addWidget(self.textEdit)

        self.Logolabel = QLabel(self.page)
        self.Logolabel.setObjectName(u"Logolabel")
        self.Logolabel.setEnabled(False)
        self.Logolabel.setPixmap(QPixmap(u":/resources/bbarolo.png"))

        self.verticalLayout_7.addWidget(self.Logolabel, 0, Qt.AlignRight|Qt.AlignBottom)

        self.stackedWidget.addWidget(self.page)
        self.Galfitpage = QWidget()
        self.Galfitpage.setObjectName(u"Galfitpage")
        self.gridLayout_3 = QGridLayout(self.Galfitpage)
        self.gridLayout_3.setSpacing(6)
        self.gridLayout_3.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_3.setObjectName(u"gridLayout_3")
        self.GalfitgroupBox = QGroupBox(self.Galfitpage)
        self.GalfitgroupBox.setObjectName(u"GalfitgroupBox")
        self.GalfitgroupBox.setEnabled(False)
        self.GalfitgroupBox.setCheckable(False)
        self.GalfitgroupBox.setChecked(False)
        self.gridLayout_14 = QGridLayout(self.GalfitgroupBox)
        self.gridLayout_14.setSpacing(6)
        self.gridLayout_14.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_14.setObjectName(u"gridLayout_14")
        self.gridLayout_19 = QGridLayout()
        self.gridLayout_19.setSpacing(6)
        self.gridLayout_19.setObjectName(u"gridLayout_19")
        self.PacheckBox = QCheckBox(self.GalfitgroupBox)
        self.PacheckBox.setObjectName(u"PacheckBox")
        self.PacheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.PacheckBox, 9, 0, 1, 1)

        self.RadseppushButton = QPushButton(self.GalfitgroupBox)
        self.RadseppushButton.setObjectName(u"RadseppushButton")

        self.gridLayout_19.addWidget(self.RadseppushButton, 0, 5, 2, 1)

        self.DenscheckBox = QCheckBox(self.GalfitgroupBox)
        self.DenscheckBox.setObjectName(u"DenscheckBox")

        self.gridLayout_19.addWidget(self.DenscheckBox, 11, 0, 1, 1)

        self.YpospushButton = QPushButton(self.GalfitgroupBox)
        self.YpospushButton.setObjectName(u"YpospushButton")

        self.gridLayout_19.addWidget(self.YpospushButton, 3, 5, 1, 1)

        self.PaPMlabel = QLabel(self.GalfitgroupBox)
        self.PaPMlabel.setObjectName(u"PaPMlabel")
        self.PaPMlabel.setAlignment(Qt.AlignCenter)

        self.gridLayout_19.addWidget(self.PaPMlabel, 9, 3, 1, 1, Qt.AlignHCenter)

        self.VrotpushButton = QPushButton(self.GalfitgroupBox)
        self.VrotpushButton.setObjectName(u"VrotpushButton")

        self.gridLayout_19.addWidget(self.VrotpushButton, 5, 5, 1, 1)

        self.wcscomboBox = QComboBox(self.GalfitgroupBox)
        self.wcscomboBox.addItem("")
        self.wcscomboBox.addItem("")
        self.wcscomboBox.addItem("")
        self.wcscomboBox.setObjectName(u"wcscomboBox")

        self.gridLayout_19.addWidget(self.wcscomboBox, 2, 3, 2, 2)

        self.PapushButton = QPushButton(self.GalfitgroupBox)
        self.PapushButton.setObjectName(u"PapushButton")

        self.gridLayout_19.addWidget(self.PapushButton, 9, 5, 1, 1)

        self.VsyspushButton = QPushButton(self.GalfitgroupBox)
        self.VsyspushButton.setObjectName(u"VsyspushButton")

        self.gridLayout_19.addWidget(self.VsyspushButton, 4, 5, 1, 1)

        self.VrotSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.VrotSpinBox.setObjectName(u"VrotSpinBox")
        self.VrotSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.VrotSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.VrotSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.VrotSpinBox.setMaximum(10000.000000000000000)

        self.gridLayout_19.addWidget(self.VrotSpinBox, 5, 2, 1, 1)

        self.IncFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.IncFilespinBox.setObjectName(u"IncFilespinBox")
        self.IncFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.IncFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.IncFilespinBox.setMinimum(1)
        self.IncFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.IncFilespinBox, 8, 7, 1, 1)

        self.label_25 = QLabel(self.GalfitgroupBox)
        self.label_25.setObjectName(u"label_25")

        self.gridLayout_19.addWidget(self.label_25, 5, 3, 1, 1)

        self.PaFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.PaFilelineEdit.setObjectName(u"PaFilelineEdit")

        self.gridLayout_19.addWidget(self.PaFilelineEdit, 9, 6, 1, 1)

        self.NringsspinBox = QSpinBox(self.GalfitgroupBox)
        self.NringsspinBox.setObjectName(u"NringsspinBox")
        self.NringsspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.NringsspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.NringsspinBox.setMaximum(500)

        self.gridLayout_19.addWidget(self.NringsspinBox, 0, 2, 1, 1)

        self.Z0checkBox = QCheckBox(self.GalfitgroupBox)
        self.Z0checkBox.setObjectName(u"Z0checkBox")
        self.Z0checkBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.Z0checkBox, 10, 0, 1, 1)

        self.YposFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.YposFilespinBox.setObjectName(u"YposFilespinBox")
        self.YposFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.YposFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.YposFilespinBox.setMinimum(1)
        self.YposFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.YposFilespinBox, 3, 7, 1, 1)

        self.label_18 = QLabel(self.GalfitgroupBox)
        self.label_18.setObjectName(u"label_18")
        sizePolicy4.setHeightForWidth(self.label_18.sizePolicy().hasHeightForWidth())
        self.label_18.setSizePolicy(sizePolicy4)

        self.gridLayout_19.addWidget(self.label_18, 11, 4, 1, 1, Qt.AlignLeft)

        self.VradpushButton = QPushButton(self.GalfitgroupBox)
        self.VradpushButton.setObjectName(u"VradpushButton")

        self.gridLayout_19.addWidget(self.VradpushButton, 7, 5, 1, 1)

        self.PaSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.PaSpinBox.setObjectName(u"PaSpinBox")
        self.PaSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.PaSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.PaSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.PaSpinBox.setMinimum(-360.000000000000000)
        self.PaSpinBox.setMaximum(360.000000000000000)

        self.gridLayout_19.addWidget(self.PaSpinBox, 9, 2, 1, 1)

        self.InccheckBox = QCheckBox(self.GalfitgroupBox)
        self.InccheckBox.setObjectName(u"InccheckBox")
        self.InccheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.InccheckBox, 8, 0, 1, 1)

        self.Z0SpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.Z0SpinBox.setObjectName(u"Z0SpinBox")
        self.Z0SpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.Z0SpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.Z0SpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.Z0SpinBox.setMaximum(10000.000000000000000)

        self.gridLayout_19.addWidget(self.Z0SpinBox, 10, 2, 1, 1)

        self.VradFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.VradFilelineEdit.setObjectName(u"VradFilelineEdit")

        self.gridLayout_19.addWidget(self.VradFilelineEdit, 7, 6, 1, 1)

        self.PaDSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.PaDSpinBox.setObjectName(u"PaDSpinBox")
        self.PaDSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.PaDSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.PaDSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.PaDSpinBox.setMaximum(360.000000000000000)
        self.PaDSpinBox.setValue(15.000000000000000)

        self.gridLayout_19.addWidget(self.PaDSpinBox, 9, 4, 1, 1)

        self.VradSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.VradSpinBox.setObjectName(u"VradSpinBox")
        self.VradSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.VradSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)

        self.gridLayout_19.addWidget(self.VradSpinBox, 7, 2, 1, 1)

        self.RadsepSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.RadsepSpinBox.setObjectName(u"RadsepSpinBox")
        self.RadsepSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.RadsepSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.RadsepSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.RadsepSpinBox.setMaximum(10000.000000000000000)

        self.gridLayout_19.addWidget(self.RadsepSpinBox, 1, 2, 1, 1)

        self.IncSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.IncSpinBox.setObjectName(u"IncSpinBox")
        self.IncSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.IncSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.IncSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.IncSpinBox.setMinimum(-90.000000000000000)
        self.IncSpinBox.setMaximum(90.000000000000000)

        self.gridLayout_19.addWidget(self.IncSpinBox, 8, 2, 1, 1)

        self.XpospushButton = QPushButton(self.GalfitgroupBox)
        self.XpospushButton.setObjectName(u"XpospushButton")

        self.gridLayout_19.addWidget(self.XpospushButton, 2, 5, 1, 1)

        self.label_13 = QLabel(self.GalfitgroupBox)
        self.label_13.setObjectName(u"label_13")

        self.gridLayout_19.addWidget(self.label_13, 1, 3, 1, 1)

        self.label_12 = QLabel(self.GalfitgroupBox)
        self.label_12.setObjectName(u"label_12")
        sizePolicy4.setHeightForWidth(self.label_12.sizePolicy().hasHeightForWidth())
        self.label_12.setSizePolicy(sizePolicy4)

        self.gridLayout_19.addWidget(self.label_12, 11, 3, 1, 1)

        self.Z0FilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.Z0FilelineEdit.setObjectName(u"Z0FilelineEdit")

        self.gridLayout_19.addWidget(self.Z0FilelineEdit, 10, 6, 1, 1)

        self.VdispSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.VdispSpinBox.setObjectName(u"VdispSpinBox")
        self.VdispSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.VdispSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.VdispSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.VdispSpinBox.setMaximum(10000.000000000000000)

        self.gridLayout_19.addWidget(self.VdispSpinBox, 6, 2, 1, 1)

        self.horizontalSpacer_18 = QSpacerItem(40, 20, QSizePolicy.Fixed, QSizePolicy.Minimum)

        self.gridLayout_19.addItem(self.horizontalSpacer_18, 0, 1, 1, 1)

        self.RadsepFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.RadsepFilelineEdit.setObjectName(u"RadsepFilelineEdit")

        self.gridLayout_19.addWidget(self.RadsepFilelineEdit, 0, 6, 2, 1)

        self.YposcheckBox = QCheckBox(self.GalfitgroupBox)
        self.YposcheckBox.setObjectName(u"YposcheckBox")
        self.YposcheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.YposcheckBox, 3, 0, 1, 1)

        self.VrotFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.VrotFilelineEdit.setObjectName(u"VrotFilelineEdit")

        self.gridLayout_19.addWidget(self.VrotFilelineEdit, 5, 6, 1, 1)

        self.IncpushButton = QPushButton(self.GalfitgroupBox)
        self.IncpushButton.setObjectName(u"IncpushButton")

        self.gridLayout_19.addWidget(self.IncpushButton, 8, 5, 1, 1)

        self.VrotFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.VrotFilespinBox.setObjectName(u"VrotFilespinBox")
        self.VrotFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.VrotFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.VrotFilespinBox.setMinimum(1)
        self.VrotFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.VrotFilespinBox, 5, 7, 1, 1)

        self.VdispcheckBox = QCheckBox(self.GalfitgroupBox)
        self.VdispcheckBox.setObjectName(u"VdispcheckBox")
        self.VdispcheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.VdispcheckBox, 6, 0, 1, 1)

        self.YposFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.YposFilelineEdit.setObjectName(u"YposFilelineEdit")

        self.gridLayout_19.addWidget(self.YposFilelineEdit, 3, 6, 1, 1)

        self.VrotcheckBox = QCheckBox(self.GalfitgroupBox)
        self.VrotcheckBox.setObjectName(u"VrotcheckBox")
        self.VrotcheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.VrotcheckBox, 5, 0, 1, 1)

        self.VdisppushButton = QPushButton(self.GalfitgroupBox)
        self.VdisppushButton.setObjectName(u"VdisppushButton")

        self.gridLayout_19.addWidget(self.VdisppushButton, 6, 5, 1, 1)

        self.IncPMlabel = QLabel(self.GalfitgroupBox)
        self.IncPMlabel.setObjectName(u"IncPMlabel")
        self.IncPMlabel.setAlignment(Qt.AlignCenter)

        self.gridLayout_19.addWidget(self.IncPMlabel, 8, 3, 1, 1, Qt.AlignHCenter)

        self.VsyscheckBox = QCheckBox(self.GalfitgroupBox)
        self.VsyscheckBox.setObjectName(u"VsyscheckBox")
        self.VsyscheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.VsyscheckBox, 4, 0, 1, 1)

        self.DensFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.DensFilelineEdit.setObjectName(u"DensFilelineEdit")

        self.gridLayout_19.addWidget(self.DensFilelineEdit, 11, 6, 1, 1)

        self.RadsepFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.RadsepFilespinBox.setObjectName(u"RadsepFilespinBox")
        self.RadsepFilespinBox.setLayoutDirection(Qt.LeftToRight)
        self.RadsepFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.RadsepFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.RadsepFilespinBox.setMinimum(1)
        self.RadsepFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.RadsepFilespinBox, 0, 7, 2, 1)

        self.label_24 = QLabel(self.GalfitgroupBox)
        self.label_24.setObjectName(u"label_24")
        sizePolicy4.setHeightForWidth(self.label_24.sizePolicy().hasHeightForWidth())
        self.label_24.setSizePolicy(sizePolicy4)

        self.gridLayout_19.addWidget(self.label_24, 6, 3, 1, 1)

        self.VradcheckBox = QCheckBox(self.GalfitgroupBox)
        self.VradcheckBox.setObjectName(u"VradcheckBox")

        self.gridLayout_19.addWidget(self.VradcheckBox, 7, 0, 1, 1)

        self.label_14 = QLabel(self.GalfitgroupBox)
        self.label_14.setObjectName(u"label_14")

        self.gridLayout_19.addWidget(self.label_14, 10, 3, 1, 1)

        self.RadsepcheckBox = QCheckBox(self.GalfitgroupBox)
        self.RadsepcheckBox.setObjectName(u"RadsepcheckBox")
        self.RadsepcheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.RadsepcheckBox, 1, 0, 1, 1)

        self.DensFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.DensFilespinBox.setObjectName(u"DensFilespinBox")
        self.DensFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.DensFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.DensFilespinBox.setMinimum(1)
        self.DensFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.DensFilespinBox, 11, 7, 1, 1)

        self.IncFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.IncFilelineEdit.setObjectName(u"IncFilelineEdit")

        self.gridLayout_19.addWidget(self.IncFilelineEdit, 8, 6, 1, 1)

        self.horizontalSpacer_5 = QSpacerItem(10, 5, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_19.addItem(self.horizontalSpacer_5, 12, 2, 1, 1)

        self.VsysSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.VsysSpinBox.setObjectName(u"VsysSpinBox")
        self.VsysSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.VsysSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.VsysSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.VsysSpinBox.setMinimum(-100000.000000000000000)
        self.VsysSpinBox.setMaximum(100000.000000000000000)

        self.gridLayout_19.addWidget(self.VsysSpinBox, 4, 2, 1, 1)

        self.Z0pushButton = QPushButton(self.GalfitgroupBox)
        self.Z0pushButton.setObjectName(u"Z0pushButton")

        self.gridLayout_19.addWidget(self.Z0pushButton, 10, 5, 1, 1)

        self.NringscheckBox = QCheckBox(self.GalfitgroupBox)
        self.NringscheckBox.setObjectName(u"NringscheckBox")
        self.NringscheckBox.setMouseTracking(True)
        self.NringscheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.NringscheckBox, 0, 0, 1, 1)

        self.VsysFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.VsysFilelineEdit.setObjectName(u"VsysFilelineEdit")

        self.gridLayout_19.addWidget(self.VsysFilelineEdit, 4, 6, 1, 1)

        self.label_26 = QLabel(self.GalfitgroupBox)
        self.label_26.setObjectName(u"label_26")

        self.gridLayout_19.addWidget(self.label_26, 4, 3, 1, 1)

        self.label_32 = QLabel(self.GalfitgroupBox)
        self.label_32.setObjectName(u"label_32")

        self.gridLayout_19.addWidget(self.label_32, 7, 3, 1, 1)

        self.PaFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.PaFilespinBox.setObjectName(u"PaFilespinBox")
        self.PaFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.PaFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.PaFilespinBox.setMinimum(1)
        self.PaFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.PaFilespinBox, 9, 7, 1, 1)

        self.XposFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.XposFilelineEdit.setObjectName(u"XposFilelineEdit")
        self.XposFilelineEdit.setAlignment(Qt.AlignLeading|Qt.AlignLeft|Qt.AlignVCenter)

        self.gridLayout_19.addWidget(self.XposFilelineEdit, 2, 6, 1, 1)

        self.VdispFilelineEdit = QLineEdit(self.GalfitgroupBox)
        self.VdispFilelineEdit.setObjectName(u"VdispFilelineEdit")

        self.gridLayout_19.addWidget(self.VdispFilelineEdit, 6, 6, 1, 1)

        self.IncDSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.IncDSpinBox.setObjectName(u"IncDSpinBox")
        self.IncDSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.IncDSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.IncDSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.IncDSpinBox.setCorrectionMode(QAbstractSpinBox.CorrectToPreviousValue)
        self.IncDSpinBox.setMaximum(90.000000000000000)
        self.IncDSpinBox.setValue(5.000000000000000)

        self.gridLayout_19.addWidget(self.IncDSpinBox, 8, 4, 1, 1)

        self.Z0FilespinBox = QSpinBox(self.GalfitgroupBox)
        self.Z0FilespinBox.setObjectName(u"Z0FilespinBox")
        self.Z0FilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.Z0FilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.Z0FilespinBox.setMinimum(1)
        self.Z0FilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.Z0FilespinBox, 10, 7, 1, 1)

        self.DensSpinBox = QDoubleSpinBox(self.GalfitgroupBox)
        self.DensSpinBox.setObjectName(u"DensSpinBox")
        self.DensSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.DensSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.DensSpinBox.setValue(1.000000000000000)

        self.gridLayout_19.addWidget(self.DensSpinBox, 11, 2, 1, 1)

        self.xposlineEdit = QLineEdit(self.GalfitgroupBox)
        self.xposlineEdit.setObjectName(u"xposlineEdit")
        self.xposlineEdit.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.xposlineEdit.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.gridLayout_19.addWidget(self.xposlineEdit, 2, 2, 1, 1)

        self.DenspushButton = QPushButton(self.GalfitgroupBox)
        self.DenspushButton.setObjectName(u"DenspushButton")

        self.gridLayout_19.addWidget(self.DenspushButton, 11, 5, 1, 1)

        self.XposFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.XposFilespinBox.setObjectName(u"XposFilespinBox")
        self.XposFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.XposFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.XposFilespinBox.setMinimum(1)
        self.XposFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.XposFilespinBox, 2, 7, 1, 1)

        self.yposlineEdit = QLineEdit(self.GalfitgroupBox)
        self.yposlineEdit.setObjectName(u"yposlineEdit")
        self.yposlineEdit.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.yposlineEdit.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.gridLayout_19.addWidget(self.yposlineEdit, 3, 2, 1, 1)

        self.VradFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.VradFilespinBox.setObjectName(u"VradFilespinBox")
        self.VradFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.VradFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.VradFilespinBox.setMinimum(1)

        self.gridLayout_19.addWidget(self.VradFilespinBox, 7, 7, 1, 1)

        self.VdispFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.VdispFilespinBox.setObjectName(u"VdispFilespinBox")
        self.VdispFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.VdispFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.VdispFilespinBox.setMinimum(1)
        self.VdispFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.VdispFilespinBox, 6, 7, 1, 1)

        self.XposcheckBox = QCheckBox(self.GalfitgroupBox)
        self.XposcheckBox.setObjectName(u"XposcheckBox")
        self.XposcheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_19.addWidget(self.XposcheckBox, 2, 0, 1, 1)

        self.GalfittoolButton = QToolButton(self.GalfitgroupBox)
        self.GalfittoolButton.setObjectName(u"GalfittoolButton")
        self.GalfittoolButton.setLayoutDirection(Qt.LeftToRight)
        self.GalfittoolButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.GalfittoolButton.setArrowType(Qt.RightArrow)

        self.gridLayout_19.addWidget(self.GalfittoolButton, 13, 0, 1, 8, Qt.AlignRight)

        self.VsysFilespinBox = QSpinBox(self.GalfitgroupBox)
        self.VsysFilespinBox.setObjectName(u"VsysFilespinBox")
        self.VsysFilespinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.VsysFilespinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.VsysFilespinBox.setMinimum(1)
        self.VsysFilespinBox.setValue(1)

        self.gridLayout_19.addWidget(self.VsysFilespinBox, 4, 7, 1, 1)


        self.gridLayout_14.addLayout(self.gridLayout_19, 0, 0, 1, 1)


        self.gridLayout_3.addWidget(self.GalfitgroupBox, 0, 1, 1, 1)

        self.stackedWidget.addWidget(self.Galfitpage)
        self.page_4 = QWidget()
        self.page_4.setObjectName(u"page_4")
        self.gridLayout_10 = QGridLayout(self.page_4)
        self.gridLayout_10.setSpacing(6)
        self.gridLayout_10.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_10.setObjectName(u"gridLayout_10")
        self.SearchgroupBox = QGroupBox(self.page_4)
        self.SearchgroupBox.setObjectName(u"SearchgroupBox")
        self.SearchgroupBox.setEnabled(False)
        self.SearchgroupBox.setCheckable(False)
        self.SearchgroupBox.setChecked(False)
        self.gridLayout_5 = QGridLayout(self.SearchgroupBox)
        self.gridLayout_5.setSpacing(6)
        self.gridLayout_5.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_5.setObjectName(u"gridLayout_5")
        self.gridLayout_22 = QGridLayout()
        self.gridLayout_22.setSpacing(6)
        self.gridLayout_22.setObjectName(u"gridLayout_22")
        self.gridLayout_22.setVerticalSpacing(8)
        self.gridLayout_22.setContentsMargins(6, -1, -1, -1)
        self.label_41 = QLabel(self.SearchgroupBox)
        self.label_41.setObjectName(u"label_41")

        self.gridLayout_22.addWidget(self.label_41, 1, 0, 1, 1)

        self.CuttypecomboBox = QComboBox(self.SearchgroupBox)
        self.CuttypecomboBox.addItem("")
        self.CuttypecomboBox.addItem("")
        self.CuttypecomboBox.setObjectName(u"CuttypecomboBox")

        self.gridLayout_22.addWidget(self.CuttypecomboBox, 1, 2, 1, 1)

        self.GrowthcheckBox = QCheckBox(self.SearchgroupBox)
        self.GrowthcheckBox.setObjectName(u"GrowthcheckBox")
        sizePolicy7 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Preferred)
        sizePolicy7.setHorizontalStretch(0)
        sizePolicy7.setVerticalStretch(0)
        sizePolicy7.setHeightForWidth(self.GrowthcheckBox.sizePolicy().hasHeightForWidth())
        self.GrowthcheckBox.setSizePolicy(sizePolicy7)
        self.GrowthcheckBox.setChecked(True)

        self.gridLayout_22.addWidget(self.GrowthcheckBox, 2, 0, 1, 1)

        self.horizontalSpacer_8 = QSpacerItem(40, 20, QSizePolicy.Fixed, QSizePolicy.Minimum)

        self.gridLayout_22.addItem(self.horizontalSpacer_8, 1, 1, 1, 1)

        self.label_42 = QLabel(self.SearchgroupBox)
        self.label_42.setObjectName(u"label_42")
        sizePolicy8 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
        sizePolicy8.setHorizontalStretch(0)
        sizePolicy8.setVerticalStretch(0)
        sizePolicy8.setHeightForWidth(self.label_42.sizePolicy().hasHeightForWidth())
        self.label_42.setSizePolicy(sizePolicy8)

        self.gridLayout_22.addWidget(self.label_42, 0, 0, 1, 1)

        self.SearchtypecomboBox = QComboBox(self.SearchgroupBox)
        self.SearchtypecomboBox.addItem("")
        self.SearchtypecomboBox.addItem("")
        self.SearchtypecomboBox.setObjectName(u"SearchtypecomboBox")

        self.gridLayout_22.addWidget(self.SearchtypecomboBox, 0, 2, 1, 1)

        self.SearchtoolButton = QToolButton(self.SearchgroupBox)
        self.SearchtoolButton.setObjectName(u"SearchtoolButton")
        self.SearchtoolButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.SearchtoolButton.setArrowType(Qt.RightArrow)

        self.gridLayout_22.addWidget(self.SearchtoolButton, 4, 0, 1, 5, Qt.AlignRight)

        self.primaryCutlineEdit = QLineEdit(self.SearchgroupBox)
        self.primaryCutlineEdit.setObjectName(u"primaryCutlineEdit")
        sizePolicy9 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.Fixed)
        sizePolicy9.setHorizontalStretch(0)
        sizePolicy9.setVerticalStretch(0)
        sizePolicy9.setHeightForWidth(self.primaryCutlineEdit.sizePolicy().hasHeightForWidth())
        self.primaryCutlineEdit.setSizePolicy(sizePolicy9)
        self.primaryCutlineEdit.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.gridLayout_22.addWidget(self.primaryCutlineEdit, 1, 3, 1, 1)

        self.Cuttype2comboBox = QComboBox(self.SearchgroupBox)
        self.Cuttype2comboBox.addItem("")
        self.Cuttype2comboBox.addItem("")
        self.Cuttype2comboBox.setObjectName(u"Cuttype2comboBox")

        self.gridLayout_22.addWidget(self.Cuttype2comboBox, 2, 2, 1, 1)

        self.SecondarycutlineEdit = QLineEdit(self.SearchgroupBox)
        self.SecondarycutlineEdit.setObjectName(u"SecondarycutlineEdit")
        sizePolicy9.setHeightForWidth(self.SecondarycutlineEdit.sizePolicy().hasHeightForWidth())
        self.SecondarycutlineEdit.setSizePolicy(sizePolicy9)
        self.SecondarycutlineEdit.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.gridLayout_22.addWidget(self.SecondarycutlineEdit, 2, 3, 1, 1)

        self.horizontalSpacer_20 = QSpacerItem(40, 10, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_22.addItem(self.horizontalSpacer_20, 3, 4, 1, 1)


        self.gridLayout_5.addLayout(self.gridLayout_22, 6, 0, 1, 1)


        self.gridLayout_10.addWidget(self.SearchgroupBox, 0, 0, 1, 1, Qt.AlignTop)

        self.stackedWidget.addWidget(self.page_4)
        self.page_5 = QWidget()
        self.page_5.setObjectName(u"page_5")
        self.gridLayout_11 = QGridLayout(self.page_5)
        self.gridLayout_11.setSpacing(6)
        self.gridLayout_11.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_11.setObjectName(u"gridLayout_11")
        self.SmoothgroupBox = QGroupBox(self.page_5)
        self.SmoothgroupBox.setObjectName(u"SmoothgroupBox")
        self.SmoothgroupBox.setEnabled(False)
        self.SmoothgroupBox.setCheckable(False)
        self.SmoothgroupBox.setChecked(False)
        self.gridLayout_6 = QGridLayout(self.SmoothgroupBox)
        self.gridLayout_6.setSpacing(6)
        self.gridLayout_6.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_6.setObjectName(u"gridLayout_6")
        self.SpatialSmoothgroupBox = QGroupBox(self.SmoothgroupBox)
        self.SpatialSmoothgroupBox.setObjectName(u"SpatialSmoothgroupBox")
        self.SpatialSmoothgroupBox.setCheckable(True)
        self.SpatialSmoothgroupBox.setChecked(False)
        self.gridLayout_24 = QGridLayout(self.SpatialSmoothgroupBox)
        self.gridLayout_24.setSpacing(6)
        self.gridLayout_24.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_24.setObjectName(u"gridLayout_24")
        self.ReducecheckBox = QCheckBox(self.SpatialSmoothgroupBox)
        self.ReducecheckBox.setObjectName(u"ReducecheckBox")

        self.gridLayout_24.addWidget(self.ReducecheckBox, 2, 0, 1, 1)

        self.SmoothOutlineEdit = QLineEdit(self.SpatialSmoothgroupBox)
        self.SmoothOutlineEdit.setObjectName(u"SmoothOutlineEdit")

        self.gridLayout_24.addWidget(self.SmoothOutlineEdit, 4, 1, 1, 1)

        self.FactorcheckBox = QCheckBox(self.SpatialSmoothgroupBox)
        self.FactorcheckBox.setObjectName(u"FactorcheckBox")
        self.FactorcheckBox.setChecked(True)

        self.gridLayout_24.addWidget(self.FactorcheckBox, 1, 0, 1, 1)

        self.label_39 = QLabel(self.SpatialSmoothgroupBox)
        self.label_39.setObjectName(u"label_39")

        self.gridLayout_24.addWidget(self.label_39, 4, 0, 1, 1)

        self.FFTcheckBox = QCheckBox(self.SpatialSmoothgroupBox)
        self.FFTcheckBox.setObjectName(u"FFTcheckBox")
        self.FFTcheckBox.setChecked(True)

        self.gridLayout_24.addWidget(self.FFTcheckBox, 3, 0, 1, 1)

        self.FactordoubleSpinBox = QDoubleSpinBox(self.SpatialSmoothgroupBox)
        self.FactordoubleSpinBox.setObjectName(u"FactordoubleSpinBox")
        sizePolicy5.setHeightForWidth(self.FactordoubleSpinBox.sizePolicy().hasHeightForWidth())
        self.FactordoubleSpinBox.setSizePolicy(sizePolicy5)
        self.FactordoubleSpinBox.setValue(2.000000000000000)

        self.gridLayout_24.addWidget(self.FactordoubleSpinBox, 1, 1, 1, 1)

        self.SmoothOutpushButton = QPushButton(self.SpatialSmoothgroupBox)
        self.SmoothOutpushButton.setObjectName(u"SmoothOutpushButton")

        self.gridLayout_24.addWidget(self.SmoothOutpushButton, 4, 2, 1, 1)

        self.BeamgroupBox = QGroupBox(self.SpatialSmoothgroupBox)
        self.BeamgroupBox.setObjectName(u"BeamgroupBox")
        sizePolicy1.setHeightForWidth(self.BeamgroupBox.sizePolicy().hasHeightForWidth())
        self.BeamgroupBox.setSizePolicy(sizePolicy1)
        self.BeamgroupBox.setCheckable(True)
        self.BeamgroupBox.setChecked(False)
        self.gridLayout_8 = QGridLayout(self.BeamgroupBox)
        self.gridLayout_8.setSpacing(6)
        self.gridLayout_8.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_8.setObjectName(u"gridLayout_8")
        self.label_7 = QLabel(self.BeamgroupBox)
        self.label_7.setObjectName(u"label_7")

        self.gridLayout_8.addWidget(self.label_7, 3, 0, 1, 1)

        self.label_2 = QLabel(self.BeamgroupBox)
        self.label_2.setObjectName(u"label_2")
        sizePolicy8.setHeightForWidth(self.label_2.sizePolicy().hasHeightForWidth())
        self.label_2.setSizePolicy(sizePolicy8)
        self.label_2.setAlignment(Qt.AlignCenter)

        self.gridLayout_8.addWidget(self.label_2, 0, 4, 1, 1)

        self.label_4 = QLabel(self.BeamgroupBox)
        self.label_4.setObjectName(u"label_4")
        sizePolicy4.setHeightForWidth(self.label_4.sizePolicy().hasHeightForWidth())
        self.label_4.setSizePolicy(sizePolicy4)
        self.label_4.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_8.addWidget(self.label_4, 1, 0, 1, 1)

        self.label_6 = QLabel(self.BeamgroupBox)
        self.label_6.setObjectName(u"label_6")

        self.gridLayout_8.addWidget(self.label_6, 2, 0, 1, 1)

        self.ObmajSpinBox = QDoubleSpinBox(self.BeamgroupBox)
        self.ObmajSpinBox.setObjectName(u"ObmajSpinBox")
        self.ObmajSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.ObmajSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ObmajSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.ObmajSpinBox.setMinimum(-1.000000000000000)
        self.ObmajSpinBox.setMaximum(10000.000000000000000)
        self.ObmajSpinBox.setValue(-1.000000000000000)

        self.gridLayout_8.addWidget(self.ObmajSpinBox, 1, 2, 1, 1)

        self.NbmajSpinBox = QDoubleSpinBox(self.BeamgroupBox)
        self.NbmajSpinBox.setObjectName(u"NbmajSpinBox")
        self.NbmajSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.NbmajSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.NbmajSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.NbmajSpinBox.setMinimum(-1.000000000000000)
        self.NbmajSpinBox.setMaximum(10000.000000000000000)
        self.NbmajSpinBox.setValue(-1.000000000000000)

        self.gridLayout_8.addWidget(self.NbmajSpinBox, 1, 4, 1, 1)

        self.label_11 = QLabel(self.BeamgroupBox)
        self.label_11.setObjectName(u"label_11")
        sizePolicy10 = QSizePolicy(QSizePolicy.Maximum, QSizePolicy.Preferred)
        sizePolicy10.setHorizontalStretch(0)
        sizePolicy10.setVerticalStretch(0)
        sizePolicy10.setHeightForWidth(self.label_11.sizePolicy().hasHeightForWidth())
        self.label_11.setSizePolicy(sizePolicy10)

        self.gridLayout_8.addWidget(self.label_11, 1, 5, 1, 1)

        self.ObminSpinBox = QDoubleSpinBox(self.BeamgroupBox)
        self.ObminSpinBox.setObjectName(u"ObminSpinBox")
        sizePolicy11 = QSizePolicy(QSizePolicy.Ignored, QSizePolicy.Preferred)
        sizePolicy11.setHorizontalStretch(0)
        sizePolicy11.setVerticalStretch(0)
        sizePolicy11.setHeightForWidth(self.ObminSpinBox.sizePolicy().hasHeightForWidth())
        self.ObminSpinBox.setSizePolicy(sizePolicy11)
        self.ObminSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.ObminSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ObminSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.ObminSpinBox.setMinimum(-1.000000000000000)
        self.ObminSpinBox.setMaximum(10000.000000000000000)
        self.ObminSpinBox.setValue(-1.000000000000000)

        self.gridLayout_8.addWidget(self.ObminSpinBox, 2, 2, 1, 1)

        self.ObpaSpinBox = QDoubleSpinBox(self.BeamgroupBox)
        self.ObpaSpinBox.setObjectName(u"ObpaSpinBox")
        self.ObpaSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.ObpaSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ObpaSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.ObpaSpinBox.setMinimum(-10000.000000000000000)
        self.ObpaSpinBox.setMaximum(10000.000000000000000)
        self.ObpaSpinBox.setValue(0.000000000000000)

        self.gridLayout_8.addWidget(self.ObpaSpinBox, 3, 2, 1, 1)

        self.label_8 = QLabel(self.BeamgroupBox)
        self.label_8.setObjectName(u"label_8")

        self.gridLayout_8.addWidget(self.label_8, 2, 3, 1, 1)

        self.label_9 = QLabel(self.BeamgroupBox)
        self.label_9.setObjectName(u"label_9")

        self.gridLayout_8.addWidget(self.label_9, 3, 3, 1, 1)

        self.NbminSpinBox = QDoubleSpinBox(self.BeamgroupBox)
        self.NbminSpinBox.setObjectName(u"NbminSpinBox")
        self.NbminSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.NbminSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.NbminSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.NbminSpinBox.setMinimum(-1.000000000000000)
        self.NbminSpinBox.setMaximum(10000.000000000000000)
        self.NbminSpinBox.setValue(-1.000000000000000)

        self.gridLayout_8.addWidget(self.NbminSpinBox, 2, 4, 1, 1)

        self.NbpaSpinBox = QDoubleSpinBox(self.BeamgroupBox)
        self.NbpaSpinBox.setObjectName(u"NbpaSpinBox")
        self.NbpaSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.NbpaSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.NbpaSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.NbpaSpinBox.setMinimum(-10000.000000000000000)
        self.NbpaSpinBox.setMaximum(10000.000000000000000)
        self.NbpaSpinBox.setValue(0.000000000000000)

        self.gridLayout_8.addWidget(self.NbpaSpinBox, 3, 4, 1, 1)

        self.label_3 = QLabel(self.BeamgroupBox)
        self.label_3.setObjectName(u"label_3")
        sizePolicy10.setHeightForWidth(self.label_3.sizePolicy().hasHeightForWidth())
        self.label_3.setSizePolicy(sizePolicy10)

        self.gridLayout_8.addWidget(self.label_3, 3, 5, 1, 1)

        self.label_10 = QLabel(self.BeamgroupBox)
        self.label_10.setObjectName(u"label_10")
        sizePolicy10.setHeightForWidth(self.label_10.sizePolicy().hasHeightForWidth())
        self.label_10.setSizePolicy(sizePolicy10)

        self.gridLayout_8.addWidget(self.label_10, 2, 5, 1, 1)

        self.label = QLabel(self.BeamgroupBox)
        self.label.setObjectName(u"label")
        sizePolicy12 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Preferred)
        sizePolicy12.setHorizontalStretch(0)
        sizePolicy12.setVerticalStretch(0)
        sizePolicy12.setHeightForWidth(self.label.sizePolicy().hasHeightForWidth())
        self.label.setSizePolicy(sizePolicy12)
        self.label.setAlignment(Qt.AlignCenter)

        self.gridLayout_8.addWidget(self.label, 0, 2, 1, 1)

        self.label_5 = QLabel(self.BeamgroupBox)
        self.label_5.setObjectName(u"label_5")
        sizePolicy4.setHeightForWidth(self.label_5.sizePolicy().hasHeightForWidth())
        self.label_5.setSizePolicy(sizePolicy4)

        self.gridLayout_8.addWidget(self.label_5, 1, 3, 1, 1)

        self.horizontalSpacer_15 = QSpacerItem(40, 20, QSizePolicy.Fixed, QSizePolicy.Minimum)

        self.gridLayout_8.addItem(self.horizontalSpacer_15, 1, 1, 1, 1)


        self.gridLayout_24.addWidget(self.BeamgroupBox, 0, 0, 1, 3)


        self.gridLayout_6.addWidget(self.SpatialSmoothgroupBox, 1, 0, 1, 3)

        self.horizontalSpacer_21 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_6.addItem(self.horizontalSpacer_21, 2, 2, 1, 1)

        self.HanninggroupBox = QGroupBox(self.SmoothgroupBox)
        self.HanninggroupBox.setObjectName(u"HanninggroupBox")
        self.HanninggroupBox.setCheckable(True)
        self.HanninggroupBox.setChecked(False)
        self.horizontalLayout_2 = QHBoxLayout(self.HanninggroupBox)
        self.horizontalLayout_2.setSpacing(6)
        self.horizontalLayout_2.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_2.setObjectName(u"horizontalLayout_2")
        self.label_40 = QLabel(self.HanninggroupBox)
        self.label_40.setObjectName(u"label_40")

        self.horizontalLayout_2.addWidget(self.label_40)

        self.HanningspinBox = QSpinBox(self.HanninggroupBox)
        self.HanningspinBox.setObjectName(u"HanningspinBox")
        self.HanningspinBox.setMinimumSize(QSize(90, 0))
        self.HanningspinBox.setLayoutDirection(Qt.RightToLeft)
        self.HanningspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.HanningspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.HanningspinBox.setMinimum(1)
        self.HanningspinBox.setSingleStep(2)

        self.horizontalLayout_2.addWidget(self.HanningspinBox)

        self.horizontalSpacer_22 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.horizontalLayout_2.addItem(self.horizontalSpacer_22)


        self.gridLayout_6.addWidget(self.HanninggroupBox, 3, 0, 1, 3)


        self.gridLayout_11.addWidget(self.SmoothgroupBox, 1, 0, 1, 1)

        self.stackedWidget.addWidget(self.page_5)
        self.page_6 = QWidget()
        self.page_6.setObjectName(u"page_6")
        self.gridLayout_13 = QGridLayout(self.page_6)
        self.gridLayout_13.setSpacing(6)
        self.gridLayout_13.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_13.setObjectName(u"gridLayout_13")
        self.MaskgroupBox = QGroupBox(self.page_6)
        self.MaskgroupBox.setObjectName(u"MaskgroupBox")
        self.MaskgroupBox.setEnabled(False)
        self.MaskgroupBox.setCheckable(True)
        self.MaskgroupBox.setChecked(False)
        self.gridLayout_9 = QGridLayout(self.MaskgroupBox)
        self.gridLayout_9.setSpacing(6)
        self.gridLayout_9.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_9.setObjectName(u"gridLayout_9")
        self.BlankcutSpinBox = QDoubleSpinBox(self.MaskgroupBox)
        self.BlankcutSpinBox.setObjectName(u"BlankcutSpinBox")
        self.BlankcutSpinBox.setLocale(QLocale(QLocale.English, QLocale.UnitedStates))
        self.BlankcutSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.BlankcutSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.BlankcutSpinBox.setMaximum(100000.000000000000000)
        self.BlankcutSpinBox.setValue(3.000000000000000)

        self.gridLayout_9.addWidget(self.BlankcutSpinBox, 2, 2, 1, 1)

        self.label_54 = QLabel(self.MaskgroupBox)
        self.label_54.setObjectName(u"label_54")

        self.gridLayout_9.addWidget(self.label_54, 2, 0, 1, 1)

        self.label_50 = QLabel(self.MaskgroupBox)
        self.label_50.setObjectName(u"label_50")

        self.gridLayout_9.addWidget(self.label_50, 1, 0, 1, 1)

        self.BlankfactordoubleSpinBox = QDoubleSpinBox(self.MaskgroupBox)
        self.BlankfactordoubleSpinBox.setObjectName(u"BlankfactordoubleSpinBox")
        self.BlankfactordoubleSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.BlankfactordoubleSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.BlankfactordoubleSpinBox.setValue(2.000000000000000)

        self.gridLayout_9.addWidget(self.BlankfactordoubleSpinBox, 1, 2, 1, 1)

        self.horizontalSpacer_12 = QSpacerItem(30, 20, QSizePolicy.Fixed, QSizePolicy.Minimum)

        self.gridLayout_9.addItem(self.horizontalSpacer_12, 1, 1, 1, 1)

        self.horizontalSpacer_11 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_9.addItem(self.horizontalSpacer_11, 1, 3, 1, 1)


        self.gridLayout_13.addWidget(self.MaskgroupBox, 6, 0, 1, 1)

        self.MomentmapsgroupBox = QGroupBox(self.page_6)
        self.MomentmapsgroupBox.setObjectName(u"MomentmapsgroupBox")
        self.MomentmapsgroupBox.setEnabled(False)
        self.MomentmapsgroupBox.setCheckable(True)
        self.MomentmapsgroupBox.setChecked(False)
        self.gridLayout_18 = QGridLayout(self.MomentmapsgroupBox)
        self.gridLayout_18.setSpacing(6)
        self.gridLayout_18.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_18.setObjectName(u"gridLayout_18")
        self.horizontalSpacer_10 = QSpacerItem(30, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_18.addItem(self.horizontalSpacer_10, 1, 3, 1, 2)

        self.label_35 = QLabel(self.MomentmapsgroupBox)
        self.label_35.setObjectName(u"label_35")

        self.gridLayout_18.addWidget(self.label_35, 3, 3, 3, 1)

        self.label_44 = QLabel(self.MomentmapsgroupBox)
        self.label_44.setObjectName(u"label_44")

        self.gridLayout_18.addWidget(self.label_44, 3, 0, 1, 1)

        self.TotalmapcheckBox = QCheckBox(self.MomentmapsgroupBox)
        self.TotalmapcheckBox.setObjectName(u"TotalmapcheckBox")
        self.TotalmapcheckBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_18.addWidget(self.TotalmapcheckBox, 3, 2, 1, 1)

        self.label_46 = QLabel(self.MomentmapsgroupBox)
        self.label_46.setObjectName(u"label_46")
        sizePolicy12.setHeightForWidth(self.label_46.sizePolicy().hasHeightForWidth())
        self.label_46.setSizePolicy(sizePolicy12)

        self.gridLayout_18.addWidget(self.label_46, 5, 0, 1, 1)

        self.MaptypecomboBox = QComboBox(self.MomentmapsgroupBox)
        self.MaptypecomboBox.addItem("")
        self.MaptypecomboBox.addItem("")
        self.MaptypecomboBox.setObjectName(u"MaptypecomboBox")

        self.gridLayout_18.addWidget(self.MaptypecomboBox, 3, 4, 3, 1)

        self.VfieldcheckBox = QCheckBox(self.MomentmapsgroupBox)
        self.VfieldcheckBox.setObjectName(u"VfieldcheckBox")
        self.VfieldcheckBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_18.addWidget(self.VfieldcheckBox, 4, 2, 1, 1)

        self.DispmapcheckBox = QCheckBox(self.MomentmapsgroupBox)
        self.DispmapcheckBox.setObjectName(u"DispmapcheckBox")
        self.DispmapcheckBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_18.addWidget(self.DispmapcheckBox, 5, 2, 1, 1)

        self.horizontalSpacer_13 = QSpacerItem(30, 20, QSizePolicy.Fixed, QSizePolicy.Minimum)

        self.gridLayout_18.addItem(self.horizontalSpacer_13, 4, 1, 1, 1)

        self.ProfilecheckBox = QCheckBox(self.MomentmapsgroupBox)
        self.ProfilecheckBox.setObjectName(u"ProfilecheckBox")
        self.ProfilecheckBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_18.addWidget(self.ProfilecheckBox, 1, 2, 1, 1)

        self.label_45 = QLabel(self.MomentmapsgroupBox)
        self.label_45.setObjectName(u"label_45")

        self.gridLayout_18.addWidget(self.label_45, 4, 0, 1, 1)

        self.label_47 = QLabel(self.MomentmapsgroupBox)
        self.label_47.setObjectName(u"label_47")

        self.gridLayout_18.addWidget(self.label_47, 1, 0, 1, 1)

        self.label_38 = QLabel(self.MomentmapsgroupBox)
        self.label_38.setObjectName(u"label_38")

        self.gridLayout_18.addWidget(self.label_38, 2, 0, 1, 1)

        self.rmsmapcheckBox = QCheckBox(self.MomentmapsgroupBox)
        self.rmsmapcheckBox.setObjectName(u"rmsmapcheckBox")
        self.rmsmapcheckBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_18.addWidget(self.rmsmapcheckBox, 2, 2, 1, 1)


        self.gridLayout_13.addWidget(self.MomentmapsgroupBox, 2, 0, 1, 1)

        self.horizontalSpacer_14 = QSpacerItem(20, 6, QSizePolicy.Minimum, QSizePolicy.Fixed)

        self.gridLayout_13.addItem(self.horizontalSpacer_14, 5, 0, 1, 1)

        self.PVgroupBox = QGroupBox(self.page_6)
        self.PVgroupBox.setObjectName(u"PVgroupBox")
        self.PVgroupBox.setEnabled(False)
        self.PVgroupBox.setCheckable(True)
        self.PVgroupBox.setChecked(False)
        self.gridLayout_7 = QGridLayout(self.PVgroupBox)
        self.gridLayout_7.setSpacing(6)
        self.gridLayout_7.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_7.setObjectName(u"gridLayout_7")
        self.PVYposSpinBox = QDoubleSpinBox(self.PVgroupBox)
        self.PVYposSpinBox.setObjectName(u"PVYposSpinBox")
        sizePolicy5.setHeightForWidth(self.PVYposSpinBox.sizePolicy().hasHeightForWidth())
        self.PVYposSpinBox.setSizePolicy(sizePolicy5)
        self.PVYposSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.PVYposSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.PVYposSpinBox.setMaximum(9999.000000000000000)

        self.gridLayout_7.addWidget(self.PVYposSpinBox, 0, 3, 1, 1)

        self.horizontalSpacer_19 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_7.addItem(self.horizontalSpacer_19, 0, 5, 1, 1)

        self.label_20 = QLabel(self.PVgroupBox)
        self.label_20.setObjectName(u"label_20")
        sizePolicy4.setHeightForWidth(self.label_20.sizePolicy().hasHeightForWidth())
        self.label_20.setSizePolicy(sizePolicy4)

        self.gridLayout_7.addWidget(self.label_20, 0, 0, 1, 1)

        self.label_22 = QLabel(self.PVgroupBox)
        self.label_22.setObjectName(u"label_22")
        sizePolicy4.setHeightForWidth(self.label_22.sizePolicy().hasHeightForWidth())
        self.label_22.setSizePolicy(sizePolicy4)

        self.gridLayout_7.addWidget(self.label_22, 0, 4, 1, 1)

        self.PVXposSpinBox = QDoubleSpinBox(self.PVgroupBox)
        self.PVXposSpinBox.setObjectName(u"PVXposSpinBox")
        sizePolicy5.setHeightForWidth(self.PVXposSpinBox.sizePolicy().hasHeightForWidth())
        self.PVXposSpinBox.setSizePolicy(sizePolicy5)
        self.PVXposSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.PVXposSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.PVXposSpinBox.setMaximum(9999.000000000000000)

        self.gridLayout_7.addWidget(self.PVXposSpinBox, 0, 1, 1, 1)

        self.label_21 = QLabel(self.PVgroupBox)
        self.label_21.setObjectName(u"label_21")
        sizePolicy10.setHeightForWidth(self.label_21.sizePolicy().hasHeightForWidth())
        self.label_21.setSizePolicy(sizePolicy10)

        self.gridLayout_7.addWidget(self.label_21, 0, 2, 1, 1)

        self.label_27 = QLabel(self.PVgroupBox)
        self.label_27.setObjectName(u"label_27")

        self.gridLayout_7.addWidget(self.label_27, 1, 0, 1, 1)

        self.PVPaSpinBox = QDoubleSpinBox(self.PVgroupBox)
        self.PVPaSpinBox.setObjectName(u"PVPaSpinBox")
        self.PVPaSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.PVPaSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.PVPaSpinBox.setMaximum(360.000000000000000)

        self.gridLayout_7.addWidget(self.PVPaSpinBox, 1, 1, 1, 1)


        self.gridLayout_13.addWidget(self.PVgroupBox, 4, 0, 1, 1)

        self.verticalSpacer_2 = QSpacerItem(20, 6, QSizePolicy.Minimum, QSizePolicy.Fixed)

        self.gridLayout_13.addItem(self.verticalSpacer_2, 3, 0, 1, 1)

        self.stackedWidget.addWidget(self.page_6)
        self.page_7 = QWidget()
        self.page_7.setObjectName(u"page_7")
        self.gridLayout_16 = QGridLayout(self.page_7)
        self.gridLayout_16.setSpacing(6)
        self.gridLayout_16.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_16.setObjectName(u"gridLayout_16")
        self.LogtextEdit = QTextEdit(self.page_7)
        self.LogtextEdit.setObjectName(u"LogtextEdit")
        self.LogtextEdit.setEnabled(True)
        palette1 = QPalette()
        brush2 = QBrush(QColor(243, 243, 243, 255))
        brush2.setStyle(Qt.SolidPattern)
        palette1.setBrush(QPalette.Active, QPalette.Base, brush2)
        palette1.setBrush(QPalette.Inactive, QPalette.Base, brush2)
        brush3 = QBrush(QColor(241, 243, 240, 255))
        brush3.setStyle(Qt.SolidPattern)
        palette1.setBrush(QPalette.Disabled, QPalette.Base, brush3)
        self.LogtextEdit.setPalette(palette1)
        font1 = QFont()
        font1.setFamilies([u"Monaco"])
        font1.setPointSize(12)
        self.LogtextEdit.setFont(font1)
        self.LogtextEdit.setReadOnly(True)

        self.gridLayout_16.addWidget(self.LogtextEdit, 0, 0, 1, 1)

        self.stackedWidget.addWidget(self.page_7)
        self.page_3 = QWidget()
        self.page_3.setObjectName(u"page_3")
        self.gridLayout_17 = QGridLayout(self.page_3)
        self.gridLayout_17.setSpacing(6)
        self.gridLayout_17.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_17.setObjectName(u"gridLayout_17")
        self.Plot2graphicsView = QGraphicsView(self.page_3)
        self.Plot2graphicsView.setObjectName(u"Plot2graphicsView")

        self.gridLayout_17.addWidget(self.Plot2graphicsView, 1, 1, 1, 1)

        self.Plot3comboBox = QComboBox(self.page_3)
        self.Plot3comboBox.addItem("")
        self.Plot3comboBox.addItem("")
        self.Plot3comboBox.addItem("")
        self.Plot3comboBox.addItem("")
        self.Plot3comboBox.addItem("")
        self.Plot3comboBox.addItem("")
        self.Plot3comboBox.addItem("")
        self.Plot3comboBox.addItem("")
        self.Plot3comboBox.setObjectName(u"Plot3comboBox")
        sizePolicy6.setHeightForWidth(self.Plot3comboBox.sizePolicy().hasHeightForWidth())
        self.Plot3comboBox.setSizePolicy(sizePolicy6)
        self.Plot3comboBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_17.addWidget(self.Plot3comboBox, 2, 0, 1, 1)

        self.Plot2comboBox = QComboBox(self.page_3)
        self.Plot2comboBox.addItem("")
        self.Plot2comboBox.addItem("")
        self.Plot2comboBox.addItem("")
        self.Plot2comboBox.addItem("")
        self.Plot2comboBox.addItem("")
        self.Plot2comboBox.addItem("")
        self.Plot2comboBox.addItem("")
        self.Plot2comboBox.addItem("")
        self.Plot2comboBox.setObjectName(u"Plot2comboBox")
        sizePolicy6.setHeightForWidth(self.Plot2comboBox.sizePolicy().hasHeightForWidth())
        self.Plot2comboBox.setSizePolicy(sizePolicy6)
        self.Plot2comboBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_17.addWidget(self.Plot2comboBox, 0, 1, 1, 1)

        self.Plot1graphicsView = QGraphicsView(self.page_3)
        self.Plot1graphicsView.setObjectName(u"Plot1graphicsView")

        self.gridLayout_17.addWidget(self.Plot1graphicsView, 1, 0, 1, 1)

        self.Plot1comboBox = QComboBox(self.page_3)
        self.Plot1comboBox.addItem("")
        self.Plot1comboBox.addItem("")
        self.Plot1comboBox.addItem("")
        self.Plot1comboBox.addItem("")
        self.Plot1comboBox.addItem("")
        self.Plot1comboBox.addItem("")
        self.Plot1comboBox.addItem("")
        self.Plot1comboBox.addItem("")
        self.Plot1comboBox.setObjectName(u"Plot1comboBox")
        sizePolicy6.setHeightForWidth(self.Plot1comboBox.sizePolicy().hasHeightForWidth())
        self.Plot1comboBox.setSizePolicy(sizePolicy6)
        self.Plot1comboBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_17.addWidget(self.Plot1comboBox, 0, 0, 1, 1)

        self.Plot3graphicsView = QGraphicsView(self.page_3)
        self.Plot3graphicsView.setObjectName(u"Plot3graphicsView")

        self.gridLayout_17.addWidget(self.Plot3graphicsView, 3, 0, 1, 1)

        self.Plot4graphicsView = QGraphicsView(self.page_3)
        self.Plot4graphicsView.setObjectName(u"Plot4graphicsView")

        self.gridLayout_17.addWidget(self.Plot4graphicsView, 3, 1, 1, 1)

        self.Plot4comboBox = QComboBox(self.page_3)
        self.Plot4comboBox.addItem("")
        self.Plot4comboBox.addItem("")
        self.Plot4comboBox.addItem("")
        self.Plot4comboBox.addItem("")
        self.Plot4comboBox.addItem("")
        self.Plot4comboBox.addItem("")
        self.Plot4comboBox.addItem("")
        self.Plot4comboBox.addItem("")
        self.Plot4comboBox.setObjectName(u"Plot4comboBox")
        sizePolicy6.setHeightForWidth(self.Plot4comboBox.sizePolicy().hasHeightForWidth())
        self.Plot4comboBox.setSizePolicy(sizePolicy6)
        self.Plot4comboBox.setLayoutDirection(Qt.RightToLeft)

        self.gridLayout_17.addWidget(self.Plot4comboBox, 2, 1, 1, 1)

        self.stackedWidget.addWidget(self.page_3)
        self.page_2 = QWidget()
        self.page_2.setObjectName(u"page_2")
        self.gridLayout_2 = QGridLayout(self.page_2)
        self.gridLayout_2.setSpacing(6)
        self.gridLayout_2.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_2.setObjectName(u"gridLayout_2")
        self.AdvancedgroupBox = QGroupBox(self.page_2)
        self.AdvancedgroupBox.setObjectName(u"AdvancedgroupBox")
        self.AdvancedgroupBox.setEnabled(False)
        self.AdvancedgroupBox.setFlat(True)
        self.AdvancedgroupBox.setCheckable(True)
        self.AdvancedgroupBox.setChecked(False)
        self.gridLayout_12 = QGridLayout(self.AdvancedgroupBox)
        self.gridLayout_12.setSpacing(6)
        self.gridLayout_12.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_12.setObjectName(u"gridLayout_12")
        self.horizontalSpacer_17 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_12.addItem(self.horizontalSpacer_17, 5, 0, 1, 5)

        self.WfunccomboBox = QComboBox(self.AdvancedgroupBox)
        self.WfunccomboBox.addItem("")
        self.WfunccomboBox.addItem("")
        self.WfunccomboBox.addItem("")
        self.WfunccomboBox.setObjectName(u"WfunccomboBox")

        self.gridLayout_12.addWidget(self.WfunccomboBox, 1, 4, 1, 1)

        self.FtypecomboBox = QComboBox(self.AdvancedgroupBox)
        self.FtypecomboBox.addItem("")
        self.FtypecomboBox.addItem("")
        self.FtypecomboBox.addItem("")
        self.FtypecomboBox.setObjectName(u"FtypecomboBox")

        self.gridLayout_12.addWidget(self.FtypecomboBox, 0, 4, 1, 1)

        self.PolynspinBox = QSpinBox(self.AdvancedgroupBox)
        self.PolynspinBox.setObjectName(u"PolynspinBox")
        sizePolicy8.setHeightForWidth(self.PolynspinBox.sizePolicy().hasHeightForWidth())
        self.PolynspinBox.setSizePolicy(sizePolicy8)
        self.PolynspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.PolynspinBox.setMinimum(-1)
        self.PolynspinBox.setMaximum(10)
        self.PolynspinBox.setValue(-1)

        self.gridLayout_12.addWidget(self.PolynspinBox, 2, 4, 1, 1)

        self.SecondstagecheckBox = QCheckBox(self.AdvancedgroupBox)
        self.SecondstagecheckBox.setObjectName(u"SecondstagecheckBox")
        sizePolicy13 = QSizePolicy(QSizePolicy.Minimum, QSizePolicy.MinimumExpanding)
        sizePolicy13.setHorizontalStretch(0)
        sizePolicy13.setVerticalStretch(0)
        sizePolicy13.setHeightForWidth(self.SecondstagecheckBox.sizePolicy().hasHeightForWidth())
        self.SecondstagecheckBox.setSizePolicy(sizePolicy13)
        self.SecondstagecheckBox.setChecked(True)

        self.gridLayout_12.addWidget(self.SecondstagecheckBox, 2, 3, 1, 1)

        self.TollineEdit = QLineEdit(self.AdvancedgroupBox)
        self.TollineEdit.setObjectName(u"TollineEdit")
        sizePolicy8.setHeightForWidth(self.TollineEdit.sizePolicy().hasHeightForWidth())
        self.TollineEdit.setSizePolicy(sizePolicy8)
        self.TollineEdit.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.gridLayout_12.addWidget(self.TollineEdit, 3, 4, 1, 1)

        self.Ftypelabel = QLabel(self.AdvancedgroupBox)
        self.Ftypelabel.setObjectName(u"Ftypelabel")
        sizePolicy14 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.MinimumExpanding)
        sizePolicy14.setHorizontalStretch(0)
        sizePolicy14.setVerticalStretch(0)
        sizePolicy14.setHeightForWidth(self.Ftypelabel.sizePolicy().hasHeightForWidth())
        self.Ftypelabel.setSizePolicy(sizePolicy14)

        self.gridLayout_12.addWidget(self.Ftypelabel, 0, 3, 1, 1)

        self.label_28 = QLabel(self.AdvancedgroupBox)
        self.label_28.setObjectName(u"label_28")

        self.gridLayout_12.addWidget(self.label_28, 6, 0, 1, 1)

        self.Tollabel = QLabel(self.AdvancedgroupBox)
        self.Tollabel.setObjectName(u"Tollabel")
        sizePolicy14.setHeightForWidth(self.Tollabel.sizePolicy().hasHeightForWidth())
        self.Tollabel.setSizePolicy(sizePolicy14)

        self.gridLayout_12.addWidget(self.Tollabel, 3, 3, 1, 1)

        self.MasktoolButton = QToolButton(self.AdvancedgroupBox)
        self.MasktoolButton.setObjectName(u"MasktoolButton")

        self.gridLayout_12.addWidget(self.MasktoolButton, 6, 3, 1, 1)

        self.Wfunclabel = QLabel(self.AdvancedgroupBox)
        self.Wfunclabel.setObjectName(u"Wfunclabel")
        sizePolicy14.setHeightForWidth(self.Wfunclabel.sizePolicy().hasHeightForWidth())
        self.Wfunclabel.setSizePolicy(sizePolicy14)

        self.gridLayout_12.addWidget(self.Wfunclabel, 1, 3, 1, 1)

        self.MaskThreshSpinBox = QDoubleSpinBox(self.AdvancedgroupBox)
        self.MaskThreshSpinBox.setObjectName(u"MaskThreshSpinBox")
        self.MaskThreshSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.MaskThreshSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.MaskThreshSpinBox.setDecimals(8)
        self.MaskThreshSpinBox.setMaximum(99999.000000000000000)

        self.gridLayout_12.addWidget(self.MaskThreshSpinBox, 6, 4, 1, 1)

        self.ErrorscheckBox = QCheckBox(self.AdvancedgroupBox)
        self.ErrorscheckBox.setObjectName(u"ErrorscheckBox")
        self.ErrorscheckBox.setLayoutDirection(Qt.LeftToRight)

        self.gridLayout_12.addWidget(self.ErrorscheckBox, 4, 0, 1, 1)

        self.label_37 = QLabel(self.AdvancedgroupBox)
        self.label_37.setObjectName(u"label_37")
        sizePolicy7.setHeightForWidth(self.label_37.sizePolicy().hasHeightForWidth())
        self.label_37.setSizePolicy(sizePolicy7)

        self.gridLayout_12.addWidget(self.label_37, 2, 0, 1, 1)

        self.label_33 = QLabel(self.AdvancedgroupBox)
        self.label_33.setObjectName(u"label_33")

        self.gridLayout_12.addWidget(self.label_33, 0, 0, 1, 1)

        self.label_34 = QLabel(self.AdvancedgroupBox)
        self.label_34.setObjectName(u"label_34")
        sizePolicy4.setHeightForWidth(self.label_34.sizePolicy().hasHeightForWidth())
        self.label_34.setSizePolicy(sizePolicy4)

        self.gridLayout_12.addWidget(self.label_34, 1, 0, 1, 1)

        self.label_36 = QLabel(self.AdvancedgroupBox)
        self.label_36.setObjectName(u"label_36")

        self.gridLayout_12.addWidget(self.label_36, 3, 0, 1, 1)

        self.line = QFrame(self.AdvancedgroupBox)
        self.line.setObjectName(u"line")
        self.line.setMinimumSize(QSize(20, 0))
        self.line.setFrameShape(QFrame.VLine)
        self.line.setFrameShadow(QFrame.Sunken)

        self.gridLayout_12.addWidget(self.line, 0, 2, 5, 1)

        self.CdensSpinBox = QDoubleSpinBox(self.AdvancedgroupBox)
        self.CdensSpinBox.setObjectName(u"CdensSpinBox")
        sizePolicy8.setHeightForWidth(self.CdensSpinBox.sizePolicy().hasHeightForWidth())
        self.CdensSpinBox.setSizePolicy(sizePolicy8)
        self.CdensSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.CdensSpinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.CdensSpinBox.setMaximum(9999.000000000000000)
        self.CdensSpinBox.setValue(10.000000000000000)

        self.gridLayout_12.addWidget(self.CdensSpinBox, 1, 1, 1, 1)

        self.NVspinBox = QSpinBox(self.AdvancedgroupBox)
        self.NVspinBox.setObjectName(u"NVspinBox")
        sizePolicy8.setHeightForWidth(self.NVspinBox.sizePolicy().hasHeightForWidth())
        self.NVspinBox.setSizePolicy(sizePolicy8)
        self.NVspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.NVspinBox.setButtonSymbols(QAbstractSpinBox.NoButtons)
        self.NVspinBox.setMinimum(-1)
        self.NVspinBox.setMaximum(9999)
        self.NVspinBox.setValue(-1)

        self.gridLayout_12.addWidget(self.NVspinBox, 2, 1, 1, 1)

        self.NormcomboBox = QComboBox(self.AdvancedgroupBox)
        self.NormcomboBox.addItem("")
        self.NormcomboBox.addItem("")
        self.NormcomboBox.addItem("")
        self.NormcomboBox.setObjectName(u"NormcomboBox")

        self.gridLayout_12.addWidget(self.NormcomboBox, 3, 1, 1, 1)

        self.MaskcomboBox = QComboBox(self.AdvancedgroupBox)
        self.MaskcomboBox.addItem("")
        self.MaskcomboBox.addItem("")
        self.MaskcomboBox.addItem("")
        self.MaskcomboBox.addItem("")
        self.MaskcomboBox.addItem("")
        self.MaskcomboBox.addItem("")
        self.MaskcomboBox.setObjectName(u"MaskcomboBox")
        sizePolicy8.setHeightForWidth(self.MaskcomboBox.sizePolicy().hasHeightForWidth())
        self.MaskcomboBox.setSizePolicy(sizePolicy8)

        self.gridLayout_12.addWidget(self.MaskcomboBox, 6, 1, 1, 1)

        self.LtypecomboBox = QComboBox(self.AdvancedgroupBox)
        self.LtypecomboBox.addItem("")
        self.LtypecomboBox.addItem("")
        self.LtypecomboBox.addItem("")
        self.LtypecomboBox.addItem("")
        self.LtypecomboBox.addItem("")
        self.LtypecomboBox.setObjectName(u"LtypecomboBox")

        self.gridLayout_12.addWidget(self.LtypecomboBox, 0, 1, 1, 1)

        self.AdriftcheckBox = QCheckBox(self.AdvancedgroupBox)
        self.AdriftcheckBox.setObjectName(u"AdriftcheckBox")
        self.AdriftcheckBox.setMinimumSize(QSize(0, 30))

        self.gridLayout_12.addWidget(self.AdriftcheckBox, 4, 1, 1, 1)


        self.gridLayout_2.addWidget(self.AdvancedgroupBox, 10, 0, 1, 1)

        self.restframe = QFrame(self.page_2)
        self.restframe.setObjectName(u"restframe")
        self.restframe.setEnabled(False)
        self.restframe.setFrameShape(QFrame.NoFrame)
        self.restframe.setFrameShadow(QFrame.Plain)
        self.restframe.setLineWidth(0)
        self.horizontalLayout_5 = QHBoxLayout(self.restframe)
        self.horizontalLayout_5.setSpacing(6)
        self.horizontalLayout_5.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_5.setObjectName(u"horizontalLayout_5")
        self.horizontalLayout_5.setContentsMargins(0, 0, 0, 0)

        self.gridLayout_2.addWidget(self.restframe, 3, 0, 1, 1)

        self.Restframe = QFrame(self.page_2)
        self.Restframe.setObjectName(u"Restframe")
        self.Restframe.setEnabled(True)
        self.Restframe.setFrameShape(QFrame.NoFrame)
        self.Restframe.setFrameShadow(QFrame.Plain)
        self.horizontalLayout_6 = QHBoxLayout(self.Restframe)
        self.horizontalLayout_6.setSpacing(6)
        self.horizontalLayout_6.setContentsMargins(11, 11, 11, 11)
        self.horizontalLayout_6.setObjectName(u"horizontalLayout_6")
        self.horizontalLayout_6.setContentsMargins(0, 0, 0, 0)
        self.redshiftlabel = QLabel(self.Restframe)
        self.redshiftlabel.setObjectName(u"redshiftlabel")
        self.redshiftlabel.setEnabled(True)

        self.horizontalLayout_6.addWidget(self.redshiftlabel)

        self.redshiftlineEdit = QLineEdit(self.Restframe)
        self.redshiftlineEdit.setObjectName(u"redshiftlineEdit")
        self.redshiftlineEdit.setEnabled(True)
        self.redshiftlineEdit.setLayoutDirection(Qt.RightToLeft)
        self.redshiftlineEdit.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.horizontalLayout_6.addWidget(self.redshiftlineEdit)

        self.restcomboBox = QComboBox(self.Restframe)
        self.restcomboBox.addItem("")
        self.restcomboBox.addItem("")
        self.restcomboBox.setObjectName(u"restcomboBox")
        self.restcomboBox.setEnabled(True)
        self.restcomboBox.setMinimumSize(QSize(120, 0))

        self.horizontalLayout_6.addWidget(self.restcomboBox)

        self.restlineEdit = QLineEdit(self.Restframe)
        self.restlineEdit.setObjectName(u"restlineEdit")
        self.restlineEdit.setEnabled(True)
        self.restlineEdit.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)

        self.horizontalLayout_6.addWidget(self.restlineEdit)


        self.gridLayout_2.addWidget(self.Restframe, 2, 0, 1, 1)

        self.FreeParametersframe = QFrame(self.page_2)
        self.FreeParametersframe.setObjectName(u"FreeParametersframe")
        self.FreeParametersframe.setEnabled(False)
        sizePolicy15 = QSizePolicy(QSizePolicy.Expanding, QSizePolicy.Maximum)
        sizePolicy15.setHorizontalStretch(0)
        sizePolicy15.setVerticalStretch(0)
        sizePolicy15.setHeightForWidth(self.FreeParametersframe.sizePolicy().hasHeightForWidth())
        self.FreeParametersframe.setSizePolicy(sizePolicy15)
        self.FreeParametersframe.setFrameShape(QFrame.StyledPanel)
        self.FreeParametersframe.setFrameShadow(QFrame.Raised)
        self.gridLayout_15 = QGridLayout(self.FreeParametersframe)
        self.gridLayout_15.setSpacing(6)
        self.gridLayout_15.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_15.setObjectName(u"gridLayout_15")
        self.yposcheckBox = QCheckBox(self.FreeParametersframe)
        self.yposcheckBox.setObjectName(u"yposcheckBox")
        sizePolicy7.setHeightForWidth(self.yposcheckBox.sizePolicy().hasHeightForWidth())
        self.yposcheckBox.setSizePolicy(sizePolicy7)

        self.gridLayout_15.addWidget(self.yposcheckBox, 3, 2, 1, 1)

        self.xposcheckBox = QCheckBox(self.FreeParametersframe)
        self.xposcheckBox.setObjectName(u"xposcheckBox")
        sizePolicy7.setHeightForWidth(self.xposcheckBox.sizePolicy().hasHeightForWidth())
        self.xposcheckBox.setSizePolicy(sizePolicy7)

        self.gridLayout_15.addWidget(self.xposcheckBox, 2, 2, 1, 1)

        self.vrotcheckBox = QCheckBox(self.FreeParametersframe)
        self.vrotcheckBox.setObjectName(u"vrotcheckBox")
        sizePolicy7.setHeightForWidth(self.vrotcheckBox.sizePolicy().hasHeightForWidth())
        self.vrotcheckBox.setSizePolicy(sizePolicy7)
        self.vrotcheckBox.setChecked(True)

        self.gridLayout_15.addWidget(self.vrotcheckBox, 2, 0, 1, 1)

        self.inccheckBox = QCheckBox(self.FreeParametersframe)
        self.inccheckBox.setObjectName(u"inccheckBox")
        sizePolicy7.setHeightForWidth(self.inccheckBox.sizePolicy().hasHeightForWidth())
        self.inccheckBox.setSizePolicy(sizePolicy7)
        self.inccheckBox.setChecked(True)

        self.gridLayout_15.addWidget(self.inccheckBox, 2, 1, 1, 1)

        self.pacheckBox = QCheckBox(self.FreeParametersframe)
        self.pacheckBox.setObjectName(u"pacheckBox")
        sizePolicy7.setHeightForWidth(self.pacheckBox.sizePolicy().hasHeightForWidth())
        self.pacheckBox.setSizePolicy(sizePolicy7)
        self.pacheckBox.setChecked(True)

        self.gridLayout_15.addWidget(self.pacheckBox, 3, 1, 1, 1)

        self.vradcheckBox = QCheckBox(self.FreeParametersframe)
        self.vradcheckBox.setObjectName(u"vradcheckBox")

        self.gridLayout_15.addWidget(self.vradcheckBox, 2, 4, 1, 1)

        self.z0checkBox = QCheckBox(self.FreeParametersframe)
        self.z0checkBox.setObjectName(u"z0checkBox")
        sizePolicy7.setHeightForWidth(self.z0checkBox.sizePolicy().hasHeightForWidth())
        self.z0checkBox.setSizePolicy(sizePolicy7)

        self.gridLayout_15.addWidget(self.z0checkBox, 3, 3, 1, 1)

        self.label_30 = QLabel(self.FreeParametersframe)
        self.label_30.setObjectName(u"label_30")
        sizePolicy12.setHeightForWidth(self.label_30.sizePolicy().hasHeightForWidth())
        self.label_30.setSizePolicy(sizePolicy12)

        self.gridLayout_15.addWidget(self.label_30, 1, 0, 1, 1)

        self.vdispcheckBox = QCheckBox(self.FreeParametersframe)
        self.vdispcheckBox.setObjectName(u"vdispcheckBox")
        sizePolicy7.setHeightForWidth(self.vdispcheckBox.sizePolicy().hasHeightForWidth())
        self.vdispcheckBox.setSizePolicy(sizePolicy7)
        self.vdispcheckBox.setChecked(True)

        self.gridLayout_15.addWidget(self.vdispcheckBox, 3, 0, 1, 1)

        self.vsyscheckBox = QCheckBox(self.FreeParametersframe)
        self.vsyscheckBox.setObjectName(u"vsyscheckBox")
        sizePolicy7.setHeightForWidth(self.vsyscheckBox.sizePolicy().hasHeightForWidth())
        self.vsyscheckBox.setSizePolicy(sizePolicy7)

        self.gridLayout_15.addWidget(self.vsyscheckBox, 2, 3, 1, 1)


        self.gridLayout_2.addWidget(self.FreeParametersframe, 1, 0, 1, 1)

        self.GalfitAdvtoolButton = QToolButton(self.page_2)
        self.GalfitAdvtoolButton.setObjectName(u"GalfitAdvtoolButton")
        self.GalfitAdvtoolButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.GalfitAdvtoolButton.setArrowType(Qt.LeftArrow)

        self.gridLayout_2.addWidget(self.GalfitAdvtoolButton, 11, 0, 1, 1, Qt.AlignRight)

        self.stackedWidget.addWidget(self.page_2)
        self.page_8 = QWidget()
        self.page_8.setObjectName(u"page_8")
        self.gridLayout_20 = QGridLayout(self.page_8)
        self.gridLayout_20.setSpacing(6)
        self.gridLayout_20.setContentsMargins(11, 11, 11, 11)
        self.gridLayout_20.setObjectName(u"gridLayout_20")
        self.horizontalSpacer_3 = QSpacerItem(40, 10, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_20.addItem(self.horizontalSpacer_3, 0, 0, 1, 1)

        self.SearchAdvgroupBox = QGroupBox(self.page_8)
        self.SearchAdvgroupBox.setObjectName(u"SearchAdvgroupBox")
        self.SearchAdvgroupBox.setEnabled(False)
        self.SearchAdvgroupBox.setFlat(True)
        self.SearchAdvgroupBox.setCheckable(True)
        self.SearchAdvgroupBox.setChecked(False)
        self.verticalLayout_11 = QVBoxLayout(self.SearchAdvgroupBox)
        self.verticalLayout_11.setSpacing(6)
        self.verticalLayout_11.setContentsMargins(11, 11, 11, 11)
        self.verticalLayout_11.setObjectName(u"verticalLayout_11")
        self.gridLayout_23 = QGridLayout()
        self.gridLayout_23.setSpacing(6)
        self.gridLayout_23.setObjectName(u"gridLayout_23")
        self.MaxAngSizeSpinBox = QDoubleSpinBox(self.SearchAdvgroupBox)
        self.MaxAngSizeSpinBox.setObjectName(u"MaxAngSizeSpinBox")
        self.MaxAngSizeSpinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.MaxAngSizeSpinBox.setMinimum(-1.000000000000000)
        self.MaxAngSizeSpinBox.setMaximum(99999.990000000005239)
        self.MaxAngSizeSpinBox.setValue(-1.000000000000000)

        self.gridLayout_23.addWidget(self.MaxAngSizeSpinBox, 10, 2, 1, 1)

        self.Minvoxellabel = QLabel(self.SearchAdvgroupBox)
        self.Minvoxellabel.setObjectName(u"Minvoxellabel")

        self.gridLayout_23.addWidget(self.Minvoxellabel, 8, 0, 1, 1)

        self.horizontalSpacer_7 = QSpacerItem(40, 20, QSizePolicy.Ignored, QSizePolicy.Minimum)

        self.gridLayout_23.addItem(self.horizontalSpacer_7, 11, 0, 1, 1)

        self.horizontalSpacer_9 = QSpacerItem(40, 20, QSizePolicy.Fixed, QSizePolicy.Minimum)

        self.gridLayout_23.addItem(self.horizontalSpacer_9, 0, 1, 1, 1)

        self.label_48 = QLabel(self.SearchAdvgroupBox)
        self.label_48.setObjectName(u"label_48")

        self.gridLayout_23.addWidget(self.label_48, 3, 0, 1, 1)

        self.MaxChanspinBox = QSpinBox(self.SearchAdvgroupBox)
        self.MaxChanspinBox.setObjectName(u"MaxChanspinBox")
        self.MaxChanspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.MaxChanspinBox.setMinimum(-1)
        self.MaxChanspinBox.setMaximum(99999)
        self.MaxChanspinBox.setValue(-1)

        self.gridLayout_23.addWidget(self.MaxChanspinBox, 9, 2, 1, 1)

        self.RejectcheckBox = QCheckBox(self.SearchAdvgroupBox)
        self.RejectcheckBox.setObjectName(u"RejectcheckBox")
        self.RejectcheckBox.setChecked(True)

        self.gridLayout_23.addWidget(self.RejectcheckBox, 12, 2, 1, 1)

        self.ThreshSpatialBox = QSpinBox(self.SearchAdvgroupBox)
        self.ThreshSpatialBox.setObjectName(u"ThreshSpatialBox")
        self.ThreshSpatialBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ThreshSpatialBox.setMinimum(-1)
        self.ThreshSpatialBox.setMaximum(10000)
        self.ThreshSpatialBox.setValue(-1)

        self.gridLayout_23.addWidget(self.ThreshSpatialBox, 2, 2, 1, 1)

        self.TwostagemergingcheckBox = QCheckBox(self.SearchAdvgroupBox)
        self.TwostagemergingcheckBox.setObjectName(u"TwostagemergingcheckBox")
        self.TwostagemergingcheckBox.setChecked(True)

        self.gridLayout_23.addWidget(self.TwostagemergingcheckBox, 12, 0, 1, 1)

        self.Mergelabel = QLabel(self.SearchAdvgroupBox)
        self.Mergelabel.setObjectName(u"Mergelabel")
        sizePolicy8.setHeightForWidth(self.Mergelabel.sizePolicy().hasHeightForWidth())
        self.Mergelabel.setSizePolicy(sizePolicy8)

        self.gridLayout_23.addWidget(self.Mergelabel, 0, 0, 1, 1)

        self.label_15 = QLabel(self.SearchAdvgroupBox)
        self.label_15.setObjectName(u"label_15")

        self.gridLayout_23.addWidget(self.label_15, 1, 0, 1, 1)

        self.AdjacentcheckBox = QCheckBox(self.SearchAdvgroupBox)
        self.AdjacentcheckBox.setObjectName(u"AdjacentcheckBox")
        self.AdjacentcheckBox.setChecked(True)

        self.gridLayout_23.addWidget(self.AdjacentcheckBox, 1, 2, 1, 1)

        self.MinChanlabel = QLabel(self.SearchAdvgroupBox)
        self.MinChanlabel.setObjectName(u"MinChanlabel")

        self.gridLayout_23.addWidget(self.MinChanlabel, 7, 0, 1, 1)

        self.MinVoxspinBox = QSpinBox(self.SearchAdvgroupBox)
        self.MinVoxspinBox.setObjectName(u"MinVoxspinBox")
        self.MinVoxspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.MinVoxspinBox.setMinimum(-1)
        self.MinVoxspinBox.setMaximum(99999)
        self.MinVoxspinBox.setValue(-1)

        self.gridLayout_23.addWidget(self.MinVoxspinBox, 8, 2, 1, 1)

        self.label_16 = QLabel(self.SearchAdvgroupBox)
        self.label_16.setObjectName(u"label_16")

        self.gridLayout_23.addWidget(self.label_16, 9, 0, 1, 1)

        self.label_49 = QLabel(self.SearchAdvgroupBox)
        self.label_49.setObjectName(u"label_49")

        self.gridLayout_23.addWidget(self.label_49, 2, 0, 1, 1)

        self.horizontalSpacer_4 = QSpacerItem(40, 20, QSizePolicy.Expanding, QSizePolicy.Minimum)

        self.gridLayout_23.addItem(self.horizontalSpacer_4, 1, 3, 1, 1)

        self.MinPixspinBox = QSpinBox(self.SearchAdvgroupBox)
        self.MinPixspinBox.setObjectName(u"MinPixspinBox")
        self.MinPixspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.MinPixspinBox.setMinimum(-1)
        self.MinPixspinBox.setMaximum(99999)
        self.MinPixspinBox.setValue(-1)

        self.gridLayout_23.addWidget(self.MinPixspinBox, 7, 2, 1, 1)

        self.horizontalSpacer_6 = QSpacerItem(40, 20, QSizePolicy.Ignored, QSizePolicy.Minimum)

        self.gridLayout_23.addItem(self.horizontalSpacer_6, 4, 0, 1, 1)

        self.MinPixlabel = QLabel(self.SearchAdvgroupBox)
        self.MinPixlabel.setObjectName(u"MinPixlabel")

        self.gridLayout_23.addWidget(self.MinPixlabel, 6, 0, 1, 1)

        self.Rejeclabel = QLabel(self.SearchAdvgroupBox)
        self.Rejeclabel.setObjectName(u"Rejeclabel")
        sizePolicy8.setHeightForWidth(self.Rejeclabel.sizePolicy().hasHeightForWidth())
        self.Rejeclabel.setSizePolicy(sizePolicy8)

        self.gridLayout_23.addWidget(self.Rejeclabel, 5, 0, 1, 1)

        self.label_17 = QLabel(self.SearchAdvgroupBox)
        self.label_17.setObjectName(u"label_17")

        self.gridLayout_23.addWidget(self.label_17, 10, 0, 1, 1)

        self.MinChanspinBox = QSpinBox(self.SearchAdvgroupBox)
        self.MinChanspinBox.setObjectName(u"MinChanspinBox")
        self.MinChanspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.MinChanspinBox.setMinimum(-1)
        self.MinChanspinBox.setMaximum(99999)
        self.MinChanspinBox.setValue(-1)

        self.gridLayout_23.addWidget(self.MinChanspinBox, 6, 2, 1, 1)

        self.ThreshVelspinBox = QSpinBox(self.SearchAdvgroupBox)
        self.ThreshVelspinBox.setObjectName(u"ThreshVelspinBox")
        self.ThreshVelspinBox.setAlignment(Qt.AlignRight|Qt.AlignTrailing|Qt.AlignVCenter)
        self.ThreshVelspinBox.setMinimum(-1)
        self.ThreshVelspinBox.setMaximum(10000)
        self.ThreshVelspinBox.setValue(-1)

        self.gridLayout_23.addWidget(self.ThreshVelspinBox, 3, 2, 1, 1)


        self.verticalLayout_11.addLayout(self.gridLayout_23)


        self.gridLayout_20.addWidget(self.SearchAdvgroupBox, 1, 0, 1, 1)

        self.SearchAdvtoolButton = QToolButton(self.page_8)
        self.SearchAdvtoolButton.setObjectName(u"SearchAdvtoolButton")
        self.SearchAdvtoolButton.setToolButtonStyle(Qt.ToolButtonTextBesideIcon)
        self.SearchAdvtoolButton.setArrowType(Qt.LeftArrow)

        self.gridLayout_20.addWidget(self.SearchAdvtoolButton, 2, 0, 1, 1, Qt.AlignRight)

        self.stackedWidget.addWidget(self.page_8)

        self.MaingridLayout.addWidget(self.stackedWidget, 1, 1, 1, 1)

        self.HidetoolButton = QToolButton(self.centralWidget)
        self.HidetoolButton.setObjectName(u"HidetoolButton")
        sizePolicy16 = QSizePolicy(QSizePolicy.Fixed, QSizePolicy.Maximum)
        sizePolicy16.setHorizontalStretch(0)
        sizePolicy16.setVerticalStretch(0)
        sizePolicy16.setHeightForWidth(self.HidetoolButton.sizePolicy().hasHeightForWidth())
        self.HidetoolButton.setSizePolicy(sizePolicy16)
        self.HidetoolButton.setArrowType(Qt.UpArrow)

        self.MaingridLayout.addWidget(self.HidetoolButton, 0, 0, 1, 1)

        self.listWidget = QListWidget(self.centralWidget)
        font2 = QFont()
        font2.setPointSize(2)
        __qlistwidgetitem = QListWidgetItem(self.listWidget)
        __qlistwidgetitem.setFont(font2);
        __qlistwidgetitem.setFlags(Qt.ItemIsDragEnabled|Qt.ItemIsUserCheckable|Qt.ItemIsEnabled);
        __qlistwidgetitem1 = QListWidgetItem(self.listWidget)
        __qlistwidgetitem1.setCheckState(Qt.Unchecked);
        __qlistwidgetitem2 = QListWidgetItem(self.listWidget)
        __qlistwidgetitem2.setCheckState(Qt.Unchecked);
        __qlistwidgetitem3 = QListWidgetItem(self.listWidget)
        __qlistwidgetitem3.setCheckState(Qt.Unchecked);
        __qlistwidgetitem4 = QListWidgetItem(self.listWidget)
        __qlistwidgetitem4.setCheckState(Qt.Unchecked);
        __qlistwidgetitem5 = QListWidgetItem(self.listWidget)
        __qlistwidgetitem5.setCheckState(Qt.Unchecked);
        __qlistwidgetitem6 = QListWidgetItem(self.listWidget)
        __qlistwidgetitem6.setFlags(Qt.ItemIsSelectable|Qt.ItemIsDragEnabled|Qt.ItemIsEnabled);
        __qlistwidgetitem7 = QListWidgetItem(self.listWidget)
        __qlistwidgetitem7.setFlags(Qt.ItemIsSelectable|Qt.ItemIsDragEnabled|Qt.ItemIsEnabled);
        self.listWidget.setObjectName(u"listWidget")
        sizePolicy17 = QSizePolicy(QSizePolicy.Preferred, QSizePolicy.Expanding)
        sizePolicy17.setHorizontalStretch(0)
        sizePolicy17.setVerticalStretch(0)
        sizePolicy17.setHeightForWidth(self.listWidget.sizePolicy().hasHeightForWidth())
        self.listWidget.setSizePolicy(sizePolicy17)
        self.listWidget.setMinimumSize(QSize(160, 200))
        self.listWidget.setFrameShape(QFrame.StyledPanel)
        self.listWidget.setLineWidth(1)
        self.listWidget.setSpacing(3)
        self.listWidget.setViewMode(QListView.ListMode)

        self.MaingridLayout.addWidget(self.listWidget, 1, 0, 1, 1, Qt.AlignTop)


        self.gridLayout.addLayout(self.MaingridLayout, 2, 1, 1, 1)

        BBaroloWindow.setCentralWidget(self.centralWidget)
        self.menuBar = QMenuBar(BBaroloWindow)
        self.menuBar.setObjectName(u"menuBar")
        self.menuBar.setGeometry(QRect(0, 0, 874, 33))
        self.menuBBaroloQt = QMenu(self.menuBar)
        self.menuBBaroloQt.setObjectName(u"menuBBaroloQt")
        BBaroloWindow.setMenuBar(self.menuBar)
        self.statusBar = QStatusBar(BBaroloWindow)
        self.statusBar.setObjectName(u"statusBar")
        BBaroloWindow.setStatusBar(self.statusBar)
        QWidget.setTabOrder(self.ParamlineEdit, self.ParampushButton)
        QWidget.setTabOrder(self.ParampushButton, self.FitslineEdit)
        QWidget.setTabOrder(self.FitslineEdit, self.FitspushButton)
        QWidget.setTabOrder(self.FitspushButton, self.BoxcheckBox)
        QWidget.setTabOrder(self.BoxcheckBox, self.XminspinBox)
        QWidget.setTabOrder(self.XminspinBox, self.XmaxspinBox)
        QWidget.setTabOrder(self.XmaxspinBox, self.YminspinBox)
        QWidget.setTabOrder(self.YminspinBox, self.YmaxspinBox)
        QWidget.setTabOrder(self.YmaxspinBox, self.ZminspinBox)
        QWidget.setTabOrder(self.ZminspinBox, self.ZmaxspinBox)
        QWidget.setTabOrder(self.ZmaxspinBox, self.ThreadspinBox)
        QWidget.setTabOrder(self.ThreadspinBox, self.ExtenspinBox)
        QWidget.setTabOrder(self.ExtenspinBox, self.AutocheckBox)
        QWidget.setTabOrder(self.AutocheckBox, self.OutfolderlineEdit)
        QWidget.setTabOrder(self.OutfolderlineEdit, self.OutfolderpushButton)
        QWidget.setTabOrder(self.OutfolderpushButton, self.AdvancedOptcheckBox)
        QWidget.setTabOrder(self.AdvancedOptcheckBox, self.AdvancedOptlineEdit)
        QWidget.setTabOrder(self.AdvancedOptlineEdit, self.KillpushButton)
        QWidget.setTabOrder(self.KillpushButton, self.listWidget)
        QWidget.setTabOrder(self.listWidget, self.textEdit)
        QWidget.setTabOrder(self.textEdit, self.wcscomboBox)
        QWidget.setTabOrder(self.wcscomboBox, self.XpospushButton)
        QWidget.setTabOrder(self.XpospushButton, self.XposFilelineEdit)
        QWidget.setTabOrder(self.XposFilelineEdit, self.XposFilespinBox)
        QWidget.setTabOrder(self.XposFilespinBox, self.YposcheckBox)
        QWidget.setTabOrder(self.YposcheckBox, self.yposlineEdit)
        QWidget.setTabOrder(self.yposlineEdit, self.YpospushButton)
        QWidget.setTabOrder(self.YpospushButton, self.YposFilelineEdit)
        QWidget.setTabOrder(self.YposFilelineEdit, self.YposFilespinBox)
        QWidget.setTabOrder(self.YposFilespinBox, self.VsyscheckBox)
        QWidget.setTabOrder(self.VsyscheckBox, self.VsysSpinBox)
        QWidget.setTabOrder(self.VsysSpinBox, self.VsyspushButton)
        QWidget.setTabOrder(self.VsyspushButton, self.VsysFilelineEdit)
        QWidget.setTabOrder(self.VsysFilelineEdit, self.VsysFilespinBox)
        QWidget.setTabOrder(self.VsysFilespinBox, self.VrotcheckBox)
        QWidget.setTabOrder(self.VrotcheckBox, self.VrotSpinBox)
        QWidget.setTabOrder(self.VrotSpinBox, self.VrotpushButton)
        QWidget.setTabOrder(self.VrotpushButton, self.VrotFilelineEdit)
        QWidget.setTabOrder(self.VrotFilelineEdit, self.VrotFilespinBox)
        QWidget.setTabOrder(self.VrotFilespinBox, self.VdispcheckBox)
        QWidget.setTabOrder(self.VdispcheckBox, self.VdispSpinBox)
        QWidget.setTabOrder(self.VdispSpinBox, self.VdisppushButton)
        QWidget.setTabOrder(self.VdisppushButton, self.VdispFilelineEdit)
        QWidget.setTabOrder(self.VdispFilelineEdit, self.VdispFilespinBox)
        QWidget.setTabOrder(self.VdispFilespinBox, self.InccheckBox)
        QWidget.setTabOrder(self.InccheckBox, self.IncSpinBox)
        QWidget.setTabOrder(self.IncSpinBox, self.IncDSpinBox)
        QWidget.setTabOrder(self.IncDSpinBox, self.IncpushButton)
        QWidget.setTabOrder(self.IncpushButton, self.IncFilelineEdit)
        QWidget.setTabOrder(self.IncFilelineEdit, self.IncFilespinBox)
        QWidget.setTabOrder(self.IncFilespinBox, self.PacheckBox)
        QWidget.setTabOrder(self.PacheckBox, self.PaSpinBox)
        QWidget.setTabOrder(self.PaSpinBox, self.PaDSpinBox)
        QWidget.setTabOrder(self.PaDSpinBox, self.PapushButton)
        QWidget.setTabOrder(self.PapushButton, self.PaFilelineEdit)
        QWidget.setTabOrder(self.PaFilelineEdit, self.PaFilespinBox)
        QWidget.setTabOrder(self.PaFilespinBox, self.Z0checkBox)
        QWidget.setTabOrder(self.Z0checkBox, self.Z0SpinBox)
        QWidget.setTabOrder(self.Z0SpinBox, self.Z0pushButton)
        QWidget.setTabOrder(self.Z0pushButton, self.Z0FilelineEdit)
        QWidget.setTabOrder(self.Z0FilelineEdit, self.Z0FilespinBox)
        QWidget.setTabOrder(self.Z0FilespinBox, self.DenscheckBox)
        QWidget.setTabOrder(self.DenscheckBox, self.DensSpinBox)
        QWidget.setTabOrder(self.DensSpinBox, self.DenspushButton)
        QWidget.setTabOrder(self.DenspushButton, self.DensFilelineEdit)
        QWidget.setTabOrder(self.DensFilelineEdit, self.DensFilespinBox)
        QWidget.setTabOrder(self.DensFilespinBox, self.SearchtypecomboBox)
        QWidget.setTabOrder(self.SearchtypecomboBox, self.CuttypecomboBox)
        QWidget.setTabOrder(self.CuttypecomboBox, self.primaryCutlineEdit)
        QWidget.setTabOrder(self.primaryCutlineEdit, self.GrowthcheckBox)
        QWidget.setTabOrder(self.GrowthcheckBox, self.Cuttype2comboBox)
        QWidget.setTabOrder(self.Cuttype2comboBox, self.SecondarycutlineEdit)
        QWidget.setTabOrder(self.SecondarycutlineEdit, self.SearchtoolButton)
        QWidget.setTabOrder(self.SearchtoolButton, self.ObmajSpinBox)
        QWidget.setTabOrder(self.ObmajSpinBox, self.NbmajSpinBox)
        QWidget.setTabOrder(self.NbmajSpinBox, self.ObminSpinBox)
        QWidget.setTabOrder(self.ObminSpinBox, self.NbminSpinBox)
        QWidget.setTabOrder(self.NbminSpinBox, self.ObpaSpinBox)
        QWidget.setTabOrder(self.ObpaSpinBox, self.NbpaSpinBox)
        QWidget.setTabOrder(self.NbpaSpinBox, self.MomentmapsgroupBox)
        QWidget.setTabOrder(self.MomentmapsgroupBox, self.ProfilecheckBox)
        QWidget.setTabOrder(self.ProfilecheckBox, self.TotalmapcheckBox)
        QWidget.setTabOrder(self.TotalmapcheckBox, self.VfieldcheckBox)
        QWidget.setTabOrder(self.VfieldcheckBox, self.DispmapcheckBox)
        QWidget.setTabOrder(self.DispmapcheckBox, self.PVgroupBox)
        QWidget.setTabOrder(self.PVgroupBox, self.PVXposSpinBox)
        QWidget.setTabOrder(self.PVXposSpinBox, self.PVYposSpinBox)
        QWidget.setTabOrder(self.PVYposSpinBox, self.PVPaSpinBox)
        QWidget.setTabOrder(self.PVPaSpinBox, self.MaskgroupBox)
        QWidget.setTabOrder(self.MaskgroupBox, self.BlankfactordoubleSpinBox)
        QWidget.setTabOrder(self.BlankfactordoubleSpinBox, self.BlankcutSpinBox)
        QWidget.setTabOrder(self.BlankcutSpinBox, self.LogtextEdit)
        QWidget.setTabOrder(self.LogtextEdit, self.vrotcheckBox)
        QWidget.setTabOrder(self.vrotcheckBox, self.vdispcheckBox)
        QWidget.setTabOrder(self.vdispcheckBox, self.inccheckBox)
        QWidget.setTabOrder(self.inccheckBox, self.pacheckBox)
        QWidget.setTabOrder(self.pacheckBox, self.xposcheckBox)
        QWidget.setTabOrder(self.xposcheckBox, self.yposcheckBox)
        QWidget.setTabOrder(self.yposcheckBox, self.vsyscheckBox)
        QWidget.setTabOrder(self.vsyscheckBox, self.z0checkBox)
        QWidget.setTabOrder(self.z0checkBox, self.AdvancedgroupBox)
        QWidget.setTabOrder(self.AdvancedgroupBox, self.GalfitAdvtoolButton)
        QWidget.setTabOrder(self.GalfitAdvtoolButton, self.SearchAdvgroupBox)
        QWidget.setTabOrder(self.SearchAdvgroupBox, self.AdjacentcheckBox)
        QWidget.setTabOrder(self.AdjacentcheckBox, self.ThreshSpatialBox)
        QWidget.setTabOrder(self.ThreshSpatialBox, self.ThreshVelspinBox)
        QWidget.setTabOrder(self.ThreshVelspinBox, self.MinChanspinBox)
        QWidget.setTabOrder(self.MinChanspinBox, self.MinPixspinBox)
        QWidget.setTabOrder(self.MinPixspinBox, self.MinVoxspinBox)
        QWidget.setTabOrder(self.MinVoxspinBox, self.MaxChanspinBox)
        QWidget.setTabOrder(self.MaxChanspinBox, self.MaxAngSizeSpinBox)
        QWidget.setTabOrder(self.MaxAngSizeSpinBox, self.TwostagemergingcheckBox)
        QWidget.setTabOrder(self.TwostagemergingcheckBox, self.RejectcheckBox)
        QWidget.setTabOrder(self.RejectcheckBox, self.SearchAdvtoolButton)
        QWidget.setTabOrder(self.SearchAdvtoolButton, self.HidetoolButton)
        QWidget.setTabOrder(self.HidetoolButton, self.xposlineEdit)
        QWidget.setTabOrder(self.xposlineEdit, self.GalfittoolButton)
        QWidget.setTabOrder(self.GalfittoolButton, self.RadseppushButton)
        QWidget.setTabOrder(self.RadseppushButton, self.RadsepFilelineEdit)
        QWidget.setTabOrder(self.RadsepFilelineEdit, self.RadsepSpinBox)
        QWidget.setTabOrder(self.RadsepSpinBox, self.VradpushButton)
        QWidget.setTabOrder(self.VradpushButton, self.VradFilelineEdit)
        QWidget.setTabOrder(self.VradFilelineEdit, self.VradSpinBox)
        QWidget.setTabOrder(self.VradSpinBox, self.VradcheckBox)
        QWidget.setTabOrder(self.VradcheckBox, self.VradFilespinBox)
        QWidget.setTabOrder(self.VradFilespinBox, self.SpatialSmoothgroupBox)
        QWidget.setTabOrder(self.SpatialSmoothgroupBox, self.ReducecheckBox)
        QWidget.setTabOrder(self.ReducecheckBox, self.SmoothOutlineEdit)
        QWidget.setTabOrder(self.SmoothOutlineEdit, self.FactorcheckBox)
        QWidget.setTabOrder(self.FactorcheckBox, self.FFTcheckBox)
        QWidget.setTabOrder(self.FFTcheckBox, self.FactordoubleSpinBox)
        QWidget.setTabOrder(self.FactordoubleSpinBox, self.SmoothOutpushButton)
        QWidget.setTabOrder(self.SmoothOutpushButton, self.BeamgroupBox)
        QWidget.setTabOrder(self.BeamgroupBox, self.HanninggroupBox)
        QWidget.setTabOrder(self.HanninggroupBox, self.HanningspinBox)
        QWidget.setTabOrder(self.HanningspinBox, self.MaptypecomboBox)
        QWidget.setTabOrder(self.MaptypecomboBox, self.rmsmapcheckBox)
        QWidget.setTabOrder(self.rmsmapcheckBox, self.redshiftlineEdit)
        QWidget.setTabOrder(self.redshiftlineEdit, self.restcomboBox)
        QWidget.setTabOrder(self.restcomboBox, self.restlineEdit)
        QWidget.setTabOrder(self.restlineEdit, self.vradcheckBox)
        QWidget.setTabOrder(self.vradcheckBox, self.ResetpushButton)
        QWidget.setTabOrder(self.ResetpushButton, self.XposcheckBox)
        QWidget.setTabOrder(self.XposcheckBox, self.NringscheckBox)
        QWidget.setTabOrder(self.NringscheckBox, self.NringsspinBox)
        QWidget.setTabOrder(self.NringsspinBox, self.RunpushButton)
        QWidget.setTabOrder(self.RunpushButton, self.RadsepcheckBox)
        QWidget.setTabOrder(self.RadsepcheckBox, self.RadsepFilespinBox)

        self.menuBar.addAction(self.menuBBaroloQt.menuAction())
        self.menuBBaroloQt.addSeparator()
        self.menuBBaroloQt.addAction(self.actionOpen_FITS_file)
        self.menuBBaroloQt.addAction(self.actionOpen_parameter_file)
        self.menuBBaroloQt.addAction(self.actionExport_parameter_file)
        self.menuBBaroloQt.addSeparator()
        self.menuBBaroloQt.addAction(self.actionReset_GUI)

        self.retranslateUi(BBaroloWindow)

        self.RunpushButton.setDefault(True)
        self.stackedWidget.setCurrentIndex(0)
        self.CuttypecomboBox.setCurrentIndex(1)
        self.SearchtypecomboBox.setCurrentIndex(1)
        self.Cuttype2comboBox.setCurrentIndex(1)
        self.Plot3comboBox.setCurrentIndex(2)
        self.Plot2comboBox.setCurrentIndex(1)
        self.Plot4comboBox.setCurrentIndex(3)
        self.WfunccomboBox.setCurrentIndex(1)
        self.FtypecomboBox.setCurrentIndex(1)
        self.NormcomboBox.setCurrentIndex(1)
        self.listWidget.setCurrentRow(-1)


        QMetaObject.connectSlotsByName(BBaroloWindow)
    # setupUi

    def retranslateUi(self, BBaroloWindow):
        BBaroloWindow.setWindowTitle(QCoreApplication.translate("BBaroloWindow", u"BBaroloGUI", None))
        self.actionOpen_FITS_file.setText(QCoreApplication.translate("BBaroloWindow", u"Open FITS file", None))
#if QT_CONFIG(shortcut)
        self.actionOpen_FITS_file.setShortcut(QCoreApplication.translate("BBaroloWindow", u"Ctrl+O", None))
#endif // QT_CONFIG(shortcut)
        self.actionOpen_parameter_file.setText(QCoreApplication.translate("BBaroloWindow", u"Open parameter file", None))
#if QT_CONFIG(shortcut)
        self.actionOpen_parameter_file.setShortcut(QCoreApplication.translate("BBaroloWindow", u"Ctrl+P", None))
#endif // QT_CONFIG(shortcut)
        self.actionExport_parameter_file.setText(QCoreApplication.translate("BBaroloWindow", u"Export parameter file", None))
#if QT_CONFIG(shortcut)
        self.actionExport_parameter_file.setShortcut(QCoreApplication.translate("BBaroloWindow", u"Ctrl+E", None))
#endif // QT_CONFIG(shortcut)
        self.actionOpen_list_of_FITS_files.setText(QCoreApplication.translate("BBaroloWindow", u"Open list of FITS files", None))
#if QT_CONFIG(shortcut)
        self.actionOpen_list_of_FITS_files.setShortcut(QCoreApplication.translate("BBaroloWindow", u"Ctrl+L", None))
#endif // QT_CONFIG(shortcut)
        self.actionReset_GUI.setText(QCoreApplication.translate("BBaroloWindow", u"Reset GUI", None))
        self.ParampushButton.setText(QCoreApplication.translate("BBaroloWindow", u"Open", None))
        self.BoxcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Box    ", None))
        self.Inputlabel.setText(QCoreApplication.translate("BBaroloWindow", u"Input file", None))
        self.FitspushButton.setText(QCoreApplication.translate("BBaroloWindow", u"Open", None))
        self.RunpushButton.setText(QCoreApplication.translate("BBaroloWindow", u"Run BBarolo!", None))
        self.AutocheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Auto Mode", None))
        self.Fitslabel.setText(QCoreApplication.translate("BBaroloWindow", u"FITS file  ", None))
        self.Xboxlabel.setText(QCoreApplication.translate("BBaroloWindow", u"x =", None))
        self.label_19.setText(QCoreApplication.translate("BBaroloWindow", u":", None))
        self.Yboxlabel.setText(QCoreApplication.translate("BBaroloWindow", u"  y =", None))
        self.label_23.setText(QCoreApplication.translate("BBaroloWindow", u":", None))
        self.Zboxlabel.setText(QCoreApplication.translate("BBaroloWindow", u"  z =", None))
        self.label_29.setText(QCoreApplication.translate("BBaroloWindow", u":", None))
        self.ResetpushButton.setText(QCoreApplication.translate("BBaroloWindow", u"Reset", None))
        self.Extenlabel.setText(QCoreApplication.translate("BBaroloWindow", u"Ext. #", None))
        self.Threadslabel.setText(QCoreApplication.translate("BBaroloWindow", u"Threads", None))
#if QT_CONFIG(tooltip)
        self.AdvancedOptcheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Add any option as a semi-colon separated list, i.e. PARAM1=VAL1; PARAM2=VAL2</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.AdvancedOptcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Advanced", None))
        self.KillpushButton.setText(QCoreApplication.translate("BBaroloWindow", u"Kill BBarolo", None))
        self.OutfolderlineEdit.setText("")
        self.OutfolderpushButton.setText(QCoreApplication.translate("BBaroloWindow", u"Open", None))
        self.Outfolderlabel.setText(QCoreApplication.translate("BBaroloWindow", u" Output folder ", None))
        self.textEdit.setHtml(QCoreApplication.translate("BBaroloWindow", u"<!DOCTYPE HTML PUBLIC \"-//W3C//DTD HTML 4.0//EN\" \"http://www.w3.org/TR/REC-html40/strict.dtd\">\n"
"<html><head><meta name=\"qrichtext\" content=\"1\" /><meta charset=\"utf-8\" /><style type=\"text/css\">\n"
"p, li { white-space: pre-wrap; }\n"
"hr { height: 1px; border-width: 0; }\n"
"li.unchecked::marker { content: \"\\2610\"; }\n"
"li.checked::marker { content: \"\\2612\"; }\n"
"</style></head><body style=\" font-family:'.AppleSystemUIFont'; font-size:13pt; font-weight:400; font-style:normal;\">\n"
"<p align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:'.Helvetica Neue DeskInterface';\">Welcome to </span><span style=\" font-family:'.Helvetica Neue DeskInterface'; font-weight:600;\">BBarolo</span><span style=\" font-family:'.Helvetica Neue DeskInterface';\">, a 3D-fitting tool to derive the kinematics of galaxies from emission-line observations. </span></p>\n"
"<p align=\"justify\" style=\"-qt-p"
                        "aragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:'.Helvetica Neue DeskInterface';\"><br /></p>\n"
"<p align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:'.Helvetica Neue DeskInterface';\">To start, insert a parameter file or a FITS file in the above fields and set the tasks in the left bar.</span></p>\n"
"<p align=\"justify\" style=\"-qt-paragraph-type:empty; margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px; font-family:'.Helvetica Neue DeskInterface';\"><br /></p>\n"
"<p align=\"justify\" style=\" margin-top:0px; margin-bottom:0px; margin-left:0px; margin-right:0px; -qt-block-indent:0; text-indent:0px;\"><span style=\" font-family:'.Helvetica Neue DeskInterface';\">When everything is set and ready, click on &quot;Run BBarolo!&quot; and enjoy. </span></p>"
                        "</body></html>", None))
        self.Logolabel.setText("")
        self.GalfitgroupBox.setTitle("")
#if QT_CONFIG(tooltip)
        self.PacheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>PA: Position angle of the rings. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.PacheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Position angle", None))
        self.RadseppushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
#if QT_CONFIG(tooltip)
        self.DenscheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>DENS: Gas column density. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.DenscheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Column density", None))
        self.YpospushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
        self.PaPMlabel.setText(QCoreApplication.translate("BBaroloWindow", u"+-", None))
        self.VrotpushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
        self.wcscomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"pixels", None))
        self.wcscomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"degrees", None))
        self.wcscomboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"sexagesimal", None))

        self.PapushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
        self.VsyspushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
#if QT_CONFIG(tooltip)
        self.IncFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_25.setText(QCoreApplication.translate("BBaroloWindow", u"km/s", None))
#if QT_CONFIG(tooltip)
        self.Z0checkBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Z0: Scale-height of the rings. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.Z0checkBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Scale height", None))
#if QT_CONFIG(tooltip)
        self.YposFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_18.setText(QCoreApplication.translate("BBaroloWindow", u"1/cm2", None))
        self.VradpushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
#if QT_CONFIG(tooltip)
        self.InccheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>INC: Inclination of the rings. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.InccheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Inclination", None))
        self.XpospushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
        self.label_13.setText(QCoreApplication.translate("BBaroloWindow", u"\"", None))
        self.label_12.setText(QCoreApplication.translate("BBaroloWindow", u"E20", None))
#if QT_CONFIG(tooltip)
        self.YposcheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>YPOS: Y-center of the rings, in pixels (starting from 0) or WCS coordinates. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.YposcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Y center", None))
        self.IncpushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
#if QT_CONFIG(tooltip)
        self.VrotFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(tooltip)
        self.VdispcheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>VDISP: Gas velocity dispersion of the rings. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.VdispcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Velocity dispersion", None))
#if QT_CONFIG(tooltip)
        self.VrotcheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>VROT: Rotation velocity of the rings. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.VrotcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Rotation velocity", None))
        self.VdisppushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
        self.IncPMlabel.setText(QCoreApplication.translate("BBaroloWindow", u"+-", None))
#if QT_CONFIG(tooltip)
        self.VsyscheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>VSYS: Systemic velocity of the galaxy. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.VsyscheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Systemic velocity", None))
#if QT_CONFIG(tooltip)
        self.RadsepFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_24.setText(QCoreApplication.translate("BBaroloWindow", u"km/s", None))
#if QT_CONFIG(tooltip)
        self.VradcheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>VRAD: Radial velocity. Default is 0.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.VradcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Radial velocity", None))
        self.label_14.setText(QCoreApplication.translate("BBaroloWindow", u"\"", None))
#if QT_CONFIG(tooltip)
        self.RadsepcheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>RADSEP: Separation between the rings. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.RadsepcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Ring width", None))
#if QT_CONFIG(tooltip)
        self.DensFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.Z0pushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
#if QT_CONFIG(tooltip)
        self.NringscheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>NRINGS: Number of rings of the 3D model. Check the box if you want the code to estimate it.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.NringscheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  Number of rings", None))
        self.label_26.setText(QCoreApplication.translate("BBaroloWindow", u"km/s", None))
        self.label_32.setText(QCoreApplication.translate("BBaroloWindow", u"km/s", None))
#if QT_CONFIG(tooltip)
        self.PaFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(tooltip)
        self.Z0FilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.xposlineEdit.setText(QCoreApplication.translate("BBaroloWindow", u"0.00", None))
        self.DenspushButton.setText(QCoreApplication.translate("BBaroloWindow", u"File", None))
#if QT_CONFIG(tooltip)
        self.XposFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p><p><br/></p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.yposlineEdit.setText(QCoreApplication.translate("BBaroloWindow", u"0.00", None))
#if QT_CONFIG(tooltip)
        self.VdispFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(tooltip)
        self.XposcheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>XPOS: X-center of the rings, in pixels (starting from 0) or WCS coordinates. Check the box if you want the code to estimate it. </p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.XposcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"  X center", None))
        self.GalfittoolButton.setText(QCoreApplication.translate("BBaroloWindow", u"Advanced", None))
#if QT_CONFIG(tooltip)
        self.VsysFilespinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Column (first=1)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.SearchgroupBox.setTitle("")
        self.label_41.setText(QCoreApplication.translate("BBaroloWindow", u"Primary Cut", None))
        self.CuttypecomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"Flux threshold", None))
        self.CuttypecomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"S/N cut", None))

        self.GrowthcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Secondary cut", None))
        self.label_42.setText(QCoreApplication.translate("BBaroloWindow", u"Type of search", None))
        self.SearchtypecomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"spectral", None))
        self.SearchtypecomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"spatial", None))

        self.SearchtoolButton.setText(QCoreApplication.translate("BBaroloWindow", u"Advanced", None))
        self.primaryCutlineEdit.setText(QCoreApplication.translate("BBaroloWindow", u"5.00", None))
        self.Cuttype2comboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"Flux threshold", None))
        self.Cuttype2comboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"S/N cut", None))

        self.SecondarycutlineEdit.setText(QCoreApplication.translate("BBaroloWindow", u"2.50", None))
        self.SmoothgroupBox.setTitle("")
        self.SpatialSmoothgroupBox.setTitle(QCoreApplication.translate("BBaroloWindow", u"Spatial Smoothing", None))
        self.ReducecheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Repixeling", None))
        self.FactorcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Factor:", None))
        self.label_39.setText(QCoreApplication.translate("BBaroloWindow", u"Output", None))
        self.FFTcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"FFT", None))
        self.SmoothOutpushButton.setText(QCoreApplication.translate("BBaroloWindow", u"Browse", None))
#if QT_CONFIG(tooltip)
        self.BeamgroupBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Manually insert beam parameters</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
#if QT_CONFIG(whatsthis)
        self.BeamgroupBox.setWhatsThis(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Whether you wish to manually insert beam parameters</p></body></html>", None))
#endif // QT_CONFIG(whatsthis)
        self.BeamgroupBox.setTitle(QCoreApplication.translate("BBaroloWindow", u"Manual Beam", None))
        self.label_7.setText(QCoreApplication.translate("BBaroloWindow", u"BPA", None))
        self.label_2.setText(QCoreApplication.translate("BBaroloWindow", u"New", None))
        self.label_4.setText(QCoreApplication.translate("BBaroloWindow", u"BMAJ", None))
        self.label_6.setText(QCoreApplication.translate("BBaroloWindow", u"BMIN", None))
        self.label_11.setText(QCoreApplication.translate("BBaroloWindow", u"' '", None))
        self.label_8.setText(QCoreApplication.translate("BBaroloWindow", u"-->", None))
        self.label_9.setText(QCoreApplication.translate("BBaroloWindow", u"-->", None))
        self.label_3.setText(QCoreApplication.translate("BBaroloWindow", u"deg", None))
        self.label_10.setText(QCoreApplication.translate("BBaroloWindow", u"' '", None))
        self.label.setText(QCoreApplication.translate("BBaroloWindow", u"Old", None))
        self.label_5.setText(QCoreApplication.translate("BBaroloWindow", u"-->", None))
        self.HanninggroupBox.setTitle(QCoreApplication.translate("BBaroloWindow", u"Hanning Smoothing", None))
        self.label_40.setText(QCoreApplication.translate("BBaroloWindow", u"Hanning window size (channels)", None))
#if QT_CONFIG(tooltip)
        self.HanningspinBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Must be an odd number!</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.MaskgroupBox.setTitle(QCoreApplication.translate("BBaroloWindow", u"Mask", None))
        self.label_54.setText(QCoreApplication.translate("BBaroloWindow", u"Signal-to-noise cut", None))
        self.label_50.setText(QCoreApplication.translate("BBaroloWindow", u"Smoothing factor", None))
        self.MomentmapsgroupBox.setTitle(QCoreApplication.translate("BBaroloWindow", u"Global spectrum and maps", None))
        self.label_35.setText(QCoreApplication.translate("BBaroloWindow", u"          Map type", None))
        self.label_44.setText(QCoreApplication.translate("BBaroloWindow", u"Column density map  ", None))
        self.TotalmapcheckBox.setText("")
        self.label_46.setText(QCoreApplication.translate("BBaroloWindow", u"Dispersion velosity map", None))
        self.MaptypecomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"MOMENT", None))
        self.MaptypecomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"GAUSSIAN", None))

        self.VfieldcheckBox.setText("")
        self.DispmapcheckBox.setText("")
        self.ProfilecheckBox.setText("")
        self.label_45.setText(QCoreApplication.translate("BBaroloWindow", u"Velocity map", None))
        self.label_47.setText(QCoreApplication.translate("BBaroloWindow", u"Global profile", None))
        self.label_38.setText(QCoreApplication.translate("BBaroloWindow", u"RMS noise map", None))
        self.rmsmapcheckBox.setText("")
        self.PVgroupBox.setTitle(QCoreApplication.translate("BBaroloWindow", u"Position-velocity diagrams", None))
        self.label_20.setText(QCoreApplication.translate("BBaroloWindow", u"Center coordinates (x,y) = (", None))
        self.label_22.setText(QCoreApplication.translate("BBaroloWindow", u") pixels", None))
        self.label_21.setText(QCoreApplication.translate("BBaroloWindow", u",", None))
        self.label_27.setText(QCoreApplication.translate("BBaroloWindow", u"Angle (degrees N->W)", None))
        self.Plot3comboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"VROT", None))
        self.Plot3comboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"VDISP", None))
        self.Plot3comboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"INC", None))
        self.Plot3comboBox.setItemText(3, QCoreApplication.translate("BBaroloWindow", u"PA", None))
        self.Plot3comboBox.setItemText(4, QCoreApplication.translate("BBaroloWindow", u"Z0", None))
        self.Plot3comboBox.setItemText(5, QCoreApplication.translate("BBaroloWindow", u"XPOS", None))
        self.Plot3comboBox.setItemText(6, QCoreApplication.translate("BBaroloWindow", u"YPOS", None))
        self.Plot3comboBox.setItemText(7, QCoreApplication.translate("BBaroloWindow", u"VSYS", None))

        self.Plot2comboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"VROT", None))
        self.Plot2comboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"VDISP", None))
        self.Plot2comboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"INC", None))
        self.Plot2comboBox.setItemText(3, QCoreApplication.translate("BBaroloWindow", u"PA", None))
        self.Plot2comboBox.setItemText(4, QCoreApplication.translate("BBaroloWindow", u"Z0", None))
        self.Plot2comboBox.setItemText(5, QCoreApplication.translate("BBaroloWindow", u"XPOS", None))
        self.Plot2comboBox.setItemText(6, QCoreApplication.translate("BBaroloWindow", u"YPOS", None))
        self.Plot2comboBox.setItemText(7, QCoreApplication.translate("BBaroloWindow", u"VSYS", None))

        self.Plot1comboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"VROT", None))
        self.Plot1comboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"VDISP", None))
        self.Plot1comboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"INC", None))
        self.Plot1comboBox.setItemText(3, QCoreApplication.translate("BBaroloWindow", u"PA", None))
        self.Plot1comboBox.setItemText(4, QCoreApplication.translate("BBaroloWindow", u"Z0", None))
        self.Plot1comboBox.setItemText(5, QCoreApplication.translate("BBaroloWindow", u"XPOS", None))
        self.Plot1comboBox.setItemText(6, QCoreApplication.translate("BBaroloWindow", u"YPOS", None))
        self.Plot1comboBox.setItemText(7, QCoreApplication.translate("BBaroloWindow", u"VSYS", None))

        self.Plot4comboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"VROT", None))
        self.Plot4comboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"VDISP", None))
        self.Plot4comboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"INC", None))
        self.Plot4comboBox.setItemText(3, QCoreApplication.translate("BBaroloWindow", u"PA", None))
        self.Plot4comboBox.setItemText(4, QCoreApplication.translate("BBaroloWindow", u"Z0", None))
        self.Plot4comboBox.setItemText(5, QCoreApplication.translate("BBaroloWindow", u"XPOS", None))
        self.Plot4comboBox.setItemText(6, QCoreApplication.translate("BBaroloWindow", u"YPOS", None))
        self.Plot4comboBox.setItemText(7, QCoreApplication.translate("BBaroloWindow", u"VSYS", None))

        self.AdvancedgroupBox.setTitle(QCoreApplication.translate("BBaroloWindow", u"Advanced", None))
        self.WfunccomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"uniform", None))
        self.WfunccomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"cosine", None))
        self.WfunccomboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"cosine^2", None))

        self.FtypecomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"chi-squared", None))
        self.FtypecomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"|m-o|", None))
        self.FtypecomboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"|m-o|/(m+o)", None))

#if QT_CONFIG(tooltip)
        self.SecondstagecheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>TWOSTAGE: Parameter regularization and second fitting stage.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.SecondstagecheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Second stage", None))
        self.TollineEdit.setText(QCoreApplication.translate("BBaroloWindow", u"1E-03", None))
#if QT_CONFIG(tooltip)
        self.Ftypelabel.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>FTYPE: Function to be minimized.</p><p>           m=model, o=observations</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.Ftypelabel.setText(QCoreApplication.translate("BBaroloWindow", u"Fitting function", None))
#if QT_CONFIG(tooltip)
        self.label_28.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>MASK: type of mask. SMOOTH (smoothing and cutting), SEARCH (source finder) or THRESHOLD (cutting) </p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_28.setText(QCoreApplication.translate("BBaroloWindow", u"Mask", None))
#if QT_CONFIG(tooltip)
        self.Tollabel.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>TOL: Tolerance of the fit.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.Tollabel.setText(QCoreApplication.translate("BBaroloWindow", u"Tolerance", None))
        self.MasktoolButton.setText(QCoreApplication.translate("BBaroloWindow", u"Mask Settings", None))
#if QT_CONFIG(tooltip)
        self.Wfunclabel.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>WFUNC: Weighting function.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.Wfunclabel.setText(QCoreApplication.translate("BBaroloWindow", u"Weighting function", None))
#if QT_CONFIG(tooltip)
        self.ErrorscheckBox.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>FLAGERRORS: Estimating the errors. </p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.ErrorscheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Errors", None))
#if QT_CONFIG(tooltip)
        self.label_37.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>NV: Number of subclouds.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_37.setText(QCoreApplication.translate("BBaroloWindow", u"Num.", None))
#if QT_CONFIG(tooltip)
        self.label_33.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>LTYPE: Layer type of scale-height.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_33.setText(QCoreApplication.translate("BBaroloWindow", u"Layer type", None))
#if QT_CONFIG(tooltip)
        self.label_34.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>CDENS: Cloud column density.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_34.setText(QCoreApplication.translate("BBaroloWindow", u"Cloud CD (1E20)", None))
#if QT_CONFIG(tooltip)
        self.label_36.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>Type of normalization, pixel by pixel (LOCAL) or azimuthally averaged (AZIMUTHAL)</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_36.setText(QCoreApplication.translate("BBaroloWindow", u"Normalization     ", None))
        self.NormcomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"LOCAL", None))
        self.NormcomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"AZIMUTHAL", None))
        self.NormcomboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"NONE", None))

        self.MaskcomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"SMOOTH", None))
        self.MaskcomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"SEARCH", None))
        self.MaskcomboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"SMOOTH&SEARCH", None))
        self.MaskcomboBox.setItemText(3, QCoreApplication.translate("BBaroloWindow", u"THRESHOLD", None))
        self.MaskcomboBox.setItemText(4, QCoreApplication.translate("BBaroloWindow", u"NEGATIVE", None))
        self.MaskcomboBox.setItemText(5, QCoreApplication.translate("BBaroloWindow", u"NONE", None))

        self.LtypecomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"Gaussian", None))
        self.LtypecomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"Sech2", None))
        self.LtypecomboBox.setItemText(2, QCoreApplication.translate("BBaroloWindow", u"Exponential", None))
        self.LtypecomboBox.setItemText(3, QCoreApplication.translate("BBaroloWindow", u"Lorentzian", None))
        self.LtypecomboBox.setItemText(4, QCoreApplication.translate("BBaroloWindow", u"Box", None))

        self.AdriftcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Asymmetric drift correction", None))
        self.redshiftlabel.setText(QCoreApplication.translate("BBaroloWindow", u"  Redshift  ", None))
        self.redshiftlineEdit.setText(QCoreApplication.translate("BBaroloWindow", u"0", None))
        self.restcomboBox.setItemText(0, QCoreApplication.translate("BBaroloWindow", u"Rest wavelength", None))
        self.restcomboBox.setItemText(1, QCoreApplication.translate("BBaroloWindow", u"Rest frequency", None))

        self.restlineEdit.setText(QCoreApplication.translate("BBaroloWindow", u"0", None))
        self.yposcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"ycenter", None))
        self.xposcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"xcenter", None))
        self.vrotcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"vrot", None))
        self.inccheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"inc", None))
        self.pacheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"pa", None))
        self.vradcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"vrad", None))
        self.z0checkBox.setText(QCoreApplication.translate("BBaroloWindow", u"scale height", None))
#if QT_CONFIG(tooltip)
        self.label_30.setToolTip(QCoreApplication.translate("BBaroloWindow", u"<html><head/><body><p>FREE: Parameters to fit.</p></body></html>", None))
#endif // QT_CONFIG(tooltip)
        self.label_30.setText(QCoreApplication.translate("BBaroloWindow", u" Free parameters", None))
        self.vdispcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"disp", None))
        self.vsyscheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"vsys", None))
        self.GalfitAdvtoolButton.setText(QCoreApplication.translate("BBaroloWindow", u"Back", None))
        self.SearchAdvgroupBox.setTitle(QCoreApplication.translate("BBaroloWindow", u"Advanced", None))
        self.Minvoxellabel.setText(QCoreApplication.translate("BBaroloWindow", u"    Minimum # of voxels", None))
        self.label_48.setText(QCoreApplication.translate("BBaroloWindow", u"    Max channel separation", None))
        self.RejectcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Reject before merging", None))
        self.TwostagemergingcheckBox.setText(QCoreApplication.translate("BBaroloWindow", u"Two Stage Merging", None))
        self.Mergelabel.setText(QCoreApplication.translate("BBaroloWindow", u"Merging criteria", None))
        self.label_15.setText(QCoreApplication.translate("BBaroloWindow", u"    Merge only adjacent", None))
        self.AdjacentcheckBox.setText("")
        self.MinChanlabel.setText(QCoreApplication.translate("BBaroloWindow", u"    Minimum # of channels   ", None))
        self.label_16.setText(QCoreApplication.translate("BBaroloWindow", u"    Maximum # of channels", None))
        self.label_49.setText(QCoreApplication.translate("BBaroloWindow", u"    Max pixel separation", None))
        self.MinPixlabel.setText(QCoreApplication.translate("BBaroloWindow", u"    Minimum # of pixels", None))
        self.Rejeclabel.setText(QCoreApplication.translate("BBaroloWindow", u"Rejection criteria", None))
        self.label_17.setText(QCoreApplication.translate("BBaroloWindow", u"    Maximum angular size ( ' )", None))
        self.SearchAdvtoolButton.setText(QCoreApplication.translate("BBaroloWindow", u"Back", None))
        self.HidetoolButton.setText(QCoreApplication.translate("BBaroloWindow", u"Tasks", None))

        __sortingEnabled = self.listWidget.isSortingEnabled()
        self.listWidget.setSortingEnabled(False)
        ___qlistwidgetitem = self.listWidget.item(1)
        ___qlistwidgetitem.setText(QCoreApplication.translate("BBaroloWindow", u"3D Fit", None));
        ___qlistwidgetitem1 = self.listWidget.item(2)
        ___qlistwidgetitem1.setText(QCoreApplication.translate("BBaroloWindow", u"3D Model", None));
        ___qlistwidgetitem2 = self.listWidget.item(3)
        ___qlistwidgetitem2.setText(QCoreApplication.translate("BBaroloWindow", u"Search", None));
        ___qlistwidgetitem3 = self.listWidget.item(4)
        ___qlistwidgetitem3.setText(QCoreApplication.translate("BBaroloWindow", u"Smooth", None));
        ___qlistwidgetitem4 = self.listWidget.item(5)
        ___qlistwidgetitem4.setText(QCoreApplication.translate("BBaroloWindow", u"Maps", None));
        ___qlistwidgetitem5 = self.listWidget.item(6)
        ___qlistwidgetitem5.setText(QCoreApplication.translate("BBaroloWindow", u"Log", None));
        ___qlistwidgetitem6 = self.listWidget.item(7)
        ___qlistwidgetitem6.setText(QCoreApplication.translate("BBaroloWindow", u"Output", None));
        self.listWidget.setSortingEnabled(__sortingEnabled)

        self.menuBBaroloQt.setTitle(QCoreApplication.translate("BBaroloWindow", u"File", None))
    # retranslateUi

