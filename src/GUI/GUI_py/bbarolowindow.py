"""
This module defines the main window for the GUI of BBarolo.
"""

########################################################################
# Copyright (C) 2023 Enrico Di Teodoro
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
########################################################################

import multiprocessing

try:
    from astropy.io import fits
    import matplotlib.pyplot as plt
    import matplotlib as mpl
except ImportError:
    raise ImportError(f"The GUI requires astropy")


from imports import *
from moc.ui_bbarolowindow import Ui_BBaroloWindow
from bbarolowindow_io import *


class MainWindow(QMainWindow, Ui_BBaroloWindow):

    def __init__(self,parent=None,args=None):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self) 
        
        self.BBexe = QDir.currentPath()+"/BBarolo"
        self.proc = None
        self.obj = ""
        self.out_path = ""
        
        self.initialize_GUI()
        self.defineSlots()
        self.shrinkWindow()
        self.plotParameters()

    
    def initialize_GUI(self):
        
        # Validators for some lineEdits
        self.xposlineEdit.setValidator(QDoubleValidator(0,10000,4,self))
        self.yposlineEdit.setValidator(QDoubleValidator(0,10000,4,self))
        self.primaryCutlineEdit.setValidator(QDoubleValidator(0,10000,10,self))
        self.SecondarycutlineEdit.setValidator(QDoubleValidator(0,10000,10,self))
        self.restlineEdit.setValidator(QRegularExpressionValidator(QRegularExpression("[0-9]+.[0-9]+e[0-9]+")))
        self.redshiftlineEdit.setValidator(QRegularExpressionValidator(QRegularExpression("[0-9]+.[0-9]+")))
        
        # Work around of current path for Mac OS apps
        execPath = QApplication.applicationDirPath()
        currentPath = QDir.currentPath()
        if currentPath=='/' and 'Contents/MacOS' in execPath:
            self.BBexe = execPath+"/BBarolo"
            currentPath = execPath[:execPath.find('/BBaroloGUI.app/Contents/MacOS')]
            QDir.setCurrent(currentPath)
            
        self.out_path = QDir.currentPath()+"/output/"
        self.OutfolderlineEdit.setText(self.out_path)
        
        self.Hide_All_3DFit_file(True)
        self.MaskThreshSpinBox.setHidden(True)
        self.KillpushButton.setEnabled(False)
        
        self.Logolabel.setPixmap(QPixmap(":/resources/bbarolo.png"))
        
        # Setting number of threads
        self.ThreadspinBox.setValue(multiprocessing.cpu_count())
    
        # Setting monospace font for output console
        fontfamily = 'Monospace'
        fontsize = 10
        if 'darwin' in sys.platform: 
            fontfamily = 'Monaco'
            fontsize = 12
        self.LogtextEdit.setFont(QFont(fontfamily,fontsize))
                
        '''
        self.plot1.setInteraction(QCP::iRangeDrag, True)
        self.plot2.setInteraction(QCP::iRangeDrag, True)
        self.plot3.setInteraction(QCP::iRangeDrag, True)
        self.plot4.setInteraction(QCP::iRangeDrag, True)
        self.plot1.setInteraction(QCP::iRangeZoom, True)
        self.plot2.setInteraction(QCP::iRangeZoom, True)
        self.plot3.setInteraction(QCP::iRangeZoom, True)
        self.plot4.setInteraction(QCP::iRangeZoom, True)
        '''
        

    def defineSlots(self):
        
        ##################################################################
        # General GUI slots
        ##################################################################
        self.AutocheckBox.stateChanged.connect(self.AutocheckBox_stateChanged)
        self.BoxcheckBox.stateChanged[int].connect(self.BoxcheckBox_stateChanged)
        self.HidetoolButton.clicked.connect(self.HidetoolButton_clicked)
        self.ResetpushButton.clicked.connect(lambda: self.resetGUI(True))
        self.AdvancedOptcheckBox.stateChanged[int].connect(self.AdvancedOptcheckBox_stateChanged)
        self.AdvancedOptlineEdit.editingFinished.connect(self.AdvancedOptlineEdit_editingFinished)
        self.FitspushButton.clicked.connect(self.FitspushButton_clicked)
        self.FitslineEdit.editingFinished.connect(self.FitslineEdit_editingFinished)
        self.ParampushButton.clicked.connect(self.ParampushButton_clicked)
        self.ParamlineEdit.editingFinished.connect(self.ParamlineEdit_editingFinished)
        self.OutfolderpushButton.clicked.connect(self.OutfolderpushButton_clicked)
        self.KillpushButton.clicked.connect(self.KillpushButton_clicked)
        self.RunpushButton.clicked.connect(self.RunpushButton_clicked)
        
        # listWidget slots
        self.listWidget.currentRowChanged[int].connect(self.listWidget_currentRowChanged)
        self.listWidget.itemChanged.connect(self.listWidget_itemChanged)
        self.listWidget.itemDoubleClicked.connect(self.listWidget_itemDoubleClicked)
        
        # Menu bar slots
        self.actionReset_GUI.triggered.connect(lambda: self.resetGUI(True))
        self.actionOpen_FITS_file.triggered.connect(self.FitspushButton_clicked)
        self.actionOpen_parameter_file.triggered.connect(self.ParampushButton_clicked)
        self.actionExport_parameter_file.triggered.connect(self.actionExport_parameter_file_triggered)
        
        #// Combobox plotting slots
        self.Plot1comboBox.activated[int].connect(self.plotParameters)
        self.Plot2comboBox.activated[int].connect(self.plotParameters)
        self.Plot3comboBox.activated[int].connect(self.plotParameters)
        self.Plot4comboBox.activated[int].connect(self.plotParameters)
                
        
        ##################################################################
        # 3DFIT slots (regular options)
        ##################################################################
        def sChanged(cb, qw, qw2=None):
            qw.setDisabled(cb.isChecked())
            if qw2: qw2.setDisabled(cb.isChecked())
            
        self.NringscheckBox.stateChanged.connect(lambda i: sChanged(self.NringscheckBox,self.NringsspinBox))
        self.RadsepcheckBox.stateChanged.connect(lambda i: sChanged(self.RadsepcheckBox,self.RadsepSpinBox))
        self.XposcheckBox.stateChanged.connect(lambda i: sChanged(self.XposcheckBox,self.xposlineEdit))
        self.YposcheckBox.stateChanged.connect(lambda i: sChanged(self.YposcheckBox,self.yposlineEdit))
        self.VsyscheckBox.stateChanged.connect(lambda i: sChanged(self.VsyscheckBox,self.VsysSpinBox))
        self.VrotcheckBox.stateChanged.connect(lambda i: sChanged(self.VrotcheckBox,self.VrotSpinBox))
        self.VdispcheckBox.stateChanged.connect(lambda i: sChanged(self.VdispcheckBox,self.VdispSpinBox))
        self.VradcheckBox.stateChanged.connect(lambda i: sChanged(self.VradcheckBox,self.VradSpinBox))
        self.InccheckBox.stateChanged.connect(lambda i: sChanged(self.InccheckBox,self.IncSpinBox,self.IncDSpinBox))
        self.PacheckBox.stateChanged.connect(lambda i: sChanged(self.PacheckBox,self.PaSpinBox,self.PaDSpinBox))
        self.Z0checkBox.stateChanged.connect(lambda i: sChanged(self.Z0checkBox,self.Z0SpinBox))
        
        self.RadseppushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.RadsepFilelineEdit,self.RadsepFilespinBox,\
                                                        self.RadsepSpinBox,self.NringsspinBox))
        self.XpospushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.XposFilelineEdit,self.XposFilespinBox,self.xposlineEdit))
        self.YpospushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.YposFilelineEdit,self.YposFilespinBox,self.yposlineEdit))
        self.VsyspushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.VsysFilelineEdit,self.VsysFilespinBox,self.VsysSpinBox))
        self.VrotpushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.VrotFilelineEdit,self.VrotFilespinBox,self.VrotSpinBox))
        self.VdisppushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.VdispFilelineEdit,self.VdispFilespinBox,self.VdispSpinBox))
        self.VradpushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.VradFilelineEdit,self.VradFilespinBox,self.VradSpinBox))
        self.IncpushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.IncFilelineEdit,self.IncFilespinBox,self.IncSpinBox))
        self.PapushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.PaFilelineEdit,self.PaFilespinBox,self.PaSpinBox))
        self.Z0pushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.Z0FilelineEdit,self.Z0FilespinBox,self.Z0SpinBox))
        self.DenspushButton.clicked.connect(lambda i: self.Selected_3DFit_file(self.DensFilelineEdit,self.DensFilespinBox,self.DensSpinBox))
        
        self.RadsepFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.RadsepFilelineEdit,self.RadsepFilespinBox,\
                                                              self.RadsepSpinBox,self.RadsepFilelineEdit.text()=='',self.NringsspinBox))
        self.XposFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.XposFilelineEdit,self.XposFilespinBox,\
                                                              self.xposlineEdit,self.XposFilelineEdit.text()==''))
        self.YposFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.YposFilelineEdit,self.YposFilespinBox,\
                                                              self.yposlineEdit,self.YposFilelineEdit.text()==''))
        self.VsysFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.VsysFilelineEdit,self.VsysFilespinBox,\
                                                              self.VsysSpinBox,self.VsysFilelineEdit.text()==''))
        self.VrotFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.VrotFilelineEdit,self.VrotFilespinBox,\
                                                              self.VrotSpinBox,self.VrotFilelineEdit.text()==''))
        self.VdispFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.VdispFilelineEdit,self.VdispFilespinBox,\
                                                              self.VdispSpinBox,self.VdispFilelineEdit.text()==''))
        self.VradFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.VradFilelineEdit,self.VradFilespinBox,\
                                                              self.VradSpinBox,self.VradFilelineEdit.text()==''))
        self.IncFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.IncFilelineEdit,self.IncFilespinBox,\
                                                              self.IncSpinBox,self.IncFilelineEdit.text()==''))
        self.PaFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.PaFilelineEdit,self.PaFilespinBox,\
                                                              self.PaSpinBox,self.PaFilelineEdit.text()==''))
        self.Z0FilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.Z0FilelineEdit,self.Z0FilespinBox,\
                                                              self.Z0SpinBox,self.Z0FilelineEdit.text()==''))
        self.DensFilelineEdit.editingFinished.connect(lambda: Hide_3DFit_file(self.DensFilelineEdit,self.DensFilespinBox,\
                                                              self.DensSpinBox,self.DensFilelineEdit.text()==''))
        self.wcscomboBox.activated[int].connect(self.wcscomboBox_activated)

        ##################################################################
        # 3DFIT slots (advanced options)
        ##################################################################
        self.GalfittoolButton.clicked.connect(lambda i: self.stackedWidget.setCurrentIndex(self.stackedWidget.count()-2))
        self.GalfitAdvtoolButton.clicked.connect(lambda i: self.stackedWidget.setCurrentIndex(1))
        self.SecondstagecheckBox.stateChanged.connect(lambda i: sChanged(self.SecondstagecheckBox,self.PolynspinBox))
        self.MasktoolButton.clicked.connect(self.MasktoolButton_clicked)
        self.MaskcomboBox.currentIndexChanged.connect(self.MaskcomboBox_currentIndexChanged)
        
        ##################################################################
        # SEARCH slots 
        ##################################################################
        self.SearchgroupBox.clicked.connect(self.SearchgroupBox_clicked)
        self.GrowthcheckBox.stateChanged.connect(self.GrowthcheckBox_stateChanged)
        self.SearchAdvgroupBox.clicked.connect(self.SearchAdvgroupBox_clicked)
        self.AdjacentcheckBox.stateChanged.connect(self.AdjacentcheckBox_stateChanged)
        
        self.SearchtoolButton.clicked.connect(lambda: self.stackedWidget.setCurrentIndex(self.stackedWidget.count()-1))
        self.SearchAdvtoolButton.clicked.connect(lambda: self.stackedWidget.setCurrentIndex(2))
        
        ##################################################################
        # SMOOTH slots
        ##################################################################
        self.BeamgroupBox.clicked.connect(lambda: self.FactorcheckBox.setChecked(not self.BeamgroupBox.isChecked()))
        self.FactorcheckBox.stateChanged[int].connect(self.FactorcheckBox_stateChanged)
        self.HanningspinBox.valueChanged.connect(self.HanningspinBox_valueChanged)
        self.SmoothOutpushButton.clicked.connect(self.SmoothOutpushButton_clicked)
        self.FactordoubleSpinBox.valueChanged.connect(self.FactordoubleSpinBox_valueChanged)
        self.NbmajSpinBox.valueChanged.connect(self.update_SmoothOutFile)
        self.ReducecheckBox.clicked.connect(self.update_SmoothOutFile)
        
    
    ############################################################
    #### General GUI slots & functions
    ############################################################
    
    def shrinkWindow(self):
        # Resizing the main window
        for i in range(5):
            QApplication.processEvents()
        self.resize(self.minimumSize())


    @Slot()
    def resetGUI(self, resetFITSfile : bool = True):
        self.enable_All()
        self.ParamlineEdit.clear()
        if resetFITSfile:
            self.FitslineEdit.clear()
        self.OutfolderlineEdit.clear()
        self.ThreadspinBox.setValue(1)
        self.ExtenspinBox.setValue(0)
        self.AdvancedOptlineEdit.clear()
        self.AdvancedOptcheckBox.setChecked(False)

        self.BoxcheckBox.setChecked(False)
        self.XminspinBox.setValue(0)
        self.XmaxspinBox.setValue(0)
        self.YminspinBox.setValue(0)
        self.YmaxspinBox.setValue(0)
        self.ZminspinBox.setValue(0)
        self.ZmaxspinBox.setValue(0)
        
        for i in [1,2,3,4,5]:
            self.setTaskFlag(i,Qt.CheckState.Unchecked)
            
        self.NringscheckBox.setChecked(False)
        self.NringsspinBox.setValue(0)
        self.RadsepcheckBox.setChecked(False)
        self.RadsepSpinBox.setValue(0.00)
        self.XposcheckBox.setChecked(False)
        self.xposlineEdit.clear()
        self.YposcheckBox.setChecked(False)
        self.yposlineEdit.clear()
        self.wcscomboBox.setCurrentIndex(0)
        self.wcscomboBox_activated(0)
        self.VsyscheckBox.setChecked(False)
        self.VsysSpinBox.setValue(0.00)
        self.VrotcheckBox.setChecked(False)
        self.VrotSpinBox.setValue(0.00)
        self.VdispcheckBox.setChecked(False)
        self.VdispSpinBox.setValue(0.00)
        self.VradcheckBox.setChecked(False)
        self.VradSpinBox.setValue(0.00)
        self.InccheckBox.setChecked(False)
        self.IncSpinBox.setValue(0.00)
        self.IncDSpinBox.setValue(10.00)
        self.PacheckBox.setChecked(False)
        self.PaSpinBox.setValue(0.00)
        self.PaDSpinBox.setValue(15.00)
        self.Z0checkBox.setChecked(False)
        self.Z0SpinBox.setValue(0.00)
        self.Z0SpinBox.setValue(0.00)
        self.Hide_All_3DFit_file(True)

        self.FreeParametersframe.setDisabled(True)
        self.vrotcheckBox.setChecked(True)
        self.vdispcheckBox.setChecked(True)
        self.inccheckBox.setChecked(True)
        self.pacheckBox.setChecked(True)
        self.xposcheckBox.setChecked(False)
        self.yposcheckBox.setChecked(False)
        self.vsyscheckBox.setChecked(False)
        self.z0checkBox.setChecked(False)
        self.vradcheckBox.setChecked(False)

        self.Restframe.setDisabled(True)
        self.redshiftlineEdit.setText("0")
        self.restcomboBox.setCurrentIndex(0)
        self.restlineEdit.setText("0")
        
        self.AdvancedgroupBox.setDisabled(True)
        self.LtypecomboBox.setCurrentIndex(0)
        self.CdensSpinBox.setValue(10)
        self.NVspinBox.setValue(-1)
        self.FtypecomboBox.setCurrentIndex(1)
        self.WfunccomboBox.setCurrentIndex(1)
        self.SecondstagecheckBox.setChecked(True)
        self.PolynspinBox.setValue(-1)
        self.NormcomboBox.setCurrentIndex(0)
        self.TollineEdit.setText("1E-03")
        self.ErrorscheckBox.setChecked(False)

        self.SearchgroupBox.setChecked(False)
        self.SearchtypecomboBox.setCurrentIndex(1)
        self.CuttypecomboBox.setCurrentIndex(1)
        self.primaryCutlineEdit.setText("5.00")
        self.Cuttype2comboBox.setCurrentIndex(1)
        self.SecondarycutlineEdit.setText("2.50")
        self.SearchAdvgroupBox.setDisabled(True)
        self.SearchAdvgroupBox.setChecked(False)
        self.AdjacentcheckBox.setChecked(True)
        self.MinPixspinBox.setValue(-1)
        self.MinChanspinBox.setValue(-1)
        self.MinVoxspinBox.setValue(-1)
        self.MaxChanspinBox.setValue(-1)
        self.MaxAngSizeSpinBox.setValue(-1)
        self.RejectcheckBox.setChecked(True)
        self.TwostagemergingcheckBox.setChecked(True)

        self.SmoothOutlineEdit.clear()
        self.FactorcheckBox.setChecked(True)
        self.FactordoubleSpinBox.setValue(2)
        self.BeamgroupBox.setChecked(False)
        self.ObmajSpinBox.setValue(-1)
        self.ObminSpinBox.setValue(-1)
        self.ObpaSpinBox.setValue(0)
        self.NbmajSpinBox.setValue(-1)
        self.NbminSpinBox.setValue(-1)
        self.NbpaSpinBox.setValue(0)
        self.FFTcheckBox.setChecked(True)
        self.ReducecheckBox.setChecked(False)

        self.MomentmapsgroupBox.setChecked(False)
        self.ProfilecheckBox.setChecked(False)
        self.TotalmapcheckBox.setChecked(False)
        self.VfieldcheckBox.setChecked(False)
        self.DispmapcheckBox.setChecked(False)

        self.MaskgroupBox.setChecked(False)
        self.BlankfactordoubleSpinBox.setValue(2.00)
        self.BlankcutSpinBox.setValue(3.00)

        self.out_path = QDir.currentPath()+"/output/"
        self.OutfolderlineEdit.setText(self.out_path)
        self.stackedWidget.setCurrentIndex(0)
        self.listWidget.setCurrentRow(0)
        
        self.initialize_GUI()


    def enable_All(self): 
        self.ParamlineEdit.setEnabled(True)
        self.FitslineEdit.setEnabled(True)
        self.ParampushButton.setEnabled(True)
        self.ThreadspinBox.setEnabled(True)
        self.Threadslabel.setEnabled(True)
        self.ExtenspinBox.setEnabled(True)
        self.Extenlabel.setEnabled(True)
        self.FitspushButton.setEnabled(True)
        self.ResetpushButton.setEnabled(True)
        self.BoxcheckBox.setEnabled(True)
        self.OutfolderlineEdit.setEnabled(True)
        self.OutfolderpushButton.setEnabled(True)
        self.Outfolderlabel.setEnabled(True)
        self.AdvancedOptcheckBox.setEnabled(True)
        self.AdvancedOptlineEdit.setEnabled(False)
        if self.AdvancedOptcheckBox.isChecked():
            self.AdvancedOptlineEdit.setEnabled(True)
            
        if self.getTaskFlag(1):
            self.GalfitgroupBox.setEnabled(True)
            self.FreeParametersframe.setEnabled(True)
            self.AdvancedgroupBox.setEnabled(True)
            self.Restframe.setEnabled(True)
        
        if self.getTaskFlag(3):
            self.SearchgroupBox.setEnabled(True)
            self.SearchAdvgroupBox.setEnabled(True)
        
        if self.getTaskFlag(4): 
            self.SmoothgroupBox.setEnabled(True)
        
        self.MomentmapsgroupBox.setEnabled(True)
        self.MaskgroupBox.setEnabled(True)
        self.AutocheckBox.setEnabled(True)
    
    
    def disable_All(self): 
        self.FitslineEdit.setDisabled(True)
        self.ParampushButton.setDisabled(True)
        self.ResetpushButton.setDisabled(True)
        self.ThreadspinBox.setDisabled(True)
        self.Threadslabel.setDisabled(True)
        self.ExtenspinBox.setDisabled(True)
        self.Extenlabel.setDisabled(True)
        self.OutfolderpushButton.setDisabled(True)
        self.OutfolderlineEdit.setDisabled(True)
        self.Outfolderlabel.setDisabled(True)
        self.ParamlineEdit.setDisabled(True)
        self.FitspushButton.setDisabled(True)
        self.BoxcheckBox.setDisabled(True)
        self.GalfitgroupBox.setDisabled(True)
        self.AdvancedgroupBox.setDisabled(True)
        self.FreeParametersframe.setDisabled(True)
        self.Restframe.setDisabled(True)
        self.SmoothgroupBox.setDisabled(True)
        self.SearchgroupBox.setDisabled(True)
        self.SearchAdvgroupBox.setDisabled(True)
        self.MomentmapsgroupBox.setDisabled(True)
        self.MaskgroupBox.setDisabled(True)
        self.AutocheckBox.setDisabled(True)
        self.AdvancedOptcheckBox.setDisabled(True)
        self.AdvancedOptlineEdit.setDisabled(True)


    @Slot()
    def AutocheckBox_stateChanged(self):
        if self.AutocheckBox.isChecked():
            self.resetGUI(False)
            self.listWidget.setEnabled(False)
            self.listWidget.setCurrentRow(1)
            self.setTaskFlag(1,Qt.CheckState.Checked)
            self.disable_All()
            self.ExtenspinBox.setEnabled(True)
            self.Extenlabel.setEnabled(True)
            self.BoxcheckBox.setEnabled(True)
        else: 
            self.enable_All()
            self.listWidget.setEnabled(True)
            
        self.FitslineEdit.setEnabled(True)
        self.FitspushButton.setEnabled(True)
        self.AutocheckBox.setEnabled(True)
    
    
    @Slot()
    def HidetoolButton_clicked(self):
        if self.listWidget.isHidden():
            self.listWidget.setHidden(False)
            self.HidetoolButton.setArrowType(Qt.ArrowType.UpArrow)
            self.HidetoolButton.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonIconOnly)
            self.MaingridLayout.addWidget(self.HidetoolButton,0,0)
        else:
            self.listWidget.setHidden(True)
            self.HidetoolButton.setArrowType(Qt.ArrowType.DownArrow)
            self.HidetoolButton.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
            if self.listWidget.currentRow()>=0:
                self.HidetoolButton.setText("   "+self.listWidget.currentItem().text())
            self.MaingridLayout.addWidget(self.HidetoolButton,0,1)
    
    
    @Slot(int)
    def BoxcheckBox_stateChanged(self, state: int):
        toenable = False if state==0 else True
        self.XmaxspinBox.setEnabled(toenable)
        self.XminspinBox.setEnabled(toenable)
        self.YmaxspinBox.setEnabled(toenable)
        self.YminspinBox.setEnabled(toenable)
        self.ZmaxspinBox.setEnabled(toenable)
        self.ZminspinBox.setEnabled(toenable)
        self.Xboxlabel.setEnabled(toenable)
        self.Yboxlabel.setEnabled(toenable)
        self.Zboxlabel.setEnabled(toenable)
     
    
    @Slot(int)
    def AdvancedOptcheckBox_stateChanged(self, state: int):
        toenable = False if state==0 else True
        self.AdvancedOptlineEdit.setEnabled(toenable)


    @Slot()
    def AdvancedOptlineEdit_editingFinished(self):
        
        try:
            pp = readAdvancedOptions(self)
            s = ''
            for k in pp: s += f'{k.upper()} = {pp[k]} \n' 
            msgBox = QMessageBox(icon=QMessageBox.Icon.Information)
            msgBox.setText('Parameters to be added:')
            msgBox.setInformativeText(s)
            ret = msgBox.exec()
        except ValueError:
            msgBox = QMessageBox(icon=QMessageBox.Icon.Critical)
            msgBox.setText("WARNING: Format of advanced is wrong")
            msgBox.setInformativeText("Please semi-colon separated options, i.e. PARAM1=VAL1; PARAM2=VAL2")
            ret = msgBox.exec()
            self.AdvancedOptlineEdit.clear()


    @Slot()
    def actionExport_parameter_file_triggered(self):
        filename, _ = QFileDialog.getSaveFileName(self,"Save Parameter file",\
                                                  QDir.currentPath(),"All files (*.*)")
        self.writeParamFile(filename)
    
    
    @Slot()
    def FitspushButton_clicked(self):
        filename, _ = QFileDialog.getOpenFileName(self,"Open FITSFILE",QDir.currentPath(),\
                            "FITS files (*.fit *.fits *.FITS *.FIT) All files (*.*)")
        if len(filename) and os.path.isfile(filename):
            self.FitslineEdit.setText(filename)
            self.FitslineEdit_editingFinished()
   
    
    @Slot()
    def FitslineEdit_editingFinished(self):
        
        def FITS_Warning(fname):
            msgBox = QMessageBox(icon=QMessageBox.Icon.Critical)
            msgBox.setText("WARNING: "+fname+" is not a readable FITS file.")
            msgBox.setInformativeText("Do you want to select a valid FITS file?")
            msgBox.setStandardButtons(QMessageBox.StandardButton.Open | \
                                      QMessageBox.StandardButton.Cancel)
            msgBox.setDefaultButton(QMessageBox.StandardButton.Open)
            ret = msgBox.exec()
            
            self.FitslineEdit.clear()
            if ret==QMessageBox.StandardButton.Open:
                self.FitspushButton_clicked()


        filename = self.FitslineEdit.text().strip()
        self.lastFilename = ''

        if not filename:
            lastFilename = filename 
            return
        
        if os.path.isfile(filename):
            if filename.find('./')==0:
                filename = filename.replace('./','')
                filename = QDir.currentPath()+filename
        
            self.FitslineEdit.setText(filename)
            
            try:
                f = fits.open(filename)
            except:
                FITS_Warning(filename)
                return
            
            h = f[self.ExtenspinBox.value()].header
            
            if self.lastFilename!=filename: 
                self.XminspinBox.setRange(1,h['NAXIS1'])
                self.XminspinBox.setValue(1)
                self.XmaxspinBox.setRange(1,h['NAXIS1'])
                self.XmaxspinBox.setValue(h['NAXIS1'])
                self.YminspinBox.setRange(1,h['NAXIS2'])
                self.YminspinBox.setValue(1)
                self.YmaxspinBox.setRange(1,h['NAXIS2'])
                self.YmaxspinBox.setValue(h['NAXIS2'])
                self.ZminspinBox.setRange(1,h['NAXIS3'])
                self.ZminspinBox.setValue(1)
                self.ZmaxspinBox.setRange(1,h['NAXIS3'])
                self.ZmaxspinBox.setValue(h['NAXIS3'])
            
            if 'BMAJ' in h: self.ObmajSpinBox.setValue(h['BMAJ']*3600.)
            if 'BMIN' in h: self.ObminSpinBox.setValue(h['BMIN']*3600.)
            if 'BPA'  in h: self.ObpaSpinBox.setValue(h['BPA']*3600.)
            
            self.NbmajSpinBox.setValue(self.FactordoubleSpinBox.value()*self.ObmajSpinBox.value())
            self.NbminSpinBox.setValue(self.FactordoubleSpinBox.value()*self.ObminSpinBox.value())
            self.NbpaSpinBox.setValue(self.ObpaSpinBox.value())
            
            self.obj = 'NONE'
            if 'OBJECT' in h: 
                self.obj = h['OBJECT']
            self.out_path = QDir.currentPath()+"/output/"+self.obj+"/"
            self.OutfolderlineEdit.setText(self.out_path)
            
            self.update_SmoothOutFile()

            if self.lastFilename!=filename:
                if self.RadsepSpinBox.value()==0 and 'BMAJ' in h: 
                    self.RadsepSpinBox.setValue(h['BMAJ']*3600.)
                self.xposlineEdit.setText(str(h['CRPIX1']-1))
                self.yposlineEdit.setText(str(h['CRPIX2']-1))
            
            f.close()
            
        else:
            FITS_Warning(filename)
        
        self.lastFilename = filename
        
    
    @Slot()
    def ParampushButton_clicked(self):
    
        filename, _ = QFileDialog.getOpenFileName(self,"Open parameter file",QDir.currentPath(),\
                                                  "Parameter files (*.par *.txt *.dat)")
        if len(filename) and os.path.isfile(filename):
            self.ParamlineEdit.setText(filename)
            self.ParamlineEdit_editingFinished()


    @Slot()
    def ParamlineEdit_editingFinished(self):
        
        def Param_Warning(fname):
            msgBox = QMessageBox(icon=QMessageBox.Icon.Critical)
            msgBox.setText(f"WARNING: Could not read in parameter file {fname}")
            msgBox.setInformativeText("Do you want to select another file?")
            msgBox.setStandardButtons(QMessageBox.StandardButton.Open | \
                                      QMessageBox.StandardButton.Cancel)
            msgBox.setDefaultButton(QMessageBox.StandardButton.Open)
            ret = msgBox.exec()
        
            self.ParamlineEdit.clear()
            if ret==QMessageBox.StandardButton.Open:
                self.ParampushButton_clicked()


        filename = self.ParamlineEdit.text().strip()
        if not filename: return
        
        if os.path.isfile(filename):
            pp = readParamFromFile(filename)
            try:
                pp = readParamFromFile(filename)
            except:
                Param_Warning(filename)
                return
            self.resetGUI()
            self.ParamlineEdit.setText(filename)
            setGUIParam(self,pp)
        else:
            Param_Warning(filename)
            return


    @Slot()
    def OutfolderpushButton_clicked(self):
    
        outdir = QFileDialog.getExistingDirectory(self,"Open Directory",QDir.currentPath(),\
                            QFileDialog.Option.ShowDirsOnly | QFileDialog.Option.DontResolveSymlinks)
        if len(outdir) and os.path.isdir(outdir):
            outdir += '/'
            self.OutfolderlineEdit.setText(outdir)
            self.out_path = outdir


    ############################################################
    #### Listwidget slots & functions
    ############################################################

    def switch_Galmod_3Dfit(self,toHide: bool):
        self.FreeParametersframe.setHidden(toHide)
        self.IncDSpinBox.setHidden(toHide)
        self.PaDSpinBox.setHidden(toHide)
        self.IncPMlabel.setHidden(toHide)
        self.PaPMlabel.setHidden(toHide)
        self.WfunccomboBox.setHidden(toHide)
        self.FtypecomboBox.setHidden(toHide)
        self.Ftypelabel.setHidden(toHide)
        self.WfunccomboBox.setHidden(toHide)
        self.Wfunclabel.setHidden(toHide)
        self.SecondstagecheckBox.setHidden(toHide)
        self.PolynspinBox.setHidden(toHide)
        self.TollineEdit.setHidden(toHide)
        self.Tollabel.setHidden(toHide)
        self.ErrorscheckBox.setHidden(toHide)
        self.ErrorscheckBox.setHidden(toHide)
        self.AdriftcheckBox.setHidden(toHide)
    
    
    def getCurrentRow(self, item: QListWidgetItem):
        currentRow=0
        itemLab = item.text()
        rows = self.listWidget.model().rowCount()

        for i in range (rows):
            rowLab = self.listWidget.item(i).text()
            if itemLab==rowLab: currentRow=i
        
        return currentRow
    
    
    @Slot(int)
    def listWidget_currentRowChanged(self,currentRow):

        if currentRow==1: self.switch_Galmod_3Dfit(False)
        else:
            self.switch_Galmod_3Dfit(True)
            currentRow -= 1

        self.stackedWidget.setCurrentIndex(currentRow)
    
    
    @Slot('QListWidgetItem*')
    def listWidget_itemChanged(self,item):
        
        currentRow = self.getCurrentRow(item)
        isChecked = False if item.checkState()==Qt.CheckState.Unchecked else True

        if self.RunpushButton.isEnabled():
            if currentRow==1:
                self.setTaskFlag(currentRow,item.checkState())
                self.GalfitgroupBox.setEnabled(isChecked)
                self.FreeParametersframe.setEnabled(isChecked)
                self.Restframe.setEnabled(isChecked)
                self.AdvancedgroupBox.setEnabled(isChecked)
            elif currentRow==2:
                self.setTaskFlag(currentRow,item.checkState())
                self.GalfitgroupBox.setEnabled(isChecked)
                self.AdvancedgroupBox.setEnabled(isChecked)
                self.Restframe.setEnabled(isChecked)
            elif currentRow==3:
                self.SearchgroupBox.setEnabled(isChecked)
                self.SearchAdvgroupBox.setEnabled(isChecked)
            elif currentRow==4:
                self.SmoothgroupBox.setEnabled(isChecked)
            elif currentRow==5:
                self.MomentmapsgroupBox.setEnabled(isChecked)
                self.PVgroupBox.setEnabled(isChecked)
                self.MaskgroupBox.setEnabled(isChecked)
                
        self.listWidget.setCurrentRow(currentRow)


    @Slot('QListWidgetItem*')
    def listWidget_itemDoubleClicked(self,item):
    
        currentRow = self.getCurrentRow(item)

        if currentRow<self.listWidget.model().rowCount()-2 and \
            currentRow>0 and self.RunpushButton.isEnabled(): 
                if item.checkState()==Qt.CheckState.Checked: 
                    item.setCheckState(Qt.CheckState.Unchecked)
                else: 
                    item.setCheckState(Qt.CheckState.Checked)


    def getTaskFlag(self,item: int) -> Qt.CheckState:
        # item = 1: 3DFIT
        # item = 2: GALMOD
        # item = 3: SEARCH
        # item = 4: SMOOTH
        # item = 5: MAPS
        return self.listWidget.item(item).checkState()
    
    
    def setTaskFlag(self,item: int, c: Qt.CheckState):
        
        self.listWidget.item(item).setCheckState(c)
        
        if item==1: #3DFIT
            if c==Qt.CheckState.Checked and self.listWidget.item(2).checkState()==Qt.CheckState.Checked:
                self.listWidget.item(2).setCheckState(Qt.CheckState.Unchecked)
        elif item==2: #GALMOD
            if c==Qt.CheckState.Checked and self.listWidget.item(1).checkState()==Qt.CheckState.Checked:
                self.listWidget.item(1).setCheckState(Qt.CheckState.Unchecked)
    

    ############################################################
    #### 3DFIT panel slots & functions
    ############################################################
    
    @Slot(int)
    def wcscomboBox_activated(self,index):
        
        if index==0:
            self.xposlineEdit.setInputMask("")
            self.yposlineEdit.setInputMask("")
            self.xposlineEdit.setText("0.00")
            self.yposlineEdit.setText("0.00")
        elif index==1:
            self.xposlineEdit.setInputMask("000.0000")
            self.yposlineEdit.setInputMask("#00.00000")
            self.xposlineEdit.setText("000.0000")
            self.yposlineEdit.setText("+00.00000")
        elif index==2:
            self.xposlineEdit.setInputMask("00:00:00.00")
            self.yposlineEdit.setInputMask("#00:00:00.000")
            self.xposlineEdit.setText("000000.00")
            self.yposlineEdit.setText("+000000.000")
        
        self.xposlineEdit.setCursorPosition(1)
        self.yposlineEdit.setCursorPosition(1)
    

    @Slot()
    def MasktoolButton_clicked(self):
        idx = self.MaskcomboBox.currentIndex()
        if idx==0: #SMOOTH
            self.setTaskFlag(5,Qt.CheckState.Checked)
            self.MaskgroupBox.setChecked(True)
            self.listWidget.setCurrentRow(5)
        elif idx==1: #SEARCH
            self.setTaskFlag(3,Qt.CheckState.Checked)
            self.SearchgroupBox.setEnabled(True)
            self.listWidget.setCurrentRow(3)
        elif idx==2: #SMOOTH&SEARCH
            self.setTaskFlag(3,Qt.CheckState.Checked)
            self.setTaskFlag(5,Qt.CheckState.Checked)
            self.MaskgroupBox.setChecked(True)
            self.SearchgroupBox.setEnabled(True)
            self.listWidget.setCurrentRow(5)

    
    @Slot(int)
    def MaskcomboBox_currentIndexChanged(self,index):
        if index==0 or index==1 or index==2:    #SMOOTH, SEARCH, SMOOTH&SEARCH
            self.MasktoolButton.setHidden(False)
            self.MaskThreshSpinBox.setHidden(True)
        elif index==3:  #THRESHOLD
            self.MasktoolButton.setHidden(True)
            self.MaskThreshSpinBox.setHidden(False)
        else:   #NEGATIVE, NONE
            self.MaskThreshSpinBox.setHidden(True)
            self.MasktoolButton.setHidden(True)
    
    
    def Selected_3DFit_file(self, le: QLineEdit, sb: QSpinBox, qw: QWidget, qw2=None):

        filename, _ = QFileDialog.getOpenFileName(self, "Open file", QDir.currentPath(), \
                      "All files (*.*)", "")

        if len(filename) and os.path.isfile(filename):
            le.setText(filename)
            Hide_3DFit_file(le,sb,qw,False)
            if qw2:
                Hide_3DFit_file(le,sb,qw2,False)
            

    def Hide_All_3DFit_file(self, Hide: bool):
        Hide_3DFit_file(self.RadsepFilelineEdit,self.RadsepFilespinBox,self.RadsepSpinBox,Hide)
        Hide_3DFit_file(self.XposFilelineEdit,self.XposFilespinBox,self.xposlineEdit,Hide)
        Hide_3DFit_file(self.YposFilelineEdit,self.YposFilespinBox,self.yposlineEdit,Hide)
        Hide_3DFit_file(self.VsysFilelineEdit,self.VsysFilespinBox,self.VsysSpinBox,Hide)
        Hide_3DFit_file(self.VrotFilelineEdit,self.VrotFilespinBox,self.VrotSpinBox,Hide)
        Hide_3DFit_file(self.VdispFilelineEdit,self.VdispFilespinBox,self.VdispSpinBox,Hide)
        Hide_3DFit_file(self.VradFilelineEdit,self.VradFilespinBox,self.VradSpinBox,Hide)
        Hide_3DFit_file(self.IncFilelineEdit,self.IncFilespinBox,self.IncSpinBox,Hide)
        Hide_3DFit_file(self.PaFilelineEdit,self.PaFilespinBox,self.PaSpinBox,Hide)
        Hide_3DFit_file(self.Z0FilelineEdit,self.Z0FilespinBox,self.Z0SpinBox,Hide)
        Hide_3DFit_file(self.DensFilelineEdit,self.DensFilespinBox,self.DensSpinBox,Hide)
    
    
    ############################################################
    #### SEARCH panel slots & functions
    ############################################################
    
    @Slot()
    def SearchgroupBox_clicked(self):
        if self.SearchgroupBox.isChecked(): 
            self.SearchAdvgroupBox.setEnabled(True)
        else:
            self.SearchAdvgroupBox.setDisabled(True)


    @Slot()
    def GrowthcheckBox_stateChanged(self):
        if self.GrowthcheckBox.isChecked():
            self.SecondarycutlineEdit.setEnabled(True)
            self.Cuttype2comboBox.setEnabled(True)
        else:
            self.SecondarycutlineEdit.setDisabled(True)
            self.Cuttype2comboBox.setDisabled(True)
    

    @Slot()
    def SearchAdvgroupBox_clicked(self):
        if self.SearchAdvgroupBox.isChecked():
            self.ThreshSpatialBox.setDisabled(True)
            self.RejectcheckBox.setChecked(True)
        else:
            self.AdjacentcheckBox.setChecked(True)
            self.RejectcheckBox.setChecked(True)
            self.TwostagemergingcheckBox.setChecked(True)


    @Slot(int)
    def AdjacentcheckBox_stateChanged(self, arg1: int):
        if (arg1==2):
            self.ThreshSpatialBox.setValue(-1)
            self.ThreshSpatialBox.setDisabled(True)
        else:
            self.ThreshSpatialBox.setEnabled(True)
            self.ThreshVelspinBox.setEnabled(True)
    

    ############################################################
    #### SMOOTH panel slots & functions
    ############################################################

    @Slot()
    def FactorcheckBox_stateChanged(self):
        if self.FactorcheckBox.isChecked():
            self.FactordoubleSpinBox.setEnabled(True)
            self.BeamgroupBox.setChecked(False)
            self.ReducecheckBox.setEnabled(True)
            self.FactordoubleSpinBox_valueChanged(self.FactordoubleSpinBox.value())
        else:
            self.FactordoubleSpinBox.setDisabled(True)
            self.BeamgroupBox.setChecked(True)
            self.ReducecheckBox.setChecked(False)
            self.ReducecheckBox.setDisabled(True)


    @Slot(int)
    def HanningspinBox_valueChanged(self, arg1: int):
        if arg1 % 2 == 0: 
            b = QMessageBox(icon=QMessageBox.Icon.Critical)
            b.setText("Only odd values allowed!")
            b.exec()
            self.HanningspinBox.setValue(arg1+1)
    
    
    @Slot()
    def SmoothOutpushButton_clicked(self):
        filename, _ = QFileDialog.getSaveFileName(self,"Save FITS file",QDir.currentPath(),\
                                                  "FITS files (*.fits *.FITS *.fit *.FIT)")
        if len(filename): 
            self.SmoothOutlineEdit.setText(filename)
    
    
    @Slot(float)
    def FactordoubleSpinBox_valueChanged(self, arg1: float):
        self.NbmajSpinBox.setValue(arg1*self.ObmajSpinBox.value())
        self.NbminSpinBox.setValue(arg1*self.ObminSpinBox.value())
        self.NbpaSpinBox.setValue(self.ObpaSpinBox.value())
    

    @Slot()
    def update_SmoothOutFile(self):
        if self.obj!='':
            smooth_name = f'{self.out_path}{self.obj}/{self.obj}_s{int(round(self.NbmajSpinBox.value()))}'
            if self.ReducecheckBox.isChecked(): 
                smooth_name+="red"
            smooth_name+=".fits"
            self.SmoothOutlineEdit.setText(smooth_name)


    ############################################################
    #### Execution slots & functions
    ############################################################
    
    @Slot()
    def KillpushButton_clicked(self):
        if self.proc and self.proc.state()==2:
            self.proc.terminate()
            QTimer.singleShot(2000, self.proc, self.proc.kill())
            self.LogtextEdit.append("\n\n ---------. Killed by the user. <----------\n")
            self.RunpushButton.setEnabled(True)
            self.KillpushButton.setEnabled(False)
    
    
    @Slot()
    def updateText(self):
        text = bytes(self.proc.readAll()).decode()
        self.LogtextEdit.insertPlainText(text)
        self.LogtextEdit.moveCursor(QTextCursor.MoveOperation.End)
        if self.getTaskFlag(1): self.plotParameters()

    
    @Slot()
    def updateExit(self):
    
        self.RunpushButton.setEnabled(True)
        self.KillpushButton.setEnabled(False)
        self.AutocheckBox.setCheckState(Qt.CheckState.Unchecked)
        self.enable_All()
        
        if self.proc.exitStatus()==QProcess.ExitStatus.NormalExit:
            if self.proc.exitCode()==0:
                self.statusBar.showMessage("  BBarolo successfully terminated")
            else:
                self.statusBar.showMessage("  BBarolo exited with errors")
        else:
            self.statusBar.showMessage("  BBarolo has been killed or crashed")
        
        outfile = QFile(self.out_path+"output.log")
        
        if outfile.open(QIODevice.OpenModeFlag.ReadWrite | \
                        QIODevice.OpenModeFlag.Truncate | \
                        QIODevice.OpenModeFlag.Text):
            out = QTextStream(outfile)
            out << self.LogtextEdit.toPlainText()
            self.proc.deleteLater()
    
    
    @Slot()
    def RunpushButton_clicked(self):
        
        self.listWidget.setCurrentRow(self.listWidget.model().rowCount()-2)
        if self.FitslineEdit.text()=='':
            self.LogtextEdit.append("ERROR: You must select at least a FITS file!")
            return
        
        self.RunpushButton.setDisabled(True)
        self.KillpushButton.setEnabled(True)
        self.disable_All()

        #self.plot1.replot()
        outfolder = self.OutfolderlineEdit.text().strip()
        if outfolder.find('~')==0:
            outfolder = os.path.expanduser('~')+outfolder[1:] 
        if outfolder[-1]!='/':
            outfolder += '/'
        
        self.out_path = outfolder
        if not os.path.exists(outfolder):
            os.makedirs(outfolder)
        
        for i in (1,2):
            if os.path.isfile(outfolder+f"rings_final{i}.txt"):
                os.remove(outfolder+f"rings_final{i}.txt")
        
        argv = []
        if self.AutocheckBox.isChecked():
            argv.append('-f')
            argv.append(getFilenameString(self))
        else:
            paramfile = outfolder+self.obj+".par"
            self.writeParamFile(paramfile)
            argv.append('-p') 
            argv.append(paramfile)
        
        self.LogtextEdit.clear()
        self.statusBar.showMessage("  BBarolo is running ...")

        self.proc = QProcess(self)
        self.proc.readyRead.connect(self.updateText)
        self.proc.readyReadStandardError.connect(self.updateText)
        self.proc.finished.connect(self.updateExit)
        self.proc.start(self.BBexe,argv)
    
    
    def writeParamFile(self, fileout: str):
        pp = readGUIParam(self)
        
        # Writing to text file
        s = "##### Input parameters for BBarolo #####\n"
        for k in pp:
            s += f'{k.upper():18s} {pp[k]} \n' 
        
        with open(fileout,'w') as f:
            f.write(s)

    
    def plotParameters(self):
        
        ...
        '''
        series = QLineSeries()
        
        for i in range(0,100):
            series.append(i, i)

        chart = QChart()
        chartView = QChartView(chart)
        

        
        self.Plot3graphicsView.setHidden(True)
        self.verticalLayout_4.addWidget(chartView,0)
        
        chart.addSeries(series)
        chart.createDefaultAxes();
        chart.setTitle("Simple Line Chart");
        ###################################
        '''
        
        '''    
        scene = QGraphicsScene()
        self.Plot3graphicsView.setScene(scene)

        fig, ax = plt.subplots(figsize=(2.5,1.5))

        ax.set_xlabel("Radius (kpc)", fontsize=8)
        ax.set_ylabel("Velocity (km/s)", fontsize=8)
        ax.grid()
        xx = [1,2,3,4,5]
        ax.plot(xx,xx)
        fig.tight_layout()
        canvas = mpl.backends.backend_qtagg.FigureCanvasQTAgg(fig)
        proxy_widget = scene.addWidget(canvas)
        #canvas.draw()
        #'''


    def plotParameters(self):
        
        return
        nsubplot = 2 if self.SecondstagecheckBox.isChecked() else 1
        labels  = ['Radius (arcs)','Vrot (km/s)','Disp (km/s)','Inc (deg)','P.A. (deg)',\
                   'Z0 (arcs)','Xpos (pix)','Ypos (pix)','Vsys (km/s)']
        
        for j in range(nsubplot):
            
            ringfile = f'{self.out_path}rings_final{j+1}.txt'
            
            try:
                data = [[float(x) for x in l.split()] for l in open(ringfile).readlines()[1:]]
                data_tuples = list(zip(*data))
                data = [list(sublist) for sublist in data_tuples]
            except:
                continue
            
            nr = len(data[0])
            
            plots = [self.plot1,self.plot2,self.plot3,self.plot4]
            index = [self.Plot1comboBox.currentIndex()+1, self.Plot2comboBox.currentIndex()+1,\
                     self.Plot3comboBox.currentIndex()+1, self.Plot4comboBox.currentIndex()+1]
            
            pen = QPen()
            pen.setWidth(2)
            if j==1 or nsubplot==1: 
                pen.setColor(QColor(178,34,34))
            else: 
                pen.setColor(QColor(137,127,127))

            for i in range(4): 
                plots[i].addGraph()
                plots[i].graph(j).setData(data[1],data[index[i]])
                plots[i].graph(j).setPen(pen)
                plots[i].graph(j).setLineStyle(QCPGraph.lsNone)
                plots[i].graph(j).setScatterStyle(QCPScatterStyle(QCPScatterStyle.ssPlusCircle, 4))
                plots[i].xAxis.setLabel(labels[0])
                plots[i].yAxis.setLabel(labels[index[i]])
                plots[i].xAxis.setRange(0, 1.10*max(data[1]))
                plots[i].yAxis.setRange(0.95*min(data[index[i]]), 1.05*max(data[index[i]]))

        self.plot1.replot()
        self.plot2.replot()
        self.plot3.replot()
        self.plot4.replot()
