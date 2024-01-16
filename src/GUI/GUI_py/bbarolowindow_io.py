"""
Utility functions to read/write from/to the GUI
"""
from imports import *


readFlag = lambda flag: True if 'true' in flag.lower() else False


def isNumber(string):
    try:
        float (string)
        return True
    except ValueError:
        return False


def getFilenameString(ui):
    
    fname = ui.FitslineEdit.text().strip()+f'[{int(ui.ExtenspinBox.value())}]'
    if ui.BoxcheckBox.isChecked():
        boxstr = f'[{int(ui.XminspinBox.value())}:{int(ui.XmaxspinBox.value())},'
        boxstr += f'{int(ui.YminspinBox.value())}:{int(ui.YmaxspinBox.value())},'
        boxstr += f'{int(ui.ZminspinBox.value())}:{int(ui.ZmaxspinBox.value())}]'
        fname += boxstr
    return fname


def getFileString(le, sb):
    return "file("+le.text()+","+str(int(sb.value()))+")"


def readFileString(s: str) -> int:
    if "file" not in s.lower(): 
        return (None, -1)
    else:
        str_list = s.split(',')
        str_list[0] = str_list[0][str_list[0].find('(')+1:]
        str_list = [ll.replace(')', '') for ll in str_list]
        fname = str_list[0]
        ncol = int(str_list[1]) if len(str_list)>1 else 1
        return (fname, ncol)


def Hide_3DFit_file(le: QLineEdit, sb: QSpinBox, qw: QWidget, Hide: bool, qw2: QWidget = None):
    le.setHidden(Hide)
    sb.setHidden(Hide)
    qw.setEnabled(Hide)
    if qw2: qw2.setEnabled(Hide)


def setModelParam(p: dict, pname: str, cb: QCheckBox, sb: QWidget, fle: QLineEdit, fsb: QSpinBox, wcscb: QComboBox = None):
    
    if pname in p:           # Parameter is given
        cb.setChecked(False)
        if wcscb: wcscb.setCurrentIndex(0)
        s = p[pname][0].strip()
        fname, col = readFileString(s)
        if fname is None:
            # Special case for XPOS and YPOS (check for WCS coordinates)
            if 'POS' in pname:  
                found = s.find('d')
                if found!=-1:
                    s = s[:found]
                    wcscb.setCurrentIndex(1)
                else:
                    if s.find(':')!=-1: 
                        wcscb.setCurrentIndex(2)
                sb.setText(s)
            else:
                # All other parameters, just set the value
                sb.setValue(float(s))
        else:
            Hide_3DFit_file(fle,fsb,sb,False)
            fle.setText(fname)
            fsb.setValue(col)
    else:                           # Parameter not given, ie. to estimate
        cb.setChecked(True)
    

def readAdvancedOptions(ui, p: dict = None, raiseError: bool = True) -> dict:
    
    if p is None: p = dict()
    
    advopt = ui.AdvancedOptlineEdit.text().strip()
    if ui.AdvancedOptcheckBox.isChecked() and advopt!='':
        advopt = advopt.split(';')
        for i in range(len(advopt)):
            try:
                par, value = advopt[i].strip().split('=')
                p[par.upper().strip()] = value.strip()
            except ValueError:
                if raiseError:
                    raise ValueError
                else:
                    pass
    return p
    

def readGUIModelParam(ui, p: dict = None) -> dict:
    
    if p is None: 
        p = dict()
    
    if not ui.RadsepFilelineEdit.isHidden():
        p['RADII'] = getFileString(ui.RadsepFilelineEdit,ui.RadsepFilespinBox)
    else:
        if not ui.NringscheckBox.isChecked():
            p['NRADII'] = ui.NringsspinBox.value()
        if not ui.RadsepcheckBox.isChecked():
            p['RADSEP'] = ui.RadsepSpinBox.value() 
        
    if not ui.XposcheckBox.isChecked():
        if ui.BoxcheckBox.isChecked() and ui.wcscomboBox.currentIndex()==0:
            xp = float(ui.xposlineEdit.text())-float(ui.XminspinBox.value())+1
            p['XPOS'] = xp
        else: 
            p['XPOS'] = ui.xposlineEdit.text()
        if ui.wcscomboBox.currentIndex()==1: p['XPOS'] += 'd'
    
    if not ui.YposcheckBox.isChecked():
        if ui.BoxcheckBox.isChecked() and ui.wcscomboBox.currentIndex()==0:
            yp = float(ui.yposlineEdit.text())-float(ui.YminspinBox.value())+1
            p['YPOS'] = yp
        else:
            p['YPOS'] = ui.yposlineEdit.text()
        if ui.wcscomboBox.currentIndex()==1: p['YPOS'] += 'd'
    
    if not ui.VsyscheckBox.isChecked():
        p['VSYS'] = ui.VsysSpinBox.value() if ui.VsysFilelineEdit.isHidden()\
                    else getFileString(ui.VsysFilelineEdit,ui.VsysFilespinBox)

    if not ui.VrotcheckBox.isChecked():
        p['VROT'] = ui.VrotSpinBox.value() if ui.VrotFilelineEdit.isHidden()\
                    else getFileString(ui.VrotFilelineEdit,ui.VrotFilespinBox)

    if not ui.VdispcheckBox.isChecked():
        p['VDISP'] = ui.VdispSpinBox.value() if ui.VdispFilelineEdit.isHidden()\
                     else getFileString(ui.VdispFilelineEdit,ui.VdispFilespinBox)
    
    if not ui.VradcheckBox.isChecked():
        p['VRAD'] = ui.VradSpinBox.value() if ui.VradFilelineEdit.isHidden()\
                    else getFileString(ui.VradFilelineEdit,ui.VradFilespinBox)

    if not ui.InccheckBox.isChecked():
        p['INC'] = ui.IncSpinBox.value() if ui.IncFilelineEdit.isHidden()\
                   else getFileString(ui.IncFilelineEdit,ui.IncFilespinBox)
        p['DELTAINC'] = ui.IncDSpinBox.value()
    
    if not ui.PacheckBox.isChecked():
        p['PA'] = ui.PaSpinBox.value() if ui.PaFilelineEdit.isHidden()\
                  else getFileString(ui.PaFilelineEdit,ui.PaFilespinBox)
        p['DELTAPA'] = ui.PaDSpinBox.value()

    if not ui.Z0checkBox.isChecked():
        p['Z0'] = ui.Z0SpinBox.value() if ui.Z0FilelineEdit.isHidden()\
                  else getFileString(ui.Z0FilelineEdit,ui.Z0FilespinBox)

    if not ui.DenscheckBox.isChecked():
        p['DENS'] = ui.DensSpinBox.value() if ui.DensSpinBox.value()\
                    else getFileString(ui.DensFilelineEdit,ui.DensFilespinBox)
    
    if float(ui.redshiftlineEdit.text())!=0:
        p['REDSHIFT'] = ui.redshiftlineEdit.text()

    if float(ui.restlineEdit.text())!=0:
        if ui.restcomboBox.currentIndex()==0:
            p['RESTWAVE'] = ui.restlineEdit.text()
        if ui.restcomboBox.currentIndex()==1: 
            p['RESTFREQ'] = ui.restlineEdit.text()

    if ui.AdvancedgroupBox.isChecked():
        p['LTYPE'] = ui.LtypecomboBox.currentIndex()+1
        p['CDENS'] = ui.CdensSpinBox.value()
        if ui.NVspinBox.value()!=-1: p['NV'] = ui.NVspinBox.value()
        
        p['MASK'] = ui.MaskcomboBox.currentText()
        if 'THRESH' in p['MASK']:
            p['THRESHOLD'] = ui.MaskThreshSpinBox.value()
        
        p['NORM'] = ui.NormcomboBox.currentText()

    if ui.AdriftcheckBox.isChecked():
        p['ADRIFT'] = True
    
    return p
    

def readGUIParam(ui, p: dict = None) -> dict:
    
    if p is None: 
        p = dict()

    p['FITSFILE'] = getFilenameString(ui)
    
    
    p['SHOWBAR'] = False
    p['THREADS'] = int(ui.ThreadspinBox.value())
    outdir = ui.OutfolderlineEdit.text().strip()
    if len(outdir)!=0: p['OUTFOLDER'] = outdir 

    # Reading GUI parameters for 3DFIT and GALMOD
    if ui.getTaskFlag(1)==Qt.CheckState.Checked:
                
        p['3DFIT'] = True
        p = readGUIModelParam(ui,p)
        
        free = ""
        if ui.vrotcheckBox.isChecked(): free += "VROT "
        if ui.vdispcheckBox.isChecked():free += "VDISP "
        if ui.inccheckBox.isChecked():  free += "INC "
        if ui.pacheckBox.isChecked():   free += "PA "
        if ui.xposcheckBox.isChecked(): free += "XPOS "
        if ui.yposcheckBox.isChecked(): free += "YPOS "
        if ui.vsyscheckBox.isChecked(): free += "VSYS "
        if ui.z0checkBox.isChecked():   free += "Z0 "
        if ui.vradcheckBox.isChecked(): free += "VRAD "
        p['FREE'] = free

        if ui.AdvancedgroupBox.isChecked():
            
            p['FTYPE'] = ui.FtypecomboBox.currentIndex()+1
            p['WFUNC'] = ui.WfunccomboBox.currentIndex()
            p['TOL'] = ui.TollineEdit.text()
            
            p['TWOSTAGE'] = ui.SecondstagecheckBox.isChecked()
            if ui.SecondstagecheckBox.isChecked():
                p['REGTYPE'] = 'bezier' if ui.PolynspinBox.value()==-1 \
                               else ui.PolynspinBox.value()
                               
            p['FLAGERRORS'] = ui.ErrorsradioButton.isChecked()

    else:
        if ui.getTaskFlag(2)==Qt.CheckState.Checked: 
            p['GALMOD'] =  True
            p = readGUIModelParam(ui,p)


    # Reading GUI parameters for SEARCH
    if ui.getTaskFlag(3)==Qt.CheckState.Checked:
        p['SEARCH'] = True
        p['SEARCHTYPE'] = ui.SearchtypecomboBox.currentText()

        cutkey = 'THRESHOLD' if ui.CuttypecomboBox.currentIndex()==0 else 'SNRCUT'
        p[cutkey] = ui.primaryCutlineEdit.text()

        if ui.GrowthcheckBox.isChecked():
            p['FLAGGROWTH'] = True
            cutkey = 'GROWTHTHRESHOLD' if ui.Cuttype2comboBox.currentIndex()==0 \
                     else 'GROWTHCUT'
            p[cutkey] = ui.SecondarycutlineEdit.text()

        if ui.SearchAdvgroupBox.isChecked():
            p['FLAGADJACENT'] = ui.AdjacentcheckBox.isChecked()
            if not p['FLAGADJACENT']:
                if ui.ThreshSpatialBox.value()!=-1:
                    p['THRESHSPATIAL'] = ui.ThreshSpatialBox.value()
                if ui.ThreshVelspinBox.value()!=-1:
                    p['THRESHVELOCITY'] = ui.ThreshVelspinBox.value()

        if ui.MinPixspinBox.value()!=-1:
            p['MINPIX'] = ui.MinPixspinBox.value()
        if ui.MinChanspinBox.value()!=-1:
            p['MINCHANNELS'] = ui.MinChanspinBox.value()
        if ui.MinVoxspinBox.value()!=-1:
            p['MINVOXELS'] = ui.MinVoxspinBox.value()
        if ui.MaxChanspinBox.value()!=-1:
            p['MAXCHANNELS'] = ui.MaxChanspinBox.value()
        if ui.MaxAngSizeSpinBox.value()!=-1:
            p['MAXANGSIZE'] = ui.MaxAngSizeSpinBox.value()
            
        p['REJECTBEFOREMERGE'] = ui.RejectcheckBox.isChecked()
        p['TWOSTAGEMERGING'] = ui.TwostagemergingcheckBox.isChecked()


    # Reading GUI parameters for SMOOTH
    if ui.getTaskFlag(4)==Qt.CheckState.Checked: 
        if ui.SpatialSmoothgroupBox.isChecked():
            p['SMOOTH'] = True
            if ui.FactorcheckBox.isChecked():
                p['FACTOR'] = ui.FactordoubleSpinBox.value()
            else:
                p['BMAJ']  =  ui.NbmajSpinBox.value() 
                p['BMIN']  = ui.NbminSpinBox.value() 
                p['BPA']   = ui.NbpaSpinBox.value() 
                p['OBMAJ'] =  ui.ObmajSpinBox.value()
                p['OBMIN'] =  ui.ObminSpinBox.value()
                p['OBPA']  =  ui.ObpaSpinBox.value()
            
            p['REDUCE'] = ui.ReducecheckBox.isChecked()
            p['FFT'] = ui.FFTcheckBox.isChecked()
            if ui.SmoothOutlineEdit.text()!='':
                p['SMOOTHOUTPUT'] = ui.SmoothOutlineEdit.text()
    
        if ui.HanninggroupBox.isChecked():
            p['SMOOTHSPEC'] = True
            p['WINDOW_TYPE'] = 'HANNING'
            p['WINDOW_SIZE'] = ui.HanningspinBox.value()


    # Reading Maps and PV parameters
    if ui.MomentmapsgroupBox.isChecked(): 
        p['GLOBALPROFILE'] = ui.ProfilecheckBox.isChecked()
        p['TOTALMAP']      = ui.TotalmapcheckBox.isChecked()
        p['VELOCITYMAP']   = ui.VfieldcheckBox.isChecked()
        p['DISPERSIONMAP'] = ui.DispmapcheckBox.isChecked()
        p['RMSMAP']        = ui.rmsmapcheckBox.isChecked()

    if ui.PVgroupBox.isChecked():
        p['FLAGPV']  = True
        p['XPOS_PV'] = ui.PVXposSpinBox.value()
        p['YPOS_PV'] = ui.PVYposSpinBox.value()
        p['PA_PV']   = ui.PVPaSpinBox.value()

    if ui.MaskgroupBox.isChecked():
        p['MASK'] = 'SMOOTH'
        p['FACTOR'] = ui.BlankfactordoubleSpinBox.value()
        p['BLANKCUT'] = ui.BlankcutSpinBox.value()
    else:
        p['MASK'] = 'NONE'
    
    # Reading in additional advanced parameters if given
    p = readAdvancedOptions(ui,p,False)
    
    return p
    

def readParamFromFile(filein: str) -> dict:
    
    # Reading in parameter file
    par = dict()
    for line in open(filein):
        l = line.strip()
        if not l: continue
        
        pname, *pvalues = l.split()
        pvalue = ' '.join(pvalues)
        
        # Getting rid of lines commented out or invalid
        if any(x in pname[0] for x in ['#','/','%']) or len(pvalues)==0:
            continue
            
        # Taking care of inline comments
        if pvalue.find('//')!=-1:
            pvalue = pvalue[:pvalue.find('//')]
        if pvalue.find('#')!=-1:
            pvalue = pvalue[:pvalue.find('#')]
        par[pname.upper()] = pvalue.split()
    
    return par


def setGUIModelParam(ui, p: dict):
        
    ui.FreeParametersframe.setEnabled(True)
    
    if 'RADII' in p:
        fname, col = readFileString(s)
        if fname:
            Hide_3DFit_file(ui.RadsepFilelineEdit,ui.RadsepFilespinBox,ui.RadsepSpinBox,False)
            ui.RadsepFilelineEdit.setText(fname)
            ui.RadsepFilespinBox.setValue(col)
    else:
        if 'NRADII' in p: 
            ui.NringscheckBox.setChecked(False)
            ui.NringsspinBox.setValue(int(p['NRADII'][0]))
        else:
            ui.NringscheckBox.setChecked(True)
            
        if 'RADSEP' in p:
            ui.RadsepcheckBox.setChecked(False)
            ui.RadsepSpinBox.setValue(float(p['RADSEP'][0]))
        else:
            ui.RadsepcheckBox.setChecked(True)

    setModelParam(p,'XPOS',ui.XposcheckBox,ui.xposlineEdit,ui.XposFilelineEdit,ui.XposFilespinBox,ui.wcscomboBox)
    setModelParam(p,'YPOS',ui.YposcheckBox,ui.yposlineEdit,ui.YposFilelineEdit,ui.YposFilespinBox,ui.wcscomboBox)
    setModelParam(p,'VSYS',ui.VsyscheckBox,ui.VsysSpinBox,ui.VsysFilelineEdit,ui.VsysFilespinBox)
    setModelParam(p,'VROT',ui.VrotcheckBox,ui.VrotSpinBox,ui.VrotFilelineEdit,ui.VrotFilespinBox)
    setModelParam(p,'VDISP',ui.VdispcheckBox,ui.VdispSpinBox,ui.VdispFilelineEdit,ui.VdispFilespinBox)
    setModelParam(p,'VRAD',ui.VradcheckBox,ui.VradSpinBox,ui.VradFilelineEdit,ui.VradFilespinBox)
    setModelParam(p,'INC',ui.InccheckBox,ui.IncSpinBox,ui.IncFilelineEdit,ui.IncFilespinBox)
    setModelParam(p,'PA',ui.PacheckBox,ui.PaSpinBox,ui.PaFilelineEdit,ui.PaFilespinBox)
    setModelParam(p,'Z0',ui.Z0checkBox,ui.Z0SpinBox,ui.Z0FilelineEdit,ui.Z0FilespinBox)
    setModelParam(p,'DENS',ui.DenscheckBox,ui.DensSpinBox,ui.DensFilelineEdit,ui.DensFilespinBox)
    
    ui.Restframe.setEnabled(True)
    if 'REDSHIFT' in p: 
        ui.redshiftlineEdit.setText(p['REDSHIFT'][0])
    if 'RESTWAVE' in p:
        ui.restcomboBox.setCurrentIndex(0)
        ui.restlineEdit.setText(p['RESTWAVE'][0])
    if 'RESTFREQ' in p:
        ui.restcomboBox.setCurrentIndex(1)
        ui.restlineEdit.setText(p['RESTFREQ'][0])

    ui.AdvancedgroupBox.setChecked(True)
    ui.AdvancedgroupBox.setEnabled(True)
    if 'LTYPE' in p:
        ui.LtypecomboBox.setCurrentIndex(int(p['LTYPE'][0])-1)
    if 'CDENS' in p:
        ui.CdensSpinBox.setValue(float(p['CDENS'][0]))
    if 'NV' in p:
        ui.NVspinBox.setValue(float(p['NV'][0]))

    if 'MASK' in p:
        mask = p['MASK'][0].strip().upper()
        if mask=='SMOOTH':
            ui.MaskcomboBox.setCurrentIndex(0)
            ui.MaskgroupBox.setChecked(True)
            if 'FACTOR' in p:
                ui.BlankfactordoubleSpinBox.setValue(float(p['FACTOR'][0]))
            if 'SNRCUT' in p:
                ui.BlankcutSpinBox.setValue(float(p['SNRCUT'][0]))

        elif 'SEARCH' in mask: 
            ui.MaskcomboBox.setCurrentIndex(1)
        elif 'SMOOTH&SEARCH' in mask:
            mask: ui.MaskcomboBox.setCurrentIndex(2)
        elif 'THRESHOLD' in mask: 
            ui.MaskcomboBox.setCurrentIndex(3)
            if 'THRESHOLD' in p:
                ui.MaskThreshSpinBox.setValue(float(p['THRESHOLD'][0]))
        elif 'NEGATIVE' in mask: 
            ui.MaskcomboBox.setCurrentIndex(4)
        else: 
            ui.MaskcomboBox.setCurrentIndex(5)
    
    if 'NORM' in p:
        norm = p['NORM'][0].strip().upper()
        if 'LOCAL' in norm:
            ui.NormcomboBox.setCurrentIndex(0)
        elif 'AZIM' in norm:
            ui.NormcomboBox.setCurrentIndex(1)
        else: 
            ui.NormcomboBox.setCurrentIndex(2)

    if 'ADRIFT' in p:
        ui.AdriftcheckBox.setChecked(readFlag(p['ADRIFT'][0]))
    

def setGUIParam(ui, p: dict):
    
    # Now reading all parameters and set the GUI elemets
    
    if 'THREADS' in p:
        ui.ThreadspinBox.setValue(int(p['THREADS'][0]))

    if 'FITSFILE' in p:
        filename = p['FITSFILE'][0]
        a, b = filename.find("[")!=-1, filename.find("]")!=-1
        ui.FitslineEdit.setText("")
        ui.FitslineEdit_editingFinished()
        
        if a and b: # There is a BOX string to decipher
            ui.BoxcheckBox.setChecked(True)
            ran = filename[filename.find("[")+1:filename.find("]")]
            ran = ran.split(',')
            filename = filename[:filename.find("[")]
            ui.FitslineEdit.setText(filename)
            ui.FitslineEdit_editingFinished()
            if not '*' in ran[0]:
                xran = ran[0].split(":")
                ui.XminspinBox.setValue(int(xran[0]))
                ui.XmaxspinBox.setValue(int(xran[1]))
            if not '*' in ran[1]:
                yran = ran[1].split(":")
                ui.YminspinBox.setValue(int(yran[0]))
                ui.YmaxspinBox.setValue(int(yran[1]))
            if not '*' in ran[2] and len(ran)>2:
                zran = ran[2].split(":")
                ui.ZminspinBox.setValue(int(zran[0]))
                ui.ZmaxspinBox.setValue(int(zran[1]))
        else:
            ui.FitslineEdit.setText(filename)
            ui.FitslineEdit_editingFinished()
    
    if 'OUTFOLDER' in p:
        ui.OutfolderlineEdit.setText(p['OUTFOLDER'][0])
    
    # Setting GUI elements for 3DFIT and GALMOD tasks
    if ('3DFIT' in p or 'GALFIT' in p) and readFlag(p['3DFIT'][0]):
        
        ui.setTaskFlag(1,Qt.CheckState.Checked)
        setGUIModelParam(ui,p)
        
        if 'DELTAINC' in p:
            ui.IncDSpinBox.setValue(float(p['DELTAINC'][0]))
        if 'DELTAPA' in p or 'DELTAPHI' in p:
            ui.PaDSpinBox.setValue(float(p['DELTAINC'][0]))

        if 'FREE' in p:
            f = ' '.join(p['FREE']).lower()
            
            ui.vrotcheckBox.setChecked(False)
            ui.vdispcheckBox.setChecked(False)
            ui.inccheckBox.setChecked(False)
            ui.pacheckBox.setChecked(False)
            ui.xposcheckBox.setChecked(False)
            ui.yposcheckBox.setChecked(False)
            ui.vsyscheckBox.setChecked(False)
            ui.z0checkBox.setChecked(False)
            ui.vradcheckBox.setChecked(False)
            
            if 'vrot' in f: ui.vrotcheckBox.setChecked(True)
            if 'disp' in f: ui.vdispcheckBox.setChecked(True)
            if 'inc'  in f: ui.inccheckBox.setChecked(True)
            if 'pa' in f or 'phi' in f: ui.pacheckBox.setChecked(True)
            if 'xpos' in f: ui.xposcheckBox.setChecked(True)
            if 'ypos' in f: ui.yposcheckBox.setChecked(True)
            if 'vsys' in f: ui.vsyscheckBox.setChecked(True)
            if 'z0'   in f: ui.z0checkBox.setChecked(True)
            if 'vrad' in f: ui.vradcheckBox.setChecked(True)
        
        if 'FTYPE' in p:
            ui.FtypecomboBox.setCurrentIndex(int(p['FTYPE'][0])-1)
        if 'WFUNC' in p:
            ui.WfunccomboBox.setCurrentIndex(int(p['WFUNC'][0]))
                
        if 'TWOSTAGE' in p:
            ui.SecondstagecheckBox.setChecked(readFlag(p['TWOSTAGE'][0]))
            if ui.SecondstagecheckBox.isChecked() and 'REGTYPE' in p:
                if p['REGTYPE'][0].lower()=='bezier':
                    ui.PolynspinBox.setValue(-1)
                else:
                    ui.PolynspinBox.setValue(int(p['REGTYPE'][0]))
        
        if 'FLAGERRORS' in p:
            ui.ErrorsradioButton.setChecked(readFlag(p['FLAGERRORS'][0]))

    else:
        ui.setTaskFlag(1,Qt.CheckState.Unchecked)
        
        if 'GALMOD' in p and readFlag(p['GALMOD'][0]):
            ui.setTaskFlag(2,Qt.CheckState.Checked)
            ui.setGUIModelParam(p)
    
    
    # Setting GUI elements for SEARCH task
    if 'SEARCH' in p and readFlag(p['SEARCH'][0]):
        ui.setTaskFlag(3,Qt.CheckState.Checked)
        
        if 'SEARCHTYPE' in p:
            idx = 0 if 'spec' in p['SEARCHTYPE'][0].lower() else 1
            ui.SearchtypecomboBox.setCurrentIndex(idx)
        
        if 'SNRCUT' in p:
            ui.CuttypecomboBox.setCurrentIndex(1)
            ui.primaryCutlineEdit.setText(p['SNRCUT'][0])
        elif 'THRESHOLD' in p:
            ui.CuttypecomboBox.setCurrentIndex(0)
            ui.primaryCutlineEdit.setText(p['THRESHOLD'][0])
    
        if 'FLAGGROWTH' in p:
            ui.GrowthcheckBox.setChecked(readFlag(p['FLAGGROWTH'][0]))
            if ui.GrowthcheckBox.isChecked():
                if 'GROWTHCUT' in p:
                    ui.Cuttype2comboBox.setCurrentIndex(1)
                    ui.SecondarycutlineEdit.setText(p['GROWTHCUT'][0])
                elif 'GROWTHTHRESHOLD' in p:
                    ui.Cuttype2comboBox.setCurrentIndex(0)
                    ui.SecondarycutlineEdit.setText(p['GROWTHTHRESHOLD'][0])
        
        ui.SearchAdvgroupBox.setEnabled(True)
        ui.SearchAdvgroupBox.setChecked(True)
        
        if 'FLAGADJACENT' in p:
            ui.AdjacentcheckBox.setChecked(readFlag(p['FLAGADJACENT'][0]))
            if ui.AdjacentcheckBox.isChecked():
                ui.ThreshSpatialBox.setValue(-1)
                ui.ThreshSpatialBox.setDisabled(True)
            else:
                ui.ThreshSpatialBox.setEnabled(True)
                if 'THRESHSPATIAL' in p:
                    ui.ThreshSpatialBox.setValue(int(p['THRESHSPATIAL'][0]))
        
        if 'THRESHVELOCITY' in p:
            ui.ThreshVelspinBox.setValue(int(p['THRESHVELOCITY'][0]))
        if 'MINCHAN' in p:
            ui.MinChanspinBox.setValue(int(p['MINCHANNELS'][0]))
        if 'MINPIX' in p:
            ui.MinPixspinBox.setValue(int(p['MINPIX'][0]))
        if 'MINVOX' in p:
            ui.MinVoxspinBox.setValue(int(p['MINVOXELS'][0]))
        if 'MAXCHANNELS' in p:
            ui.MaxChanspinBox.setValue(int(p['MAXCHANNELS'][0]))
        if 'MAXANGSIZE' in p:
            ui.MaxAngSizeSpinBox.setValue(float(p['MAXANGSIZE'][0]))
        if 'REJECTBEFOREMERGE' in p:
            ui.RejectcheckBox.setChecked(readFlag(p['REJECTBEFOREMERGE'][0]))
        if 'TWOSTAGEMERGING' in p:
            ui.TwostagemergingcheckBox.setChecked(readFlag(p['TWOSTAGEMERGING'][0]))
        
    else:
        ui.setTaskFlag(3,Qt.CheckState.Unchecked)
    

    # Setting GUI elements for SMOOTH task
    if 'SMOOTH' in p and readFlag(p['SMOOTH'][0]):
        ui.setTaskFlag(4,Qt.CheckState.Checked)
        ui.SpatialSmoothgroupBox.setChecked(True)
        if 'FACTOR' in p:
            ui.FactorcheckBox.setChecked(True)
            ui.FactordoubleSpinBox.setValue(float(p['FACTOR'][0]))
        else:
            beamflag = any(x in p for x in ['BMAJ','BMIN','BPA'])
            ui.BeamgroupBox.setChecked(beamflag)
            if 'BMAJ' in p: ui.NbmajSpinBox.setValue(float(p['BMAJ'][0]))
            if 'BMIN' in p: ui.NbminSpinBox.setValue(float(p['BMIN'][0]))
            if 'BPA' in p: ui.NbpaSpinBox.setValue(float(p['BPA'][0]))
            if 'OBMAJ' in p: ui.ObmajSpinBox.setValue(float(p['OBMAJ'][0]))
            if 'OBMIN' in p: ui.ObminSpinBox.setValue(float(p['OBMIN'][0]))
            if 'OBPA' in p: ui.ObpaSpinBox.setValue(readFlag(p['OBPA'][0]))

        if 'FFT' in p: ui.FFTcheckBox.setChecked(readFlag(p['FFT'][0]))
        if 'REDUCE' in p: ui.ReducecheckBox.setChecked(readFlag(p['REDUCE'][0]))
    
    else:
        ui.setTaskFlag(4,Qt.CheckState.Unchecked)
        
    # Setting GUI elements for SMOOTHSPEC task
    if 'SMOOTHSPEC' in p and readFlag(p['SMOOTHSPEC'][0]):
        ui.setTaskFlag(4,Qt.CheckState.Checked)
        if 'WINDOW_TYPE' in p and 'HANN' in p['WINDOW_TYPE'][0].upper():  
            ui.HanninggroupBox.setChecked(True)
        if 'WINDOW_SIZE' in p:
            hanning_window = int(p['WINDOW_SIZE'][0])
            if hanning_window%2==0: hanning_window+=1
            ui.HanningspinBox.setValue(hanning_window)
    

    # Setting GUI elements for MAPS task
    if 'GLOBALPROFILE' in p:
        ui.ProfilecheckBox.setChecked(readFlag(p['GLOBALPROFILE'][0]))
    if 'TOTALMAP' in p:
        ui.TotalmapcheckBox.setChecked(readFlag(p['TOTALMAP'][0]))
    if 'VELOCITYMAP' in p:
        ui.VfieldcheckBox.setChecked(readFlag(p['VELOCITYMAP'][0]))
    if 'DISPERSIONMAP' in p:
        ui.DispmapcheckBox.setChecked(readFlag(p['DISPERSIONMAP'][0]))
    if 'RMSMAP' in p:
        ui.rmsmapcheckBox.setChecked(readFlag(p['RMSMAP'][0]))
    if 'MAPTYPE' in p:
        if 'GAUS' in p['MAPTYPE'][0]: ui.MaptypecomboBox.setCurrentIndex(1)
        else: ui.MaptypecomboBox.setCurrentIndex(0)
    
    mapflag = any([ui.ProfilecheckBox.isChecked(),ui.TotalmapcheckBox.isChecked(),\
                   ui.VfieldcheckBox.isChecked(),ui.DispmapcheckBox.isChecked(),\
                   ui.rmsmapcheckBox.isChecked()])
    ui.MomentmapsgroupBox.setChecked(mapflag)
    
    if 'PVSLICE' in p and readFlag(p['PVSLICE'][0]):
        ui.PVgroupBox.setChecked(True)
        if 'XPOS_PV' in p:
            posval = float(p['XPOS_PV'][0]) if isNumber(p['XPOS_PV'][0]) else 0
            ui.PVXposSpinBox.setValue(posval)
        if 'YPOS_PV' in p:
            posval = float(p['YPOS_PV'][0]) if isNumber(p['YPOS_PV'][0]) else 0
            ui.PVYposSpinBox.setValue(posval)
        if 'PA_PV' in p:
            posval = float(p['PA_PV'][0]) if isNumber(p['PA_PV'][0]) else 0
            ui.PVPaSpinBox.setValue(posval)
    else:
        ui.PVgroupBox.setChecked(False)
        
    
    if 'MASK' in p and 'NONE' not in p['MASK'][0].upper():
        ui.MaskgroupBox.setChecked(True)
        if'BLANKCUT' in p:
            ui.BlankcutSpinBox.setValue(float(p['BLANKCUT'][0]))
        if 'FACTOR' in p:
            ui.BlankfactordoubleSpinBox.setValue(float(p['FACTOR'][0]))
    else:
        ui.MaskgroupBox.setChecked(False)

    flag = any([ui.MomentmapsgroupBox.isChecked(),ui.PVgroupBox.isChecked(),ui.MaskgroupBox])
    if flag:
        ui.setTaskFlag(5,Qt.CheckState.Checked)
        
