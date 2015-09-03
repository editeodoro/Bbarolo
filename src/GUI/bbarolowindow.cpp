/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 Bbarp;p is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@unibo.it
-----------------------------------------------------------------------*/

#include "bbarolowindow.h"
#include "ui_bbarolowindow.h"
#include <QtGui>
#include <QListWidgetItem>
// For redirecting stdout
//#include "consolestream.h"
//#include "q_streamdebug.h"
//ConsoleStream cout(&std::cout,ui->LogtextEdit); <--- I prefer this one
//new Q_DebugStream(std::cout, ui->LogtextEdit);

BBaroloWindow::BBaroloWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::BBaroloWindow)
{
    ui->setupUi(this);
    ui->xposlineEdit->setValidator(new QDoubleValidator(0,10000,4,this));
    ui->yposlineEdit->setValidator(new QDoubleValidator(0,10000,4,this));
    ui->primaryCutlineEdit->setValidator(new QDoubleValidator(0,10000,10,this));
    ui->SecondarycutlineEdit->setValidator(new QDoubleValidator(0,10000,10,this));

    obj = "";

    // Work around of current path for MAC OSX>=10.9 bug
    QString execPath = QApplication::applicationDirPath();
    QDir currentPath = QDir::currentPath();
    if (currentPath.absolutePath()=="/" && execPath.indexOf("Contents/MacOS")!=-1) {
        currentPath = QApplication::applicationDirPath();
        currentPath.cdUp();
        currentPath.cdUp();
        currentPath.cdUp();
        QDir::setCurrent(currentPath.absolutePath());
    }

    out_path = QDir::currentPath()+"/output/";
    ui->OutfolderlineEdit->setText(out_path);

    ui->plot1->setInteraction(QCP::iRangeDrag, true);
    ui->plot2->setInteraction(QCP::iRangeDrag, true);
    ui->plot3->setInteraction(QCP::iRangeDrag, true);
    ui->plot4->setInteraction(QCP::iRangeDrag, true);
    ui->plot1->setInteraction(QCP::iRangeZoom, true);
    ui->plot2->setInteraction(QCP::iRangeZoom, true);
    ui->plot3->setInteraction(QCP::iRangeZoom, true);
    ui->plot4->setInteraction(QCP::iRangeZoom, true);

    Hide_All_3DFit_file(true);
    ui->Logolabel->setPixmap(QPixmap(currentPath.path()+"/BBaroloGUI.app/Contents/Resources/Bbarolo.tiff"));
    ui->AutocheckBox->setHidden(true);
    ui->MaskThreshSpinBox->setHidden(true);

    QFont log_font;
#ifdef MACOSX
    log_font.setFamily(QString::fromUtf8("Monaco"));
    log_font.setPointSize(12);
#else
    log_font.setFamily(QString::fromUtf8("Monospace"));
    log_font.setPointSize(10);
#endif
    ui->LogtextEdit->setFont(log_font);

}

BBaroloWindow::~BBaroloWindow()
{
    delete ui;
}

void BBaroloWindow::resetGUI() {
    enable_All();
    ui->ParamlineEdit->clear();
    ui->FitslineEdit->clear();
    ui->OutfolderlineEdit->clear();

    ui->BoxcheckBox->setChecked(false);
    ui->XminspinBox->setValue(0);
    ui->XmaxspinBox->setValue(0);
    ui->YminspinBox->setValue(0);
    ui->YmaxspinBox->setValue(0);
    ui->ZminspinBox->setValue(0);
    ui->ZmaxspinBox->setValue(0);

    set3DFitFlag(Qt::Unchecked);
    ui->NringscheckBox->setChecked(false);
    ui->NringsspinBox->setValue(0);
    ui->RadsepcheckBox->setChecked(false);
    ui->RadsepSpinBox->setValue(0.00);
    ui->XposcheckBox->setChecked(false);
    ui->xposlineEdit->clear();
    ui->YposcheckBox->setChecked(false);
    ui->yposlineEdit->clear();
    ui->wcscomboBox->setCurrentIndex(0);
    on_wcscomboBox_activated(0);
    ui->VsyscheckBox->setChecked(false);
    ui->VsysSpinBox->setValue(0.00);
    ui->VrotcheckBox->setChecked(false);
    ui->VrotSpinBox->setValue(0.00);
    ui->VdispcheckBox->setChecked(false);
    ui->VdispSpinBox->setValue(0.00);
    ui->InccheckBox->setChecked(false);
    ui->IncSpinBox->setValue(0.00);
    ui->IncDSpinBox->setValue(10.00);
    ui->PacheckBox->setChecked(false);
    ui->PaSpinBox->setValue(0.00);
    ui->PaDSpinBox->setValue(15.00);
    ui->Z0checkBox->setChecked(false);
    ui->Z0SpinBox->setValue(0.00);
    Hide_All_3DFit_file(true);

    ui->FreeParametersframe->setDisabled(true);
    ui->vrotcheckBox->setChecked(true);
    ui->vdispcheckBox->setChecked(true);
    ui->inccheckBox->setChecked(true);
    ui->pacheckBox->setChecked(true);
    ui->xposcheckBox->setChecked(false);
    ui->yposcheckBox->setChecked(false);
    ui->vsyscheckBox->setChecked(false);
    ui->z0checkBox->setChecked(false);

    ui->AdvancedgroupBox->setDisabled(true);
    ui->LtypecomboBox->setCurrentIndex(0);
    ui->CdensSpinBox->setValue(10);
    ui->NVspinBox->setValue(-1);
    ui->FtypecomboBox->setCurrentIndex(1);
    ui->WfunccomboBox->setCurrentIndex(1);
    ui->SecondstagecheckBox->setChecked(true);
    ui->PolynspinBox->setValue(-1);
    ui->NormcomboBox->setCurrentIndex(0);
    ui->TollineEdit->setText("1E-03");
    ui->ErrorsradioButton->setChecked(false);

    setSearchFlag(Qt::Unchecked);
    ui->SearchgroupBox->setChecked(false);
    ui->SearchtypecomboBox->setCurrentIndex(1);
    ui->CuttypecomboBox->setCurrentIndex(1);
    ui->primaryCutlineEdit->setText("5.00");
    ui->Cuttype2comboBox->setCurrentIndex(1);
    ui->SecondarycutlineEdit->setText("2.50");
    ui->SearchAdvgroupBox->setDisabled(true);
    ui->SearchAdvgroupBox->setChecked(false);
    ui->AdjacentcheckBox->setChecked(true);
    ui->MinPixspinBox->setValue(-1);
    ui->MinChanspinBox->setValue(-1);
    ui->MinVoxspinBox->setValue(-1);
    ui->MaxChanspinBox->setValue(-1);
    ui->MaxAngSizeSpinBox->setValue(-1);
    ui->RejectcheckBox->setChecked(true);
    ui->TwostagemergingcheckBox->setChecked(true);

    setSmoothFlag(Qt::Unchecked);
    ui->SmoothOutlineEdit->clear();
    ui->FactorcheckBox->setChecked(true);
    ui->FactordoubleSpinBox->setValue(2);
    ui->BeamgroupBox->setChecked(false);
    ui->ObmajSpinBox->setValue(-1);
    ui->ObminSpinBox->setValue(-1);
    ui->ObpaSpinBox->setValue(0);
    ui->NbmajSpinBox->setValue(-1);
    ui->NbminSpinBox->setValue(-1);
    ui->NbpaSpinBox->setValue(0);
    ui->FFTcheckBox->setChecked(true);
    ui->ReducecheckBox->setChecked(false);

    ui->MomentmapsgroupBox->setChecked(false);
    ui->ProfilecheckBox->setChecked(false);
    ui->TotalmapcheckBox->setChecked(false);
    ui->VfieldcheckBox->setChecked(false);
    ui->DispmapcheckBox->setChecked(false);

    ui->MaskgroupBox->setChecked(false);
    ui->BlankfactordoubleSpinBox->setValue(2.00);
    ui->BlankcutSpinBox->setValue(3.00);

    out_path = QDir::currentPath()+"/output/";
    ui->stackedWidget->setCurrentIndex(0);
    ui->listWidget->setCurrentRow(0);
}

void BBaroloWindow::enable_All() {
    ui->ParamlineEdit->setEnabled(true);
    ui->FitslineEdit->setEnabled(true);
    ui->ParampushButton->setEnabled(true);
    ui->FitspushButton->setEnabled(true);
    ui->ResetpushButton->setEnabled(true);
    ui->BoxcheckBox->setEnabled(true);
    ui->OutfolderlineEdit->setEnabled(true);
    ui->OutfolderpushButton->setEnabled(true);

    if (get3DFitFlag()) {
        ui->GalfitgroupBox->setEnabled(true);
        ui->FreeParametersframe->setEnabled(true);
        ui->AdvancedgroupBox->setEnabled(true);
    }
    if (getSmoothFlag()) ui->SmoothgroupBox->setEnabled(true);
    if (getSearchFlag()) {
        ui->SearchgroupBox->setEnabled(true);
        ui->SearchAdvgroupBox->setEnabled(true);
    }
    ui->MomentmapsgroupBox->setEnabled(true);
    ui->MaskgroupBox->setEnabled(true);
    ui->AutocheckBox->setEnabled(true);
}

void BBaroloWindow::disable_All() {
    ui->FitslineEdit->setDisabled(true);
    ui->ParampushButton->setDisabled(true);
    ui->ResetpushButton->setDisabled(true);
    ui->OutfolderpushButton->setDisabled(true);
    ui->OutfolderlineEdit->setDisabled(true);
    ui->ParamlineEdit->setDisabled(true);
    ui->FitspushButton->setDisabled(true);
    ui->BoxcheckBox->setDisabled(true);
    ui->GalfitgroupBox->setDisabled(true);
    ui->AdvancedgroupBox->setDisabled(true);
    ui->FreeParametersframe->setDisabled(true);
    ui->SmoothgroupBox->setDisabled(true);
    ui->SearchgroupBox->setDisabled(true);
    ui->SearchAdvgroupBox->setDisabled(true);
    ui->MomentmapsgroupBox->setDisabled(true);
    ui->MaskgroupBox->setDisabled(true);
    ui->AutocheckBox->setDisabled(true);
}

void BBaroloWindow::on_AutocheckBox_stateChanged()
{
    if(ui->AutocheckBox->isChecked()) disable_All();
    else enable_All();
    ui->FitslineEdit->setEnabled(true);
    ui->FitspushButton->setEnabled(true);
    ui->AutocheckBox->setEnabled(true);
}

void BBaroloWindow::on_BoxcheckBox_stateChanged()
{
    if (ui->BoxcheckBox->isChecked()) {
        ui->XmaxspinBox->setEnabled(true);
        ui->XminspinBox->setEnabled(true);
        ui->YmaxspinBox->setEnabled(true);
        ui->YminspinBox->setEnabled(true);
        ui->ZmaxspinBox->setEnabled(true);
        ui->ZminspinBox->setEnabled(true);
        ui->Xboxlabel->setEnabled(true);
        ui->Yboxlabel->setEnabled(true);
        ui->Zboxlabel->setEnabled(true);
    }
    else {
        ui->XmaxspinBox->setDisabled(true);
        ui->XminspinBox->setDisabled(true);
        ui->YmaxspinBox->setDisabled(true);
        ui->YminspinBox->setDisabled(true);
        ui->ZmaxspinBox->setDisabled(true);
        ui->ZminspinBox->setDisabled(true);
        ui->Xboxlabel->setDisabled(true);
        ui->Yboxlabel->setDisabled(true);
        ui->Zboxlabel->setDisabled(true);
    }
}


//////////////////////
/// 3Dfit panel slots
//////////////////////

void BBaroloWindow::on_NringscheckBox_stateChanged()
{
    if (ui->NringscheckBox->isChecked()) ui->NringsspinBox->setDisabled(true);
    else ui->NringsspinBox->setEnabled(true);
}

void BBaroloWindow::on_RadsepcheckBox_stateChanged()
{
    if (ui->RadsepcheckBox->isChecked()) ui->RadsepSpinBox->setDisabled(true);
    else ui->RadsepSpinBox->setEnabled(true);
}

void BBaroloWindow::on_XposcheckBox_stateChanged()
{
    if (ui->XposcheckBox->isChecked()) ui->xposlineEdit->setDisabled(true);
    else ui->xposlineEdit->setEnabled(true);
}

void BBaroloWindow::on_YposcheckBox_stateChanged()
{
    if (ui->YposcheckBox->isChecked()) ui->yposlineEdit->setDisabled(true);
    else ui->yposlineEdit->setEnabled(true);
}

void BBaroloWindow::on_wcscomboBox_activated(int index)
{
    switch (index) {
    case 0:
        ui->xposlineEdit->setInputMask("");
        ui->yposlineEdit->setInputMask("");
        ui->xposlineEdit->setText(tr("0.00"));
        ui->yposlineEdit->setText(tr("0.00"));
        break;
    case 1:
        ui->xposlineEdit->setInputMask("000.0000;_");
        ui->yposlineEdit->setInputMask("#00.0000;_");
        ui->xposlineEdit->setText(tr("000.0000"));
        ui->yposlineEdit->setText(tr("+00.0000"));
        break;
    case 2:
        ui->xposlineEdit->setInputMask("00:00:00.00;_");
        ui->yposlineEdit->setInputMask("#00:00:00.00;_");
        ui->xposlineEdit->setText(tr("000000.00"));
        ui->yposlineEdit->setText(tr("+000000.00"));
        break;
    }

    ui->xposlineEdit->setCursorPosition(1);
    ui->yposlineEdit->setCursorPosition(1);
}

void BBaroloWindow::on_VsyscheckBox_stateChanged()
{
    if (ui->VsyscheckBox->isChecked()) ui->VsysSpinBox->setDisabled(true);
    else ui->VsysSpinBox->setEnabled(true);
}

void BBaroloWindow::on_VrotcheckBox_stateChanged()
{
    if (ui->VrotcheckBox->isChecked()) ui->VrotSpinBox->setDisabled(true);
    else ui->VrotSpinBox->setEnabled(true);
}

void BBaroloWindow::on_VdispcheckBox_stateChanged()
{
    if (ui->VdispcheckBox->isChecked()) ui->VdispSpinBox->setDisabled(true);
    else ui->VdispSpinBox->setEnabled(true);
}

void BBaroloWindow::on_InccheckBox_stateChanged()
{
    if (ui->InccheckBox->isChecked()) {
        ui->IncSpinBox->setDisabled(true);
        ui->IncDSpinBox->setDisabled(true);
    }
    else {
        ui->IncSpinBox->setEnabled(true);
        ui->IncDSpinBox->setEnabled(true);
    }
}

void BBaroloWindow::on_PacheckBox_stateChanged()
{
    if (ui->PacheckBox->isChecked()) {
        ui->PaSpinBox->setDisabled(true);
        ui->PaDSpinBox->setDisabled(true);
    }
    else {
        ui->PaSpinBox->setEnabled(true);
        ui->PaDSpinBox->setEnabled(true);
    }
}

void BBaroloWindow::on_Z0checkBox_stateChanged()
{
    if (ui->Z0checkBox->isChecked()) ui->Z0SpinBox->setDisabled(true);
    else ui->Z0SpinBox->setEnabled(true);
}

void BBaroloWindow::on_SecondstagecheckBox_stateChanged()
{
    if(ui->SecondstagecheckBox->isChecked()) ui->PolynspinBox->setEnabled(true);
    else ui->PolynspinBox->setDisabled(true);

}

void BBaroloWindow::on_MasktoolButton_clicked()
{

    if (ui->MaskcomboBox->currentIndex()==0) {
        ui->listWidget->setCurrentRow(5);
        ui->MaskgroupBox->setChecked(true);

    }
    else if (ui->MaskcomboBox->currentIndex()==1) {
        ui->listWidget->setCurrentRow(3);
        ui->SearchgroupBox->setEnabled(true);
    }

}

void BBaroloWindow::on_MaskcomboBox_currentIndexChanged(int index)
{
    if (index==0) {
        ui->MasktoolButton->setHidden(false);
        ui->MaskThreshSpinBox->setHidden(true);
    }
    else if (index==1) {
        ui->MasktoolButton->setHidden(false);
        ui->MaskThreshSpinBox->setHidden(true);
    }
    else if (index==2) {
        ui->MaskThreshSpinBox->setHidden(false);
        ui->MasktoolButton->setHidden(true);
    }
    else {
        ui->MaskThreshSpinBox->setHidden(true);
        ui->MasktoolButton->setHidden(true);
    }

}

void BBaroloWindow::on_RadseppushButton_clicked(){Selected_3DFit_file(ui->RadsepFilelineEdit,ui->RadsepFilespinBox,ui->RadsepSpinBox); Hide_3DFit_file(ui->RadsepFilelineEdit,ui->RadsepFilespinBox,ui->NringsspinBox, false);}
void BBaroloWindow::on_XpospushButton_clicked(){Selected_3DFit_file(ui->XposFilelineEdit,ui->XposFilespinBox,ui->xposlineEdit);}
void BBaroloWindow::on_YpospushButton_clicked(){Selected_3DFit_file(ui->YposFilelineEdit,ui->YposFilespinBox,ui->yposlineEdit);}
void BBaroloWindow::on_VsyspushButton_clicked(){Selected_3DFit_file(ui->VsysFilelineEdit,ui->VsysFilespinBox,ui->VsysSpinBox);}
void BBaroloWindow::on_VrotpushButton_clicked(){Selected_3DFit_file(ui->VrotFilelineEdit,ui->VrotFilespinBox,ui->VrotSpinBox);}
void BBaroloWindow::on_VdisppushButton_clicked(){Selected_3DFit_file(ui->VdispFilelineEdit,ui->VdispFilespinBox,ui->VdispSpinBox);}
void BBaroloWindow::on_IncpushButton_clicked(){Selected_3DFit_file(ui->IncFilelineEdit,ui->IncFilespinBox,ui->IncSpinBox);}
void BBaroloWindow::on_PapushButton_clicked(){Selected_3DFit_file(ui->PaFilelineEdit,ui->PaFilespinBox,ui->PaSpinBox);}
void BBaroloWindow::on_Z0pushButton_clicked(){Selected_3DFit_file(ui->Z0FilelineEdit,ui->Z0FilespinBox,ui->Z0SpinBox);}
void BBaroloWindow::on_DenspushButton_clicked(){Selected_3DFit_file(ui->DensFilelineEdit,ui->DensFilespinBox,ui->DensSpinBox);}

void BBaroloWindow::on_RadsepFilelineEdit_editingFinished(){Hide_3DFit_file(ui->RadsepFilelineEdit,ui->RadsepFilespinBox,ui->RadsepSpinBox,ui->RadsepFilelineEdit->text().isEmpty());Hide_3DFit_file(ui->RadsepFilelineEdit,ui->RadsepFilespinBox,ui->NringsspinBox,ui->RadsepFilelineEdit->text().isEmpty());}
void BBaroloWindow::on_XposFilespinBox_editingFinished(){Hide_3DFit_file(ui->XposFilelineEdit,ui->XposFilespinBox,ui->xposlineEdit,ui->XposFilelineEdit->text().isEmpty());}
void BBaroloWindow::on_YposFilelineEdit_editingFinished(){Hide_3DFit_file(ui->YposFilelineEdit,ui->YposFilespinBox,ui->yposlineEdit,ui->YposFilelineEdit->text().isEmpty());}
void BBaroloWindow::on_VsysFilelineEdit_editingFinished(){Hide_3DFit_file(ui->VsysFilelineEdit,ui->VsysFilespinBox,ui->VsysSpinBox,ui->VsysFilelineEdit->text().isEmpty());}
void BBaroloWindow::on_VrotFilelineEdit_editingFinished() {Hide_3DFit_file(ui->VrotFilelineEdit,ui->VrotFilespinBox,ui->VrotSpinBox,ui->VrotFilelineEdit->text().isEmpty());}
void BBaroloWindow::on_VdispFilelineEdit_editingFinished(){Hide_3DFit_file(ui->VdispFilelineEdit,ui->VdispFilespinBox,ui->VdispSpinBox,ui->VdispFilelineEdit->text().isEmpty());}
void BBaroloWindow::on_IncFilelineEdit_editingFinished(){Hide_3DFit_file(ui->IncFilelineEdit,ui->IncFilespinBox,ui->IncSpinBox,ui->IncFilelineEdit->text().isEmpty());}
void BBaroloWindow::on_PaFilelineEdit_editingFinished(){Hide_3DFit_file(ui->PaFilelineEdit,ui->PaFilespinBox,ui->PaSpinBox,ui->PaFilelineEdit->text().isEmpty());}
void BBaroloWindow::on_Z0FilelineEdit_editingFinished(){Hide_3DFit_file(ui->Z0FilelineEdit,ui->Z0FilespinBox,ui->Z0SpinBox,ui->Z0FilelineEdit->text().isEmpty());}
void BBaroloWindow::on_DensFilelineEdit_editingFinished(){Hide_3DFit_file(ui->DensFilelineEdit,ui->DensFilespinBox,ui->DensSpinBox,ui->DensFilelineEdit->text().isEmpty());}

void BBaroloWindow::on_GalfittoolButton_clicked(){ui->stackedWidget->setCurrentIndex(ui->stackedWidget->count()-2);}
void BBaroloWindow::on_GalfitAdvtoolButton_clicked(){ui->stackedWidget->setCurrentIndex(1);}

void BBaroloWindow::Hide_3DFit_file(QLineEdit *le, QSpinBox *sb, QWidget *qw, bool Hide) {
        le->setHidden(Hide);
        sb->setHidden(Hide);
        qw->setEnabled(Hide);
}

void BBaroloWindow::Hide_All_3DFit_file(bool Hide) {
    Hide_3DFit_file(ui->RadsepFilelineEdit,ui->RadsepFilespinBox,ui->RadsepSpinBox,Hide);
    Hide_3DFit_file(ui->XposFilelineEdit,ui->XposFilespinBox,ui->xposlineEdit,Hide);
    Hide_3DFit_file(ui->YposFilelineEdit,ui->YposFilespinBox,ui->yposlineEdit,Hide);
    Hide_3DFit_file(ui->VsysFilelineEdit,ui->VsysFilespinBox,ui->VsysSpinBox,Hide);
    Hide_3DFit_file(ui->VrotFilelineEdit,ui->VrotFilespinBox,ui->VrotSpinBox,Hide);
    Hide_3DFit_file(ui->VdispFilelineEdit,ui->VdispFilespinBox,ui->VdispSpinBox,Hide);
    Hide_3DFit_file(ui->IncFilelineEdit,ui->IncFilespinBox,ui->IncSpinBox,Hide);
    Hide_3DFit_file(ui->PaFilelineEdit,ui->PaFilespinBox,ui->PaSpinBox,Hide);
    Hide_3DFit_file(ui->Z0FilelineEdit,ui->Z0FilespinBox,ui->Z0SpinBox,Hide);
    Hide_3DFit_file(ui->DensFilelineEdit,ui->DensFilespinBox,ui->DensSpinBox,Hide);
}

void BBaroloWindow::Selected_3DFit_file(QLineEdit *le, QSpinBox *sb, QWidget *qw) {

    QString filename = QFileDialog::getOpenFileName(this, tr("Open file"), QDir::currentPath(), tr("All files (*.*)"), 0, QFileDialog::DontUseNativeDialog);
    if (!filename.isEmpty()) {
        le->setText(filename);
        Hide_3DFit_file(le,sb,qw, false);
    }
}


//////////////////////
/// Search panel slots
//////////////////////
void BBaroloWindow::on_SearchgroupBox_clicked()
{
    if (ui->SearchgroupBox->isChecked()) {
        ui->SearchAdvgroupBox->setEnabled(true);
    }
    else {
        ui->SearchAdvgroupBox->setDisabled(true);
    }
}

void BBaroloWindow::on_GrowthcheckBox_stateChanged()
{
    if (ui->GrowthcheckBox->isChecked()) {
        ui->SecondarycutlineEdit->setEnabled(true);
        ui->Cuttype2comboBox->setEnabled(true);
    }
    else {
        ui->SecondarycutlineEdit->setDisabled(true);
        ui->Cuttype2comboBox->setDisabled(true);
    }
}

void BBaroloWindow::on_SearchAdvgroupBox_clicked()
{
    if (ui->SearchAdvgroupBox->isChecked()){
        ui->ThreshSpatialBox->setDisabled(true);
        ui->RejectcheckBox->setChecked(true);
    }
    else {
        ui->AdjacentcheckBox->setChecked(true);
        ui->RejectcheckBox->setChecked(true);
        ui->TwostagemergingcheckBox->setChecked(true);
    }
}


void BBaroloWindow::on_AdjacentcheckBox_stateChanged(int arg1)
{
    if (arg1==Qt::Checked) {
        ui->ThreshSpatialBox->setValue(-1);
        ui->ThreshSpatialBox->setDisabled(true);
    }
    else {
        ui->ThreshSpatialBox->setEnabled(true);
        ui->ThreshVelspinBox->setEnabled(true);
    }
}


void BBaroloWindow::on_SearchtoolButton_clicked(){ui->stackedWidget->setCurrentIndex(ui->stackedWidget->count()-1);}
void BBaroloWindow::on_SearchAdvtoolButton_clicked(){ui->stackedWidget->setCurrentIndex(2);}

//////////////////////
/// Smooth panel slots
//////////////////////
void BBaroloWindow::on_BeamgroupBox_clicked()
{
    if (ui->BeamgroupBox->isChecked()) {
        ui->FactorcheckBox->setChecked(false);
    }
    else ui->FactorcheckBox->setChecked(true);
}

void BBaroloWindow::on_FactorcheckBox_stateChanged()
{
    if(ui->FactorcheckBox->isChecked()) {
        ui->FactordoubleSpinBox->setEnabled(true);
        ui->BeamgroupBox->setChecked(false);
        ui->ReducecheckBox->setEnabled(true);
        on_FactordoubleSpinBox_valueChanged(ui->FactordoubleSpinBox->value());
    }
    else {
        ui->FactordoubleSpinBox->setDisabled(true);
        ui->BeamgroupBox->setChecked(true);
        ui->ReducecheckBox->setChecked(false);
        ui->ReducecheckBox->setDisabled(true);
    }
}

void BBaroloWindow::on_NbmajSpinBox_valueChanged(double arg1)
{
    if (obj!=""){
        QString smooth_name = "./output/"+obj+"/"+obj+"_s"+QString::number(lround(arg1));
        if (ui->ReducecheckBox->isChecked()) smooth_name+="red";
        smooth_name+=".fits";
        ui->SmoothOutlineEdit->setText(smooth_name);
    }
}

void BBaroloWindow::on_SmoothOutpushButton_clicked()
{
    QString filename = QFileDialog::getSaveFileName(this,tr("Save FITS file"),QDir::currentPath(),tr("FITS files (*.fits *.FITS *.fit *.FIT)"));
    if (!filename.isNull()) ui->SmoothOutlineEdit->setText(filename);
}

void BBaroloWindow::on_FactordoubleSpinBox_valueChanged(double arg1)
{
    ui->NbmajSpinBox->setValue(arg1*ui->ObmajSpinBox->value());
    ui->NbminSpinBox->setValue(arg1*ui->ObminSpinBox->value());
    ui->NbpaSpinBox->setValue(ui->ObpaSpinBox->value());
}

void BBaroloWindow::on_ReducecheckBox_clicked(bool checked)
{
    if (obj!=""){
        QString smooth_name = "./output/"+obj+"/"+obj+"_s"+QString::number(lround(ui->NbmajSpinBox->value()));
        if (checked) smooth_name+="red";
        smooth_name+=".fits";
        ui->SmoothOutlineEdit->setText(smooth_name);
    }
}


void BBaroloWindow::on_actionExport_parameter_file_triggered()
{
    QString filename = QFileDialog::getSaveFileName(this,tr("Save Parameter file"),QDir::currentPath(),tr("All files (*.*)"));
    writeParamFile(filename);
}


void BBaroloWindow::switch_Galmod_3Dfit(bool toHide) {
    ui->FreeParametersframe->setHidden(toHide);
    ui->IncDSpinBox->setHidden(toHide);
    ui->PaDSpinBox->setHidden(toHide);
    ui->IncPMlabel->setHidden(toHide);
    ui->PaPMlabel->setHidden(toHide);
    ui->WfunccomboBox->setHidden(toHide);
    ui->FtypecomboBox->setHidden(toHide);
    ui->Ftypelabel->setHidden(toHide);
    ui->WfunccomboBox->setHidden(toHide);
    ui->Wfunclabel->setHidden(toHide);
    ui->SecondstagecheckBox->setHidden(toHide);
    ui->PolynspinBox->setHidden(toHide);
    ui->Polynlabel->setHidden(toHide);
    ui->TollineEdit->setHidden(toHide);
    ui->Tollabel->setHidden(toHide);
    ui->ErrorsradioButton->setHidden(toHide);
}

/////////////////////////////////
/// List widget slots & Functions
/////////////////////////////////
int BBaroloWindow::getCurrentRow(QListWidgetItem *item) {

    int currentRow=0;
    QString itemLab = item->text();
    int rows = ui->listWidget->model()->rowCount();

    for (int i=0; i<rows; i++) {
        QString rowLab = ui->listWidget->item(i)->text();
        if (itemLab==rowLab) currentRow=i;
    }
    return currentRow;
}

void BBaroloWindow::on_listWidget_currentRowChanged(int currentRow)
{
    if (currentRow==1) switch_Galmod_3Dfit(false);
    else {
        switch_Galmod_3Dfit(true);
        currentRow--;
    }
    ui->stackedWidget->setCurrentIndex(currentRow);
}


void BBaroloWindow::on_listWidget_itemChanged(QListWidgetItem *item)
{
    QString itemLab = item->text();
    int currentRow=getCurrentRow(item);

    if (ui->RunpushButton->isEnabled()) {
        switch (currentRow) {
        case 1:
            set3DFitFlag(item->checkState());
            ui->GalfitgroupBox->setEnabled(item->checkState());
            ui->FreeParametersframe->setEnabled(item->checkState());
            ui->AdvancedgroupBox->setEnabled(item->checkState());
            break;
        case 2:
            set3DModelFlag(item->checkState());
            ui->GalfitgroupBox->setEnabled(item->checkState());
            ui->AdvancedgroupBox->setEnabled(item->checkState());
        case 3:
            ui->SearchgroupBox->setEnabled(item->checkState());
            ui->SearchAdvgroupBox->setEnabled(item->checkState());
            break;
        case 4:
            ui->SmoothgroupBox->setEnabled(item->checkState());
        default:
            break;
        }
    }

    ui->listWidget->setCurrentRow(currentRow);

}

void BBaroloWindow::on_listWidget_itemDoubleClicked(QListWidgetItem *item)
{
    int currentRow=getCurrentRow(item);

    if (currentRow<ui->listWidget->model()->rowCount()-2 &&
        currentRow>0 && ui->RunpushButton->isEnabled()) {
        if (item->checkState()==Qt::Checked) item->setCheckState(Qt::Unchecked);
        else item->setCheckState(Qt::Checked);
    }
}

void BBaroloWindow::on_HidetoolButton_clicked()
{
    if (ui->listWidget->isHidden()) {
        ui->listWidget->setHidden(false);
        ui->HidetoolButton->setArrowType(Qt::UpArrow);
        ui->HidetoolButton->setToolButtonStyle(Qt::ToolButtonIconOnly);
        ui->MaingridLayout->addWidget(ui->HidetoolButton, 0, 0);
    }
    else {
        ui->listWidget->setHidden(true);
        ui->HidetoolButton->setArrowType(Qt::DownArrow);
        ui->HidetoolButton->setToolButtonStyle(Qt::ToolButtonTextBesideIcon);
        if(ui->listWidget->currentRow()>=0)
            ui->HidetoolButton->setText("   "+ui->listWidget->currentItem()->text());
        ui->MaingridLayout->addWidget(ui->HidetoolButton, 0, 1);

    }
}


void BBaroloWindow::set3DFitFlag(Qt::CheckState c) {
    ui->listWidget->item(1)->setCheckState(c);
    if (c==Qt::Checked && ui->listWidget->item(2)->checkState()==Qt::Checked)
        ui->listWidget->item(2)->setCheckState(Qt::Unchecked);
}

bool BBaroloWindow::get3DFitFlag() {return ui->listWidget->item(1)->checkState();};

void BBaroloWindow::set3DModelFlag(Qt::CheckState c) {
    ui->listWidget->item(2)->setCheckState(c);
    if (c==Qt::Checked && ui->listWidget->item(1)->checkState()==Qt::Checked)
        ui->listWidget->item(1)->setCheckState(Qt::Unchecked);
}

bool BBaroloWindow::get3DModelFlag() {return ui->listWidget->item(2)->checkState();};
void BBaroloWindow::setSearchFlag(Qt::CheckState c) {ui->listWidget->item(3)->setCheckState(c);};
bool BBaroloWindow::getSearchFlag() {return ui->listWidget->item(3)->checkState();};
void BBaroloWindow::setSmoothFlag(Qt::CheckState c) {ui->listWidget->item(4)->setCheckState(c);};
bool BBaroloWindow::getSmoothFlag() {return ui->listWidget->item(4)->checkState();};
void BBaroloWindow::setMapsFlag(Qt::CheckState c) {ui->listWidget->item(5)->setCheckState(c);};
bool BBaroloWindow::getMapsFlag() {return ui->listWidget->item(5)->checkState();};

std::string BBaroloWindow::getFileString(QLineEdit *le, QSpinBox *sb) {
    QString to_return = "file("+le->text()+","+QString::number(sb->value())+")";
    return to_return.toStdString();
}

int BBaroloWindow::readFileString(QString &s) {
    if (!s.contains("ile")) return -1;
    else {
        QStringList list = s.split(",");
        list[0].remove(0,list[0].indexOf("(")+1);
        for (int i=1; i<list.size();i++) {
            list[i].remove(")");
        }
        s = list[0];
        return list[1].toInt();

    }
}

