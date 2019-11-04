/*-----------------------------------------------------------------------
 This program is free software; you can redistribute it and/or modify it
 under the terms of the GNU General Public License as published by the
 Free Software Foundation; either version 2 of the License, or (at your
 option) any later version.

 BBarolo is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 for more details.

 You should have received a copy of the GNU General Public License
 along with BBarolo; if not, write to the Free Software Foundation,
 Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA

 Correspondence concerning BBarolo may be directed to:
    Internet email: enrico.diteodoro@gmail.com
-----------------------------------------------------------------------*/

#include<iostream>
#include<fstream>
#include <iomanip>
#include <sys/stat.h>
#include <QtGui>
#include "bbarolowindow.h"
#include "ui_bbarolowindow.h"
#include "qcustomplot.h"
#include "../Arrays/param.hh"
#include "../Arrays/header.hh"
#include "../Utilities/utils.hh"


void BBaroloWindow::on_FitslineEdit_editingFinished()
{

    QString filename = ui->FitslineEdit->text().trimmed();
    static QString lastFilename = "";

    if (filename=="") {lastFilename = filename; return;}
    if(fexists(filename.toStdString())) {
        
        if (filename.toStdString().find("./")==0) {
            filename.remove(0,1);
            filename = QDir::currentPath()+filename;
        }
        
        ui->FitslineEdit->setText(filename);

        Header *h = new Header;
        h->header_read(filename.toStdString());

        if (lastFilename!=filename) {
            ui->XminspinBox->setRange(1,h->DimAx(0));
            ui->XminspinBox->setValue(1);
            ui->XmaxspinBox->setRange(1,h->DimAx(0));
            ui->XmaxspinBox->setValue(h->DimAx(0));
            ui->YminspinBox->setRange(1,h->DimAx(1));
            ui->YminspinBox->setValue(1);
            ui->YmaxspinBox->setRange(1,h->DimAx(1));
            ui->YmaxspinBox->setValue(h->DimAx(1));
            ui->ZminspinBox->setRange(1,h->DimAx(2));
            ui->ZminspinBox->setValue(1);
            ui->ZmaxspinBox->setRange(1,h->DimAx(2));
            ui->ZmaxspinBox->setValue(h->DimAx(2));
        }

        ui->ObmajSpinBox->setValue(h->Bmaj()*arcsconv(h->Cunit(0)));
        ui->ObminSpinBox->setValue(h->Bmin()*arcsconv(h->Cunit(1)));
        ui->ObpaSpinBox->setValue(h->Bpa());
        ui->NbmajSpinBox->setValue(ui->FactordoubleSpinBox->value()*ui->ObmajSpinBox->value());
        ui->NbminSpinBox->setValue(ui->FactordoubleSpinBox->value()*ui->ObminSpinBox->value());
        ui->NbpaSpinBox->setValue(ui->ObpaSpinBox->value());

        obj = QString::fromStdString(h->Name());
        out_path = QDir::currentPath()+"/output/"+obj+"/";
        ui->OutfolderlineEdit->setText(out_path);

        QString smo_out = out_path+obj+"_s";
        if (ui->FactorcheckBox->isChecked())
            smo_out += QString::number(lround(ui->FactordoubleSpinBox->value()*ui->ObmajSpinBox->value()));
        else smo_out += (QString::number(lround(ui->NbmajSpinBox->value()))+".fits");
        if (ui->ReducecheckBox->isChecked()) smo_out+="red";
        smo_out+=".fits";
        ui->SmoothOutlineEdit->setText(smo_out);

        if (lastFilename!=filename) {
            if (ui->RadsepSpinBox->value()==0) ui->RadsepSpinBox->setValue(h->Bmaj()*arcsconv(h->Cunit(0)));
            ui->xposlineEdit->setText(QString::fromStdString(to_string(h->Crpix(0)-1,2)));
            ui->yposlineEdit->setText(QString::fromStdString(to_string(h->Crpix(1)-1,2)));
        }
        delete h;
    }
    else {
        QMessageBox msgBox;
        msgBox.setText("WARNING: "+filename+" does not exist.");
        msgBox.setInformativeText("Do you want to select a valid FITS file?");
        msgBox.setStandardButtons(QMessageBox::Open | QMessageBox::Cancel);
        msgBox.setDefaultButton(QMessageBox::Open);
        int ret = msgBox.exec();

        switch (ret) {
            case QMessageBox::Open:
                on_FitspushButton_clicked();
                break;
          case QMessageBox::Cancel:
                ui->FitslineEdit->clear();
                break;
        }
    }

    lastFilename = filename;
}

void BBaroloWindow::on_FitspushButton_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this, tr("Open FITSFILE"), QDir::currentPath(), tr("FITS files (*.fit *.fits *.FITS *.FIT);; All files (*.*)"), 0, QFileDialog::DontUseNativeDialog);
    if (!filename.isNull()) {
        ui->FitslineEdit->setText(filename);
        on_FitslineEdit_editingFinished();
    }
}

void BBaroloWindow::on_ParampushButton_clicked()
{
    QString filename = QFileDialog::getOpenFileName(this, tr("Open parameter file"), QDir::currentPath(), tr("Parameter files (*.par *.txt *.dat)"), 0, QFileDialog::DontUseNativeDialog);
    if( !filename.isNull()) {
        resetGUI();
        ui->ParamlineEdit->setText(filename);
        readParamFromFile(filename.toStdString());
    }
}

void BBaroloWindow::on_OutfolderpushButton_clicked()
{
    QString outdir = QFileDialog::getExistingDirectory(this, tr("Open Directory"),QDir::currentPath(),QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);
    if(!outdir.isNull()) {
        outdir.append("/");
        ui->OutfolderlineEdit->setText(outdir);
        out_path = outdir;
    }
}

void BBaroloWindow::updateText()
{
    QString appendText(proc->readAll());
    ui->LogtextEdit->insertPlainText (appendText);
    ui->LogtextEdit->moveCursor (QTextCursor::End);

    if (get3DFitFlag()) plotParameters();
}

void BBaroloWindow::updateExit()
{
    ui->RunpushButton->setEnabled(true);
    enable_All();
    
    if (proc->exitCode()==0 && proc->exitStatus()==0) ui->statusBar->showMessage("  BBarolo successfully terminated");
    else {
        if (proc->exitStatus()==0)  ui->statusBar->showMessage("  BBarolo exited with errors");
        else ui->statusBar->showMessage("  BBarolo has been killed");
    }

    QFile outfile(out_path+"output.log");
    outfile.open(QIODevice::ReadWrite |QIODevice::Truncate | QIODevice::Text);
    QTextStream out(&outfile);
    out << ui->LogtextEdit->toPlainText() << endl;
    proc->deleteLater();

}

void BBaroloWindow::on_KillpushButton_clicked()
{
    if (proc->state()==2) {
        proc->terminate();
        QTimer::singleShot(2000, proc, SLOT(kill()));
        ui->LogtextEdit->append("\n\n ----------> Killed by the user. <----------\n");
        ui->RunpushButton->setEnabled(true);
    }
}

void BBaroloWindow::on_RunpushButton_clicked() {

    ui->listWidget->setCurrentRow(ui->listWidget->model()->rowCount()-2);
    if (ui->FitslineEdit->text().isEmpty()) {
        ui->LogtextEdit->append("ERROR: You must select at least a FITS file!");
        return;
    }
    ui->RunpushButton->setDisabled(true);
    disable_All();

    ui->plot1->replot();
    std::string outfolder = ui->OutfolderlineEdit->text().trimmed().toStdString();
    checkHome(outfolder);
    if (outfolder[outfolder.size()-1]!='/') outfolder.append("/");
    out_path = QString::fromStdString(outfolder);
    mkdirp (out_path.toStdString().c_str());

    QFile::remove(out_path+"ringlog1.txt");
    QFile::remove(out_path+"ringlog2.txt");

    QString paramfile = out_path+obj+".par";
    writeParamFile(paramfile);


    QString argv1="-p", argv2=paramfile;
    if (ui->AutocheckBox->isChecked()) {
        argv1 = "-f";
        argv2 = ui->FitslineEdit->text();
    }

    ui->LogtextEdit->clear();
    ui->statusBar->showMessage("  BBarolo is running ...");

    proc = new QProcess(this);
    QString cmd = QApplication::applicationDirPath()+"/BBarolo "+argv1+" "+argv2;
    connect(proc, SIGNAL(readyRead()), this, SLOT(updateText()));
    connect(proc, SIGNAL(readyReadStandardError()), this, SLOT(updateText()));
    connect(proc, SIGNAL(finished(int)), this, SLOT(updateExit()));
    proc->start(cmd);
    
}

void BBaroloWindow::readParamFromFile(std::string filein) {

    Param *par = new Param;
    par->readParams(filein);

    QString filename = QString::fromStdString(par->getImageFile());
    bool isBox = filename.contains("[") && filename.contains("]");
    ui->FitslineEdit->setText("");
    on_FitslineEdit_editingFinished();
    ui->ThreadspinBox->setValue(par->getThreads());
    if (isBox) {
        ui->BoxcheckBox->setChecked(true);
        int firstp = filename.indexOf("[");
        int lastp = filename.indexOf("]");
        QStringList range = filename.mid(firstp+1,lastp-firstp-1).split(",",QString::SkipEmptyParts);
        filename = filename.mid(0,firstp);
        ui->FitslineEdit->setText(filename);
        on_FitslineEdit_editingFinished();
        if (!range[0].contains("*")) {
            QStringList xrange = range[0].split(":",QString::SkipEmptyParts);
            ui->XminspinBox->setValue(xrange[0].toInt());
            ui->XmaxspinBox->setValue(xrange[1].toInt());
        }
        if (!range[1].contains("*")) {
            QStringList yrange = range[1].split(":",QString::SkipEmptyParts);
            ui->YminspinBox->setValue(yrange[0].toInt());
            ui->YmaxspinBox->setValue(yrange[1].toInt());
        }
        if (!range[2].contains("*")) {
            QStringList zrange = range[2].split(":",QString::SkipEmptyParts);
            ui->ZminspinBox->setValue(zrange[0].toInt());
            ui->ZmaxspinBox->setValue(zrange[1].toInt());
        }
    }
    else {
        ui->FitslineEdit->setText(filename);
        on_FitslineEdit_editingFinished();
    }

    QString outdir = QString::fromStdString(par->getOutfolder());
    if (outdir!="") ui->OutfolderlineEdit->setText(outdir);

    if (par->getflagGalFit()) {
        set3DFitFlag(Qt::Checked);
        readModelParam(par);
        GALFIT_PAR &p = par->getParGF();
        
        ui->IncDSpinBox->setValue(p.DELTAINC);
        ui->PaDSpinBox->setValue(p.DELTAPHI);

        std::string free = p.FREE;
        free = makelower(free);
        int found = free.find("vrot");
        if (found>=0) ui->vrotcheckBox->setChecked(true);
        else ui->vrotcheckBox->setChecked(false);
        found = free.find("vdisp");
        if (found>=0) ui->vdispcheckBox->setChecked(true);
        else ui->vdispcheckBox->setChecked(false);
        found = free.find("inc");
        if (found>=0) ui->inccheckBox->setChecked(true);
        else ui->inccheckBox->setChecked(false);
        found = free.find("pa");
        if (found>=0) ui->pacheckBox->setChecked(true);
        else ui->pacheckBox->setChecked(false);
        found = free.find("xpos");
        if (found>=0) ui->xposcheckBox->setChecked(true);
        else ui->xposcheckBox->setChecked(false);
        found = free.find("ypos");
        if (found>=0) ui->yposcheckBox->setChecked(true);
        else ui->yposcheckBox->setChecked(false);
        found = free.find("vsys");
        if (found>=0) ui->vsyscheckBox->setChecked(true);
        else ui->vsyscheckBox->setChecked(false);
        found = free.find("z0");
        if (found>=0) ui->z0checkBox->setChecked(true);
        else ui->z0checkBox->setChecked(false);
        found = free.find("vrad");
        if (found>=0) ui->vradcheckBox->setChecked(true);
        else ui->vradcheckBox->setChecked(false);

        ui->FtypecomboBox->setCurrentIndex(p.FTYPE-1);
        ui->WfunccomboBox->setCurrentIndex(p.WFUNC);

        ui->SecondstagecheckBox->setChecked(p.TWOSTAGE);
        if (p.TWOSTAGE) {
            if (makelower(p.POLYN)=="bezier") ui->PolynspinBox->setValue(-1);
            else ui->PolynspinBox->setValue(atoi(p.POLYN.c_str()));
        }

        if (p.flagERRORS) ui->ErrorsradioButton->setChecked(true);
        else ui->ErrorsradioButton->setChecked(false);
    }
    else {
        set3DFitFlag(Qt::Unchecked);
        if (par->getflagGalMod()) {
            set3DModelFlag(Qt::Checked);
            readModelParam(par);
        }
    }

    if (par->getflagSearch()) setSearchFlag(Qt::Checked);
    else setSearchFlag(Qt::Unchecked);
    if (getSearchFlag()) {
        SEARCH_PAR p = par->getParSE();
        if (p.searchType=="spectral") ui->SearchtypecomboBox->setCurrentIndex(0);
        else ui->SearchtypecomboBox->setCurrentIndex(1);
        if (p.UserThreshold) {
            ui->CuttypecomboBox->setCurrentIndex(0);
            ui->primaryCutlineEdit->setText(QString::number(p.threshold,'g',10));
        }
        else {
            ui->CuttypecomboBox->setCurrentIndex(1);
            ui->primaryCutlineEdit->setText(QString::number(p.snrCut,'g',10));
        }
        if (p.flagGrowth) {
            ui->GrowthcheckBox->setChecked(true);
            if (p.flagUserGrowthT) {
                ui->Cuttype2comboBox->setCurrentIndex(0);
                ui->SecondarycutlineEdit->setText(QString::number(p.growthThreshold,'g',10));
            }
            else {
                ui->Cuttype2comboBox->setCurrentIndex(1);
                ui->SecondarycutlineEdit->setText(QString::number(p.growthCut,'g',10));
            }
        }
        else ui->GrowthcheckBox->setChecked(false);

        ui->SearchAdvgroupBox->setEnabled(true);
        ui->SearchAdvgroupBox->setChecked(true);
        if (p.flagAdjacent) {
            ui->ThreshSpatialBox->setValue(-1);
            ui->ThreshSpatialBox->setDisabled(true);
            ui->AdjacentcheckBox->setChecked(true);
        }
        else {
            ui->AdjacentcheckBox->setChecked(false);
            ui->ThreshSpatialBox->setEnabled(true);
            ui->ThreshSpatialBox->setValue(p.threshSpatial);
        }

        ui->ThreshVelspinBox->setValue(p.threshVelocity);
        ui->MinChanspinBox->setValue(p.minChannels);
        ui->MinPixspinBox->setValue(p.minPix);
        ui->MinVoxspinBox->setValue(p.minVoxels);
        ui->RejectcheckBox->setChecked(p.RejectBeforeMerge);
        ui->MaxChanspinBox->setValue(p.maxChannels);
        ui->MaxAngSizeSpinBox->setValue(p.maxAngSize);
        ui->TwostagemergingcheckBox->setChecked(p.TwoStageMerging);
    }

    setSmoothFlag(Qt::Unchecked);
    if (par->getflagSmooth()) {
        setSmoothFlag(Qt::Checked);
        ui->SpatialSmoothgroupBox->setChecked(true);
        if (par->getFactor()!=-1) {
            ui->FactorcheckBox->setChecked(true);
            ui->FactordoubleSpinBox->setValue(par->getFactor());
        }
        else {
            ui->BeamgroupBox->setChecked(true);
            ui->NbmajSpinBox->setValue(par->getBmaj());
            ui->NbminSpinBox->setValue(par->getBmin());
            ui->NbpaSpinBox->setValue(par->getBpa());
            if (par->getOBmaj()!=-1) ui->ObmajSpinBox->setValue(par->getOBmaj());
            if (par->getOBmin()!=-1) ui->ObminSpinBox->setValue(par->getOBmin());
            if (par->getOBpa()!=0) ui->ObpaSpinBox->setValue(par->getOBmin());
        }
        ui->FFTcheckBox->setChecked(par->getflagFFT());
        ui->ReducecheckBox->setChecked(par->getflagReduce());
    } 
    if (par->getflagHanning()) {
        setSmoothFlag(Qt::Checked);
        ui->HanninggroupBox->setChecked(true);
        size_t hanning_window = par->getHanningWindow();
        if (hanning_window%2==0) hanning_window+=1;
        ui->HanningspinBox->setValue(hanning_window);
    }
     
    ui->MomentmapsgroupBox->setChecked(par->getMaps());
    if (ui->MomentmapsgroupBox->isChecked()) {
        ui->ProfilecheckBox->setChecked(par->getGlobProf());
        ui->TotalmapcheckBox->setChecked(par->getTotalMap());
        ui->VfieldcheckBox->setChecked(par->getVelMap());
        ui->DispmapcheckBox->setChecked(par->getDispMap());
        ui->rmsmapcheckBox->setChecked(par->getRMSMap());
        if (makelower(par->getMapType())=="gaussian") ui->MaptypecomboBox->setCurrentIndex(1);
        else ui->MaptypecomboBox->setCurrentIndex(0);
    }

    ui->PVgroupBox->setChecked(par->getFlagPV());
    if (ui->PVgroupBox->isChecked()) {
        ui->PVXposSpinBox->setValue(par->getXPOS_PV());
        ui->PVYposSpinBox->setValue(par->getYPOS_PV());
        ui->PVPaSpinBox->setValue(par->getPA_PV());
    }

    bool masking = par->getMASK()=="NONE" ? false : true;
    ui->MaskgroupBox->setChecked(masking);
    if (ui->MaskgroupBox) {
        ui->BlankcutSpinBox->setValue(par->getBlankCut());
        if (par->getFactor()!=-1) ui->BlankfactordoubleSpinBox->setValue(par->getFactor());
    }
}

void BBaroloWindow::writeParamFile(QString file) {

    using namespace std;

    int n=20;
    std::ofstream out(file.toStdString().c_str());
    out << left;

    QString filename =  ui->FitslineEdit->text().trimmed();
    if (ui->BoxcheckBox->isChecked()) {
        filename.append("["+QString::number(ui->XminspinBox->value())+":"+QString::number(ui->XmaxspinBox->value())+","
                           +QString::number(ui->YminspinBox->value())+":"+QString::number(ui->YmaxspinBox->value())+","
                           +QString::number(ui->ZminspinBox->value())+":"+QString::number(ui->ZmaxspinBox->value())+"]");
    }
    out << setw(n) << "FITSFILE" << filename.toStdString() << endl;

    QString outdir = ui->OutfolderlineEdit->text().trimmed();
    if (!outdir.isNull()) out << setw(n) << "OUTFOLDER" << outdir.toStdString() << endl;
    out << setw(n) << "THREADS" << ui->ThreadspinBox->value() << endl;
    
    if (get3DFitFlag()) {
        out << setw(n) << "3DFIT" << "true" << endl;
        writeModelParam(out);

        string FREE = "";
        if (ui->vrotcheckBox->isChecked()) FREE += "VROT ";
        if (ui->vdispcheckBox->isChecked()) FREE += "VDISP ";
        if (ui->inccheckBox->isChecked()) FREE += "INC ";
        if (ui->pacheckBox->isChecked()) FREE += "PA ";
        if (ui->xposcheckBox->isChecked()) FREE += "XPOS ";
        if (ui->yposcheckBox->isChecked()) FREE += "YPOS ";
        if (ui->vsyscheckBox->isChecked()) FREE += "VSYS ";
        if (ui->z0checkBox->isChecked()) FREE += "Z0 ";
        if (ui->vradcheckBox->isChecked()) FREE += "VRAD ";
        out << setw(n) << "FREE" << FREE << endl;

        if (ui->AdvancedgroupBox->isChecked()) {

            out << setw(n) << "FTYPE" << ui->FtypecomboBox->currentIndex()+1 << endl
                << setw(n) << "WFUNC" << ui->WfunccomboBox->currentIndex() << endl
                << setw(n) << "TOL" << ui->TollineEdit->text().toStdString() << endl;
            if (ui->SecondstagecheckBox->isChecked()) {
                out << setw(n) << "TWOSTAGE" << "true" << endl;
                if (ui->PolynspinBox->value()==-1)
                    out << setw(n) << "POLYN" << "bezier" << endl;
                else out << setw(n) << "POLYN" << ui->PolynspinBox->value() << endl;
            }
            else out << setw(n) << "TWOSTAGE" << "false" << endl;
            if (ui->ErrorsradioButton->isChecked())
                out << setw(n) << "flagErrors" << "true" << endl;
        }
    }
    else {
        out << setw(n) << "3DFIT" << "false" << endl;
        if (get3DModelFlag()) {
            out << setw(n) << "GALMOD" << "true" << endl;
            writeModelParam(out);
        }
    }

    if (getSearchFlag()) {
        out << setw(n) << "SEARCH" << "true" << endl;
        out << setw(n) << "searchType";
        if (ui->SearchtypecomboBox->currentIndex()==0) out << "spectral" << endl;
        else out << "spatial" << endl;
        if (ui->CuttypecomboBox->currentIndex()==0) out << setw(n) << "threshold";
        else out << setw(n) << "snrCut";
        out << ui->primaryCutlineEdit->text().toStdString() << endl;
        if (ui->GrowthcheckBox->isChecked()) {
            out << setw(n) << "flagGrowth" << "true" << endl;
            if (ui->Cuttype2comboBox->currentIndex()==0) out << setw(n) << "growthThreshold";
            else out << setw(n) << "growthCut";
            out << ui->SecondarycutlineEdit->text().toStdString() << endl;
        }
        if (ui->SearchAdvgroupBox->isChecked()) {
            if (ui->AdjacentcheckBox->isChecked()) {
                out << setw(n) << "flagAdjacent" << "true" << endl;
            }
             else {
                out << setw(n) << "flagAdjacent" << "false" << endl;
                if (ui->ThreshSpatialBox->value()!=-1)
                    out << setw(n) << "threshSpatial" << ui->ThreshSpatialBox->value() << endl;
                if (ui->ThreshVelspinBox->value()!=-1)
                    out << setw(n) << "threshVelocity"<< ui->ThreshVelspinBox->value() << endl;
            }
        }
        if (ui->MinPixspinBox->value()!=-1)
            out << setw(n) << "minPix" << ui->MinPixspinBox->value() << endl;
        if (ui->MinChanspinBox->value()!=-1)
            out << setw(n) << "minChannels" << ui->MinChanspinBox->value() << endl;
        if (ui->MinVoxspinBox->value()!=-1)
            out << setw(n) << "minVoxels" << ui->MinVoxspinBox->value() << endl;
        if (ui->MaxChanspinBox->value()!=-1)
            out << setw(n) << "maxChannels" << ui->MaxChanspinBox->value() << endl;
        if (ui->MaxAngSizeSpinBox->value()!=-1)
            out << setw(n) << "maxAngsize" << ui->MaxAngSizeSpinBox->value() << endl;
        if (ui->RejectcheckBox->isChecked())
            out << setw(n) << "RejectBeforeMerge" << "true" << endl;
        else out << setw(n) << "RejectBeforeMerge" << "false" << endl;
        if (ui->TwostagemergingcheckBox->isChecked())
            out << setw(n) << "TwoStageMerging" << "true" << endl;
        else out << setw(n) << "TwoStageMerging" << "false" << endl;
    }
    else out << setw(n) << "SEARCH" << "false" << endl;

    if (getSmoothFlag()) {
        if (ui->SpatialSmoothgroupBox->isChecked()) {
            out << setw(n) << "SMOOTH" << "true" << endl;
            if (ui->FactorcheckBox->isChecked()) {
                out << setw(n) << "Factor" << ui->FactordoubleSpinBox->value() << endl;
                if (ui->ReducecheckBox->isChecked())
                    out << setw(n) << "Reduce" << "true" << endl;
                else out << setw(n) << "Reduce" << "false" << endl;
            }
            else {
                out << setw(n) << "BMAJ" << ui->NbmajSpinBox->value() << endl
                    << setw(n) << "BMIN" << ui->NbminSpinBox->value() << endl
                    << setw(n) << "BPA"  << ui->NbpaSpinBox->value() << endl
                    << setw(n) << "OBMAJ"<< ui->ObmajSpinBox->value() << endl
                    << setw(n) << "OBMIN"<< ui->ObminSpinBox->value() << endl
                    << setw(n) << "OBPA" << ui->ObpaSpinBox->value() << endl;
            }
            if (ui->FFTcheckBox->isChecked())
                out << setw(n) << "FFT" << "true" << endl;
            else out << setw(n) << "FFT" << "false" << endl;
            if (ui->SmoothOutlineEdit->text()!="")
                out << setw(n) << "SmoothOutput" << ui->SmoothOutlineEdit->text().toStdString() << endl;
            }
        else out << setw(n) << "SMOOTH" << "false" << endl;
        
        if (ui->HanninggroupBox->isChecked()) {
            out << setw(n) << "HANNING" << "true" << endl;
            out << setw(n) << "HANNING_SIZE" << ui->HanningspinBox->value();
        }
        else out << setw(n) << "HANNING" << "false" << endl;
    }


    if (ui->MomentmapsgroupBox->isChecked()) {
       if (ui->ProfilecheckBox->isChecked())
            out << setw(n) << "globalProfile" << "true" << endl;
       if (ui->TotalmapcheckBox->isChecked())
            out << setw(n) << "totalMap" << "true" << endl;
       if (ui->VfieldcheckBox->isChecked())
            out << setw(n) << "velocityMap" << "true" << endl;
       if (ui->DispmapcheckBox->isChecked())
            out << setw(n) << "dispersionMap" << "true" << endl;
       if (ui->rmsmapcheckBox->isChecked())
            out << setw(n) << "rmsMap" << "true" << endl;
    }

    if (ui->PVgroupBox->isChecked()) {
        out << setw(n) << "flagPV" << "true" << endl;
        out << setw(n) << "XPOS_PV" << ui->PVXposSpinBox->value() << endl;
        out << setw(n) << "YPOS_PV" << ui->PVYposSpinBox->value() << endl;
        out << setw(n) << "PA_PV"   << ui->PVPaSpinBox->value() << endl;
    }

    if (ui->MaskgroupBox->isChecked()) {
        out << setw(n) << "blankCube" << "true" << endl
            << setw(n) << "factor" << ui->BlankfactordoubleSpinBox->value() << endl
            << setw(n) << "blankCut" << ui->BlankcutSpinBox->value() << endl;
    }
    else out << setw(n) << "blankCube" << "false" << endl;

    out << setw(n) << "showBar" << "false" << endl;
    out.close();
}

void BBaroloWindow::readModelParam(Param *p) {

    ui->FreeParametersframe->setEnabled(true);
    GALFIT_PAR &par = p->getParGF();
    QString s;
    int col = -1;
    if (par.RADII!="-1") {
        s = QString::fromStdString(par.RADII);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->RadsepFilelineEdit,ui->RadsepFilespinBox,ui->RadsepSpinBox,false);
            ui->RadsepFilelineEdit->setText(s);
            ui->RadsepFilespinBox->setValue(col);
        }
    }
    else {
        if (par.NRADII==-1) ui->NringscheckBox->setChecked(true);
        else {
            ui->NringscheckBox->setChecked(false);
            ui->NringsspinBox->setValue(par.NRADII);
        }
        if (par.RADSEP==-1) ui->RadsepcheckBox->setChecked(true);
        else {
            ui->RadsepcheckBox->setChecked(false);
            ui->RadsepSpinBox->setValue(par.RADSEP);
        }
    }
    if (par.XPOS=="-1") ui->XposcheckBox->setChecked(true);
    else {
        ui->XposcheckBox->setChecked(false);
        ui->wcscomboBox->setCurrentIndex(0);
        s = QString::fromStdString(par.XPOS);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->XposFilelineEdit,ui->XposFilespinBox,ui->yposlineEdit,false);
            ui->XposFilelineEdit->setText(s);
            ui->XposFilespinBox->setValue(col);
        }
        else {
            std::string xpos = par.XPOS;
            int found = xpos.find("d");
            if (found!=-1) {
                xpos.erase(found,xpos.size()-1);
                ui->wcscomboBox->setCurrentIndex(1);
            }
            else {
                found = xpos.find(':');
                if (found!=-1) ui->wcscomboBox->setCurrentIndex(2);

            }
            ui->xposlineEdit->setText(QString::fromStdString(xpos));
        }
    }
    if (par.YPOS=="-1") ui->YposcheckBox->setChecked(true);
    else {
        ui->YposcheckBox->setChecked(false);
        ui->wcscomboBox->setCurrentIndex(0);
        s = QString::fromStdString(par.YPOS);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->YposFilelineEdit,ui->YposFilespinBox,ui->yposlineEdit,false);
            ui->YposFilelineEdit->setText(s);
            ui->YposFilespinBox->setValue(col);
        }
        else {
            std::string ypos = par.YPOS;
            int found = ypos.find('d');
            if (found!=-1) {
                ypos.erase(found,ypos.size()-1);
                ui->wcscomboBox->setCurrentIndex(1);
            }
            else {
                found = ypos.find(':');
                if (found!=-1) ui->wcscomboBox->setCurrentIndex(2);
            }
            ui->yposlineEdit->setText(QString::fromStdString(ypos));
        }
    }


    if (par.VSYS=="-1") ui->VsyscheckBox->setChecked(true);
    else {
        ui->VsyscheckBox->setChecked(false);
        s = QString::fromStdString(par.VSYS);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->VsysFilelineEdit,ui->VsysFilespinBox,ui->VsysSpinBox,false);
            ui->VsysFilelineEdit->setText(s);
            ui->VsysFilespinBox->setValue(col);
        }
        else ui->VsysSpinBox->setValue(atof(par.VSYS.c_str()));
    }
    if (par.VROT=="-1") ui->VrotcheckBox->setChecked(true);
    else {
        ui->VrotcheckBox->setChecked(false);
        s = QString::fromStdString(par.VROT);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->VrotFilelineEdit,ui->VrotFilespinBox,ui->VrotSpinBox,false);
            ui->VrotFilelineEdit->setText(s);
            ui->VrotFilespinBox->setValue(col);
        }
        else ui->VrotSpinBox->setValue(atof(par.VROT.c_str()));
    }
    if (par.VDISP=="-1") ui->VdispcheckBox->setChecked(true);
    else {
        ui->VdispcheckBox->setChecked(false);
        s = QString::fromStdString(par.VDISP);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->VdispFilelineEdit,ui->VdispFilespinBox,ui->VdispSpinBox,false);
            ui->VdispFilelineEdit->setText(s);
            ui->VdispFilespinBox->setValue(col);
        }
        else ui->VdispSpinBox->setValue(atof(par.VDISP.c_str()));
    }
    if (par.VRAD=="-1") ui->VradcheckBox->setChecked(true);
    else {
        ui->VradcheckBox->setChecked(false);
        s = QString::fromStdString(par.VRAD);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->VradFilelineEdit,ui->VradFilespinBox,ui->VradSpinBox,false);
            ui->VradFilelineEdit->setText(s);
            ui->VradFilespinBox->setValue(col);
        }
        else ui->VradSpinBox->setValue(atof(par.VRAD.c_str()));
    }
    if (par.INC=="-1") ui->InccheckBox->setChecked(true);
    else {
        ui->InccheckBox->setChecked(false);
        s = QString::fromStdString(par.INC);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->IncFilelineEdit,ui->IncFilespinBox,ui->IncSpinBox,false);
            ui->IncFilelineEdit->setText(s);
            ui->IncFilespinBox->setValue(col);
        }
        else ui->IncSpinBox->setValue(atof(par.INC.c_str()));
    }
    if (par.PHI=="-1") ui->PacheckBox->setChecked(true);
    else {
        ui->PacheckBox->setChecked(false);
        s = QString::fromStdString(par.PHI);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->PaFilelineEdit,ui->PaFilespinBox,ui->PaSpinBox,false);
            ui->PaFilelineEdit->setText(s);
            ui->PaFilespinBox->setValue(col);
        }
        else ui->PaSpinBox->setValue(atof(par.PHI.c_str()));
    }
    if (par.Z0=="-1") ui->Z0checkBox->setChecked(true);
    else {
        ui->Z0checkBox->setChecked(false);
        s = QString::fromStdString(par.Z0);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->Z0FilelineEdit,ui->Z0FilespinBox,ui->Z0SpinBox,false);
            ui->Z0FilelineEdit->setText(s);
            ui->Z0FilespinBox->setValue(col);
        }
        else ui->Z0SpinBox->setValue(atof(par.Z0.c_str()));
    }
    if (par.DENS=="-1") ui->DenscheckBox->setChecked(true);
    else {
        ui->DenscheckBox->setChecked(false);
        s = QString::fromStdString(par.DENS);
        col = readFileString(s);
        if (col != -1) {
            Hide_3DFit_file(ui->DensFilelineEdit,ui->DensFilespinBox,ui->DensSpinBox,false);
            ui->DensFilelineEdit->setText(s);
            ui->DensFilespinBox->setValue(col);
        }
        else ui->DensSpinBox->setValue(atof(par.DENS.c_str()));
    }
    
    ui->Restframe->setEnabled(true);
    if (par.REDSHIFT!=0) ui->redshiftlineEdit->setText(QString::number(par.REDSHIFT));
    if (par.RESTWAVE[0]!=-1) {
        ui->restcomboBox->setCurrentIndex(0);
        ui->restlineEdit->setText(QString::number(par.RESTWAVE[0]));
    }
    else if (par.RESTFREQ[0]!=-1) {
        ui->restcomboBox->setCurrentIndex(1);
        ui->restlineEdit->setText(QString::number(par.RESTFREQ[0]));
    }
    
    ui->AdvancedgroupBox->setChecked(true);
    ui->AdvancedgroupBox->setEnabled(true);
    ui->LtypecomboBox->setCurrentIndex(par.LTYPE-1);
    ui->CdensSpinBox->setValue(par.CDENS);
    ui->NVspinBox->setValue(par.NV);
    if (p->getMASK()=="SMOOTH") {
        ui->MaskcomboBox->setCurrentIndex(0);
        ui->MaskgroupBox->setChecked(true);
        ui->BlankfactordoubleSpinBox->setValue(p->getFactor());
        ui->BlankcutSpinBox->setValue(p->getParSE().snrCut);
    }
    else if (p->getMASK()=="SEARCH") ui->MaskcomboBox->setCurrentIndex(1);
    else if (p->getMASK()=="SMOOTH&SEARCH") ui->MaskcomboBox->setCurrentIndex(2);
    else if (p->getMASK()=="THRESHOLD") {
        ui->MaskcomboBox->setCurrentIndex(3);
        ui->MaskThreshSpinBox->setValue(p->getParSE().threshold);
    }
    else if (p->getMASK()=="NEGATIVE") ui->MaskcomboBox->setCurrentIndex(4);
    else ui->MaskcomboBox->setCurrentIndex(5);

    if (par.NORM=="LOCAL") ui->NormcomboBox->setCurrentIndex(0);
    else if (par.NORM=="AZIM") ui->NormcomboBox->setCurrentIndex(1);
    else ui->NormcomboBox->setCurrentIndex(2);
    
    if (par.flagADRIFT) ui->AdriftcheckBox->setChecked(true);
}

void BBaroloWindow::writeModelParam(std::ofstream &out) {

    int n=20;
    if (!ui->RadsepFilelineEdit->isHidden()) {
        out << setw(n) << "RADII" << getFileString(ui->RadsepFilelineEdit,ui->RadsepFilespinBox) << endl;
    }
    else {
        if (!ui->NringscheckBox->isChecked())
            out << setw(n) << "NRADII" << ui->NringsspinBox->value() << endl;
        if (!ui->RadsepcheckBox->isChecked())
            out << setw(n) << "RADSEP" << ui->RadsepSpinBox->value() << endl;
    }
    if (!ui->XposcheckBox->isChecked()) {
        if (ui->BoxcheckBox->isChecked() && ui->wcscomboBox->currentIndex()==0) {
            float xp = ui->xposlineEdit->text().toFloat()-ui->XminspinBox->value()+1;
            out << setw(n)  << "XPOS" << xp;
        }
        else out << setw(n) << "XPOS" << ui->xposlineEdit->text().toStdString();
        if (ui->wcscomboBox->currentIndex()==1) out << "d";
        out << endl;
    }
    if (!ui->YposcheckBox->isChecked()) {
        if (ui->BoxcheckBox->isChecked() && ui->wcscomboBox->currentIndex()==0) {
            float yp = ui->yposlineEdit->text().toFloat()-ui->YminspinBox->value()+1;
            out << setw(n)  << "YPOS" << yp;
        }
        else out << setw(n) << "YPOS" << ui->yposlineEdit->text().toStdString();
        if (ui->wcscomboBox->currentIndex()==1) out << "d";
        out << endl;
    }
    if (!ui->VsyscheckBox->isChecked()) {
        out << setw(n) << "VSYS";
        if (!ui->VsysFilelineEdit->isHidden()) out << getFileString(ui->VsysFilelineEdit,ui->VsysFilespinBox);
        else out << ui->VsysSpinBox->value();
        out << endl;
    }
    if (!ui->VrotcheckBox->isChecked()) {
        out << setw(n) << "VROT";
        if (!ui->VrotFilelineEdit->isHidden()) out << getFileString(ui->VrotFilelineEdit,ui->VrotFilespinBox);
        else out << ui->VrotSpinBox->value();
        out << endl;
    }
    if (!ui->VdispcheckBox->isChecked()) {
        out << setw(n) << "VDISP";
        if (!ui->VdispFilelineEdit->isHidden()) out << getFileString(ui->VdispFilelineEdit,ui->VdispFilespinBox);
        else out << ui->VdispSpinBox->value();
        out << endl;
    }
    if (!ui->VradcheckBox->isChecked()) {
        out << setw(n) << "VRAD";
        if (!ui->VradFilelineEdit->isHidden()) out << getFileString(ui->VradFilelineEdit,ui->VradFilespinBox);
        else out << ui->VradSpinBox->value();
        out << endl;
    }
    if (!ui->InccheckBox->isChecked()) {
        out << setw(n) << "INC";
        if (!ui->IncFilelineEdit->isHidden()) out << getFileString(ui->IncFilelineEdit,ui->IncFilespinBox);
        else out << ui->IncSpinBox->value();
        out << endl;
        out << setw(n) << "DELTAINC" << ui->IncDSpinBox->value() << endl;
    }
    if (!ui->PacheckBox->isChecked()) {
        out << setw(n) << "PA";
        if (!ui->PaFilelineEdit->isHidden()) out << getFileString(ui->PaFilelineEdit,ui->PaFilespinBox);
        else out << ui->PaSpinBox->value();
        out << endl;
        out << setw(n) << "DELTAPA"  << ui->PaDSpinBox->value() << endl;
    }

    if (!ui->Z0checkBox->isChecked()) {
        out << setw(n) << "Z0";
        if (!ui->Z0FilelineEdit->isHidden()) out << getFileString(ui->Z0FilelineEdit,ui->Z0FilespinBox);
        else out << ui->Z0SpinBox->value();
        out << endl;
    }

    if (!ui->DenscheckBox->isChecked()) {
        out << setw(n) << "DENS";
        if (!ui->DensFilelineEdit->isHidden()) out << getFileString(ui->DensFilelineEdit,ui->DensFilespinBox);
        else out << ui->DensSpinBox->value();
        out << endl;
    }

    if (ui->redshiftlineEdit->text().toFloat()!=0) 
        out << setw(n) << "REDSHIFT" << ui->redshiftlineEdit->text().toStdString() << endl;
    
    if (ui->restlineEdit->text().toFloat()!=0) {
        if (ui->restcomboBox->currentIndex()==0) out << setw(n) << "RESTWAVE" << ui->restlineEdit->text().toStdString() << endl;
        if (ui->restcomboBox->currentIndex()==1) out << setw(n) << "RESTFREQ" << ui->restlineEdit->text().toStdString() << endl;
    } 
    

    if (ui->AdvancedgroupBox->isChecked()) {
        out << setw(n) << "LTYPE" << ui->LtypecomboBox->currentIndex()+1 << endl
            << setw(n) << "CDENS" << ui->CdensSpinBox->value() << endl;
        if (ui->NVspinBox->value()!=-1)
            out << setw(n) << "NV" << ui->NVspinBox->value() << endl;

        out << setw(n) << "MASK";
        if (ui->MaskcomboBox->currentIndex()==0) out << "SMOOTH" << endl;
        else if (ui->MaskcomboBox->currentIndex()==1) out << "SEARCH" << endl;
        else if (ui->MaskcomboBox->currentIndex()==2) out << "SMOOTH&SEARCH" << endl;
        else if (ui->MaskcomboBox->currentIndex()==3) {
            out << "THRESHOLD" << endl;
            out << setw(n) << "threshold" << ui->MaskThreshSpinBox->value() << endl;
        }
        else if (ui->MaskcomboBox->currentIndex()==4) out << "NEGATIVE" << endl;
        else out << "NONE" << endl;

        out << setw(n) << "NORM";
        if (ui->NormcomboBox->currentIndex()==0) out << "LOCAL";
        else if (ui->NormcomboBox->currentIndex()==1) out << "AZIM";
        else out << "NONE";
        out << endl;
    }

    if (ui->AdriftcheckBox->isChecked()) out << setw(n) << "ADRIFT" << "true";
}

void BBaroloWindow::plotParameters(){

    const unsigned nsubplot = ui->SecondstagecheckBox->isChecked() ? 2 : 1;
    const unsigned MAXPAR = 9;

    for (unsigned j=0; j<nsubplot; j++) {
        std::vector<std::vector<double> > allData;
        std::string ringlog = out_path.toStdString()+"ringlog"+(QString::number(j+1)).toStdString()+".txt";
        if (!getData(allData,ringlog,false)) continue;

        int nr = allData.size();
        std::vector<QVector<double> > data(9);
        float max_val[MAXPAR] = {0,0,0,0,0,0,0,0,-300000};
        float min_val[MAXPAR] = {0,0,0,90,360,10000,100000,100000,300000};
        QString labels[MAXPAR] = {"Radius (arcs)","Vrot (km/s)","Disp. (km/s)",
                                  "Inc. (deg)","P.A. (deg)","Z0 (arcs)","Xpos (pix)","Ypos (pix)","Vsys (km/s)"};
        QCustomPlot *plt[4] = {ui->plot1,ui->plot2,ui->plot3,ui->plot4};
        int index[4] = {ui->Plot1comboBox->currentIndex()+1, ui->Plot2comboBox->currentIndex()+1,
                        ui->Plot3comboBox->currentIndex()+1, ui->Plot4comboBox->currentIndex()+1};

        for (int i=0; i<nr-1;i++) {
            data[0].append(allData[i+1][1]);     // Radius
            data[1].append(allData[i+1][2]);     // VROT
            data[2].append(allData[i+1][3]);     // VDISP
            data[3].append(allData[i+1][4]);     // INC
            data[4].append(allData[i+1][5]);     // PA
            data[5].append(allData[i+1][6]);     // Z0
            data[6].append(allData[i+1][9]);     // XPOS
            data[7].append(allData[i+1][10]);    // YPOS
            data[8].append(allData[i+1][11]);    // VSYS
            if (data[0][i]>max_val[0])      max_val[0]=data[0][i];
            if (data[1][i]>max_val[1])      max_val[1]=data[1][i];
            if (data[2][i]>max_val[2])      max_val[2]=data[2][i];
            if (data[3][i]>max_val[3])      max_val[3]=data[3][i];
            else if (data[3][i]<min_val[3]) min_val[3]=data[3][i];
            if (data[4][i]>max_val[4])      max_val[4]=data[4][i];
            else if (data[4][i]<min_val[4]) min_val[4]=data[4][i];
            if (data[5][i]>max_val[5])      max_val[5]=data[5][i];
            else if (data[5][i]<min_val[5]) min_val[5]=data[5][i];
            if (data[6][i]>max_val[6])      max_val[6]=data[6][i];
            else if (data[6][i]<min_val[6]) min_val[6]=data[6][i];
            if (data[7][i]>max_val[7])      max_val[7]=data[7][i];
            else if (data[7][i]<min_val[7]) min_val[7]=data[7][i];
            if (data[8][i]>max_val[8])      max_val[8]=data[8][i];
            else if (data[8][i]<min_val[8]) min_val[8]=data[8][i];
        }

        if (!(ui->RadsepcheckBox->isChecked() && ui->NringscheckBox->isChecked()))
            max_val[0] = ui->RadsepSpinBox->value()*ui->NringsspinBox->value();

        QPen pen;
        pen.setWidth(2);
        if (j==1 || nsubplot==1) pen.setColor(QColor(178,34,34));
        else pen.setColor(QColor(137,127,127));

        for (int i=0; i<4; i++) {
            plt[i]->addGraph();
            plt[i]->graph(j)->setData(data[0],data[index[i]]);
            plt[i]->graph(j)->setPen(pen);
            plt[i]->graph(j)->setLineStyle(QCPGraph::lsNone);
            plt[i]->graph(j)->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssPlusCircle, 4));
            plt[i]->xAxis->setLabel(labels[0]);
            plt[i]->yAxis->setLabel(labels[index[i]]);
            plt[i]->xAxis->setRange(0, 1.10*max_val[0]);
            plt[i]->yAxis->setRange(0.95*min_val[index[i]], 1.05*max_val[index[i]]);
        }
    }
    ui->plot1->replot();
    ui->plot2->replot();
    ui->plot3->replot();
    ui->plot4->replot();

}
