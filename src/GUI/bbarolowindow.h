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

#ifndef BBAROLOWINDOW_H
#define BBAROLOWINDOW_H

#include <QMainWindow>
#include <QProcess>
#include <QtGui>
#include <QSpinBox>
#include <QListWidgetItem>

#include "../Arrays/param.hh"

namespace Ui {
class BBaroloWindow;
}

class BBaroloWindow : public QMainWindow
{
    Q_OBJECT
    
public:
    explicit BBaroloWindow(QWidget *parent = 0);
    ~BBaroloWindow();

private slots:
    void on_ParampushButton_clicked();
    void on_FitspushButton_clicked();
    void on_FitslineEdit_editingFinished();
    void on_BoxcheckBox_stateChanged();
    void on_AutocheckBox_stateChanged();
    void on_ResetpushButton_clicked() {resetGUI();};


    // Slots for 3Dfit panel
    void on_NringscheckBox_stateChanged();
    void on_RadsepcheckBox_stateChanged();
    void on_XposcheckBox_stateChanged();
    void on_YposcheckBox_stateChanged();
    void on_wcscomboBox_activated(int index);
    void on_VsyscheckBox_stateChanged();
    void on_VrotcheckBox_stateChanged();
    void on_VdispcheckBox_stateChanged();
    void on_InccheckBox_stateChanged();
    void on_PacheckBox_stateChanged();
    void on_Z0checkBox_stateChanged();
    void on_SecondstagecheckBox_stateChanged();
    void on_MasktoolButton_clicked();
    void on_MaskcomboBox_currentIndexChanged(int index);
    void on_RadseppushButton_clicked();
    void on_RadsepFilelineEdit_editingFinished();
    void on_XpospushButton_clicked();
    void on_XposFilespinBox_editingFinished();
    void on_YposFilelineEdit_editingFinished();
    void on_YpospushButton_clicked();
    void on_VsyspushButton_clicked();
    void on_VsysFilelineEdit_editingFinished();
    void on_VrotpushButton_clicked();
    void on_VrotFilelineEdit_editingFinished();
    void on_VdisppushButton_clicked();
    void on_VdispFilelineEdit_editingFinished();
    void on_IncpushButton_clicked();
    void on_IncFilelineEdit_editingFinished();
    void on_PapushButton_clicked();
    void on_PaFilelineEdit_editingFinished();
    void on_Z0pushButton_clicked();
    void on_Z0FilelineEdit_editingFinished();
    void on_DenspushButton_clicked();
    void on_DensFilelineEdit_editingFinished();
    void on_GalfittoolButton_clicked();
    void on_GalfitAdvtoolButton_clicked();
    void Hide_3DFit_file(QLineEdit *le, QSpinBox *sp, QWidget *qw, bool Hide);
    void Hide_All_3DFit_file(bool Hide);
    void Selected_3DFit_file (QLineEdit *le, QSpinBox *sb, QWidget *qw);


    // Slots for Search panel
    void on_SearchgroupBox_clicked();
    void on_GrowthcheckBox_stateChanged();
    void on_SearchAdvgroupBox_clicked();
    void on_AdjacentcheckBox_stateChanged(int arg1);
    void on_SearchtoolButton_clicked();
    void on_SearchAdvtoolButton_clicked();


    // Slots for Smooth panel
    void on_BeamgroupBox_clicked();
    void on_FactorcheckBox_stateChanged();
    void on_NbmajSpinBox_valueChanged(double arg1);
    void on_SmoothOutpushButton_clicked();
    void on_FactordoubleSpinBox_valueChanged(double arg1);
    void on_ReducecheckBox_clicked(bool checked);

    // Menu bar slots
    void on_actionOpen_FITS_file_triggered(){on_FitspushButton_clicked();}
    void on_actionOpen_parameter_file_triggered() {on_ParampushButton_clicked();}
    void on_actionExport_parameter_file_triggered();
    void on_actionReset_GUI_triggered() {resetGUI();};


    // Running slots
    void on_RunpushButton_clicked();
    void on_KillpushButton_clicked();
    void updateText();
    void updateExit();

    // Combobox plotting slots
    void on_Plot1comboBox_activated(int) {plotParameters();};
    void on_Plot2comboBox_activated(int) {plotParameters();};
    void on_Plot3comboBox_activated(int) {plotParameters();};
    void on_Plot4comboBox_activated(int) {plotParameters();};



    ///////////////////
    void switch_Galmod_3Dfit(bool toHide);
    int getCurrentRow(QListWidgetItem *item);
    void on_listWidget_currentRowChanged(int currentRow);

    void on_listWidget_itemChanged(QListWidgetItem *item);

    void on_listWidget_itemDoubleClicked(QListWidgetItem *item);


    void on_HidetoolButton_clicked();

    void set3DFitFlag(Qt::CheckState);
    bool get3DFitFlag();
    void set3DModelFlag(Qt::CheckState);
    bool get3DModelFlag();
    void setSearchFlag(Qt::CheckState);
    bool getSearchFlag();
    void setSmoothFlag(Qt::CheckState);
    bool getSmoothFlag();
    void setMapsFlag(Qt::CheckState);
    bool getMapsFlag();




    void on_OutfolderpushButton_clicked();

private:
    void readParamFromFile(std::string filein);
    void readModelParam(Param *par);
    void writeParamFile(QString file="param.par");
    void writeModelParam(std::ofstream &out);
    std::string getFileString(QLineEdit *le, QSpinBox *sb);
    int readFileString(QString &s);


    void resetGUI();
    void disable_All();
    void enable_All();

    void plotParameters();

    Ui::BBaroloWindow *ui;
    QProcess* proc;
    QString obj;
    QString out_path;



};

#endif // BBAROLOWINDOW_H
