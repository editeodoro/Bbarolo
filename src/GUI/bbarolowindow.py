# QT-related modules
#import moc.icons_rc
from PyQt5 import QtGui,QtWidgets,QtCore
from ui_bbarolowindow import Ui_BBaroloWindow


class MainWindow(QtWidgets.QMainWindow, Ui_BBaroloWindow):
    
    def __init__(self, parent=None,args=None):
        super(MainWindow, self).__init__(parent)
        self.setupUi(self)
''''  
        # Parameters for login
        self.auth  = {'username': None, 'password': None, 'valid': False}
        self.subID = {'ID': None, 'valid': False}
        self.req = requests.Session()
        self.isrunning = False
                
        # Setting Icons
        self.icon_success = QtGui.QPixmap(":icons/success.png") 
        self.icon_failure = QtGui.QPixmap(":icons/unsuccess.png")
        #icon = self.style().standardIcon(QtWidgets.QStyle.SP_DialogApplyButton);
        #pixmap = icon.pixmap(QtCore.QSize(64, 64));
        #self.LoginImage.setPixmap(pixmap)
        
        # The multiprocess that handle the submission
        self.thread = arXivSubmission()
        self.thread.finished.connect(self.ExitSuccess)
         
        # A error message dialog
        self.errMsg = QtWidgets.QMessageBox(self)
        self.errMsg.setIcon(QtWidgets.QMessageBox.Critical)
        self.errMsg.setWindowTitle("ArXiv Submitter Error")        
        self.errMsg.finished.connect(self.ResetErrMsg)
        
        # Message Dialog when all went well
        self.successMsg = QtWidgets.QMessageBox(self)
        self.successMsg.setStandardButtons(QtWidgets.QMessageBox.Close)
        self.successMsg.finished.connect(self.CloseSuccessMsg)
        
        # Get the local time of next submission and set it
        self.nextsubtime = QtCore.QDateTime.currentDateTime()
        self.dateTimeEdit.setMinimumDate(self.nextsubtime.date())
        self.setNextSubmissionTime()
        self.dateTimeEdit.setDateTime(self.nextsubtime)
        
        # A countdown to the submission and clock
        self.timer = QtCore.QTimer()
        self.timer.timeout.connect(self.countdown)
        self.clock = QtCore.QTimer()
        self.clock.timeout.connect(self.currentTime)
        self.clock.start(1000)
        
        # Set some fixed size objects 
        #self.setFixedSize(500,500)
        self.MinspinBox.setFixedSize(30,24)
        self.LoginImage.setFixedSize(20,20)
        #self.textEdit.setFixedSize(100,50)
        self.IDImage.setFixedSize(20,20)
        # Set approriate font size for some GUI elements
        font = self.LatencytoolButton.font()
        font.setPointSize(font.pointSize()-1)
        self.LatencytoolButton.setFont(font)
        self.ClocktoolButton.setFont(font)
        self.Latencylabel.setFont(font)
        self.Clocklabel.setFont(font)
        self.IDgroupBox.setStyleSheet("QGroupBox { font-weight: bold; } ")
        self.SubgroupBox.setStyleSheet("QGroupBox { font-weight: bold; } ")
        self.PrefgroupBox.setStyleSheet("QGroupBox { font-weight: bold; } ")
        
        # Windows for managing browsers and drivers
        self.windri  = DriverWindow(self)
        self.setBrowsers()
        
        # My signature
        sign = QtWidgets.QLabel("Designed by EdT - 2017   ")
        sign = QtWidgets.QLabel("<a href='mailto:enrico.diteodoro@gmail.com' style='color:#808080; text-decoration:none' >Designed by EdT - 2017   </a>")
        sign.setOpenExternalLinks(True)
        sign.setAlignment(QtCore.Qt.AlignRight)
        sign.setStyleSheet("font: italic 11pt")
        sign.setToolTip("Mail to EdT")
        self.statusBar.addWidget(sign)
        
        # Slots
        self.PWcheckBox.clicked[bool].connect(self.ShowHidePassword)
        self.PWlineEdit.editingFinished.connect(self.checkCredentials)
        self.UserlineEdit.editingFinished.connect(self.checkCredentials)
        self.SubIDspinBox.editingFinished.connect(self.checkSubID)
        self.SubmitpushButton.clicked.connect(self.SubmitButtonClicked)        
        self.DriverstoolButton.clicked.connect(self.OpenDriverWindow)
        self.BrowsercomboBox.currentIndexChanged[str].connect(self.BrowserChanged)
        self.LatencytoolButton.clicked.connect(self.testLatency)
        self.ClocktoolButton.clicked.connect(self.testClock)
        
        # Hide preference box and set up action to unhide it
        self.showSettings(False)
        self.actionShow_advanced_settings.triggered[bool].connect(self.showSettings)
        
        # Load stuff from the command line, if present
        if args:
            if args.user: self.UserlineEdit.setText(args.user)
            if args.password: self.PWlineEdit.setText(args.password)
            if args.id: self.SubIDspinBox.setValue(args.id)
            #if args.time: ...
            #if args.browser: ...
            if args.nmin: self.MinspinBox.setValue(args.nmin)
            #if args.test: self.TestcheckBox.setChecked(args.test)
            if args.nodelaycorrection: self.CorrectDelaycheckBox.setChecked(False)
            
        
    
    def countdown(self):
        now = QtCore.QDateTime.currentDateTime()
        delta = now.secsTo(self.subtime)
        if delta >= 0:
            h,m,s = int(delta/3600.), int(delta/60.)%60, delta%60.
            self.Countdownlabel.setText('<font color="#B22222"> Time to submission: %02d:%02d:%02d</font>'%(h,m,s)) 
        else: self.timer.stop()
 
        
    def currentTime(self):
        now = QtCore.QDateTime.currentDateTime()
        self.currentTimelabel.setText('<font color="#038d49"> Current time: '+now.toString('hh:mm:ss (dd MMM yyyy)')) 
        if now > self.nextsubtime:
            self.setNextSubmissionTime()
            
    
    def setBrowsers(self):
        """
        This function updates the BrowsercomboBox based on the available drivers
        registered in self.windri.manager class
        """
        self.BrowsercomboBox.clear()
        drivers = self.windri.manager.drivers
        for key, value in drivers.items():
            if value[1]: self.BrowsercomboBox.addItem(key)
         
        # Defaulting to PhantomJS
        if drivers["Phantom JS"][1]: 
            self.BrowsercomboBox.setCurrentIndex(self.BrowsercomboBox.count()-1)


    def testLatency(self):
        """ 
        This function tests the network latency and print the median
        latency and error (3sigma). The value denotes the connection time
        """
        self.Latencylabel.setText("Benchmarking arxiv.org...")
        self.Latencylabel.repaint()
        QtWidgets.QApplication.processEvents()
        try:
            av,med,var,mad,latency = getLatency("http://arxiv.org",10)
        except:
            self.Latencylabel.setText("Something went wrong. Try again")
        else:
            self.Latencylabel.setText("<b> (%.0f +/- %.0f) msec"%(med*1000,3*1000*1.486*mad))
        
    
    def testClock(self):
        """ 
        This function calculates print the delay between the local clock and the ntp.server
        """
        self.Clocklabel.setText("") 
        try:
            delay = getDelay()
        except:
            self.Clocklabel.setText("Something went wrong. Try again")
        else:
            self.Clocklabel.setText("<b> Delay ~ %.3f sec"%(delay.total_seconds()))
                
    
    def EnableDisableAll(self):
        """ Enabling or disabling GUI object according to what is happening """
        if self.isrunning:
            self.isrunning=False
            # Enable back the disabled elements
            self.IDgroupBox.setEnabled(True)
            self.SubgroupBox.setEnabled(True)
            self.SubgroupBox.setEnabled(True)
            self.PrefgroupBox.setEnabled(True)
            self.TestcheckBox.setEnabled(True)
            self.Countdownlabel.clear()
            self.Newsublabel.setEnabled(True)
            self.SubmitpushButton.setText("Submit!")
            self.textEdit.clear()
            self.LatencytoolButton.setEnabled(True)
            self.Latencylabel.setEnabled(True)
            self.ClocktoolButton.setEnabled(True)
            self.Clocklabel.setEnabled(True)
        else:    
            self.isrunning = True
            self.IDgroupBox.setDisabled(True)
            self.SubgroupBox.setDisabled(True)
            self.SubgroupBox.setDisabled(True)
            self.PrefgroupBox.setDisabled(True)
            self.TestcheckBox.setDisabled(True)
            self.Newsublabel.setDisabled(True)
            self.SubmitpushButton.setText("Stop")
            self.LatencytoolButton.setDisabled(True)
            self.Latencylabel.setDisabled(True)
            self.ClocktoolButton.setDisabled(True)
            self.Clocklabel.setDisabled(True)
            
            
    def ResetErrMsg(self):
        """ 
        Just reset the error dialog box 
        """
        self.errMsg.setIcon(QtWidgets.QMessageBox.Critical)
        self.errMsg.setInformativeText("")
        self.errMsg.setText("")
      
        
    def ExitSuccess(self):
        """ 
        Submission went ok, so let's just pop up a message box and reset 
        """
        self.successMsg.setWindowFlags(self.successMsg.windowFlags() & (~QtCore.Qt.WindowCloseButtonHint))
        
        # Add ok button to message box
        okbutton = self.successMsg.button(QtWidgets.QMessageBox.Close)
        # Waiting for 20s for the submission to be correctly processed
        okbutton.setDisabled(True)
        self.successMsg.setIcon(QtWidgets.QMessageBox.Information)
        self.successMsg.setText("Your submission is being processed")
        self.successMsg.setInformativeText("Please wait until further notice...")
        self.successMsg.show()
        
        def allGood():
            okbutton.setEnabled(True)
            self.successMsg.setIconPixmap(self.icon_success.scaled(60,60))
            self.successMsg.setWindowTitle("All done!")
            self.successMsg.setText("Congratulations!")
            self.successMsg.setInformativeText("Your paper has been submitted at %s. \n\nPlease wait for confirmation email and then press OK to close the browser."%(str(self.thread.submitted)))
        
        # Waiting 20s and then show allGood message
        QtCore.QTimer.singleShot(20000, allGood)
        
        
    
    def setNextSubmissionTime(self):
        """
        This function finds the localtime of next deadline for submission
        """
        localtime = QtCore.QDateTime.currentDateTime()
        timezone_ny = QtCore.QTimeZone(QtCore.QByteArray(b'America/New_York'))
        timenow_ny = localtime.toTimeZone(timezone_ny)

        if timenow_ny.time()<QtCore.QTime(14,0,0):
            subtime_ny  = QtCore.QDateTime(timenow_ny.date(),QtCore.QTime(14,0,0),timezone_ny)
        else: 
            subtime_ny  = QtCore.QDateTime(timenow_ny.date().addDays(1),QtCore.QTime(14,0,0),timezone_ny)
    
        # Skipping Saturdays and Sundays
        day = subtime_ny.date().dayOfWeek()
        if day==6: subtime_ny = subtime_ny.addDays(2)
        if day==7: subtime_ny = subtime_ny.addDays(1)
        
        # Converting to local time
        self.nextsubtime = subtime_ny.toTimeZone(localtime.timeZone())
        self.Newsublabel.setText('<font color="#0066CC"><b> Next deadline: %s'%(self.nextsubtime.toString("ddd dd/MM/yyyy hh:mm:ss")))
        
        
    @QtCore.pyqtSlot()
    def CloseSuccessMsg(self):
        self.thread.close()
        self.EnableDisableAll()
    
    
    @QtCore.pyqtSlot(bool)
    def ShowHidePassword(self,checked):
        """ 
        Show and hide password when checkbox clicked 
        """
        if checked: self.PWlineEdit.setEchoMode(QtWidgets.QLineEdit.Normal)
        else: self.PWlineEdit.setEchoMode(QtWidgets.QLineEdit.Password)


    @QtCore.pyqtSlot()
    def checkCredentials(self):
        """ 
        Check whether the inserted credentials are correct 
        """
        user, passw = self.UserlineEdit.text(), self.PWlineEdit.text()
        isUpdated = user!=self.auth['username'] or passw!=self.auth['password']
        if isUpdated:
            self.auth = {'username': user, 'password': passw, 'valid': False}
            self.LoginImage.clear()
            self.LoginImage.setToolTip("")
            if (user and passw):
                try:
                    # Here I send a request to arxiv to login
                    p = self.req.post("https://arxiv.org/user/login", data=self.auth)
                except: a = -2 # No internet connection
                else:   a = p.text.find("Login succeeded")

                if a==-2:
                    self.errMsg.setText("I could not contact arXiv server. \n Check your internet connection.")
                    self.errMsg.show()
                elif a==-1: # Login failed
                    self.LoginImage.setPixmap(self.icon_failure)
                    self.LoginImage.setToolTip("Login failed. Check your credentials")
                else: # Successful login
                    self.LoginImage.setPixmap(self.icon_success)
                    self.LoginImage.setToolTip("Login successful")
                    self.auth['valid'] = True
                    # Try to get self.ID automatically
                    if not self.subID["ID"] and self.SubIDspinBox.value()==0:
                        b = p.text.find('<a id="submit/')
                        c = p.text.find('" class="confirm_delete" href="https://arxiv.org/user/')
                        if b and c > 0:
                            idstr = p.text[b+14:c]
                            try: self.SubIDspinBox.setValue(float(idstr))
                            except: pass
    
    
    @QtCore.pyqtSlot()
    def checkSubID(self):
        """ 
        Check whether the inserted submission ID exists and is ready 
        """
        idsub = int(self.SubIDspinBox.value())
        isUpdated = idsub!=self.subID['ID']
        #if isUpdated or (not isUpdated and not self.subID['valid']):
        self.subID = {'ID': idsub, 'valid': False}
        self.IDImage.clear()
        self.IDImage.setToolTip("")
        if self.auth['valid'] and idsub!=0:
            try:
                # Here I send a request to arxiv to resume the submission
                p = self.req.get("https://arxiv.org/submit/"+str(idsub)+"/resume")
            except: 
                a = -2 # If no internet connection
            else:
                a = p.text.find("/"+str(idsub)+"/")
                    
                    
            if a==-2:
                self.errMsg.setText("I could not contact arXiv server. \n Check your internet connection.")
                self.errMsg.show()  
            elif a==-1: # Submission does not exist
                self.IDImage.setPixmap(self.icon_failure)
                self.IDImage.setToolTip("Submission does not exist. Check the ID")
            else: # Submission is ok
                 # Check whether the submission is ready to go (Submit button)
                b = p.text.find('form action="https://arxiv.org/submit/'+str(idsub)+'/submit"')
                c = p.text.find('div id="submit_button" style="display: none"')
                readytoSubmit = b!=-1 and c==-1
                if readytoSubmit:     
                    self.IDImage.setPixmap(self.icon_success)
                    self.IDImage.setToolTip("Submission ID is ok")
                    self.subID['valid'] = True
                else:
                    self.IDImage.setPixmap(self.icon_failure)
                    self.IDImage.setToolTip("Submission exists but is not ready")
                    self.errMsg.setText("Submission %s exists but..."%str(idsub))
                    self.errMsg.setInformativeText("... is not complete. Submit button in Preview page must be unhidden!")
                    self.errMsg.show()


    @QtCore.pyqtSlot()
    def OpenDriverWindow(self):
        self.windri.RestoreDefaults()   
        self.windri.show()
    
    
    @QtCore.pyqtSlot()
    def CloseDriverWindow(self):
        self.windri.close()
        self.windri.msg.setIcon(QtWidgets.QMessageBox.Information)
        self.windri.msg.setInformativeText("")
        self.windri.msg.setText("")
        self.setBrowsers()
        
    
    @QtCore.pyqtSlot(str)
    def BrowserChanged(self,brow):
        # Show a special message if Safari is selected.
        # Safari driver needs to be activated first
        if brow=="Safari":
            self.errMsg.setIcon(QtWidgets.QMessageBox.Warning)
            self.errMsg.setText("Make sure that automation is turned on in Safari or ArXivSubmitter will crash!")
            self.errMsg.setInformativeText("In Safari menu: Develop -> Allow Remote Automation")
            self.errMsg.show()
     
        
    @QtCore.pyqtSlot()
    def SubmitButtonClicked(self):
        """ 
        This function handles the submission procedure  
        """
        if self.isrunning:
            # The function is already running, so kill it and re-enable the GUI 
            self.timer.stop()
            # Disconnecting ExitSuccess slot and terminate thread. 
            # Not elegant, but I could not find better solution
            self.thread.finished.disconnect()
            self.thread.terminate()
            self.thread.wait()
            self.thread.close()
            self.thread.finished.connect(self.ExitSuccess)
            self.EnableDisableAll()
        else:
            self.checkSubID()
            readyToGo = self.auth['valid'] and self.subID['valid']
            if readyToGo:
                # Setting submission time
                self.subtime = self.dateTimeEdit.dateTime()
                timestr = self.subtime.toString("dd/MM/yy hh:mm:ss.zzz")
                now = QtCore.QDateTime.currentDateTime()
                # Checking if time is appropriate and starting countdown
                delta = now.secsTo(self.subtime)
                if delta<0:
                    self.errMsg.setText("Submission time is in the past:\n %s"%timestr)
                    self.errMsg.show()
                    return
                if delta<20:
                    self.errMsg.setText("Your submission time is in %d sec."%delta)
                    self.errMsg.setInformativeText("Please allow at least 20 seconds for the whole process.")
                    self.errMsg.show()
                    return
  
                # Starting countdown
                self.timer.start(1000)
                # Disable all fields
                self.EnableDisableAll()
                  
                # Now gathering all parameters needed        
                user, passw, subid = self.auth['username'], self.auth['password'], self.subID['ID']
                timestr = self.subtime.toString("dd/MM/yy hh:mm:ss.zzz")
                browser = self.BrowsercomboBox.currentText()
                minutes = self.MinspinBox.value() 
                subflag = not self.TestcheckBox.isChecked()
                delayflag = self.CorrectDelaycheckBox.isChecked()
                
                # Writing some output info
                if not subflag:
                    self.textEdit.append('<font color="#FF2828">THIS IS JUST A TEST! </font>')
                    self.textEdit.append('Paper %d will be displayed at %s'
                                         %(self.subID['ID'], self.subtime.toString("hh:mm:ss.zzz (dd/MM/yy)")))
                else:
                    self.textEdit.append("Paper %d will be submitted at %s"
                                         %(self.subID['ID'], self.subtime.toString("hh:mm:ss.zzz (dd/MM/yy)")))
                if delta>minutes*60 and not browser=="Phantom JS" :
                    self.textEdit.append("%s will open %d minutes before submission time."%(browser,minutes))
                if browser=="Phantom JS": self.textEdit.append("You can relax now, I'll take care of the rest.")
                else: self.textEdit.append("Please do not interact with the browser.")
                
                
                # Start the daemon for submission
                self.thread.setParams(user,passw,subid,timestr,browser,minutes,subflag,delayflag)
                self.thread.start()
                
            else:
                self.errMsg.setText("Please insert all the correct information before submitting.")
                self.errMsg.setInformativeText("Username, password and/or ID are not verified.")
                self.errMsg.show()
                

    @QtCore.pyqtSlot(bool)
    def showSettings(self,checked):
        if (checked): 
            self.PrefgroupBox.setHidden(False)
            self.TestcheckBox.setHidden(False)
        else: 
            self.TestcheckBox.setHidden(True)
            self.PrefgroupBox.setHidden(True)

        # Resizing the main window
        for i in range(0, 5):
              # The minimum size is not computed until some events are processed in the event loop. 
              # Just process the event loop for some iterations and then resize to minimum.
              QtWidgets.QApplication.processEvents()
        self.resize(self.minimumSizeHint())
'''