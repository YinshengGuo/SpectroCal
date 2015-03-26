import re

from guidata.qt.QtGui import QApplication,QWidget, QFileDialog,QPushButton,QMainWindow,QListWidget,QLabel, QTextBrowser,QLineEdit#, QVBoxLayout, QHBoxLayout
from guidata.qt.QtCore import SIGNAL,QRect

#---Import plot widget base class
from guiqwt.curve import CurvePlot
from guiqwt.plot import PlotManager
from guiqwt.builder import make
#from guidata.configtools import get_icon
#---
import numpy as np
from scipy.optimize import minimize
from scipy.stats import linregress


#==============================================================================
# to do list
# compile
#==============================================================================


class interface():
    def update_status_bar(self,message):
        self.statusbar.showMessage(message)
        
    def func_load_spectrum(self):
        self.spectrum_curve_plot.del_all_items()
        self.spectrum_file_path=QFileDialog.getOpenFileName() # file opening might fail
        self.spectrum_file_path=unicode(self.spectrum_file_path)
        try:
            cache=np.genfromtxt(self.spectrum_file_path) # throws exception if clicked 'cancel' in file dialog
        except IOError as err: 
            message="IOError occurred when reading calibration file: "+err.message
            self.update_status_bar(message)
            return
        except ValueError as err: # fixed: ref wavelength can go through
            message="ValueError occurred when reading calibration file: "+err.message
            self.update_status_bar(message)
            return 
        if cache.ndim==2:
            try: # throws exception if text file format is wrong
                self.xdata=cache[:,0] # 1st col is pixel
                self.ydata=cache[:,1] # 2nd col is intensity
                self.spectrum_curve_item=make.curve(self.xdata,self.ydata,color='k',linewidth=2)
                self.spectrum_range_item=make.range(self.xdata[-1]/2,self.xdata[-1]/2+20)
                self.spectrum_curve_plot.add_item(self.spectrum_curve_item)
                self.spectrum_curve_plot.add_item(self.spectrum_range_item)
                self.spectrum_curve_plot.do_autoscale()
                message=self.spectrum_file_path+" is loaded successfully"
                self.update_status_bar(message)
            except IndexError as err: # need full capture of all kinds of errors
                message="IndexError occurred when loading calibration spectrum: "+err.message
                self.update_status_bar(message)
                return
        if cache.ndim==1:                
            try: # throws exception if text file format is wrong
                self.ydata=cache[:,0]
                self.xdata=np.arange(cache.size)
                self.spectrum_curve_item=make.curve(self.xdata,self.ydata,color='k',linewidth=2)
                self.spectrum_range_item=make.range(self.xdata[-1]/2,self.xdata[-1]/2+20)
                self.spectrum_curve_plot.add_item(self.spectrum_curve_item)
                self.spectrum_curve_plot.add_item(self.spectrum_range_item)
                self.spectrum_curve_plot.do_autoscale()
                message=self.spectrum_file_path+" is loaded successfully"
                self.update_status_bar(message)
            except IndexError as err: # need full capture of all kinds of errors
                message="IndexError occurred when loading calibration spectrum: "+err.message
                self.update_status_bar(message)
                return
                
    def Gaussian(self,xdata,parameters):
        center,sigma,amp,offset=parameters
        ydata=offset+amp*np.exp((-0.5)*(xdata-center)**2/sigma**2)
        return ydata
    def errfunc(self,parameters,xdata,ydata):
        delta=self.Gaussian(xdata,parameters)-ydata
        total_error=np.sum(delta**2)
        return total_error       
    def func_fit_spectrum(self):
        if self.xdata==[]: # check if spectrum is loaded
            message="No calibration spectrum available"
            self.update_status_bar(message)
            return
        try:
            low_bound,high_bound=self.spectrum_range_item.get_range() # low and high bound might be reversed
        except AttributeError as err: # happens when clicked without any spectrum loaded
            message="AttributeError occurred when locating calibration peak: "+err.message
            self.update_status_bar(message)
            return
        # print low_bound,high_bound
        if low_bound>high_bound: # range item is 'polar'
            low_bound,high_bound=high_bound,low_bound
            
        index=(self.xdata>low_bound)*(self.xdata<high_bound)
        xdata_to_fit=self.xdata[index]
        ydata_to_fit=self.ydata[index]
        center=xdata_to_fit[ydata_to_fit.argmax()]
        sigma=1 # small width fit narrow peak better
        amp=ydata_to_fit.max()-ydata_to_fit.min()
        offset=ydata_to_fit.min()
        p0=[center,sigma,amp,offset]
        #y0=Gaussian(xdata_to_fit,p0)                        
        min_result=minimize(self.errfunc,p0[:],args=(xdata_to_fit,ydata_to_fit),method='Nelder-Mead',options={'xtol':1e-9,'maxfev':1e5})
#        print min_result
        self.fit_results.setText(unicode(min_result))
        p1=min_result.x
        y1=self.Gaussian(xdata_to_fit,p1)
        #fit0=make.curve(xdata_to_fit,y0,color='r',linewidth=2)
        fit1=make.curve(xdata_to_fit,y1,color='b',linewidth=2)
        #spectrum_curve_plot.add_item(fit0)
        self.spectrum_curve_plot.add_item(fit1)
        self.spectrum_curve_plot.replot()
        self.peak_location_label.setText("Current peak location: %3.2f"%p1[0]) # assumes 512 pixels
        
    def func_clear_spectrum(self):        
        if self.xdata==[]: # check if spectrum is loaded
            message="No calibration spectrum available"
            self.update_status_bar(message)
            return
        self.spectrum_curve_plot.del_all_items()
        self.spectrum_curve_item=make.curve(self.xdata,self.ydata,color='k',linewidth=2)
        self.spectrum_curve_plot.add_item(self.spectrum_curve_item)
        self.spectrum_range_item=make.range(self.xdata[-1]/2,self.xdata[-1]/2+20) # this item cannot be accessed by outside code, need to use self.blah
        self.spectrum_curve_plot.add_item(self.spectrum_range_item)
        self.spectrum_curve_plot.do_autoscale()

    # functions on pixel list for calibration        
    def func_add_cal_pixel(self):
        text=self.peak_location_label.text()
        regexp=r'.* (?P<num>\d*\.\d*)'
        regexp=re.compile(regexp)
        hit=regexp.match(text)
        if hit is not None:
            entry=hit.group("num")
            self.cal_pixel_list.addItem(entry)
        else:
            message="Fail to extract peak location"
            self.update_status_bar(message)
            return            
    def func_del_cal_pixel(self):
        self.cal_pixel_list.takeItem(self.cal_pixel_list.currentRow())        
    def func_clear_cal_pixel(self):
        self.cal_pixel_list.clear()
    
    def func_load_ref_wvlen(self):
        self.ref_wvlen_list.clear() # clear all previous content
        ref_file_name=QFileDialog.getOpenFileName() # error handling, path input?
        ref_file_name=unicode(ref_file_name) # if clicked "cancel" in dialog, this gives u''
        if ref_file_name==u'': # value None didn't catch the 'cancel' case
            message="No reference wavelength is loaded"
            self.update_status_bar(message)
            return
        try:
            f=open(ref_file_name)
            ref_entries=f.readlines() # get rid of return at the end
        except IOError as err:
            message="IOError occurred when loading reference wavelength: "+err.message # no message with this error?
            self.update_status_bar(message)
            return
        ref_entries=[item.rstrip() for item in ref_entries]
        self.ref_wvlen_list.addItems(ref_entries)        
    
    # functions on wavelength list for calibration
    def func_add_cal_wvlen(self):
        item=self.ref_wvlen_list.currentItem()
        if item is None:
            message='No reference wavelength is selected'
            self.update_status_bar(message)
            return
        entry=item.text()
#        print entry
        regexp=r'(?P<name>.*)\s(?P<wavelength>\d*\.\d*)'
        regexp=re.compile(regexp)
        try:
            hit=regexp.match(entry)
            entry=hit.group('wavelength') # assumes hit has something
        except TypeError as err:
            message="No wavelength entry is added:\t"+err.message
            self.update_status_bar(message)
            return
        except AttributeError as err: # occurs when hit is None
            message="No wavelength entry is added:\t"+err.message
            self.update_status_bar(message)
            return
        self.cal_wvlen_list.addItem(entry)
            
    def func_del_cal_wvlen(self):
        self.cal_wvlen_list.takeItem(self.cal_wvlen_list.currentRow())
    def func_clear_cal_wvlen(self):
        self.cal_wvlen_list.clear()
        
    def func_calibrate_button(self):
        if self.cal_pixel_list.count()==0 or self.cal_wvlen_list.count()==0:
            message="Cannot calibrate with an empty list"
            self.update_status_bar(message)
            return
        if self.cal_pixel_list.count()!=self.cal_wvlen_list.count():
            message="Pixel and wavelength lists differ in length"
            self.update_status_bar(message)
            return
            
        pixel_list=[]
        # below: from QListWidget to numpy array
        for index in range(self.cal_pixel_list.count()):
            pixel_list.append(self.cal_pixel_list.item(index))
        pixel_list=[unicode(item.text()) for item in pixel_list] # unicode() returns unicode string, item.text() returns Qstring object
        pixel_array=[float(item) for item in pixel_list] # float() takes single value, error occurs if letters are present
        pixel_array=np.array(pixel_array)       
#        print pixel_list
        
        wvlen_list=[]
        for index in range(self.cal_wvlen_list.count()):
            wvlen_list.append(self.cal_wvlen_list.item(index))
        wvlen_list=[unicode(item.text()) for item in wvlen_list]
        wvlen_array=[float(item) for item in wvlen_list]
        wvlen_array=np.array(wvlen_array)
        
        slope,intercept,r_value,p_value,std_err=linregress(pixel_array,wvlen_array)
        self.slope=slope
        self.intercept=intercept        
        
#        print 'calibration result:\n',slope,intercept,r_value,p_value,std_err
        cal_res="Calibration result:\n"+"slope: {0}\nintercept: {1}\nr value: {2}\np value: {3}\nstd err: {4}\n\n".format(slope,intercept,r_value,p_value,std_err)        
        cal_res=cal_res+"Pixel list:\n"+" ".join(pixel_list)+'\n'
        cal_res=cal_res+"Wavelength list:\n"+" ".join(wvlen_list)+'\n'        
        self.calibration_results.setText(cal_res)
        
    def func_save_cal_nanometer(self):
        # need to have a valid calibration spectrum, use this name with '.cal' extension
        try:
            calibration_file_path=self.spectrum_file_path[:-4]+'_nanometer.cal'
        except TypeError as err:# occurs when clicked without loading in spectrum
            message="TypeError:\t"+err.message
            self.update_status_bar(message)
            return        
        
        x_wavelength=self.slope*self.xdata+self.intercept                
        
        np.savetxt(calibration_file_path,x_wavelength,fmt="%.4e")            
        message="calibrated to nanometer:\t"+calibration_file_path
        self.update_status_bar(message)
            
    def func_save_cal_wavenumber(self):
        # need to have a valid calibration spectrum, use this name with '.cal' extension
        try:
            calibration_file_path=self.spectrum_file_path[:-4]+'_wavenumber.cal'
        except TypeError as err:# occurs when clicked without loading in spectrum
            message="TypeError:\t"+err.message
            self.update_status_bar(message)
            return        
        
        x_wavelength=self.slope*self.xdata+self.intercept                
        # catch ValueErrors here
        excitation_wvlen=float(self.laser_wavelength.text())
        x_wavenumber=10.0**7/excitation_wvlen-10.0**7/x_wavelength
        
        np.savetxt(calibration_file_path,x_wavenumber,fmt="%.4e")
        message="calibrated to wavenumber:\t"+calibration_file_path
        self.update_status_bar(message)
        
    def set_interface_geometry(self):
        self.main_window.setGeometry(QRect(50,50,1070,640))
        # cal spectrum and its associates
        self.spectrum_curve_plot.setGeometry(QRect(20,20,600,300))
        self.fit_button.setGeometry(QRect(20,340,100,50))
        self.clear_button.setGeometry(QRect(20,400,100,50))
        self.load_spectrum_button.setGeometry(QRect(20,460,100,50))
        self.peak_location_label.setGeometry(QRect(140,350,200,20))
        self.fit_results.setGeometry(QRect(140,380,200,130))        

        # cal pixel list and its associates
        self.cal_pixel_list.setGeometry(QRect(640,20,100,200))        
        self.add_cal_pixel_button.setGeometry(QRect(640,220,100,50))
        self.del_cal_pixel_button.setGeometry(QRect(640,270,100,50))
        self.clear_cal_pixel_button.setGeometry(QRect(640,320,100,50))        
        
        # cal wavelength list and associates
        self.cal_wvlen_list.setGeometry(QRect(760,20,100,200))
        self.add_cal_wvlen_button.setGeometry(QRect(760,220,100,50))    
        self.del_cal_wvlen_button.setGeometry(QRect(760,270,100,50))    
        self.clear_cal_wvlen_button.setGeometry(QRect(760,320,100,50))        
        
        # ref wavelength list and loading button
        self.ref_wvlen_list.setGeometry(QRect(880,20,170,450))
        self.load_ref_wvlen_button.setGeometry(QRect(880,470,170,50))        
        
        # calibrate button
        self.calibrate_button.setGeometry(QRect(640,370,220,50))
        self.save_wavenumber.setGeometry(QRect(640,420,220,50))
        self.save_nanometer.setGeometry(QRect(640,470,220,50))
        
        self.calibration_results.setGeometry(QRect(420,340,210,140))
        self.laser_wavelength.setGeometry(QRect(560,485,70,25))
        self.laser_label.setGeometry(QRect(420,485,140,25))

    def __init__(self):
        self.main_window=QMainWindow()
        
        self.spectrum_file_path=None
        self.xdata=[]
        self.ydata=[]
        self.slope=None
        self.intercept=None
        
        self.ui=QWidget()
        # spectrum display and peak fitting
        self.spectrum_curve_plot=CurvePlot()
        self.spectrum_curve_plot.setParent(self.ui)        
        self.spectrum_curve_item=make.curve(self.xdata,self.ydata,color='k',linewidth=2)
        self.spectrum_curve_plot.add_item(self.spectrum_curve_item)
        
        self.load_spectrum_button=QPushButton()
        self.load_spectrum_button.setText("Load spectrum")
        self.load_spectrum_button.setParent(self.ui)
        self.ui.connect(self.load_spectrum_button,SIGNAL("clicked()"),self.func_load_spectrum)
 
        self.fit_button=QPushButton()
        self.fit_button.setText("fit \npeak location")
        self.fit_button.setParent(self.ui)
        self.ui.connect(self.fit_button,SIGNAL("clicked()"),self.func_fit_spectrum)
        
        self.clear_button=QPushButton()
        self.clear_button.setText("reset spectrum")
        self.clear_button.setParent(self.ui)        
        self.ui.connect(self.clear_button,SIGNAL("clicked()"),self.func_clear_spectrum)
        
        self.peak_location_label=QLabel()
        self.peak_location_label.setText("Current peak location")
        self.peak_location_label.setParent(self.ui)

        self.fit_results=QTextBrowser()
        self.fit_results.setParent(self.ui)
        self.fit_results.setText("Peak fitting results")
        
        # list of located pixels for calibration, and associates
        self.cal_pixel_list=QListWidget()
        self.cal_pixel_list.setParent(self.ui)        
    
        self.add_cal_pixel_button=QPushButton()
        self.add_cal_pixel_button.setText("add to \ncal pixel")
        self.add_cal_pixel_button.setParent(self.ui)
        self.ui.connect(self.add_cal_pixel_button,SIGNAL("clicked()"),self.func_add_cal_pixel)
        
        self.del_cal_pixel_button=QPushButton()
        self.del_cal_pixel_button.setText("delete \ncal pixel")
        self.del_cal_pixel_button.setParent(self.ui)
        self.ui.connect(self.del_cal_pixel_button,SIGNAL("clicked()"),self.func_del_cal_pixel)
        
        self.clear_cal_pixel_button=QPushButton()
        self.clear_cal_pixel_button.setText("clear \ncal pixel")
        self.clear_cal_pixel_button.setParent(self.ui)
        self.ui.connect(self.clear_cal_pixel_button,SIGNAL("clicked()"),self.func_clear_cal_pixel)
        
        # list of reference wavelengths
        self.ref_wvlen_list=QListWidget()
        self.ref_wvlen_list.setParent(self.ui)
      
        self.load_ref_wvlen_button=QPushButton()
        self.load_ref_wvlen_button.setText("Load \nref wavelength")
        self.load_ref_wvlen_button.setParent(self.ui)
        self.ui.connect(self.load_ref_wvlen_button,SIGNAL("clicked()"),self.func_load_ref_wvlen)    
    
        # selected wavelength for calibration, and associates
        self.cal_wvlen_list=QListWidget()
        self.cal_wvlen_list.setParent(self.ui)   
     
        self.add_cal_wvlen_button=QPushButton()
        self.add_cal_wvlen_button.setText("add to \ncal wavelength")
        self.add_cal_wvlen_button.setParent(self.ui)
        self.ui.connect(self.add_cal_wvlen_button,SIGNAL("clicked()"),self.func_add_cal_wvlen)
      
        self.del_cal_wvlen_button=QPushButton()
        self.del_cal_wvlen_button.setText("delete \ncal wavelength")
        self.del_cal_wvlen_button.setParent(self.ui)
        self.ui.connect(self.del_cal_wvlen_button,SIGNAL("clicked()"),self.func_del_cal_wvlen)
        
        self.clear_cal_wvlen_button=QPushButton()
        self.clear_cal_wvlen_button.setText("clear \ncal wavelength")
        self.clear_cal_wvlen_button.setParent(self.ui)
        self.ui.connect(self.clear_cal_wvlen_button,SIGNAL("clicked()"),self.func_clear_cal_wvlen)

        # calibration
        self.calibrate_button=QPushButton()
        self.calibrate_button.setText("Calibrate")
        self.calibrate_button.setParent(self.ui)
        self.ui.connect(self.calibrate_button,SIGNAL("clicked()"),self.func_calibrate_button)
        
        self.calibration_results=QTextBrowser()
        self.calibration_results.setParent(self.ui)
        self.calibration_results.setText("Calibration results here")
        
        self.laser_label=QLabel()
        self.laser_label.setParent(self.ui)
        self.laser_label.setText("Laser wavelength (nm)")        
        
        self.laser_wavelength=QLineEdit()
        self.laser_wavelength.setParent(self.ui)
        self.laser_wavelength.setText("632.82")
        
        self.save_wavenumber=QPushButton()
        self.save_wavenumber.setText("Save calibration - wavenumber")
        self.save_wavenumber.setParent(self.ui)
        self.ui.connect(self.save_wavenumber,SIGNAL("clicked()"),self.func_save_cal_wavenumber)

        self.save_nanometer=QPushButton()
        self.save_nanometer.setText("Save calibration - nanometer")
        self.save_nanometer.setParent(self.ui)
        self.ui.connect(self.save_nanometer,SIGNAL("clicked()"),self.func_save_cal_nanometer)
        
        self.menubar=self.main_window.menuBar()
        self.fileMenu=self.menubar.addMenu('File')
        self.helpMenu=self.menubar.addMenu('Help')
        
        
        self.toolbar=self.main_window.addToolBar("curve toolbar")
        
        self.statusbar=self.main_window.statusBar()
        self.update_status_bar("Initiation has finished")
    
        self.manager=PlotManager(self.main_window)
        self.manager.add_plot(self.spectrum_curve_plot)
        self.manager.add_toolbar(self.toolbar)    
        self.manager.register_all_curve_tools() # make range slide, enable RMB menu

        self.main_window.setCentralWidget(self.ui) # single layer parent-children relationship
        self.set_interface_geometry()
        self.main_window.setWindowTitle("Raman Calibration Helper")
   

def run():
    app=QApplication([])
    stuff=interface()
    stuff.main_window.show()
    app.exec_()

if __name__=="__main__":
    run()