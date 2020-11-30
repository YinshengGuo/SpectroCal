# for version history and change log, see readme.md

# ----------------
from guidata.qt import QtCore, QtGui
from guidata.qt.QtGui import (QApplication, QWidget, QFileDialog, QPushButton,
                              QMainWindow, QListWidget, QLabel, QTextBrowser, 
                              QLineEdit, QSizePolicy, QFont)
from guidata.qt.QtCore import QRect  # ,Signal

# ---Import plot widget base class
from guiqwt.curve import CurvePlot
from guiqwt.plot import PlotManager
from guiqwt.builder import make
#from guidata.configtools import get_icon

# ---
import re
import numpy as np
from scipy.stats import linregress
from lmfit.models import Model, ConstantModel, GaussianModel, LorentzianModel


# to fix the issue of displaying on high resolution screen

# import PyQt5
# from PyQt5 import QtCore, QtGui, uic, QtWidgets
# from PyQt5.QtWidgets import QWidget

# if hasattr(QtCore.Qt, 'AA_EnableHighDpiScaling'):
#     QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, False)
# else:
#     print("No AA_EnableHighDpiScaling")

# if hasattr(QtCore.Qt, 'AA_UseHighDpiPixmaps'):
#     QApplication.setAttribute(QtCore.Qt.AA_UseHighDpiPixmaps, False)
# else:
#     print("No AA_UseHighDpiPixmaps")
    
# QtGui.QApplication.setAttribute(QtCore.Qt.AA_Use96Dpi)

QApplication.setAttribute(QtCore.Qt.AA_EnableHighDpiScaling, True)
# ----------------


def print_library_versions():
    import guidata
    import guiqwt
    import numpy
    import scipy
    import lmfit
    print("Libraries in use: ",
          f"guidata: {guidata.__version__}",
          f"guidata.qt: {guidata.qt.__version__}",
          f'guiqwt: {guiqwt.__version__}',
          f'numpy: {numpy.__version__}',
          f"scipy: {scipy.__version__}",
          f"lmfit: {lmfit.__version__}",
          sep='\n')
    return None


class interface():
    def update_status_bar(self, message):
        self.statusbar.showMessage(message)

    def func_load_spectrum(self):
        self.spectrum_curve_plot.del_all_items()
        self.spectrum_file_path, _ = QFileDialog.getOpenFileName()  # file opening might fail
        # self.spectrum_file_path=unicode(self.spectrum_file_path)
        self.spectrum_file_path = str(self.spectrum_file_path)
        print(f"spectrum_file_name {self.spectrum_file_path}")        
        if self.spectrum_file_path == u'':  # value None didn't catch the 'cancel' case
            message = "No reference spectrum is loaded"
            self.update_status_bar(message)
            return None
        
        try:
            # throws exception if clicked 'cancel' in file dialog
            cache = np.genfromtxt(self.spectrum_file_path)
        except IOError as err:
            # print(dir(err), err.args, err.errno, err.strerror, # err.characters_written,
            #       err.filename,err.filename2,err.winerror)
            message = "IOError occurred when reading calibration file: "+err.message
            self.update_status_bar(message)
            return None
        except ValueError as err:  # fixed: ref wavelength can go through
            message = "ValueError occurred when reading calibration file: "+err.message
            self.update_status_bar(message)
            return None
        except:
            message = "Loading spectrum failed"
            self.update_status_bar(message)
            return None
        if cache.ndim == 2:
            try:  # throws exception if text file format is wrong
                self.xdata = cache[:, 0]  # 1st col is pixel
                self.ydata = cache[:, 1]  # 2nd col is intensity
                self.spectrum_curve_item = make.curve(
                    self.xdata, self.ydata, color='k', linewidth=2)
                self.spectrum_range_item = make.range(
                    self.xdata[-1]/2, self.xdata[-1]/2+50)
                self.spectrum_curve_plot.add_item(self.spectrum_curve_item)
                self.spectrum_curve_plot.add_item(self.spectrum_range_item)
                self.spectrum_curve_plot.do_autoscale()
                message = self.spectrum_file_path+" is loaded successfully"
                self.update_status_bar(message)
            except IndexError as err:  # need full capture of all kinds of errors
                message = "IndexError occurred when loading 2-column calibration spectrum: "+err.message
                self.update_status_bar(message)
                return
        if cache.ndim == 1:
            try:  # throws exception if text file format is wrong
                self.ydata = cache
                self.xdata = np.arange(cache.size)
                self.spectrum_curve_item = make.curve(
                    self.xdata, self.ydata, color='k', linewidth=2)
                self.spectrum_range_item = make.range(
                    self.xdata[-1]/2, self.xdata[-1]/2+50)
                self.spectrum_curve_plot.add_item(self.spectrum_curve_item)
                self.spectrum_curve_plot.add_item(self.spectrum_range_item)
                self.spectrum_curve_plot.do_autoscale()
                message = self.spectrum_file_path+" is loaded successfully"
                self.update_status_bar(message)
            except IndexError as err:  # need full capture of all kinds of errors
                message = "IndexError occurred when loading 1-column calibration spectrum: "+err.message
                self.update_status_bar(message)
                return

    def Gaussian(self, xdata, parameters):
        center, sigma, amp, offset = parameters
        ydata = offset+amp*np.exp((-0.5)*(xdata-center)**2/sigma**2)
        return ydata
        
    def Lorentzian(self, xdata, parameters):
        center, sigma, amp, offset = parameters
        ydata = offset+amp*1/((xdata-center)**2+(sigma/2)**2)
        return ydata

    def errfunc(self, parameters, xdata, ydata):
        #delta = self.Gaussian(xdata, parameters)-ydata
        delta = self.Lorentzian(xdata, parameters)-ydata
        total_error = np.sum(delta**2)
        return total_error

    def func_fit_spectrum(self):
        if self.xdata is []:  # check if spectrum is loaded
            message = "No calibration spectrum available"
            self.update_status_bar(message)
            return
        try:
            # low and high bound might be reversed
            low_bound, high_bound = self.spectrum_range_item.get_range()
        except AttributeError as err:  # happens when clicked without any spectrum loaded
            # message = "AttributeError occurred when locating calibration peak: "+err.message
            message = "No calibration spectrum is loaded."
            self.update_status_bar(message)
            return
        # print low_bound,high_bound
        if low_bound > high_bound:  # range item is 'polar'
            low_bound, high_bound = high_bound, low_bound

        index = (self.xdata > low_bound)*(self.xdata < high_bound)
        xdata_to_fit = self.xdata[index]
        ydata_to_fit = self.ydata[index]

        pk_mdl = ConstantModel(prefix='const_')+GaussianModel(prefix='gauss_')
        p0 = pk_mdl.make_params()
        p0['const_c'].set(value=xdata_to_fit.min())
        p0['gauss_center'].set(value=(low_bound+high_bound)/2)
        p0['gauss_sigma'].set(value=1)
        p0['gauss_amplitude'].set(value=0)
        
        #pk_mdl = ConstantModel(prefix='const_')+LorentzianModel(prefix='lortzn_')
        #p0 = pk_mdl.make_params()
        #p0['const_c'].set(value=xdata_to_fit.min())
        #p0['lortzn_center'].set(value=(low_bound+high_bound)/2)
        #p0['lortzn_sigma'].set(value=1)
        #p0['lortzn_amplitude'].set(value=0)

        result = pk_mdl.fit(ydata_to_fit, x=xdata_to_fit, params=p0)
        y1 = result.best_fit
        p1 = result.params

        # self.fit_results.setText(unicode(result.fit_report(show_correl=False)))
        self.fit_results.setText(str(result.fit_report(show_correl=False)))

        # fit0=make.curve(xdata_to_fit,y0,color='r',linewidth=2)
        fit1 = make.curve(xdata_to_fit, y1, color='b', linewidth=2)
        # spectrum_curve_plot.add_item(fit0)
        self.spectrum_curve_plot.add_item(fit1)
        self.spectrum_curve_plot.replot()
        self.peak_location_label.setText(
            f"Current peak location: {p1['gauss_center'].value:3.3f}")  # assumes 512 pixels
            #"Current peak location: %3.2f" % p1['gauss_center'].value)  # assumes 512 pixels
            #"Current peak location: %3.2f" % p1['lortzn_center'].value)
        return

    def func_clear_spectrum(self):
        if self.xdata is []:  # check if spectrum is loaded
            message = "No calibration spectrum available"
            self.update_status_bar(message)
            return
        self.spectrum_curve_plot.del_all_items()
        self.spectrum_curve_item = make.curve(
            self.xdata, self.ydata, color='k', linewidth=2)
        self.spectrum_curve_plot.add_item(self.spectrum_curve_item)
        # this item cannot be accessed by outside code, need to use self.blah
        try:
            # range item is 100 px wide at the center of spectrum
            self.spectrum_range_item = make.range(self.xdata[-1]/2, self.xdata[-1]/2+100)
        except IndexError: # list index out of range, when no calibration spectrum is loaded
            message="No calibration spectrum is loaded."
            self.update_status_bar(message)
            return        
        self.spectrum_curve_plot.add_item(self.spectrum_range_item)
        self.spectrum_curve_plot.do_autoscale()
        
    def func_zoomout_spectrum(self):
        self.spectrum_curve_plot.do_autoscale()

    # functions on pixel list for calibration
    def func_add_cal_pixel(self):
        text = self.peak_location_label.text()
        regexp = r'.* (?P<num>\d*\.\d*)'
        regexp = re.compile(regexp)
        hit = regexp.match(text)
        if hit is not None:
            entry = hit.group("num")
            self.cal_pixel_list.addItem(entry)
        else:
            message = "Fail to extract peak location"
            self.update_status_bar(message)
            return

    def func_del_cal_pixel(self):
        self.cal_pixel_list.takeItem(self.cal_pixel_list.currentRow())

    def func_clear_cal_pixel(self):
        self.cal_pixel_list.clear()

    def func_load_ref_wvlen(self):
        self.ref_wvlen_list.clear()  # clear all previous content
        self.ref_file_name,_ = QFileDialog.getOpenFileName()  # error handling, path input?
        # ref_file_name=unicode(ref_file_name) # if clicked "cancel" in dialog, this gives u''
        # if clicked "cancel" in dialog, this gives u''
        self.ref_file_name = str(self.ref_file_name)
        print("ref_file_name", self.ref_file_name)
        # print("_", _)
        if self.ref_file_name == u'':  # value None didn't catch the 'cancel' case
            message = "No reference wavelength is loaded"
            self.update_status_bar(message)
            return        
        try:
            f = open(self.ref_file_name)
            ref_entries = f.readlines()  # get rid of return at the end
        except IOError as err:
            # print(dir(err))
            # print(err.args, err.errno, err.strerror,err. # err.characters_written,
            #       err.filename,err.filename2,err.winerror)
            message = "IOError occurred when loading reference wavelength: " # + err.msg  # no message with this error?
            self.update_status_bar(message)
            return
        ref_entries = [item.rstrip() for item in ref_entries]
        self.ref_wvlen_list.addItems(ref_entries)

    # functions on wavelength list for calibration
    def func_add_cal_wvlen(self):
        item = self.ref_wvlen_list.currentItem()
        if item is None:
            message = 'No reference wavelength is selected'
            self.update_status_bar(message)
            return
        entry = item.text()
#        print entry
        regexp = r'(?P<name>.*)\s(?P<wavelength>-{0,1}\d*\.\d*)'
        regexp = re.compile(regexp)
        try:
            hit = regexp.match(entry)
            entry = hit.group('wavelength')  # assumes hit has something
        except TypeError as err:
            message = "No wavelength entry is added:\t"+err.message
            self.update_status_bar(message)
            return
        except AttributeError as err:  # occurs when hit is None
            message = "No wavelength entry is added:\t"+err.message
            self.update_status_bar(message)
        except:
            message = "func_add_cal_wvlen(): unexpected error occurred." # err.message
            self.update_status_bar(message)
            return
        self.cal_wvlen_list.addItem(entry)

    def func_del_cal_wvlen(self):
        self.cal_wvlen_list.takeItem(self.cal_wvlen_list.currentRow())

    def func_clear_cal_wvlen(self):
        self.cal_wvlen_list.clear()

    def func_calibrate_button(self):
        if self.cal_pixel_list.count() <= 1 or self.cal_wvlen_list.count() <= 1:
            message = "Need at least two points to calibrate."
            self.update_status_bar(message)
            return None
        if self.cal_pixel_list.count() != self.cal_wvlen_list.count():
            message = "Pixel and wavelength lists differ in length"
            self.update_status_bar(message)
            return None

        pixel_list = []
        # below: from QListWidget to numpy array
        for index in range(self.cal_pixel_list.count()):
            pixel_list.append(self.cal_pixel_list.item(index))
        # pixel_list=[unicode(item.text()) for item in pixel_list] # unicode() returns unicode string, item.text() returns Qstring object
        # unicode() returns unicode string, item.text() returns Qstring object
        pixel_list = [str(item.text()) for item in pixel_list]
        # float() takes single value, error occurs if letters are present
        pixel_array = [float(item) for item in pixel_list]
        pixel_array = np.array(pixel_array)
#        print pixel_list

        wvlen_list = []
        for index in range(self.cal_wvlen_list.count()):
            wvlen_list.append(self.cal_wvlen_list.item(index))
        #wvlen_list=[unicode(item.text()) for item in wvlen_list]
        wvlen_list = [str(item.text()) for item in wvlen_list]
        wvlen_array = [float(item) for item in wvlen_list]
        wvlen_array = np.array(wvlen_array)

        slope, intercept, r_value, p_value, std_err = linregress(pixel_array, wvlen_array)
        self.slope = slope
        self.intercept = intercept

#        print 'calibration result:\n',slope,intercept,r_value,p_value,std_err
        cal_res = "Calibration result:\n" + \
            "slope: {0}\nintercept: {1}\nr value: {2}\np value: {3}\nstd err: {4}\n\n".format(
                slope, intercept, r_value, p_value, std_err)
        cal_res = cal_res+"Pixel list:\n"+" ".join(pixel_list)+'\n'
        cal_res = cal_res+"Wavelength list:\n"+" ".join(wvlen_list)+'\n'
        self.calibration_results.setText(cal_res)
        return None

    def func_calibrate_shift_button(self):
        if self.cal_pixel_list.count() <= 1 or self.cal_wvlen_list.count() <= 1:
            message = "Cannot calibrate with less than two points"
            self.update_status_bar(message)
            return
        if self.cal_pixel_list.count() != self.cal_wvlen_list.count():
            message = "Pixel and wavelength lists differ in length"
            self.update_status_bar(message)
            return

        pixel_list = []
        # below: from QListWidget to numpy array
        for index in range(self.cal_pixel_list.count()):
            pixel_list.append(self.cal_pixel_list.item(index))
        # unicode() returns unicode string, item.text() returns Qstring object
        pixel_list = [str(item.text()) for item in pixel_list]
        # float() takes single value, error occurs if letters are present
        pixel_array = [float(item) for item in pixel_list]
        pixel_array = np.array(pixel_array)
        # print pixel_list
        shift_list = []
        for index in range(self.cal_wvlen_list.count()):
            shift_list.append(self.cal_wvlen_list.item(index))
        shift_list = [str(item.text()) for item in shift_list]
        shift_array = [float(item) for item in shift_list]
        shift_array = np.array(shift_array)

        def cal_shift(x, k, b, l):
            '''convert x(px) to y(cm-1) : px -> nm -> cm-1'''
            y = 1e-2*1e9*(1/l-1/(k*x+b))
            return y
        mdl = Model(cal_shift, prefix='cal_Shift_')
        # k,b,l=-0.095,944,532 # initial guess, k&b roughly calculated from two points with assumed laser wavelength l
        p0 = mdl.make_params()
        p0['cal_Shift_k'].set(value=1)
        p0['cal_Shift_b'].set(value=1)
        # True,max=l+1e-3,min=l-1e-3)
        p0['cal_Shift_l'].set(value=float(self.laser_wavelength.text()), vary=False)
        result = mdl.fit(shift_array, x=pixel_array, params=p0)
        p1 = result.params
        self.slope = p1['cal_Shift_k']
        self.intercept = p1['cal_Shift_b']
        # print(result.fit_report())

        cal_res = result.fit_report()
        cal_res = cal_res+"\nPixel list:\n"+" ".join(pixel_list)+'\n'
        cal_res = cal_res+"Raman shift list:\n"+" ".join(shift_list)+'\n'
        self.calibration_results.setText(cal_res)
        return None

    def func_save_fitting_results(self):
        # need to have a valid calibration spectrum, use this name with '.cal' extension
        try:
            calibration_file_path = self.spectrum_file_path[:-4]+'_fit.cal'
        except TypeError as err: # occurs when spectrum_file_path is None, when clicked without loading in spectrum
            # print(dir(err), err.__doc__, err.args, sep='\n')
            message = "No calibration spectrum is loaded." # f"TypeError:\t{err.args[0]}." # err.message no longer available
            self.update_status_bar(message)
            return
        
        if self.slope is None or self.intercept is None:
            message="No calibration is completed or available to be saved."
            self.update_status_bar(message)
            return            
        
        with open(calibration_file_path, 'w') as f:
            f.write(self.calibration_results.toPlainText())
        message = "linear regression info:\t"+calibration_file_path
        self.update_status_bar(message)
        return None

    def func_save_cal_nanometer(self):
        try: # need to have a valid calibration spectrum, use this name with '.cal' extension
            calibration_file_path = self.spectrum_file_path[:-4]+'_nm.cal'
        except TypeError as err: # occurs when spectrum_file_path is None, when clicked without loading in spectrum
            message = "No calibration spectrum is loaded." # f"TypeError:\t{err.args[0]}." # err.message
            self.update_status_bar(message)
            return
        
        try: # exception occurs when slope and intercept are None, when no linear regression is done
            x_wavelength = self.slope*self.xdata+self.intercept
        except TypeError as err:
            message="No calibration is completed or available to be saved."
            self.update_status_bar(message)
            return
        
        np.savetxt(calibration_file_path, x_wavelength, fmt="%.4e")
        message = "calibrated to nanometer:\t"+calibration_file_path
        self.update_status_bar(message)
        return None

    def func_save_cal_wavenumber(self):
        # need to have a valid calibration spectrum, use this name with '.cal' extension
        try:
            calibration_file_path=self.spectrum_file_path[:-4]+'_wavenumber.cal'
        except TypeError as err:  # occurs when spectrum_file_path is None, when clicked without loading in spectrum
            message = "No calibration spectrum is loaded." # f"TypeError:\t{err.args[0]}." # err.message
            self.update_status_bar(message)
            return
        
        try: # exception occurs when slope and intercept are None, when no linear regression is done
            x_wavelength = self.slope*self.xdata+self.intercept
        except TypeError as err:
            message="No calibration is completed or available to be saved."
            self.update_status_bar(message)
            return

        # catch ValueErrors here
        excitation_wvlen = float(self.laser_wavelength.text())
        x_wavenumber = 1e7/excitation_wvlen-1e7/x_wavelength # 10.0**7/excitation_wvlen-10.0**7/x_wavelength

        np.savetxt(calibration_file_path, x_wavenumber, fmt="%.4e")
        message = "calibrated to wavenumber:\t"+calibration_file_path
        self.update_status_bar(message)
        return None

    def func_save_cal_eV(self):
        # need to have a valid calibration spectrum, use this name with '.cal' extension
        try:
            calibration_file_path = self.spectrum_file_path[:-4]+'_eV.cal'
        except TypeError as err: # occurs when spectrum_file_path is None, when clicked without loading in spectrum
            message = "No calibration spectrum is loaded." # f"TypeError:\t{err.args[0]}." # err.message
            self.update_status_bar(message)
            return
        
        try: # exception occurs when slope and intercept are None, when no linear regression is done
            x_wavelength = self.slope*self.xdata+self.intercept
        except TypeError as err: 
            message="No calibration is completed or available to be saved."
            self.update_status_bar(message)
            return
        
        h = 6.626e-34  # unit: Js
        c = 299792458e9  # unit: nm/s
        J2eV = 1.60218e-19  # unit: J/eV
        x_eV = h*c/x_wavelength/J2eV
        np.savetxt(calibration_file_path, x_eV, fmt="%.4e")
        message = "calibrated to eV:\t"+calibration_file_path
        self.update_status_bar(message)
        return None

    def setup_interface_geometry(self):
        '''
        Full range: 1000x500
        Window size: 1020x640
        QRect(upperleft_x,upperleft_y,span_x,span_y)
        '''
        self.main_window.setGeometry(QRect(
            100, 100, 1040, 605))  # upperleft x and y are arbitrary starting coord on screen
        self.main_window.setSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred)
        
        # cal spectrum and its associates
        self.spectrum_curve_plot.setGeometry(QRect(20, 20, 600, 300))
        self.fit_button.setGeometry(QRect(20, 340, 100, 50))
        self.fit_button.setSizePolicy(QSizePolicy.Preferred,QSizePolicy.Preferred)
        # self.setSizePolicy(QtWidgets.QSizePolicy(QtWidgets.QSizePolicy.MinimumExpanding,
        #                                         QtWidgets.QSizePolicy.MinimumExpanding))
        self.clear_button.setGeometry(QRect(20, 395, 100, 50))
        self.load_spectrum_button.setGeometry(QRect(20, 450, 100, 50))
        self.peak_location_label.setGeometry(QRect(140, 340, 200, 20))
        self.fit_results.setGeometry(QRect(140, 360, 240, 140))

        # cal pixel list and its associates
        self.cal_pixel_list.setGeometry(QRect(640, 20, 100, 180))
        self.add_cal_pixel_button.setGeometry(QRect(640, 200, 100, 45))
        self.del_cal_pixel_button.setGeometry(QRect(640, 245, 100, 45))
        self.clear_cal_pixel_button.setGeometry(QRect(640, 290, 100, 45))
        # cal wavelength list and associates
        self.cal_wvlen_list.setGeometry(QRect(760, 20, 100, 180))
        self.add_cal_wvlen_button.setGeometry(QRect(760, 200, 100, 45))
        self.del_cal_wvlen_button.setGeometry(QRect(760, 245, 100, 45))
        self.clear_cal_wvlen_button.setGeometry(QRect(760, 290, 100, 45))

        # ref wavelength list and loading button
        self.ref_wvlen_list.setGeometry(QRect(880, 20, 140, 370))
        self.load_ref_wvlen_button.setGeometry(QRect(880, 390, 140, 60))
        self.laser_label.setGeometry(QRect(880, 450, 140, 25))
        self.laser_wavelength.setGeometry(QRect(880, 475, 70, 25))

        # calibrate button
        self.calibrate_button.setGeometry(QRect(640, 345, 110, 50))
        self.calibrate_shift_button.setGeometry(QRect(750, 345, 110, 50))
        self.save_fitting_results.setGeometry(QRect(640, 400, 110, 50))
        self.save_wavenumber.setGeometry(QRect(750, 400, 110, 50))
        self.save_nanometer.setGeometry(QRect(640, 450, 110, 50))
        self.save_eV.setGeometry(QRect(750, 450, 110, 50))
        self.calibration_results.setGeometry(QRect(420, 345, 210, 160))
        return None
    
    def setup_ui(self):
        
        self.ui = QWidget()
        
        # explicitly set font size, to fix display issue on high res screen
        txtQFont=QFont('Arial',10)
        
        # spectrum display and peak fitting
        self.spectrum_curve_plot = CurvePlot()
        
        # print(dir(self.spectrum_curve_plot),
        #       self.spectrum_curve_plot.setAxisFont.__doc__,
        #       self.spectrum_curve_plot.AXIS_IDS,
        #       self.spectrum_curve_plot.AXIS_NAMES,
        #       sep='\n') # display documentation 
        self.spectrum_curve_plot.setAxisFont(0, QFont('Arial', 10))
        self.spectrum_curve_plot.setAxisFont(2, QFont('Arial', 10))
        
        self.spectrum_curve_plot.setParent(self.ui)
        self.spectrum_curve_item = make.curve(
            self.xdata, self.ydata, color='k', linewidth=2)
        self.spectrum_curve_plot.add_item(self.spectrum_curve_item)

        self.load_spectrum_button = QPushButton()
        self.load_spectrum_button.setText("Load spectrum")
        self.load_spectrum_button.setFont(txtQFont)
        self.load_spectrum_button.setParent(self.ui)
        # self.ui.connect(self.load_spectrum_button,Signal("clicked()"),self.func_load_spectrum)
        self.load_spectrum_button.clicked.connect(self.func_load_spectrum)

        self.fit_button = QPushButton()
        self.fit_button.setText("fit \npeak location")
        self.fit_button.setFont(txtQFont)
        self.fit_button.setParent(self.ui)
        # self.ui.connect(self.fit_button,Signal("clicked()"),self.func_fit_spectrum)
        self.fit_button.clicked.connect(self.func_fit_spectrum)

        self.clear_button = QPushButton()
        self.clear_button.setText("reset spectrum")
        self.clear_button.setFont(txtQFont)
        self.clear_button.setParent(self.ui)
        # self.ui.connect(self.clear_button,Signal("clicked()"),self.func_clear_spectrum)
        self.clear_button.clicked.connect(self.func_clear_spectrum)
        #self.clear_button.clicked.connect(self.func_zoomout_spectrum)

        self.peak_location_label = QLabel()
        self.peak_location_label.setText("Current peak location")
        self.peak_location_label.setFont(txtQFont)
        self.peak_location_label.setParent(self.ui)

        self.fit_results = QTextBrowser()
        self.fit_results.setFont(txtQFont)
        self.fit_results.setParent(self.ui)
        self.fit_results.setText("Peak fitting results")

        # list of located pixels for calibration, and associates
        self.cal_pixel_list = QListWidget()
        self.cal_pixel_list.setParent(self.ui)

        self.add_cal_pixel_button = QPushButton()
        self.add_cal_pixel_button.setText("add to \ncal pixel")
        self.add_cal_pixel_button.setFont(txtQFont)
        self.add_cal_pixel_button.setParent(self.ui)
        # self.ui.connect(self.add_cal_pixel_button,Signal("clicked()"),self.func_add_cal_pixel)
        self.add_cal_pixel_button.clicked.connect(self.func_add_cal_pixel)

        self.del_cal_pixel_button = QPushButton()
        self.del_cal_pixel_button.setText("delete \ncal pixel")
        self.del_cal_pixel_button.setFont(txtQFont)
        self.del_cal_pixel_button.setParent(self.ui)
        # self.ui.connect(self.del_cal_pixel_button,Signal("clicked()"),self.func_del_cal_pixel)
        self.del_cal_pixel_button.clicked.connect(self.func_del_cal_pixel)

        self.clear_cal_pixel_button = QPushButton()
        self.clear_cal_pixel_button.setText("clear \ncal pixel")
        self.clear_cal_pixel_button.setFont(txtQFont)
        self.clear_cal_pixel_button.setParent(self.ui)
        # self.ui.connect(self.clear_cal_pixel_button,Signal("clicked()"),self.func_clear_cal_pixel)
        self.clear_cal_pixel_button.clicked.connect(self.func_clear_cal_pixel)

        # list of reference wavelengths
        self.ref_wvlen_list = QListWidget()
        self.ref_wvlen_list.setParent(self.ui)

        self.load_ref_wvlen_button = QPushButton()
        self.load_ref_wvlen_button.setText("Load reference\n Wavelength/\nRaman shift")
        self.load_ref_wvlen_button.setFont(txtQFont)
        self.load_ref_wvlen_button.setParent(self.ui)
        # self.ui.connect(self.load_ref_wvlen_button,Signal("clicked()"),self.func_load_ref_wvlen)
        self.load_ref_wvlen_button.clicked.connect(self.func_load_ref_wvlen)

        # selected wavelength for calibration, and associates
        self.cal_wvlen_list = QListWidget()
        self.cal_wvlen_list.setParent(self.ui)

        self.add_cal_wvlen_button = QPushButton()
        self.add_cal_wvlen_button.setText("add to \ncal reference")
        self.add_cal_wvlen_button.setFont(txtQFont)
        self.add_cal_wvlen_button.setParent(self.ui)
        # self.ui.connect(self.add_cal_wvlen_button,Signal("clicked()"),self.func_add_cal_wvlen)
        self.add_cal_wvlen_button.clicked.connect(self.func_add_cal_wvlen)

        self.del_cal_wvlen_button = QPushButton()
        self.del_cal_wvlen_button.setText("delete \ncal reference")
        self.del_cal_wvlen_button.setFont(txtQFont)
        self.del_cal_wvlen_button.setParent(self.ui)
        # self.ui.connect(self.del_cal_wvlen_button,Signal("clicked()"),self.func_del_cal_wvlen)
        self.del_cal_wvlen_button.clicked.connect(self.func_del_cal_wvlen)

        self.clear_cal_wvlen_button = QPushButton()
        self.clear_cal_wvlen_button.setText("clear \ncal reference")
        self.clear_cal_wvlen_button.setFont(txtQFont)
        self.clear_cal_wvlen_button.setParent(self.ui)
        # self.ui.connect(self.clear_cal_wvlen_button,Signal("clicked()"),self.func_clear_cal_wvlen)
        self.clear_cal_wvlen_button.clicked.connect(self.func_clear_cal_wvlen)

        # calibration
        self.calibrate_button = QPushButton()
        self.calibrate_button.setText("Calibrate with\nWavelength")
        self.calibrate_button.setFont(txtQFont)
        self.calibrate_button.setParent(self.ui)
        # self.ui.connect(self.calibrate_button,Signal("clicked()"),self.func_calibrate_button)
        self.calibrate_button.clicked.connect(self.func_calibrate_button)

        self.calibrate_shift_button = QPushButton()
        self.calibrate_shift_button.setText("Calibrate with\nRaman shift")
        self.calibrate_shift_button.setFont(txtQFont)
        self.calibrate_shift_button.setParent(self.ui)
        self.calibrate_shift_button.clicked.connect(
            self.func_calibrate_shift_button)

        self.calibration_results = QTextBrowser()
        self.calibration_results.setParent(self.ui)
        self.calibration_results.setText("Calibration results here")
        self.calibration_results.setFont(txtQFont)

        self.laser_label = QLabel()
        self.laser_label.setParent(self.ui)
        self.laser_label.setText("Raman pump (nm)")
        self.laser_label.setFont(txtQFont)

        self.laser_wavelength = QLineEdit()
        self.laser_wavelength.setParent(self.ui)
        self.laser_wavelength.setText("532.00") # "632.82"
        self.laser_wavelength.setFont(txtQFont)

        self.save_fitting_results = QPushButton()
        self.save_fitting_results.setText("Save cal - fit")
        self.save_fitting_results.setFont(txtQFont)
        self.save_fitting_results.setParent(self.ui)
        # self.ui.connect(self.save_fitting_results,Signal("clicked()"),self.func_save_fitting_results)
        self.save_fitting_results.clicked.connect(
            self.func_save_fitting_results)

        self.save_wavenumber = QPushButton()
        self.save_wavenumber.setText("Save cal - cm-1")
        self.save_wavenumber.setFont(txtQFont)
        self.save_wavenumber.setParent(self.ui)
        # self.ui.connect(self.save_wavenumber,Signal("clicked()"),self.func_save_cal_wavenumber)
        self.save_wavenumber.clicked.connect(self.func_save_cal_wavenumber)

        self.save_nanometer = QPushButton()
        self.save_nanometer.setText("Save cal - nm")
        self.save_nanometer.setFont(txtQFont)
        self.save_nanometer.setParent(self.ui)
        # self.ui.connect(self.save_nanometer,Signal("clicked()"),self.func_save_cal_nanometer)
        self.save_nanometer.clicked.connect(self.func_save_cal_nanometer)

        self.save_eV = QPushButton()
        self.save_eV.setText("Save cal - eV")
        self.save_eV.setFont(txtQFont)
        self.save_eV.setParent(self.ui)
        # self.ui.connect(self.save_eV,Signal("clicked()"),self.func_save_cal_eV)
        self.save_eV.clicked.connect(self.func_save_cal_eV)
        return None

    def __init__(self):
        self.main_window = QMainWindow()

        self.spectrum_file_path = None
        self.xdata = []
        self.ydata = []
        self.slope = None
        self.intercept = None

        self.main_window.setWindowTitle("Raman Calibration Helper")
        
        self.menubar = self.main_window.menuBar()
        self.fileMenu = self.menubar.addMenu('File')
        self.helpMenu = self.menubar.addMenu('Help')
        
        # single layer parent-children relationship
        # self.ui = QWidget()
        self.setup_ui()
        self.main_window.setCentralWidget(self.ui)
        self.setup_interface_geometry()
        
        self.toolbar = self.main_window.addToolBar("curve toolbar")

        self.statusbar = self.main_window.statusBar()
        self.update_status_bar("Initiation has finished")

        self.manager = PlotManager(self.main_window)
        self.manager.add_plot(self.spectrum_curve_plot)
        self.manager.add_toolbar(self.toolbar)
        self.manager.register_all_curve_tools()  # make range slide, enable RMB menu
        
        return None


def SpectroCal():
    import sys
    app = QApplication([])
    stuff = interface()
    stuff.main_window.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    print_library_versions()
    SpectroCal()
