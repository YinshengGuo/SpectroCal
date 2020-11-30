# Spectral Calibration 
GUI built with QT via guidata and guiqwt

## Version History
----
### v0.82 20201129
0. Update from WinPython 3.4.4.1 to WinPython 3.7.7.0  
    
    Main difference and main changes: 

1. fixed issue porting code from PyQt4 to PyQt5

    In RamanCal 0.8, guidata 1.7.6 uses guidata.qt 5.14.1

    In RamanCal 0.7, guidata 1.7.5 uses guidata.qt 4.11.4

    Some calling signatures have been changed between Qt4 and Qt5, for example:  
    
    ref_file_name = QFileDialog.getOpenFileName() # v0.7 is changed to 
    ref_file_name, _ = QFileDialog.getOpenFileName() # v0.8

    QFileDialog.getOpenFileName() return format is changed to return a tuple in v0.8. 
    file name is at index 0, option is at index 1.
    
    If keep using the v0.7 call signature, 
    ref_file_name becomes a tuple which does not give a valid file path.

2. fixed PyQt5 issue of displaying on high resolution screens
    
    explicitly set QApplication attribute after library import
    
    explicitly set QFont for all text
    
    to get around display issues on high resolution screens

3. updated python syntax using f-string
    
    In RamanCal 0.7, f-string is not available in python 3.4.4.1

4. updated python syntax of error attributes
    
    In RamanCal 0.8, error no longer has error.message attribute, error.args[0] should be used instead

5. added several exception handling cases

6. added codes to fit lorentzian lineshape using lmfit

### v0.81 20200507
0. migrate from PyQt4 to PyQt5, fix library issues

### v0.8 20170308
0. parse negative/antiStokes Raman shift

### v0.7 (20160923)
0. library versions  
    python 3.4.4, guidata 1.7.5, guiqwt 3.0.2
    numpy 1.10.4, scipy 0.17.0, lmfit 0.9.2
1. migrated from py27 to py34
    - connect signal to slot

        old: parentQWidget.connect(button,SIGNAL('clicked()'),func_button)

        new: button.clicked.connect(func_button)

        ref: https://wiki.qt.io/Signals_and_Slots_in_PySide
    - unicode() -> str()

        unicode() of Python2 is equivalent to str() in Python3
        
        ref: http://stackoverflow.com/questions/6812031/how-to-make-unicode-string-with-python3

2. added button and function to calibrate with Raman shift standards

3. adjusted widgets size and layout

### v0.6 (20160808)
0. library versions  
    python 2.7.9, guidata 1.6.2, guiqwt 2.3.2
    numpy 1.9.1, scipy 0.14.0, lmfit 0.9.2

1. use lmfit package for fitting

    more robust fitting in general

    able to fit negative peaks (FRIKES)

2. adjust peak fitting output text width

3. notification bar signal from 'save-cal fit' button

## To-do list
- output peak fitting list
- compile into executable
- fix this: is [], changed elementwise == comparison to identity is