from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QStackedWidget, QWidget, QPushButton, QHBoxLayout
)

import sys
from Software.Processor import *
from Software.Predictor.PtyChoPy import * 
from PyQt5.QtCore import pyqtSignal
import numpy as np
import scipy
from Software.Predictor.PtyChoPy.ErrorPopup import *

class PtyChoPy(QMainWindow):
    iterations_complete = pyqtSignal(object, str)

    def __init__(self, cmap = "binary"):
        super().__init__()
        self.cmap = cmap
        self.setWindowTitle("PtyChoPy")
        self.setGeometry(100, 100, 1400, 800)
        self.processed_trace = None
        self.init_ui()

    def init_vars(self, data):
        self.raw_trace, self.raw_delay, self.raw_wvl = map(data.get, ['trace', 'delay', 'wvl'])
        self.start_processing()
        self.show()
        
    def start_processing(self):
        trace0 = self.step0_processing(self.raw_trace)

        # proceed to step 1: crop 
        self.raw_data = [trace0, self.raw_delay, self.raw_wvl]
        self.handle_stacked_widgets()

    def init_ui(self):
        # Main layout
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        self.layout = QVBoxLayout(main_widget)

        # Navigation buttons
        nav_layout = QHBoxLayout()
        return_button = QPushButton("<< Return")
        proceed_button = QPushButton("Proceed >>")
        nav_layout.addWidget(return_button)
        nav_layout.addWidget(proceed_button)
        self.layout.addLayout(nav_layout)
        return_button.clicked.connect(self.navigate_back)
        proceed_button.clicked.connect(self.navigate_forward)

        # Stacked widget to manage multiple screens
        self.screens_widget = QStackedWidget()
        self.layout.addWidget(self.screens_widget)
        self.init_stacked_widget()

    def init_stacked_widget(self):
        self.cropper = Crop(self, self.cmap)
        self.FFT_filter = FFT_Filter(self, cmap = self.cmap)
        self.t_smoother = Smooth_Time(self, cmap = self.cmap)
        self.f_smoother = Smooth_Freq(self, cmap = self.cmap)
        self.init_guess = InitGuess(self, cmap = self.cmap)
        self.iterator = Recovery(self, cmap = self.cmap)
        self.iterator.divergence_occuring.connect(self.handle_div_popup)
        self.init_guess.divergence_occuring.connect(self.handle_div_popup)

        self.screens_widget.addWidget(self.cropper)
        self.screens_widget.addWidget(self.FFT_filter)
        self.screens_widget.addWidget(self.t_smoother)
        self.screens_widget.addWidget(self.f_smoother)
        self.screens_widget.addWidget(self.init_guess)
        self.screens_widget.addWidget(self.iterator)

    def handle_stacked_widgets(self):
        if self.screens_widget.currentIndex() == 0:
            # always use raw trace for cropping
            self.cropper.init_vars(self.raw_data)
            self.setWindowTitle("PtyChoPy: Cropping")

        elif self.screens_widget.currentIndex() == 1:
            self.cropped_data = [self.cropper.cropped_trace, self.cropper.cropped_delay, self.cropper.cropped_wvl]
            self.FFT_filter.init_vars(*self.cropped_data)
            self.setWindowTitle("PtyChoPy: FFT Filtering")

        elif self.screens_widget.currentIndex() == 2:
            self.filtered_data = [self.FFT_filter.filtered_trace, self.cropped_data[1], self.cropped_data[2]]
            self.t_smoother.init_vars(self.filtered_data) 
            self.setWindowTitle("PtyChoPy: Time Smoothing")

        elif self.screens_widget.currentIndex() == 3:
            self.tsmooth_data = [self.t_smoother.smooth_trace, self.filtered_data[1], self.filtered_data[2]]

            # do interpolation before we proceed to smooth f
            interp_trace, interp_delay, interp_wvl, interp_f = self.step4_interpolate(self.tsmooth_data)
            self.interped_data = [interp_trace, interp_delay, interp_wvl, interp_f]
            self.f_smoother.init_vars(self.interped_data)            
            self.setWindowTitle("PtyChoPy: Frequency Smoothing")

        # recovery screens
        elif self.screens_widget.currentIndex() == 4:
            self.init_data = [self.f_smoother.smooth_trace, self.interped_data[1], self.interped_data[2], self.interped_data[3]]
            self.setWindowTitle("PtyChoPy: Initializing Guess (loading...)")            # to display while waiting for ptycho constraint
            self.init_guess.guess_ready.connect(self.handle_guess)
            self.init_guess.init_vars(*self.init_data)
            self.setWindowTitle("PtyChoPy: Initial Guess")

        elif self.screens_widget.currentIndex() == 5:
            self.it_max = self.init_guess.it_max
            self.iterator.recovery_done.connect(self.handle_recovery_completed)
            self.iterator.init_vars(self.pro_trace_int, self.guess_matrix, self.pro_delay, self.pro_time, self.pro_f, self.it_max)
            self.setWindowTitle("PtyChoPy: Recovery")

    def navigate_forward(self):
        """Navigate to the next screen."""
        current_index = self.screens_widget.currentIndex()
        if current_index < self.screens_widget.count() - 1:
            self.screens_widget.setCurrentIndex(current_index + 1)
            self.handle_stacked_widgets()

    def navigate_back(self):
        """Navigate to the previous screen."""
        current_index = self.screens_widget.currentIndex()
        if current_index > 0:
            self.screens_widget.setCurrentIndex(current_index - 1)
            self.handle_stacked_widgets()

    @staticmethod
    def step0_processing(data_matrix):
        # perform outlier removal
        threshold = 4 * np.mean(np.min(data_matrix))
        ind_outliers = np.where(data_matrix < threshold)[0]
        for ind in ind_outliers:
            if 0 < ind < len(data_matrix) - 1:
                data_matrix[ind] = (data_matrix[ind - 1] + data_matrix[ind + 1]) / 2

        # Translate trace by minimum and remove negative points
        lambda_profile_pre = np.mean(data_matrix, axis = 1) / np.max(np.mean(data_matrix, axis = 1))
        data_matrix = np.clip(data_matrix - np.min(lambda_profile_pre), 0, None)
        data_matrix /= np.max(data_matrix)

        # Perform background suppression using grey opening
        background = scipy.ndimage.grey_opening(data_matrix, size=(1000, 1000))
        data_matrix = np.clip(data_matrix - background, 0, None)
        data_matrix /= np.max(data_matrix)
        return data_matrix

    @staticmethod
    def step4_interpolate(data):
        trace, delay, wvl = data
        frequency = (3e8) / (wvl * 1e-9)       # not a constant step

        # Construct frequency vector with constant step
        f_min = np.min(frequency)
        f_max = np.max(frequency)
        Df = f_max - f_min
        df = np.abs(frequency[1:] - frequency[:-1])
        df_min = np.min(df)
        N_f = int(np.floor(Df / df_min))
        interp_f = np.linspace(f_min, f_min + (N_f*df_min), N_f + 1)

        # Construct time vector with constant step
        Ddelay = np.min(delay[1:] - delay[:-1])
        interp_delay = np.arange(delay[0], delay[-1] + (Ddelay / 2), Ddelay / 2)

        # interpolate 
        new_shape = (len(interp_f) / trace.shape[0], len(interp_delay) / trace.shape[1])
        interp_trace = scipy.ndimage.zoom(np.abs(trace), new_shape, order=3)
        interp_trace = np.clip(interp_trace, 0, None)
        interp_trace /= np.max(interp_trace)
        
        interp_wvl = np.linspace(wvl[0], wvl[-1], len(interp_f))
        return interp_trace, interp_delay, interp_wvl, interp_f

    def handle_guess(self, data):
        self.it_max = self.init_guess.it_max
        self.pro_delay, self.pro_wvl, self.pro_f, self.pro_time, self.pro_trace_int, self.guess_matrix = data.values()

    def handle_div_popup(self):
        popup = ErrorPopup(self)
        response = popup.exec_()  # Blocks until the dialog is closed

        if response == QDialog.Accepted:
            while self.screens_widget.count() > 0:
                widget = self.screens_widget.widget(0)
                self.screens_widget.removeWidget(widget)
                widget.deleteLater()

            # self.layout.removeWidget(self.screens_widget)
            self.init_stacked_widget()              # reinitialize widgets
            self.screens_widget.setCurrentIndex(0)  # Return to cropping screen and try again
            self.handle_stacked_widgets()

    def handle_recovery_completed(self, data):
        data_dict = {"trace": data[0], 
                     "pulse": data[1], 
                     "switch": np.flip(data[2]), 
                     "delay": self.pro_delay,
                     "wvl": self.pro_wvl, 
                     "freq": self.pro_f, 
                     "time": self.pro_time,
                     "pro_trace": self.pro_trace_int.T,
                     "error": data[3], 
                     }
        self.iterations_complete.emit(data_dict, "PtyChoPy")
        self.close()

if __name__ == "__main__":
    # file_path = r"C:\Users\lepar\OneDrive - INRS\Documents\INRS\MSc\thesis\PtyPy\Experimental_traces\20201201_2emeFROst.txt"
    # data = np.loadtxt(file_path)
    # wavelength = data[2:, 0]
    # delay = data[0, 2:]
    # raw_trace = data[2:, 2:]
    # raw_trace /= np.max(raw_trace)

    file_path = r"data\XP_FROSt_traces\20201201_FRost_3mm fused silica.txt"
    data = np.loadtxt(file_path)
    wavelength = data[2:, 0]
    delay = data[0, 2:]
    raw_trace = data[2:, 2:]
    raw_trace /= np.max(raw_trace)

    app = QApplication(sys.argv)
    main_window = PtyChoPy()
    main_window.show()
    main_window.init_vars([raw_trace, delay, wavelength])

    sys.exit(app.exec_())
