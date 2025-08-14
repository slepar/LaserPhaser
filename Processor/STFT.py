import sys
import numpy as np
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QHBoxLayout, QPushButton, 
    QWidget, QLabel, QLineEdit, QSpinBox, QFileDialog, QMessageBox
)
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar
)
from matplotlib.figure import Figure
from scipy.signal import stft
import pickle as pkl
from matplotlib.gridspec import GridSpec
import scipy.ndimage

class STFT(QMainWindow):
    def __init__(self, pulse, time, cmap = "binary"):
        super().__init__()
        self.setWindowTitle("STFT Analyzer")
        self.setGeometry(100, 100, 1200, 800)
        self.cmap = cmap
        self.time = time
        self.raw_pulse = pulse.numpy() if not isinstance(pulse, np.ndarray) else pulse
        self.interpolated_pulse = None

        # set up ui
        self.init_ui()
        
        # proceed to plot pulse
        self.plot_pulse(self.raw_pulse, self.time)

    def init_ui(self):
        # Central Widget
        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        
        # Layouts
        self.main_layout = QVBoxLayout(self.central_widget)
        self.input_layout = QHBoxLayout()
        self.plot_layout = QVBoxLayout()

        # Input Section
        self.spin_window_size = QSpinBox()
        self.spin_window_size.setRange(1, 1024)
        self.spin_window_size.setValue(256)
        self.spin_interp_N = QSpinBox()
        self.spin_interp_N.setRange(128, 2048)
        self.spin_interp_N.setValue(512)
        self.button_compute = QPushButton("Compute STFT")
        self.input_layout.addWidget(QLabel("Interpolation N:"))
        self.input_layout.addWidget(self.spin_interp_N)
        self.input_layout.addWidget(QLabel("Window Size:"))
        self.input_layout.addWidget(self.spin_window_size)
        self.input_layout.addWidget(self.button_compute)
        self.main_layout.addLayout(self.input_layout)

        # Plot Section
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.plot_layout.addWidget(self.toolbar)
        self.plot_layout.addWidget(self.canvas)
        self.main_layout.addLayout(self.plot_layout)

        # Connections
        self.spin_interp_N.editingFinished.connect(self.interp_pulse)
        self.button_compute.clicked.connect(self.compute_stft)
        self.canvas.mpl_connect("button_press_event", self.on_press)
        self.canvas.mpl_connect("button_release_event", self.on_release)

    def plot_pulse(self, pulse, time):
        try:
            amplitude = np.abs(pulse)/np.max(np.abs(pulse))
            phase = np.unwrap(np.angle(pulse))

            self.figure.clear()
            self.amp_ax = self.figure.add_subplot(111)
            self.amp_line, = self.amp_ax.plot(time, amplitude, color = "red", linestyle = "-", label = "Amplitude")
            self.pha_ax = self.amp_ax.twinx()
            self.pha_line, = self.pha_ax.plot(time, phase, color = "red", linestyle = "--", label = "Phase")
            self.amp_ax.set_title('Predicted Pulse')
            self.amp_ax.set_xlabel('Time [fs]')
            self.amp_ax.set_ylabel('Amplitude [a.u.]')
            self.pha_ax.set_ylabel('Phase [rad]')
            self.amp_ax.legend(handles=[self.amp_line, self.pha_line], labels=["Amplitude", "Phase"], loc=2)

            self.amp_ax.set_xlim(time.min(), time.max())
            self.pha_ax.set_xlim(time.min(), time.max())

            self.amp_ax.set_ylim(amplitude.min() - 0.1, amplitude.max() + 0.1)
            self.pha_ax.set_ylim(phase.min() - 0.5, phase.max() + 0.5)

            self.canvas.draw()
    
        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to plot pulse: {e}")

    def on_press(self, event):
        """Handle mouse button press - when pressing, show raw data"""
        # Check if the event occurred in the target subplot
        if event.inaxes == self.amp_ax or event.inaxes == self.pha_ax:
            amplitude = np.abs(self.raw_pulse)/np.max(np.abs(self.raw_pulse))
            phase = np.unwrap(np.angle(self.raw_pulse))
            self.amp_line.set_data(self.time, amplitude)
            self.pha_line.set_data(self.time, phase)
            self.amp_ax.set_xlim(self.time.min(), self.time.max())
            self.pha_ax.set_xlim(self.time.min(), self.time.max())
            self.amp_ax.set_ylim(amplitude.min() - 0.1, amplitude.max() + 0.1)
            self.pha_ax.set_ylim(phase.min() - 0.5, phase.max() + 0.5)
            self.canvas.draw()

    def on_release(self, event):
        """Handle mouse button release - when released, plot processed."""
        # Check if the event occurred in the target subplot
        if event.inaxes == self.amp_ax or event.inaxes == self.pha_ax:
            amplitude = np.abs(self.interpolated_pulse)/np.max(np.abs(self.interpolated_pulse))
            phase = np.unwrap(np.angle(self.interpolated_pulse))
            self.amp_line.set_data(self.interpolated_time, amplitude)
            self.pha_line.set_data(self.interpolated_time, phase)
            self.amp_ax.set_xlim(self.interpolated_time.min(), self.interpolated_time.max())
            self.pha_ax.set_xlim(self.interpolated_time.min(), self.interpolated_time.max())
            self.amp_ax.set_ylim(amplitude.min() - 0.1, amplitude.max() + 0.1)
            self.pha_ax.set_ylim(phase.min() - 0.5, phase.max() + 0.5)
            self.canvas.draw()

    def interp_pulse(self):
        try:
            # itnerpolate pulse first 
            interp_N = self.spin_interp_N.value()
            self.interpolated_time = np.linspace(self.time[0], self.time[-1], interp_N)
            new_shape = (interp_N / self.raw_pulse.shape[0])
            self.interpolated_pulse = scipy.ndimage.zoom(self.raw_pulse, new_shape, order=3)
            self.plot_pulse(self.interpolated_pulse, self.interpolated_time)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to interpolate: {e}")

    def compute_stft(self):
        try:
            if self.interpolated_pulse is None:       # ensyre pulse is interpolated before stft
                self.interp_pulse()

            # Parse signal and parameters
            window_size = self.spin_window_size.value()

            # Perform STFT
            fs = 1/(self.interpolated_time[1] - self.interpolated_time[0])        # self.time is in [fs]
            f, t, Zxx = stft(self.interpolated_pulse, fs = fs, nperseg = window_size, return_onesided = False)
            Zxx = np.fft.fftshift(Zxx, axes = 0)  
            f = np.fft.fftshift(f)
        
            self.plot_stft(t, f, Zxx)

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to compute STFT: {e}")

    def plot_stft(self, t, f, Zxx):
        try:
            # interpolate to make is prettier
            new_shape = (512 / Zxx.shape[0], 512 / Zxx.shape[1])
            Zxx = scipy.ndimage.zoom(np.abs(Zxx), new_shape, order=3)
            Zxx /= np.max(Zxx)
            f = np.linspace(min(f), max(f), 512)
            t = np.linspace(-t[-1]/2, t[-1]/2, 512)

            # # Plot the STFT
            self.figure.clear()
            gs = GridSpec(2, 2, width_ratios=[1, 4], height_ratios=[4, 1], figure=self.figure)

            ax1 = self.figure.add_subplot(gs[0, 0])
            ax1.plot(-np.mean(np.abs(Zxx), axis = 1), f, color='black')
            ax1.set_ylabel('Frequency, $f - f_0$ [PHz]')
            ax1.set_xticks([])
            ax2 = self.figure.add_subplot(gs[0, 1])
            pcm = ax2.pcolormesh(t, f, np.abs(Zxx), shading='gouraud', cmap=self.cmap)
            ax2.set_xticks([])
            ax2.set_yticks([])
            ax3 = self.figure.add_subplot(gs[1, 1])
            ax3.plot(t, np.mean(np.abs(Zxx), axis = 0), color='black')
            ax3.set_xlabel('Time, $t$ [fs]')
            ax3.set_yticks([])

            self.canvas.draw()

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to plot STFT: {e}")

import os
import pickle
if __name__ == "__main__":

    dir = r"Recovered_Data\raw_experimental_data.txt"
    with open(os.path.join(dir, "PCP.pkl"), 'rb') as f:
        pcp_data = pickle.load(f)
    with open(os.path.join(dir, "DFN.pkl"), 'rb') as f:
        ai_data = pickle.load(f)

    print(pcp_data.keys(), ai_data.keys())

    data = np.loadtxt(os.path.join(dir, "raw_data.txt"))
    raw_delay = data[0, 2:]
    raw_wvl = data[2:, 0]
    raw_trace = data[2:, 2:]
    raw_trace /= np.max(raw_trace)

    app = QApplication(sys.argv)
    plot = STFT(ai_data['ai_pulse'], ai_data['ai_time'])
    plot.show()

    sys.exit(app.exec_())

