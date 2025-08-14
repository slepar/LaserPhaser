import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QLabel, QHBoxLayout, QDoubleSpinBox, QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
from matplotlib.gridspec import GridSpec

class Smooth_Freq(QWidget):
    def __init__(self, parent=None, cmap = "binary"):
        super().__init__(parent)
        self.cmap = cmap
        self.filter_order = 6
        self.init_ui()

        self.setWindowTitle("Smoothing Window")
        self.setGeometry(200, 200, 1200, 800)

    def init_vars(self, data):
        self.raw_trace, self.raw_delay, self.raw_wvl, self.raw_freq = data
        self.magnitude = int(np.floor(np.log10(abs(np.max(self.raw_freq)))))

        # Spectre Non Perturbé
        self.raw_spectrum = np.mean(self.raw_trace[:, :10], axis = 1)
        self.raw_spectrum /= np.max(self.raw_spectrum)

        self.init_plot()
        self.perform_smoothing()

    def init_ui(self):
        # Add input area
        self.input_widget = QWidget()
        self.input_layout = QHBoxLayout()
        self.filter_min_spinbox = QDoubleSpinBox()
        self.filter_min_spinbox.setValue(0.2)
        self.filter_max_spinbox = QDoubleSpinBox()
        self.filter_max_spinbox.setValue(10)
        self.input_layout.addWidget(QLabel(f"Filter minimum:"))
        self.input_layout.addWidget(self.filter_min_spinbox)
        self.input_layout.addWidget(QLabel(f"Filter maximum:"))
        self.input_layout.addWidget(self.filter_max_spinbox)
        self.input_widget.setLayout(self.input_layout)
        self.filter_min_spinbox.valueChanged.connect(self.perform_smoothing)
        self.filter_max_spinbox.valueChanged.connect(self.perform_smoothing)

        # Create figure and canvas
        self.plot_widget = QWidget()
        self.plot_layout = QVBoxLayout()
        self.fig = Figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.plot_layout.addWidget(self.toolbar)
        self.plot_layout.addWidget(self.canvas)
        self.plot_widget.setLayout(self.plot_layout)
        self.canvas.mpl_connect("button_press_event", self.on_press)
        self.canvas.mpl_connect("button_release_event", self.on_release)

        # Set layout
        self.main_layout = QVBoxLayout()
        self.main_layout.addWidget(self.input_widget)
        self.main_layout.addWidget(self.plot_widget)
        self.setLayout(self.main_layout)

    def init_plot(self):
        self.fig.clear()
        gs = GridSpec(3, 1, figure = self.fig)

        self.ax0 = self.fig.add_subplot(gs[0, 0])
        self.line1, = self.ax0.plot(self.raw_freq, np.abs(self.raw_spectrum), color = "red", label = "Raw")
        self.line2, = self.ax0.plot(self.raw_freq, np.abs(self.raw_spectrum), color = "blue", label = "Filter")
        self.line3, = self.ax0.plot(self.raw_freq, np.abs(self.raw_spectrum), color = "green", label = "Filtered")
        self.ax0.set_xlabel('Frequency, $f$ [Hz]')
        self.ax0.set_ylabel('Amplitude [a.u.]')
        self.ax0.legend(loc='upper left')
        self.ax0.set_title("Spectrum")

        # Subfigure trace
        self.ax1 = self.fig.add_subplot(gs[1:, 0])  # Use the last 2 rows for this plot
        self.img1 = self.ax1.imshow(
            self.raw_trace,
            extent=[self.raw_delay[0], self.raw_delay[-1], self.raw_wvl[-1], self.raw_wvl[0]],
            aspect="auto",
            cmap=self.cmap,
        )
        self.ax1.set_title("Frequency Smoothed Trace")
        self.ax1.set_xlabel('Time Delay, $\\tau$ [fs]')
        self.ax1.set_ylabel('Wavelength, $\\lambda$ [nm]')
        self.cbar1 = self.fig.colorbar(self.img1, ax = self.ax1, orientation = 'vertical')
        self.cbar1.set_label('Intensity [a.u.]', rotation = 270, labelpad = 15)

        self.fig.tight_layout()
        self.canvas.draw_idle()
        
    def perform_smoothing(self):
        self.f_max = self.filter_max_spinbox.value() * (10**self.magnitude)
        self.f_min = self.filter_min_spinbox.value() * (10**self.magnitude)

        # Filtrage
        self.filter_vector = np.exp(-((self.raw_freq - np.abs(self.f_min + self.f_max) / 2) / (np.abs(self.f_min - self.f_max) / 2)) ** self.filter_order)
        self.filter_matrix, _ = np.meshgrid(self.filter_vector, self.raw_delay) 
        self.smooth_trace = self.raw_trace.T * self.filter_matrix
        
        # # multiply by 1/f^2
        M_f3, _ = np.meshgrid(self.raw_freq, self.raw_delay)
        M_f3[np.isnan(M_f3)] = 0
        M_f3 /= np.max(M_f3)
        self.smooth_trace /= M_f3**2
        self.smooth_trace = self.smooth_trace.T  #/ np.max(self.smooth_trace)

        # update the plot
        self.update_plot(self.smooth_trace)

    def on_press(self, event):
        if event.inaxes == self.ax1:
            self.update_plot(self.raw_trace)

    def on_release(self, event):
        if event.inaxes == self.ax1:
            self.update_plot(self.smooth_trace)

    def update_plot(self, trace):
        # trace /= np.max(trace)
        self.img1.set_data(trace)
        self.img1.set_clim(vmin = trace.min(), vmax = trace.max())
        self.cbar1.update_normal(self.img1)

        self.line1.set_data(self.raw_freq, np.abs(self.raw_spectrum))
        self.line2.set_data(self.raw_freq, np.abs(self.filter_vector))
        self.line3.set_data(self.raw_freq, np.abs(self.raw_spectrum * self.filter_vector))

        self.canvas.draw_idle()

if __name__ == "__main__":
    app = QApplication(sys.argv)

    # Example data
    file_path = r"C:\Users\lepar\OneDrive - INRS\Documents\INRS\MSc\thesis\PtyPy\Experimental_traces\20201201_2emeFROst.txt"
    data = np.loadtxt(file_path)

    wavelength = data[2:, 0]
    delay = data[0, 2:]
    raw_trace = data[2:, 2:]
    raw_trace /= np.max(raw_trace)

    main_window = Smooth_Freq()
    main_window.init_vars([raw_trace, delay, wavelength, wavelength])
    main_window.show()

    # freq = (3e8) / (wavelength * 1e-6)       # not a constant step

    # # Spectre Non Perturbé
    # spectrum = np.mean(raw_trace[:, :10], axis = 1)
    # spectrum /= np.max(spectrum)

    # # Filtrage
    # filter_bounds = [np.min(freq), np.max(freq)] #[f*1e14 for f in filter_bounds]
    # filter_order = 6
    # filter = np.exp(-((freq - (filter_bounds[0] + filter_bounds[1]) / 2) / ((filter_bounds[0] - filter_bounds[1]) / 2)) ** filter_order)
    # filter_matrix, _ = np.meshgrid(delay, filter)
    # trace_smoothed = raw_trace * filter_matrix

    # # # multiply by 1/f^2
    # M_f3, M_delais3 = np.meshgrid(delay, freq)
    # M_f3[np.isnan(M_f3)] = 0
    # M_f3 /= np.max(M_f3)
    # trace_smoothed /= M_f3**2

    # plt.figure(figsize=(6, 6))
    # plt.subplot(3, 1, 1)
    # plt.imshow(raw_trace, aspect='auto', extent=[delay.min(), delay.max(), freq.min(), freq.max()])
    # # plt.ylim([V_f2.min(), V_f2.max()])
    # plt.ylabel('f (Hz)')
    # plt.xlabel('t (fs)')
    # plt.subplot(3, 1, 2)
    # plt.plot(freq, np.abs(spectrum))
    # plt.plot(freq, np.abs(filter))
    # plt.plot(freq, np.abs(spectrum * filter))
    # plt.xlabel('f (Hz)')
    # plt.ylabel('t (fs)')
    # plt.subplot(3, 1, 3)
    # plt.imshow(trace_smoothed, aspect='auto', extent=[delay.min(), delay.max(), freq.min(), freq.max()])
    # # plt.xlim([V_f3.min(), V_f3.max()])
    # plt.xlabel('f (Hz)')
    # plt.ylabel('t (fs)')
    # plt.tight_layout()
    # plt.show()

    sys.exit(app.exec_())
