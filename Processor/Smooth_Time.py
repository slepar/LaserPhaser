import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QLabel, QSpinBox, QHBoxLayout, QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
from scipy.signal import convolve2d
from matplotlib.gridspec import GridSpec

class Smooth_Time(QWidget):
    def __init__(self, parent = None, cmap = "binary"):
        super().__init__(parent)
        self.Order_smooth = 6
        self.Size_smooth = 0.9
        self.smooth_trace = None
        self.cmap = cmap
        self.setWindowTitle("Smoothing Window")
        self.setGeometry(200, 200, 1200, 800)

        self.init_ui()

    def init_vars(self, data):
        self.raw_trace, self.raw_delay, self.raw_wvl = data

        self.perform_smoothing()
        self.init_plot()

    def init_ui(self):
        # Add input area
        self.input_widget = QWidget()
        self.input_layout = QHBoxLayout()
        self.filter_length_spinbox = QSpinBox()
        self.filter_length_spinbox.setMinimum(1)
        self.filter_length_spinbox.setValue(10)
        self.input_layout.addWidget(QLabel("Filter Length:"))
        self.input_layout.addWidget(self.filter_length_spinbox)
        self.input_widget.setLayout(self.input_layout)
        self.filter_length_spinbox.editingFinished.connect(self.perform_smoothing)

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
        gs = GridSpec(3, 1, figure = self.fig)  # 3 rows for proportions

        self.ax0 = self.fig.add_subplot(gs[0, 0])
        self.line, = self.ax0.plot(self.x_smooth, self.V_smooth, color = "blue")
        self.ax0.set_ylabel('Amplitude [a.u.]')
        self.ax0.set_title('Smoothing Vector')

        self.ax1 = self.fig.add_subplot(gs[1:, 0]) 
        self.img1 = self.ax1.imshow(
            np.abs(self.smooth_trace),
            extent = [self.raw_delay[0], self.raw_delay[-1], self.raw_wvl[-1], self.raw_wvl[0]],
            aspect = "auto",
            cmap = self.cmap,
        )
        self.ax1.set_title("Smoothed Trace")
        self.ax1.set_xlabel('Time Delay, $\\tau$ [fs]')
        self.ax1.set_ylabel('Wavelength, $\\lambda$ [nm]')
        self.cbar1 = self.fig.colorbar(self.img1, ax = self.ax1, orientation = 'vertical')
        self.cbar1.set_label('Intensity [a.u.]', rotation = 270, labelpad = 15)

        self.fig.tight_layout()
        self.canvas.draw()
        
    def perform_smoothing(self):
        self.N_smooth = self.filter_length_spinbox.value()

        # Copy first and last rows with smoothing
        self.smooth_trace = self.raw_trace.copy()
        for _ in range(self.N_smooth):        
            self.smooth_trace = np.concatenate((self.raw_trace[:, [0]], self.smooth_trace, self.raw_trace[:, [-1]]), axis=1)
        for _ in range(self.N_smooth):        
            self.smooth_trace = np.concatenate((self.smooth_trace[[0], :], self.smooth_trace, self.smooth_trace[[-1], :]), axis=0)

        # Create smoothing vector and matrix
        self.x_smooth = np.arange(self.N_smooth)
        self.V_smooth = np.exp(-((self.x_smooth - (self.N_smooth - 1) / 2) / (self.Size_smooth * (self.N_smooth - 1) / 2)) ** self.Order_smooth)
        M_smooth = self.V_smooth.reshape(1, -1)

        # Convolve smoothing matrix with trace and crop
        self.smooth_trace = convolve2d(self.smooth_trace, M_smooth, mode='same')
        self.smooth_trace = self.smooth_trace[self.N_smooth:self.smooth_trace.shape[0] - self.N_smooth, 
                                              self.N_smooth:self.smooth_trace.shape[1] - self.N_smooth]

        # update the plot
        self.init_plot()

    def on_press(self, event):
        if event.inaxes == self.ax1:
            self.update_plot(self.raw_trace)

    def on_release(self, event):
        if event.inaxes == self.ax1:
            self.update_plot(self.smooth_trace)

    def update_plot(self, trace):      
        self.img1.set_data(np.abs(trace))
        self.img1.set_clim(vmin = np.abs(trace).min(), vmax = np.abs(trace).max())
        self.cbar1.update_normal(self.img1)
        self.canvas.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)

    # Example data
    file_path = r"C:\Users\lepar\OneDrive - INRS\Documents\INRS\MSc\thesis\PtyPy\Experimental_traces\20201201_2emeFROst.txt"
    data = np.loadtxt(file_path)

    wavelength = data[2:, 0]
    delay = data[0, 2:]
    raw_trace = data[2:, 2:]
    raw_trace /= np.max(raw_trace)

    # Remove outliers    
    # trace0 = step0_processing(raw_trace)

    # crop window here 

    # FFT filter window
    main_window = Smooth_Time()
    main_window.init_vars(raw_trace, delay, wavelength)
    main_window.show()

    # print(main_window.filtered_trace)

    sys.exit(app.exec_())
