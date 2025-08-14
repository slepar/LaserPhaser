import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QLabel, QSpinBox, QHBoxLayout, QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar)
from matplotlib.figure import Figure
import pyfftw

class FFTPlot(QWidget):
    def __init__(self, data, parent=None, cmap = "binary"):
        super().__init__(parent)
        self.raw_trace, self.delay, self.wvl = data
        self.t_slice = [self.delay[len(self.delay)//4], self.delay[-len(self.delay)//4]]
        self.N_slices = [np.argmin(np.abs(self.delay - t)) for t in self.t_slice]
        self.cmap = cmap
        self.filtered_trace = None
        self.dragging_line = None  # Track which line is being dragged

        self.init_plot_space()

    def init_plot_space(self):
        plot_layout = QVBoxLayout(self)
        self.fig = Figure(figsize=(12, 6))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        plot_layout.addWidget(self.canvas)

        self.fig.clear()
        self.ax1, self.img1, self.cbar1, vlines1 = self.init_2D_plot(231, "Raw Trace")
        self.ax2, self.img2, self.cbar2, vlines2 = self.init_2D_plot(232, "Filtered Trace")
        self.ax3, self.ax3_raw_lines, self.ax3_filt_lines, _ = self.init_1D_plot(233, "Trace: Raw and Filtered")
        self.ax4, self.img4, self.cbar4, vlines4 = self.init_2D_plot(234,  "Raw FFT Trace")
        self.ax5, self.img5, self.cbar5, vlines5 = self.init_2D_plot(235,  "Filtered FFT Trace")
        self.ax6, self.ax6_raw_lines, self.ax6_filt_lines, self.filter_line = self.init_1D_plot(236, "FFT: Raw and Filtered")

        # reassemble lines for simplicity
        self.vlines1 = [vlines1[0], vlines2[0], vlines4[0], vlines5[0]]
        self.vlines2 = [vlines1[1], vlines2[1], vlines4[1], vlines5[1]]

        self.fig.tight_layout()
        self.canvas.draw()

        self.canvas.mpl_connect('button_press_event', self.on_press)
        self.canvas.mpl_connect('button_release_event', self.on_release)
        self.canvas.mpl_connect('motion_notify_event', self.on_motion)

    def init_2D_plot(self, position, title):
        ax = self.fig.add_subplot(position)
        img = ax.imshow(
            np.zeros((128, 128)),
            # extent = [np.min(x), np.max(x), np.min(y), np.max(y)],
            aspect = "auto",
            cmap = self.cmap,
        )
        ax.set_title(title)
        ax.set_xlabel('Time Delay, $\\tau$ [fs]')
        ax.set_ylabel("Frequency, $f$ [Hz]" if "FFT" in title else "Wavelength, $\\lambda$ [nm]")
        vline1 = ax.axvline(x = self.t_slice[0], color = "green", linestyle = "--", linewidth = 2)
        vline2 = ax.axvline(x = self.t_slice[1], color = "red", linestyle = "--", linewidth = 2)
        cbar = self.fig.colorbar(img, ax = ax, orientation = "vertical")
        cbar.set_label("Intensity [a.u.]", rotation = 270, labelpad = 15)
        return ax, img, cbar, [vline1, vline2]

    def init_1D_plot(self, position, title):
        ax = self.fig.add_subplot(position)
        line1_filt, = ax.plot([], [], "g-", label="Filtered")
        line2_filt, = ax.plot([], [], "r-")
        line1_raw, = ax.plot([], [], "g--", label="Raw")
        line2_raw, = ax.plot([], [], "r--")
        filter, = ax.plot([], [], "blue")
        ax.set_ylim(-0.1, 1.1)
        ax.set_title(title)
        ax.set_xlabel("Frequency, $f$ [Hz]" if "FFT" in title else "Wavelength, $\\lambda$ [nm]")
        ax.set_ylabel("Amplitude [a.u.]")
        ax.legend()
        return ax, [line1_raw, line2_raw], [line1_filt, line2_filt], filter
    
    def update_2D(self, trace, delay, wvl, position_ID):
        img, cbar = {1: (self.img1, self.cbar1), 
                2: (self.img2, self.cbar2),
                4: (self.img4, self.cbar4), 
                5: (self.img5, self.cbar5),
                }.get(position_ID, (None, None, None, None))
        img.set_data(np.abs(trace))
        img.set_clim(vmin = np.abs(trace).min(), vmax = np.abs(trace).max())
        cbar.update_normal(img)
        img.set_extent([delay[0], delay[-1], wvl[-1], wvl[0]])

    def update_1D(self, x, y_raw, y_filt, position_ID):
        ax, raw_lines, filt_lines = {3: (self.ax3, self.ax3_raw_lines, self.ax3_filt_lines), 
                                     6: (self.ax6, self.ax6_raw_lines, self.ax6_filt_lines),
                                    }.get(position_ID, (None, None))
        for i, slice_idx in enumerate(self.N_slices):
            raw_lines[i].set_data(x, y_raw[:, slice_idx])
            filt_lines[i].set_data(x, y_filt[:, slice_idx])
        ax.relim()
        ax.autoscale()
        ax.set_ylim(-0.1, 1.1)

    def update_filter_line(self, x, filter):    
        self.filter_line.set_data(x, filter)
        roi = np.where(filter != 0)[0]
        self.ax6.set_xlim(np.min(x[roi - len(filter)//6]), np.max(x[roi + len(filter)//6]))
        self.ax6.set_ylim(-0.1, 1.1)

    def update_all_plots(self, filt_trace, trace_fft, trace_fft_filt, freq, filter):
        self.filtered_trace = np.abs(filt_trace)
        self.trace_fft = np.abs(trace_fft) / np.max(np.abs(trace_fft))
        self.trace_fft_filt = np.abs(trace_fft_filt) / np.max(np.abs(trace_fft_filt))
        self.freq = freq
        self.filter = filter

        self.update_2D(self.filtered_trace, self.delay, self.wvl, 2)
        self.update_2D(self.trace_fft, self.delay, self.freq, 4)
        self.update_2D(self.trace_fft_filt, self.delay, self.freq, 5)
        self.update_1D(self.wvl, self.raw_trace, self.filtered_trace, 3)
        self.update_1D(self.freq, self.trace_fft, self.trace_fft_filt, 6)
        self.update_filter_line(self.freq, self.filter)
        self.canvas.draw()

    def on_press(self, event):
        """Handle mouse button press - when pressing, show raw data"""
        # Check if the event occurred in the target subplot
        if event.inaxes == self.ax1 or event.inaxes == self.ax2:
            if any(abs(event.xdata - vline.get_xdata()[0]) < 10 for vline in self.vlines1):
                self.dragging_line = "vmin"
            elif any(abs(event.xdata - vline.get_xdata()[0]) < 10 for vline in self.vlines2):
                self.dragging_line = "vmax"

    def on_release(self, event):
        """Handle mouse button release - when released, plot processed."""
        # Check if the event occurred in the target subplot
        if event.inaxes in (self.ax1, self.ax2, self.ax4, self.ax5) and self.dragging_line is not None:
            if self.dragging_line == "vmin":
                self.t_slice = [int(event.xdata), self.t_slice[1]]
            elif self.dragging_line == "vmax":
                self.t_slice = [self.t_slice[0], int(event.xdata)]
            self.update_vlines()
            self.dragging_line = None

    def on_motion(self, event):
        """Handle mouse motion to drag lines."""
        if (event.inaxes == self.ax1 or event.inaxes == self.ax2) and self.dragging_line is not None:
            if self.dragging_line == "vmin":
                self.t_slice = [int(event.xdata), self.t_slice[1]]
            elif self.dragging_line == "vmax":
                self.t_slice = [self.t_slice[0], int(event.xdata)]
            self.update_vlines()

    def update_vlines(self):
        self.N_slices = [np.argmin(np.abs(self.delay - t)) for t in self.t_slice]
        if self.dragging_line == "vmin":
            for vline in self.vlines1:
                vline.set_xdata([self.t_slice[0], self.t_slice[0]])
            self.update_1D_vlines("vmin")

        elif self.dragging_line == "vmax":
            for vline in self.vlines2:
                vline.set_xdata([self.t_slice[1], self.t_slice[1]])
            self.update_1D_vlines("vmax")
        self.canvas.draw()

    def update_1D_vlines(self, line):
        if line == "vmin":
            self.ax3_raw_lines[0].set_data(self.wvl, np.abs(self.raw_trace[:, self.N_slices[0]]))
            self.ax6_raw_lines[0].set_data(self.freq, np.abs(self.trace_fft[:, self.N_slices[0]]))
            self.ax3_filt_lines[0].set_data(self.wvl, np.abs(self.filtered_trace[:, self.N_slices[0]]))
            self.ax6_filt_lines[0].set_data(self.freq, np.abs(self.trace_fft_filt[:, self.N_slices[0]]))
        else:
            self.ax6_raw_lines[1].set_data(self.freq, np.abs(self.trace_fft[:, self.N_slices[1]]))
            self.ax3_raw_lines[1].set_data(self.wvl, np.abs(self.raw_trace[:, self.N_slices[1]]))
            self.ax3_filt_lines[1].set_data(self.wvl, np.abs(self.filtered_trace[:, self.N_slices[1]]))
            self.ax6_filt_lines[1].set_data(self.freq, np.abs(self.trace_fft_filt[:, self.N_slices[1]]))

class FFT_Filter(QWidget):
    def __init__(self, parent=None, cmap = "binary"):
        super().__init__(parent)
        self.filtered_trace = None
        self.cmap = cmap
        self.plot_widget = None

        self.setWindowTitle("FFT Filtering Window")
        self.setGeometry(200, 200, 1400, 800)
        self.init_ui()

    def init_vars(self, raw_trace, delay, wavelength):
        self.raw_delay = delay
        self.raw_wvl = wavelength
        self.raw_trace = raw_trace

        if self.plot_widget is not None:
            self.main_layout.removeWidget(self.plot_widget)

        self.plot_widget = FFTPlot(data = [raw_trace, delay, wavelength], cmap = self.cmap)
        self.main_layout.addWidget(self.plot_widget)
        self.plot_widget.update_2D(self.raw_trace, self.raw_delay, self.raw_wvl, 1)           # plotting raw data

        self.perform_fft_filter()

    def init_ui(self):
        # Add input area
        self.input_widget = QWidget()
        self.input_layout = QHBoxLayout()
        self.filter_length_spinbox = QSpinBox()
        self.filter_length_spinbox.setMinimum(1)
        self.filter_length_spinbox.setValue(15)
        self.input_layout.addWidget(QLabel("Filter Length:"))
        self.input_layout.addWidget(self.filter_length_spinbox)
        self.input_widget.setLayout(self.input_layout)
        self.filter_length_spinbox.editingFinished.connect(self.perform_fft_filter)

        # Set layout
        self.main_layout = QVBoxLayout()
        self.main_layout.addWidget(self.input_widget)
        self.setLayout(self.main_layout)
        
    def perform_fft_filter(self):
        self.N_filter = self.filter_length_spinbox.value()

        # Perform FFT on each row of spectrogram
        self.trace_fft = self._fft_(self.raw_trace)
        self.freq = np.arange(0, len(self.raw_wvl), 1) - (len(self.raw_wvl)) // 2

        # Create the hyper-Gaussian filter function and matrix
        self.Order_filter = 80
        self.hg_filter_fun = np.exp(- (self.freq / self.N_filter) ** self.Order_filter)
        self.hg_filter = np.meshgrid(self.raw_delay, self.hg_filter_fun)[1]  

        # Apply the filter and ifft
        self.trace_fft_filt = self.trace_fft * self.hg_filter
        self.filtered_trace = self._ifft_(self.trace_fft_filt)
        self.filtered_trace /= np.max(np.abs(self.filtered_trace))

        # update the plot
        self.plot_widget.update_all_plots(self.filtered_trace, self.trace_fft, self.trace_fft_filt, self.freq, self.hg_filter_fun)

    @staticmethod
    def _ifft_(fft_trace):
        fftw_obj = pyfftw.builders.ifft(np.zeros_like(fft_trace[:, 0], dtype = np.complex128))
        fftshift = np.fft.fftshift
        ifftshift = np.fft.ifftshift
        trace = np.zeros_like(fft_trace, dtype = np.complex128)
        for i in range(fft_trace.shape[1]):
            trace[:, i] = ifftshift(fftw_obj(fftshift(fft_trace[:, i])))
        return trace
    
    @staticmethod
    def _fft_(trace):
        fftw_obj = pyfftw.builders.fft(np.zeros_like(trace[:, 0], dtype = np.complex128))
        fftshift = np.fft.fftshift
        ifftshift = np.fft.ifftshift
        fft_trace = np.zeros_like(trace, dtype = np.complex128)
        for i in range(trace.shape[1]):
            fft_trace[:, i] = fftshift(fftw_obj(ifftshift(trace[:, i])))
        return fft_trace

if __name__ == "__main__":
    app = QApplication(sys.argv)

    # Example data
    file_path = r"C:\Users\lepar\OneDrive - INRS\Documents\INRS\MSc\thesis\PtyPy\Experimental_traces\20201201_2emeFROst.txt"
    data = np.loadtxt(file_path)

    wavelength = data[2:, 0]
    delay = data[0, 2:]
    raw_trace = data[2:, 2:]
    raw_trace /= np.max(raw_trace)

    print(delay[0], delay[-1])
    # Remove outliers    
    # trace0 = step0_processing(raw_trace)
    # crop window here 

    # FFT filter window
    main_window = FFT_Filter()
    main_window.init_vars(raw_trace, delay, wavelength)
    main_window.show()

    # print(main_window.filtered_trace)

    sys.exit(app.exec_())
