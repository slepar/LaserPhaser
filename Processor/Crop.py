import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QHBoxLayout, QPushButton, QTableWidgetItem, QTableWidget, QVBoxLayout, QWidget
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector
import os
from PyQt5.QtCore import pyqtSignal

class Crop(QWidget):
    cropping_done = pyqtSignal(object)
    def __init__(self, parent = None, cmap = "binary"):
        super().__init__(parent)
        self.cropped_trace = None
        self.cmap = cmap
        
        self.setWindowTitle("Cropping Window")
        self.resize(800, 600)
        self.init_ui()

    def init_vars(self, data):
        self.raw_trace, self.delay, self.wavelength = data
        self.init_plot()
        self._autocalc_range()
        self.crop_trace()

    def init_ui(self):
        # input widget
        self.input_widget = QWidget()
        self.input_layout = QHBoxLayout()
        self.table = QTableWidget(2, 2)  # 2 rows and 2 columns
        self.table.setHorizontalHeaderLabels(["Minimum", "Maximum"])
        self.table.setVerticalHeaderLabels(["Time Delay [fs]", "Wavelength [nm]"])
        self.table.resizeColumnsToContents()
        self.table.resizeRowsToContents()
        self.table.cellChanged.connect(self._update_selector_from_table)
        self.table.setStyleSheet("""
            QTableWidget {background: transparent; border: none;}
            QHeaderView::section {background: transparent; padding: 0px;}
        """)
        self.input_layout.addWidget(self.table)
        self.input_widget.setLayout(self.input_layout)

        # Create figure and canvas
        self.plot_widget = QWidget()
        self.plot_layout = QVBoxLayout()
        self.fig = Figure(figsize = (6, 4))
        self.canvas = FigureCanvas(self.fig)
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.plot_layout.addWidget(self.toolbar)
        self.plot_layout.addWidget(self.canvas)
        self.plot_widget.setLayout(self.plot_layout)
        
        # put widget into layout
        self.main_layout = QVBoxLayout()
        self.main_layout.addWidget(self.input_widget)
        self.main_layout.addWidget(self.plot_widget)
        self.setLayout(self.main_layout)

    def _autocalc_range(self):
        # auto position range at 1% intensity for wavelength
        wind = np.where(np.mean(self.raw_trace, axis = 1) >= 0.01)[0]
        self.wmin, self.wmax = self.wavelength[min(wind)], self.wavelength[max(wind)]
        
        # auto position range at 1/3 of the way for t
        self.tmin, self.tmax = self.delay[len(self.delay)//3], self.delay[-len(self.delay)//3]

        self._update_table()

    def init_plot(self):
        # for selecting the region
        self.ax1 = self.fig.add_subplot(121)
        self.img1 = self.ax1.imshow(self.raw_trace, aspect="auto", cmap=self.cmap,
            extent=[self.delay[0], self.delay[-1], self.wavelength[-1], self.wavelength[0]],
        )
        self.ax1.set_title("Drag to select a region")
        self.ax1.set_xlabel('Time Delay, $\\tau$ [fs]')
        self.ax1.set_ylabel('Wavelength, $\\lambda$ [nm]')
        self.cbar1 = self.fig.colorbar(self.img1, ax = self.ax1, orientation = 'vertical')
        self.cbar1.set_label('Intensity [a.u.]', rotation = 270, labelpad = 15)

        self.ax2 = self.fig.add_subplot(122)
        self.img2 = self.ax2.imshow(np.zeros_like(self.raw_trace), aspect="auto", cmap=self.cmap, )
        self.ax2.set_xlabel('Time Delay, $\\tau$ [fs]')
        self.ax2.set_ylabel('Wavelength, $\\lambda$ [nm]')
        self.cbar2 = self.fig.colorbar(self.img2, ax = self.ax2, orientation = 'vertical')
        self.cbar2.set_label('Intensity [a.u.]', rotation = 270, labelpad = 15)
        self.ax2.set_title("Selected Region")

        # Initialize RectangleSelector
        self.rect_selector = RectangleSelector(
            self.ax1, self._on_select, useblit=True,
            button=[1],  # Left mouse button
            minspanx=5,
            minspany=5,
            spancoords="data",
            interactive=True,
        )

        self.fig.tight_layout()
        self.canvas.draw()

    def _on_select(self, eclick, erelease):
        """
        Callback for rectangle selection.
        """
        # Get selected region boundaries
        self.tmin, self.tmax, self.wmin, self.wmax = self.rect_selector.extents
        self.crop_trace()

    def _update_table(self):
        # Update UI with the selected region 
        # temporarily stop signal from sending between selector and ui, otherwise is loops
        self.table.blockSignals(True)  # Temporarily block signals
        self.table.setItem(0, 0, QTableWidgetItem(f"{self.tmin:0.02f}"))
        self.table.setItem(0, 1, QTableWidgetItem(f"{self.tmax:0.02f}"))
        self.table.setItem(1, 0, QTableWidgetItem(f"{self.wmin:0.02f}"))
        self.table.setItem(1, 1, QTableWidgetItem(f"{self.wmax:0.02f}"))
        self.table.blockSignals(False)  # Re-enable signals

    def crop_trace(self):
        t_ind_min, t_ind_max = np.argmin(np.abs(self.delay - self.tmin)), np.argmin(np.abs(self.delay - self.tmax))
        w_ind_min, w_ind_max = np.argmin(np.abs(self.wavelength - self.wmin)), np.argmin(np.abs(self.wavelength - self.wmax))

        self.cropped_trace = self.raw_trace[w_ind_min:w_ind_max, t_ind_min:t_ind_max] / np.max(self.raw_trace[w_ind_min:w_ind_max, t_ind_min:t_ind_max])
        self.cropped_delay = self.delay[t_ind_min:t_ind_max]
        self.cropped_wvl = self.wavelength[w_ind_min:w_ind_max]
        self._update_fig2(self.cropped_trace, self.cropped_delay, self.cropped_wvl)
        self._update_table()

    def _update_fig2(self, trace, delay, wvl):
        self.img2.set_data(trace)
        self.img2.set_extent([delay[0], delay[-1], wvl[-1], wvl[0]])
        self.img2.set_clim(vmin = trace.min(), vmax = trace.max())
        self.cbar2.update_normal(self.img2)
        self.canvas.draw()

    def _update_selector_from_table(self):
        self.tmin = float(self.table.item(0, 0).text())
        self.tmax = float(self.table.item(0, 1).text())
        self.wmin = float(self.table.item(1, 0).text())
        self.wmax = float(self.table.item(1, 1).text())
            
        # Update RectangleSelector's extents and recrop
        self.rect_selector.extents = (self.tmin, self.tmax, self.wmin, self.wmax)
        self.crop_trace()
    
    def closeEvent(self, event):
        # Emit the signal with the data
        self.cropping_done.emit([self.cropped_trace, self.cropped_delay, self.cropped_wvl])
        super().closeEvent(event)

if __name__ == "__main__":
    app = QApplication(sys.argv)

    # Example data
    # file_path = os.path.join(os.getcwd(), 'Experimental_traces', '20201201_FRost_3mm fused silica.txt')
    file_path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\data\20201201_FRost_3mm fused silica.txt\raw_experimental_data.txt"
    data = np.loadtxt(file_path)

    wavelength = data[2:, 0]
    delay = data[0, 2:]
    raw_trace = data[2:, 2:]
    raw_trace /= np.max(raw_trace)

    # remove outliers    
    # trace0 = step0_processing(raw_trace)

    # crop window
    main_window = Crop()
    main_window.init_vars([raw_trace, delay, wavelength])
    main_window.show()
    trace1, delay1, delay1 = main_window.cropped_trace, main_window.cropped_delay, main_window.cropped_wvl

    sys.exit(app.exec_())
