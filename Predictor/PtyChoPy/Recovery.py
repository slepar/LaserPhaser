import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QLabel, QWidget, QVBoxLayout, QMessageBox, QDialog
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib.figure import Figure
from Software.Predictor.PtyChoPy.Cython_wrappers.fc_iteration_HIO import *
from Software.Predictor.PtyChoPy.Cython_wrappers.fc_iteration_ER import *
from PyQt5.QtCore import QThread, pyqtSignal
from matplotlib.gridspec import GridSpec
import math 

class IterationThread(QThread):
    """ Threading class that allows the recovery iterations to run independently of the GUI updates. 
        This way we can update the plot and display line for every n_plots during iterations. """
    iteration_done = pyqtSignal(int, object, object, object, object)
    def __init__(self, XP_trace, initguess_matrix, time, delay, freq, it_max, parent=None):
        super().__init__(parent)
        self.XP_trace, self.M_Produit_t_t0_guess, self.t, self.t0, self.freq, self.it_max = XP_trace, initguess_matrix, time, delay, freq, it_max

        self.k_max = 10
        self.it_Objt_amp_2_amp_n_phi = 100
        self.it_Objt_amp_n_phi_2_amp = 10000
        self.it_HIO_2_ER = 2

    def run(self):
        M_Produit_t_t0_it = self.M_Produit_t_t0_guess.copy()
        for it in range(self.it_max):
            M_Produit_t_t0_it, Trace_amp_omg_t0_it, Chp_t_it, Obj_t_it, Err_it = self.iteration_step(
                it, M_Produit_t_t0_it)           

            # FOR TESTING
            # Err_it = math.nan

            self.iteration_done.emit(it, Trace_amp_omg_t0_it, Chp_t_it, Obj_t_it, Err_it)
            if math.isnan(Err_it):
                break

    def iteration_step(self, it, M_Produit_t_t0_it):
        Case_objet = 1 if it < self.it_Objt_amp_2_amp_n_phi or it >= self.it_Objt_amp_n_phi_2_amp else 0
        if it < self.it_HIO_2_ER:
            Trace_amp_omg_t0_it, M_Produit_t_t0_it, Obj_t_it, Chp_t_it, Err_it = fc_iteration_HIO(
                self.XP_trace, M_Produit_t_t0_it, self.t, self.t0, self.k_max, Case_objet)
        else:
            Trace_amp_omg_t0_it, M_Produit_t_t0_it, Obj_t_it, Chp_t_it, Err_it = fc_iteration_ER(
                self.XP_trace, M_Produit_t_t0_it, self.t, self.t0, self.k_max, Case_objet)
        return M_Produit_t_t0_it, Trace_amp_omg_t0_it, Chp_t_it, Obj_t_it, Err_it
    
class Recovery(QWidget):
    recovery_done = pyqtSignal(object)
    divergence_occuring = pyqtSignal()

    def __init__(self, parent = None, cmap = "binary"):
        super().__init__(parent)
        self.cmap = cmap
        self.init_ui()
        
    def init_ui(self):
        layout = QVBoxLayout()
        self.display_line = QLabel("Starting recovery...")
        self.figure = Figure(figsize = (8, 6))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.display_line)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)      

    def init_plot(self):
        gs = GridSpec(3, 2, figure = self.figure) 
        self.ax1 = self.figure.add_subplot(gs[1, 0])
        self.img1 = self.ax1.imshow(np.abs(self.XP_trace_int).T, aspect = 'auto', cmap = self.cmap,
                                    extent = [self.t0.min(), self.t0.max(), self.freq.min(), self.freq.max()])
        self.ax1.set_ylabel('Frequency, $f$ [Hz]')
        self.ax1.set_xlabel('Time Delay, $\\tau$ [fs]')
        self.ax1.set_title('Processed FROSt')
        self.cbar1 = self.figure.colorbar(self.img1, ax = self.ax1, orientation = 'vertical')
        self.cbar1.set_label('Intensity [a.u.]', rotation = 270, labelpad = 15)

        self.ax2 = self.figure.add_subplot(gs[1, 1])
        self.img2 = self.ax2.imshow(np.zeros_like(self.XP_trace_int), aspect = 'auto', cmap = self.cmap,
                                    extent = [self.t0.min(), self.t0.max(), self.freq.min(), self.freq.max()])
        self.ax2.set_ylabel('Frequency, $f$ [Hz]')
        self.ax2.set_xlabel('Time Delay, $\\tau$ [fs]')
        self.ax2.set_title('Reconstructed FROSt')
        self.cbar2 = self.figure.colorbar(self.img2, ax = self.ax2, orientation='vertical')
        self.cbar2.set_label('Intensity [a.u.]', rotation = 270, labelpad = 15)

        self.ax3_amp = self.figure.add_subplot(gs[2, 0])
        self.ax3_line1, = self.ax3_amp.plot([], [], color = "red", linestyle = "-")
        self.ax3_amp.set_title('Reconstructed Pulse')
        self.ax3_phase = self.ax3_amp.twinx()
        self.ax3_line2, = self.ax3_phase.plot([], [], color = "red", linestyle = "--")
        self.ax3_amp.set_xlabel('Time, $t$ [fs]')
        self.ax3_amp.set_ylabel('Amplitude [a.u.]')
        self.ax3_phase.set_ylabel('Phase [rads]', rotation = 270, labelpad = 15)

        self.ax4 = self.figure.add_subplot(gs[2, 1])
        self.ax4_line, = self.ax4.plot([], [], color = "blue", linestyle = "-")
        self.ax4.set_xlabel('Time, $t$ [fs]')
        self.ax4.set_ylabel('Amplitude [a.u.]')
        self.ax4.set_title('Reconstructed Switch')
        
        self.ax5 = self.figure.add_subplot(gs[0, :2])
        self.ax5_line, = self.ax5.plot([], [], color = "black", linestyle = "-")
        self.ax5.set_xlabel('Iterations, $i$')
        self.ax5.set_ylabel('Error [%]')
        self.ax5.set_title('Recovery Error')
        self.ax5.set_xlim([0, self.it_max + 1])
        self.ax5.set_ylim([0, 100])

        self.figure.tight_layout()
        self.canvas.draw()

    def init_vars(self, XP_trace_int, initguess_matrix, delay, time, freq, it_max):
        # variables for iterations
        self.it_max = it_max
        self.n_plot = it_max // 1 # plots every iterations
        self.V_Err = np.empty(self.it_max + 1)

        # Checkers for shapes
        if freq.shape != time.shape:
            raise ValueError("Time and frequency vectors must have same number of points.")
        elif not (XP_trace_int.shape[1] == time.shape[0] and XP_trace_int.shape[0] == delay.shape[0]):
            XP_trace_int = XP_trace_int.T
            if not (XP_trace_int.shape[1] == time.shape[0] and XP_trace_int.shape[0] == delay.shape[0]):
                raise ValueError("Trace must have same number of points as time and delay vectors.")
        elif XP_trace_int.shape != initguess_matrix.shape:
            raise ValueError("Initial guess matrix must have same shape as trace.")

        self.XP_trace_int, self.M_Produit_t_t0_guess, self.t, self.t0, self.freq, self.it_max = XP_trace_int, initguess_matrix, time, delay, freq, it_max
        self.XP_trace_int /= np.max(self.XP_trace_int)

        self.init_plot()
        self.init_iterations()

    def init_iterations(self):
        self.thread = IterationThread(self.XP_trace_int, self.M_Produit_t_t0_guess, self.t, self.t0, self.freq, self.it_max)
        self.thread.iteration_done.connect(self.update_display)
        self.thread.start()

    def update_display(self, it, Trace_amp_omg_t0_it, Chp_t_it, Obj_t_it, Err_it):
        self.V_Err[it] = Err_it
        self.display_line.setText(f"Iteration {it + 1}/{self.it_max} completed with a recovery error of {Err_it*10:0.02f}%.")    
        if math.isnan(Err_it):
            self.divergence_occuring.emit()
        
        predicted_trace_int = np.abs(Trace_amp_omg_t0_it.T)**2 / np.max(np.abs(Trace_amp_omg_t0_it.T)**2)
        amplitude = np.abs(Chp_t_it[0]) / np.max(np.abs(Chp_t_it[0]))
        roi = np.where(amplitude >= 0.01)[0]

        self.img2.set_data(np.abs(predicted_trace_int))
        self.img2.set_clim(vmin = np.abs(predicted_trace_int).min(), vmax = np.abs(predicted_trace_int).max())
        self.cbar2.update_normal(self.img2)
        self.ax3_line1.set_data(self.t[roi], amplitude[roi])
        self.ax3_line2.set_data(self.t[roi], np.unwrap(np.angle(Chp_t_it[0, roi])))
        self.ax3_amp.set_title(f'Reconstructed Pulse, it = {it + 1} / {self.it_max}')
        self.ax3_phase.relim()
        self.ax3_amp.relim()
        self.ax3_amp.autoscale_view()
        self.ax3_phase.autoscale_view()
        self.ax4_line.set_data(self.t, np.flip(np.abs(Obj_t_it[0]) / np.max(np.abs(Obj_t_it[0]))))
        self.ax4.relim()
        self.ax4.autoscale_view()
        self.ax4.set_title(f'Reconstructed Switch, it = {it + 1} / {self.it_max}')
        
        self.ax5_line.set_data(np.arange(it + 1), (self.V_Err[:it + 1]*10))
        self.ax5.relim()
        self.ax5.autoscale_view()
        self.canvas.draw()
        
        # send final recovery to main window
        if it == self.it_max - 1:
            self.recovery_done.emit([predicted_trace_int, Chp_t_it[0], Obj_t_it[0], self.V_Err])

import os 

if __name__ == '__main__':

    app = QApplication(sys.argv)

    path = r"C:\Users\lepar\OneDrive - INRS\Documents\INRS\MSc\thesis\PtyPy\recovered_traces\20201201_FRost_3mm fused silica.txt"
    recovered_data = {}
    for filename in os.listdir(path):
        if filename.endswith(".txt"):
            filepath = os.path.join(path, filename)
            try:
                recovered_data[os.path.splitext(filename)[0]] = np.loadtxt(filepath)
            except ValueError:      # handles the case where pulse is complex
                recovered_data[os.path.splitext(filename)[0]] = np.loadtxt(filepath, dtype = complex)

    print(recovered_data.keys())

    widget2 = Recovery()
    widget2.show()

    widget2.init_vars(recovered_data['processed_trace_int'], recovered_data['initial_guess'], recovered_data['processed_delay'], 
                                             recovered_data['processed_time'], recovered_data['processed_freq'], 10)
    sys.exit(app.exec_())
