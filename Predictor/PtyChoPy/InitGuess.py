import numpy as np
from Software.Predictor.PtyChoPy.Cython_wrappers.fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN.contrainte_ptychographic import call_fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN
from Software.Predictor.PtyChoPy.Cython_wrappers.fc_decalage_temporelle_vecteur_t_de_t0_BEN.decalage_temporelle_vecteur import call_fc_decalage_temporelle_vecteur_t_de_t0_BEN
import sys
from PyQt5.QtWidgets import QApplication, QLabel, QSpinBox, QHBoxLayout, QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import (FigureCanvasQTAgg as FigureCanvas, NavigationToolbar2QT as NavigationToolbar,)
from matplotlib.figure import Figure
import scipy.ndimage
import pyfftw
from PyQt5.QtCore import pyqtSignal

class InitGuess(QWidget):
    guess_ready = pyqtSignal(object)
    divergence_occuring = pyqtSignal()

    def __init__(self, parent = None, cmap = "binary"):
        super().__init__(parent)
        self.it_max = 10
        self.cmap = cmap
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()

        # inputs
        self.input_widget = QWidget()
        self.input_layout = QHBoxLayout()
        self.input_layout.addWidget(QLabel("Number of iterations: "))
        self.N_its_spinbox = QSpinBox()
        self.N_its_spinbox.setValue(self.it_max)
        self.input_layout.addWidget(self.N_its_spinbox)
        self.input_widget.setLayout(self.input_layout)
        layout.addWidget(self.input_widget)
        self.N_its_spinbox.valueChanged.connect(self.update_iters)

        # figure
        self.figure = Figure(figsize=(8, 6))
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def init_vars(self, raw_trace, delay, wavelength, freq):
        # 
        # interp and pad takes about 1/2 second
        self.XP_trace_int, self.delay, self.wvl, self.freq, self.time = self.interp_and_pad(raw_trace, delay, wavelength, freq)

        # Checkers for shapes (0.0007 s)
        if self.freq.shape != self.time.shape:
            raise ValueError("Time and frequency vectors must have same number of points.")
        elif not (self.XP_trace_int.shape[1] == self.time.shape[0] and self.XP_trace_int.shape[0] == self.delay.shape[0]):
            self.XP_trace_int = self.XP_trace_int.T
            if not (self.XP_trace_int.shape[1] == self.time.shape[0] and self.XP_trace_int.shape[0] == self.delay.shape[0]):
                raise ValueError("Trace must have same number of points as time and delay vectors.")

        self.f0 = np.mean(self.freq)
        self.Trace_amp_XP_omg_t0 = np.abs(np.sqrt(self.XP_trace_int))
        self.N, self.N0 = len(self.time), len(self.delay)
        self.V_omg = 2 * np.pi * self.freq
        self.omg0 = 2 * np.pi * self.f0

        self.init_guess()

    def _ifft_(self, trace):
        fftw_obj = pyfftw.builders.ifft(np.zeros_like(trace[0, :], dtype = np.complex128))
        fftshift = np.fft.fftshift
        ifftshift = np.fft.ifftshift
        ifft_trace = np.zeros_like(trace, dtype = np.complex128)
        for i in range(self.N0):
            ifft_trace[i, :] = ifftshift(fftw_obj(fftshift(trace[i, :])))
        return ifft_trace
    
    def _fft_(self, trace):
        fftw_obj = pyfftw.builders.fft(np.zeros_like(trace[0, :], dtype = np.complex128))
        fftshift = np.fft.fftshift
        ifftshift = np.fft.ifftshift
        fft_trace = np.zeros_like(trace, dtype = np.complex128)
        for i in range(self.N0):
            fft_trace[i, :] = fftshift(fftw_obj(ifftshift(trace[i, :])))
        return fft_trace
    
    def init_guess(self):    
        self.Produit_amp_t_t0 = self._ifft_(self.Trace_amp_XP_omg_t0)

        self.Obj_t_guess, self.Chp_t_guess = call_fc_Contrainte_Ptychographic_Objet_porte_ou_non_BEN(
            self.time, self.delay, 10, 1, self.Produit_amp_t_t0.T)

        M_Obj_decalage_guess = call_fc_decalage_temporelle_vecteur_t_de_t0_BEN(self.time, self.delay, self.Obj_t_guess, 1)
        
        M_Chp_t_guess, _ = np.meshgrid(self.Chp_t_guess[0, :], self.delay)
        self.M_Produit_t_t0_guess = M_Chp_t_guess * M_Obj_decalage_guess.T
        
        self.M_Produit_amp_omg_t0_guess = self._fft_(self.M_Produit_t_t0_guess)

        self.update_plot()
        self.emit_data()

    def update_plot(self):
        try:        
            amplitude = np.abs(self.Chp_t_guess[0, :]) / np.max(np.abs(self.Chp_t_guess[0, :]))
            switch = np.flip(np.abs(self.Obj_t_guess[0, :]) / np.max(np.abs(self.Obj_t_guess[0, :])))
            roi = np.where(amplitude >= 0.05)[0]
            self.M_Produit_t_t0_guess /= np.max(np.abs(self.M_Produit_t_t0_guess))
            self.Trace_amp_XP_omg_t0 /= np.max(np.abs(self.Trace_amp_XP_omg_t0))
            self.Produit_amp_t_t0 /= np.max(np.abs(self.Produit_amp_t_t0))
            self.M_Produit_amp_omg_t0_guess /= np.max(np.abs(self.M_Produit_amp_omg_t0_guess))
        except:             # if guess is empty
            self.divergence_occuring.emit()
            pass

        # Clear the figure and plot the data
        self.figure.clear()
        ax3 = self.figure.add_subplot(321)
        ax3.plot(self.time[roi], amplitude[roi], color = "red", linestyle = "-")
        ax3.set_title(f'Initial Guess $E(t)$')
        ax_phase = ax3.twinx()
        ax_phase.plot(self.time[roi], np.unwrap(np.angle(self.Chp_t_guess[0, roi])), color = "red", linestyle = "--")
        ax3.set_xlabel('Time, $t$ [fs]')
        ax3.set_ylabel('Amplitude [a.u.]')
        ax_phase.set_ylabel('Phase [rads]', rotation = 270, labelpad = 15)

        ax4 = self.figure.add_subplot(322)
        ax4.plot(self.time, switch, color = "blue", linestyle = "-")
        ax4.set_xlabel('Time, $t$ [fs]')
        ax4.set_ylabel('Amplitude [a.u.]')
        ax4.set_title('Initial Guess $S(t)$')

        ax1 = self.figure.add_subplot(323)
        img1 = ax1.imshow(np.abs(self.Trace_amp_XP_omg_t0).T, aspect = 'auto', cmap = self.cmap,
                   extent = [self.delay[0], self.delay[-1], self.wvl[-1], self.wvl[0]])
        ax1.set_ylabel('Wavelength, $\\lambda$ [nm]')
        ax1.set_xlabel('Time Delay, $\\tau$ [fs]')
        ax1.set_title('Processed $E(\\lambda, \\tau)$')
        cbar1 = self.figure.colorbar(img1, ax = ax1, orientation = 'vertical')
        cbar1.set_label('Amplitude [a.u.]', rotation = 270, labelpad = 15)

        ax2 = self.figure.add_subplot(324)
        img2 = ax2.imshow(np.abs(self.Produit_amp_t_t0).T, aspect = 'auto', cmap = self.cmap,
                   extent = [self.delay[0], self.delay[-1], self.time[0], self.time[-1]])
        ax2.set_ylabel('Time, $t$ [fs]')
        ax2.set_xlabel('Time Delay, $\\tau$ [fs]')
        ax2.set_title('Processed $E(t, \\tau)$')
        cbar2 = self.figure.colorbar(img2, ax = ax2, orientation = 'vertical')
        cbar2.set_label('Amplitude [a.u.]', rotation = 270, labelpad = 15)
        ax2.set_ylim([self.time[roi].min(), self.time[roi].max()])

        ax5 = self.figure.add_subplot(325)
        img5 = ax5.imshow(abs(self.M_Produit_amp_omg_t0_guess).T, aspect = 'auto', cmap = self.cmap,
                           extent = [self.delay[0], self.delay[-1], self.wvl[-1], self.wvl[0]])
        ax5.set_ylabel('Wavelength, $\\lambda$ [nm]')
        ax5.set_xlabel('Time Delay, $\\tau$ [fs]')
        ax5.set_title('Initial Guess $E(\\lambda, \\tau)$')
        cbar5 = self.figure.colorbar(img5, ax = ax5, orientation = 'vertical')
        cbar5.set_label('Amplitude [a.u.]', rotation = 270, labelpad = 15)

        ax6 = self.figure.add_subplot(326)
        img6 = ax6.imshow(abs(self.M_Produit_t_t0_guess).T, aspect='auto', cmap = self.cmap,
                           extent=[self.delay.min(), self.delay.max(), self.time.min(), self.time.max()])
        ax6.set_ylabel('Time, $t$ [fs]')
        ax6.set_xlabel('Time Delay, $\\tau$ [fs]')
        ax6.set_title('Initial Guess $E(t, \\tau)$')
        cbar6 = self.figure.colorbar(img6, ax = ax6, orientation = 'vertical')
        cbar6.set_label('Amplitude [a.u.]', rotation = 270, labelpad = 15)
        ax6.set_ylim([self.time[roi].min(), self.time[roi].max()])

        self.figure.tight_layout()
        self.canvas.draw()

    def update_iters(self):
        self.it_max = self.N_its_spinbox.value()

    @staticmethod
    def interp_and_pad(trace, delay, wvl, freq):
        # def prog02_step1(self):
        dt_delais_SUR_dt = 1
        Dt_SUR_Dt_delais = 4
        Multiplier_dt_delais_par = 1

        # interpolate over the time direction
        dt_delais = delay[9] - delay[8]
        dt_new = dt_delais * Multiplier_dt_delais_par
        interp_delay = np.arange(min(delay), max(delay) + dt_new, dt_new)
        scale_factor = len(interp_delay) / len(delay)
        trace /= np.max(trace)
        interp_trace = scipy.ndimage.zoom(trace, (1, scale_factor), order = 3)

        # interpolate over the frequency direction
        dt = (interp_delay[9] - interp_delay[8]) / dt_delais_SUR_dt  
        Dt_delais = interp_delay[-1] - interp_delay[0]
        Dt = Dt_delais * Dt_SUR_Dt_delais
        time = np.arange(-Dt / 2, Dt / 2 + dt, dt)
        N_t = len(time)
        Df1 = 1 / (dt * 1e-15)  # Convert fs to s for frequency calculation

        new_f = (((np.arange(N_t) / (N_t - 1)) - 0.5) * Df1) + np.mean(freq)
        scale_factor = len(new_f) / len(freq)
        interp_trace = scipy.ndimage.zoom(interp_trace, (scale_factor, 1), order=3)
        new_wvl = np.linspace(wvl[0], wvl[-1], len(new_f))

        # def prog02_step2(self):

        # padding with mean of edges
        Dt_delais_new = (time[-1] - time[0])*1.05
        n_ajout_pts_s2 = int(np.floor((Dt_delais_new / dt_delais - Dt_delais / dt_delais) / 2))
        M_spectre_haut = np.tile(np.mean(interp_trace[:, :10], axis = 1), (n_ajout_pts_s2, 1)).T
        M_spectre_bas = np.tile(np.mean(interp_trace[:, -10:], axis = 1), (n_ajout_pts_s2, 1)).T
        new_trace = np.hstack((M_spectre_haut, interp_trace, M_spectre_bas))
        new_delay = np.concatenate([
            interp_delay[0] - np.arange(n_ajout_pts_s2, 0, -1) * dt_delais,
            interp_delay,
            interp_delay[-1] + np.arange(1, n_ajout_pts_s2 + 1) * dt_delais
        ])

        # def prog03_step1(self):
        Coupure_fond = 0.0005
        new_trace = np.nan_to_num(new_trace, nan = 0)  
        new_trace = np.clip(np.abs(new_trace) - Coupure_fond * np.max(new_trace), 0, None)
        new_trace /= np.max(new_trace)

        return new_trace, new_delay, new_wvl, new_f, time
    
    def emit_data(self):
        data = {"delay": self.delay, "wavelength": self.wvl, "frequency": self.freq, "time": self.time,
                "processed_trace": self.XP_trace_int, "guess_matrix": self.M_Produit_t_t0_guess}
        self.guess_ready.emit(data)

if __name__ == '__main__':

    app = QApplication(sys.argv)

    # Example data
    file_path = r"C:\Users\lepar\OneDrive - INRS\Documents\INRS\MSc\thesis\PtyPy\Experimental_traces\20201201_FRost_3mm fused silica.txt"
    data = np.loadtxt(file_path)

    wavelength = data[2:, 0]
    delay = data[0, 2:]
    raw_trace = np.abs(data[2:, 2:])
    raw_trace /= np.max(raw_trace)
    raw_freq = (3e8) / (wavelength * 1e-6)       # not a constant step

    main_window = InitGuess()
    main_window.init_vars(raw_trace, delay, wavelength, raw_freq)
    main_window.show()

    sys.exit(app.exec_())
