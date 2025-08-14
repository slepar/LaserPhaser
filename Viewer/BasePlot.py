from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QApplication
)
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar,
)
import numpy as np
import sys
import pickle
import scipy
import copy

class BasePlot(QWidget):
    """ This is the class that all plots will inherit from. It will hold all the common functions between them. """
    def __init__(self, parent = None, cmap = "binary"):
        super().__init__(parent)
        self.cmap = cmap
        self.ax_to_watch = None
        self.init_ui()

    def init_ui(self, figsize = (18, 6)):
        layout = QVBoxLayout()
        self.fig = Figure(figsize = figsize)
        self.canvas = FigureCanvas(self.fig)      
        self.toolbar = NavigationToolbar(self.canvas, self)
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    @staticmethod
    def add_line(ax, label="", color='red', style='-'):
        line, = ax.plot([], [], color=color, linestyle=style, label=label)
        return line
    
    def init_2D_plot(self, gs_position, title, xlabel = "Time Delay, $\\tau$ [fs]", ylabel = "Wavelength, $\\lambda$ [nm]"):
        ax = self.fig.add_subplot(gs_position)
        img = ax.imshow(np.zeros((128, 128)), aspect = 'auto', cmap = self.cmap)
        cbar = self.fig.colorbar(img, ax = ax, orientation = 'vertical') 
        cbar.set_label('Intensity [a.u.]', rotation = 270, labelpad = 15)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return ax, img, cbar

    def init_1D_plot(self, gs_position, title, xlabel = "Time, $t$ [fs]", color = "red", phase = True):
        ax_amp = self.fig.add_subplot(gs_position)
        ax_amp.set_title(title)
        ax_amp.set_xlabel(xlabel)
        ax_amp.set_ylabel("Amplitude [a.u.]")
        amp_line = self.add_line(ax_amp, '', color, '-')

        # Create a twin axis for phase
        if phase is True:
            ax_pha = ax_amp.twinx()
            ax_pha.set_ylabel("Phase, $\\phi$ [rad]", rotation=270, labelpad=15)
            pha_line = self.add_line(ax_pha, '', color, '--')
            return ax_amp, ax_pha, amp_line, pha_line
        else:
            return ax_amp, amp_line
    
    """ Calls for modulating plots. """
    @staticmethod
    def update_2D_plot(data, img, cbar, use_pro_trace = False):        
        trace = data['pro_trace'] if use_pro_trace is True else data['trace']
        img.set_data(trace)
        img.set_clim(vmin = np.min(trace), vmax = np.max(trace))
        cbar.update_normal(img)
        img.set_extent((data['delay'].min(), data['delay'].max(), data['wvl'].max(), data['wvl'].min()))             # flipped wvl is necessary

    def update_1D_plot(self, x_data, y_data, ax_amp, line_amp, ax_phase = None, line_phase = None):
        amp, phase = self.decompose_complex(y_data)
        roi = self.get_roi(amp, 0.025)
        line_amp.set_data(x_data, amp)
        ax_amp.set_xlim(np.min(x_data[roi]), np.max(x_data[roi]))
        ax_amp.set_ylim([-0.1, 1.75])

        if line_phase is not None:
            line_phase.set_data(x_data, phase)
            ax_phase.set_xlim(np.min(x_data[roi]), np.max(x_data[roi]))
            ax_phase.set_ylim([-4*np.pi, 2.5*np.pi])
        else:
            ax_amp.set_ylim([-0.1, 1.25])

        self.canvas.draw_idle()

    @staticmethod
    def add_legend(label, ax, line, loc = "upper right", ncols = 2):
        line.set_label(label)
        ax.legend(loc = loc, ncols = ncols) 

    def update_fwhm(self, time, pulse, ax, line, label = ""):
        fwhm = self.calc_fwhm(time, pulse)
        self.add_legend(f"{label}: {fwhm:0.02f} fs", ax, line, "upper left", 2)
        self.canvas.draw_idle()

    """ Misc functions used for various things. """
    def change_cmap(self, cmap_name):
        for ax in self.fig.get_axes():
            for img in ax.images:
                img.set_cmap(cmap_name)
        self.cmap = cmap_name
        self.canvas.draw_idle()

    def clear_plot(self):
        for ax in self.fig.get_axes():
            for img in ax.images:
                img.set_data(np.zeros_like(img.get_array()))
            for line in ax.lines:
                line.set_data([], [])
        self.canvas.draw_idle()

    @staticmethod
    def load_pkl(path):
        with open(path, 'rb') as f:
            data = pickle.load(f)
        return data

    @staticmethod
    def calc_energy(roi, window):
        roi_power = np.sum(roi ** 2)
        total_power = np.sum(window ** 2)
        return roi_power / total_power

    @staticmethod
    def decompose_complex(c_vec):
        return np.abs(c_vec) / np.max(np.abs(c_vec)), np.angle(c_vec)

    @staticmethod
    def get_roi(data, threshold1 = 0.1, threshold2 = None):
        if threshold2 is None:
            return np.where(data >= threshold1)[0]
        else:
            return np.where((data >= threshold1) & (data <= threshold2))[0]

    def calc_fwhm(self, time, pulse):
        if len(pulse) <= 128:
            time = self.interpolate(time, 2*len(time))
            pulse = self.interpolate(pulse, len(time))
        amplitude, _ = self.decompose_complex(pulse)
        ind = self.get_roi(amplitude, np.max(amplitude)/2)
        return np.abs(time[ind[0]] - time[ind[-1]])

    def polyfit_phase(self, wvl, spectrum):
        amp, phase_w = self.decompose_complex(spectrum)
        phase_w = np.unwrap(phase_w)
        roi = self.get_roi(amp, 0.01)
        coefficients = np.polyfit(wvl[roi], phase_w[roi], 11)
        fitted_phase_w = np.polyval(coefficients, wvl[roi])

        # calculate derivatives 
        gdd = 2 * coefficients[-3]  # Coefficient of omega^2 * 2!
        tod = 6 * coefficients[-4]  # Coefficient of omega^3 * 3!
        fitted_phase_w -= fitted_phase_w[len(fitted_phase_w)//2]
        fit_error = (np.sqrt(np.mean((phase_w[roi] - fitted_phase_w)**2))/(np.max(phase_w[roi]) - np.max(fitted_phase_w)))*100
        return wvl[roi], fitted_phase_w, gdd, tod, fit_error

    @staticmethod
    def interpolate(data, new_shape):
        if isinstance(new_shape, int):
            return scipy.ndimage.zoom(data, new_shape / len(data), order = 3)
        else:
            zoom_factors = (new_shape[0] / data.shape[0], new_shape[1] / data.shape[1])
            return scipy.ndimage.zoom(data, zoom_factors, order = 3)
        
    @staticmethod
    def normalize(data):
        return np.abs(data) / np.max(np.abs(data))

    def add_phase_fit(self, wvl, spectrum, ax):
        self.new_line = self.add_line(ax, "", 'blue', '--') if len(ax.lines) == 1 else self.new_line
        fitted_wvl, fitted_phase_w, gdd, tod, fit_error = self.polyfit_phase(wvl, spectrum)
        self.new_line.set_data(fitted_wvl, fitted_phase_w)
        self.add_legend(f'Err = {fit_error:0.02f}%\nGDD = {gdd:0.02f}, TOD = {tod:0.02f}', ax, self.new_line, loc = "upper left", ncols = 1)
        self.canvas.draw_idle()

    def on_crop_change(self, event):
        if self.ax_to_watch is None:
            pass
        
        xlim = self.ax_to_watch.get_xlim()
        if (self.previous_lim != xlim):
            self.previous_lim = xlim 
            line_to_watch = self.ax_to_watch.lines[0]
            x_data = line_to_watch.get_xdata()
            y_data = line_to_watch.get_ydata()
            ylim = self.ax_to_watch.get_ylim()
            # try:
            visible_indices = (x_data >= xlim[0]) & (x_data <= xlim[1]) & (y_data >= ylim[0]) & (y_data <= ylim[1])
            visible_x = x_data[visible_indices]
            visible_y = y_data[visible_indices]
                
            # visible_x = x_data[self.get_roi(x_data, *xlim)]
            # visible_y = y_data[self.get_roi(y_data, *ylim)]

            energy_present = self.calc_energy(visible_y, y_data)
            fwhm = self.calc_fwhm(visible_x, visible_y)
            self.add_legend(f"FWHM: {fwhm:0.02f} fs\nE = {energy_present*100:0.02f}%", 
                            self.ax_to_watch, line_to_watch, "upper left")
            self.canvas.draw_idle()
            # except:
            #     pass

    """ For calculating errors. """
    def check_dims(self, x, y, z):
        if x.shape != z.shape[1]:
            x = self.interpolate(x, z.shape[1])
        if y.shape != z.shape[0]:
            y = self.interpolate(y, z.shape[0])
        return x, y, z
    
    def match_traces(self, data1, data2, string1 = 'trace', string2 = 'trace'):
        # avoids rewriting data in the og dict
        temp_data1, temp_data2 = copy.deepcopy(data1), copy.deepcopy(data2)
        
        # Extract and process data
        y1 = temp_data1['sh_wvl'] if 'sh_wvl' in temp_data1 else temp_data1['wvl']
        y2 = temp_data2['sh_wvl'] if 'sh_wvl' in temp_data2 else temp_data2['wvl']
        x1, z1 = temp_data1['delay'], self.convert_from_tensor(temp_data1[string1])
        x2, z2 = temp_data2['delay'], self.convert_from_tensor(temp_data2[string2])

        # interpolate to same length first
        x1, y1, z1 = self.check_dims(x1, y1, z1)
        x2, y2, z2 = self.check_dims(x2, y2, z2)

        x1, x2, y1, y2, z1, z2 = self.clip_matrix_regions(x1, x2, y1, y2, z1, z2)
        new_shape = [min(z1.shape[0], z2.shape[0]), min(z1.shape[1], z2.shape[1])]

        x1, x2 = self.interpolate(x1, new_shape[1]), self.interpolate(x2, new_shape[1])
        y1, y2 = self.interpolate(y1, new_shape[0]), self.interpolate(y2, new_shape[0])
        z1, z2 = self.normalize(self.interpolate(z1, new_shape)), self.normalize(self.interpolate(z2, new_shape))

        # Update the copied datasets
        temp_data1.update({'delay': x1, 'wvl': y1, string1: z1})
        temp_data2.update({'delay': x2, 'wvl': y2, string2: z2})
        return temp_data1, temp_data2

    def clip_matrix_regions(self, x1, x2, y1, y2, z1, z2):
        # determine region for axis 1
        lower_lim, upper_lim = max(np.min(x1), np.min(x2)), min(np.max(x1), np.max(x2))
        roi_x1, roi_x2 = self.get_roi(x1, lower_lim, upper_lim), self.get_roi(x2, lower_lim, upper_lim)

        # repeat for axis 2
        lower_lim, upper_lim = max(np.min(y1), np.min(y2)), min(np.max(y1), np.max(y2))
        roi_y1, roi_y2 = self.get_roi(y1, lower_lim, upper_lim), self.get_roi(y2, lower_lim, upper_lim)
        return x1[roi_x1], x2[roi_x2], y1[roi_y1], y2[roi_y2], z1[roi_y1][:, roi_x1], z2[roi_y2][:, roi_x2]

    def clip_vector_regions(self, x1, x2, y1, y2):
        lower_lim, upper_lim = max(np.min(x1), np.min(x2)), min(np.max(x1), np.max(x2))
        roi_x1, roi_x2 = self.get_roi(x1, lower_lim, upper_lim), self.get_roi(x2, lower_lim, upper_lim)
        return x1[roi_x1], x2[roi_x2], y1[roi_x1], y2[roi_x2]
    
    @staticmethod
    def calc_matrix_err(xp, pred):        
        numerator = np.sum(xp * pred)
        denominator = np.sqrt(np.sum(xp**2) * np.sum(pred**2))
        return np.sqrt(1 - (numerator / denominator)**2)*100

    @staticmethod
    def convert_from_tensor(tensor):
        # is either an array or a tensor (from AI prediction)
        return tensor if isinstance(tensor, np.ndarray) else tensor.numpy()

    def match_vectors(self, x1, x2, y1, y2):
        y1, y2 = self.convert_from_tensor(y1), self.convert_from_tensor(y2)           
        x1, x2, y1, y2 = self.clip_vector_regions(x1, x2, y1, y2)
        new_length = max(len(x1), len(x2))

        # Interpolate both vectors to the new length
        x1, x2 = self.interpolate(x1, new_length), self.interpolate(x2, new_length)
        y1, y2 = self.normalize(self.interpolate(y1, new_length)), self.normalize(self.interpolate(y2, new_length))
        return x1, x2, y1, y2
    
    def calc_vector_err(self, x1, x2, y1, y2):
        x1, x2, y1, y2 = self.match_vectors(x1, x2, y1, y2)
        j_err = self.calc_mse(y1, y2) 
        return np.abs(j_err)*100
    
    @staticmethod
    def calc_mse(vec1, vec2):
        return np.mean((vec1 - vec2) ** 2)

    @staticmethod
    def _fft_(pulse):
        return np.fft.fftshift(np.fft.fft(np.fft.ifftshift(pulse)))

    @staticmethod
    def _ifft_(spectrum):
        return np.fft.ifftshift(np.fft.ifft(np.fft.fftshift(spectrum)))
    
if __name__ == "__main__":
    import os

    app = QApplication(sys.argv)
    
    path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\DFN.pkl"
    path2 = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\PCP.pkl"
    path3 = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\FROG.pkl"

    # frog_path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20250115-SHG-FROG_ProbeAttWedge_1030nm_10fsStepSize"
    data = np.loadtxt(r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\raw_data.txt")
    raw_delay, raw_wvl, raw_trace = data[0, 1:], data[1:, 0], data[1:, 1:] / np.max(data[1:, 1:])

    window = BasePlot()
    window.resize(1800, 800)
    
    data1 = window.load_pkl(path)
    data2 = window.load_pkl(path2)
    data3 = window.load_pkl(path3)


    window.show()

    print(data2.keys())

    sys.exit(app.exec_())
