try:
    from Software.Viewer.BasePlot import BasePlot
except:
    from BasePlot import BasePlot

class PredictionPlot(BasePlot):
    def __init__(self, parent=None, cmap="binary"):
        super().__init__(parent, cmap)
        self.previous_lim = None
        self.ax_to_watch = None
        self.init()

    def init(self):
        # Add plots
        gs = self.fig.add_gridspec(2, 3)
        self.ax1, self.img1, self.cbar1 = self.init_2D_plot(gs[0, 0], 'Experimental Trace')
        self.ax2, self.img2, self.cbar2 = self.init_2D_plot(gs[0, 1], 'Processed Trace')
        self.ax3, self.img3, self.cbar3 = self.init_2D_plot(gs[0, 2], 'Recovered Trace')

        self.t_amp_ax, self.t_pha_ax, self.t_amp_line, self.t_pha_line = self.init_1D_plot(gs[1, 0], "Predicted Pulse")
        self.w_amp_ax, self.w_pha_ax, self.w_amp_line, self.w_pha_line = self.init_1D_plot(gs[1, 1], "Predicted Spectrum", xlabel="Wavelength, $\\lambda$ [nm]")
        self.switch_ax, self.switch_line = self.init_1D_plot(gs[1, 2], "Predicted Switch", phase=False)
        self.fig.tight_layout()
        self.canvas.draw_idle()

        # Assign axis to watch for zoom adjustments
        self.previous_lim = self.t_amp_ax.get_xlim()
        self.ax_to_watch = self.t_amp_ax
        self.canvas.mpl_connect("draw_event", self.on_crop_change)

    def prep_data(self, data):
        if 'spectrum' not in data:
            data['spectrum'] = self._fft_(data['pulse'])
        if len(data['spectrum']) != len(data['wvl']):
            data['spectrum'] = self.interpolate(data['spectrum'], len(data['wvl']))
        if 'switch' in data and len(data['switch']) != data['trace'].shape[1]:
            data['switch'] = self.interpolate(data['switch'], data['trace'].shape[1])
        if len(data['pulse']) < 256:
            data['pulse'] = self.interpolate(data['pulse'], 512)
            data['time'] = self.interpolate(data['time'], 512)
        return data
    
    def plot(self, raw_data, data):
        self.switch_plot_visible('switch' in data)

        if data['trace'].shape != raw_data['trace'].shape:
            data['pro_trace'] = self.interpolate(data['pro_trace'], raw_data['trace'].shape)
            data['trace'] = self.interpolate(data['trace'], raw_data['trace'].shape)
            data['wvl'] = self.interpolate(data['wvl'], data['trace'].shape[0])
            data['delay'] = self.interpolate(data['delay'], data['trace'].shape[1])
        self.data = self.prep_data(data)

        # plot 1D graphs
        self.update_1D_plot(self.data['time'], self.data['pulse'], self.t_amp_ax, self.t_amp_line, self.t_pha_ax, self.t_pha_line)
        self.update_1D_plot(self.data['wvl'], self.data['spectrum'], self.w_amp_ax, self.w_amp_line, self.w_pha_ax, self.w_pha_line)
        self.add_phase_fit(self.data['wvl'], self.data['spectrum'], self.w_pha_ax)
        if 'switch' in self.data:
            self.update_1D_plot(self.data['delay'], self.data['switch'], self.switch_ax, self.switch_line)

        # first we check processing vs raw
        temp_raw, temp_processed = self.match_traces(raw_data, self.data, string1 = 'trace', string2 = 'pro_trace')
        self.ax2.set_title(f"Processing Error: {self.calc_matrix_err(temp_raw['trace'], temp_processed['pro_trace']):0.02f}%")
        self.update_2D_plot(temp_raw, self.img1, self.cbar1)
        self.update_2D_plot(temp_processed, self.img2, self.cbar2, True)

        # next we compare recovered to processed
        temp_processed, temp_recovered = self.match_traces(temp_processed, self.data, string1 = 'pro_trace', string2 = 'trace')
        self.ax3.set_title(f"Recovery Error: {self.calc_matrix_err(temp_processed['pro_trace'], temp_recovered['trace']):0.02f}%")
        self.update_2D_plot(temp_recovered, self.img3, self.cbar3)

    def switch_plot_visible(self, show = True):
        self.switch_ax.set_visible(show)  # Hide the axis
        self.switch_line.set_visible(show)  # Hide the line (if applicable)
        self.fig.tight_layout()  # Adjust layout to account for hidden plot
        self.canvas.draw_idle()  # Redraw the canvas

if __name__ == "__main__":
    import os
    import sys
    from PyQt5.QtWidgets import QApplication
    import numpy as np
    
    app = QApplication(sys.argv)
    
    # path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20250115-SHG-FROG_ProbeAttWedge_1030nm_HigherRes_5fsStepSize_Meas3\FROG.pkl"
    # data = np.loadtxt(r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20250115-SHG-FROG_ProbeAttWedge_1030nm_HigherRes_5fsStepSize_Meas3\raw_data.txt")

    path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\PCP.pkl"
    data = np.loadtxt(r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\raw_data.txt")

    raw_delay, raw_wvl, raw_trace = data[0, 1:], data[1:, 0], data[1:, 1:] / np.max(data[1:, 1:])

    window = PredictionPlot()
    window.resize(1400, 800)
    data = window.load_pkl(path)

    raw_data = {'delay': raw_delay, 'wvl': raw_wvl, 'trace': raw_trace}
    window.plot(raw_data, data)
    window.show()
    

    print(data.keys())

    sys.exit(app.exec_())


