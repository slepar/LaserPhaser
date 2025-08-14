try:
    from Software.Viewer.BasePlot import BasePlot
except:
    from BasePlot import BasePlot
    
class ComparisonPlot(BasePlot):
    def __init__(self, parent=None, cmap="binary"):
        super().__init__(parent, cmap)
        self.init()
        self.frog_data = None

    def init(self):
        # Add plots
        gs = self.fig.add_gridspec(2, 3)
        self.ax1, self.img1, self.cbar1 = self.init_2D_plot(gs[0, 0], 'Experimental Trace')
        self.ax2, self.img2, self.cbar2 = self.init_2D_plot(gs[0, 1], 'Recovered Trace 1')
        self.ax3, self.img3, self.cbar3 = self.init_2D_plot(gs[0, 2], 'Recovered Trace 2')

        self.t_amp_ax, self.t_pha_ax, self.t_amp_line, self.t_pha_line = self.init_1D_plot(gs[1, 0], "Predicted Pulses")
        self.w_amp_ax, self.w_pha_ax, self.w_amp_line, self.w_pha_line = self.init_1D_plot(gs[1, 1], "Predicted Spectras", xlabel="Wavelength, $\\lambda$ [nm]")
        self.switch_ax, self.switch_line = self.init_1D_plot(gs[1, 2], "Predicted Switches", phase=False)
        self.fig.tight_layout()
        self.canvas.draw_idle()

    def init_1D_plot(self, gs_position, title, xlabel = "Time, $t$ [fs]", labels = ["Data1", "Data2"], colors = ['red', 'blue'], phase = True):
        """ This overwrites the function in BasePlot so we can add more than one line. """
        ax_amp = self.fig.add_subplot(gs_position)
        ax_amp.set_title(title)
        ax_amp.set_xlabel(xlabel)
        ax_amp.set_ylabel("Amplitude [a.u.]")

        amp_lines = []
        for i in range(len(labels)):
            amp_lines.append(self.add_line(ax_amp, "", colors[i]))

        # Create a twin axis for phase
        if phase is True:
            ax_pha = ax_amp.twinx()
            ax_pha.set_ylabel("Phase, $\\phi$ [rad]", rotation=270, labelpad=15)

            pha_lines = []
            for i in range(len(labels)):
                pha_lines.append(self.add_line(ax_pha, "", colors[i], '--'))
                
            return ax_amp, ax_pha, amp_lines, pha_lines
        else:
            return ax_amp, amp_lines
    
    def prep_data(self, dataset):
        # interpolate and assure all parties are present
        if 'spectrum' not in dataset:
            dataset['spectrum'] = self._fft_(dataset['pulse'])
            dataset['wvl'] = self.interpolate(dataset['wvl'], len(dataset['spectrum']))
        if len(dataset['pulse']) != len(dataset['time']):
            dataset['pulse'] = self.interpolate(dataset['pulse'], max(256, 2*len(dataset['pulse'])))
            dataset['time'] = self.interpolate(dataset['time'], len(dataset['pulse']))
        if ('switch' in dataset) and (len(dataset['switch']) != len(dataset['delay'])):
            dataset['delay'] = self.interpolate(dataset['delay'], len(dataset['switch']))
        return dataset
    
    def calc_all_errs(self, data1, data2):
        # calculate errors on vectors
        pulse_err = self.calc_vector_err(data1['time'], data2['time'], data1['pulse'], data2['pulse'])
        spectrum_err = self.calc_vector_err(data1['wvl'], data2['wvl'], data1['spectrum'], data2['spectrum'])
        self.t_amp_ax.set_title(f"Pulse Error: {pulse_err:0.02f}%")
        self.w_amp_ax.set_title(f"Spectrum Error: {spectrum_err:0.02f}%")
        if 'switch' in data1 and 'switch' in data2:
            switch_err = self.calc_vector_err(data1['delay'], data2['delay'], data1['switch'], data2['switch'])
            self.switch_ax.set_title(f"Switch Error: {switch_err:0.02f}%")
        else:
            self.switch_ax.set_title(f"Switch Error: N/A %")

    def plot(self, raw_data, data1, data2 = None, labels = ["Data1", "Data2"]):        
        self.clear_plot()

        self.datasets = ((self.prep_data(data1), self.prep_data(data2)) if data2 is not None else [self.prep_data(data1)])
        if data2 is not None:
            self.calc_all_errs(*self.datasets)

        # update 1D plots
        for i in range(len(self.datasets)):
            self.update_1D_plot(self.datasets[i]['time'], self.datasets[i]['pulse'], self.t_amp_ax, self.t_amp_line[i], self.t_pha_ax, self.t_pha_line[i])
            self.update_1D_plot(self.datasets[i]['wvl'], self.datasets[i]['spectrum'], self.w_amp_ax, self.w_amp_line[i], self.w_pha_ax, self.w_pha_line[i])
            self.update_1D_plot(self.datasets[i]['delay'], self.datasets[i]['switch'], self.switch_ax, self.switch_line[i])
            self.update_fwhm(self.datasets[i]['time'], self.datasets[i]['pulse'], self.t_amp_ax, self.t_amp_line[i], labels[i])
            self.add_legend(f"{labels[i]}", self.w_amp_ax, self.w_amp_line[i], "upper left", 2)
            self.add_legend(f"{labels[i]}", self.switch_ax, self.switch_line[i], "upper left", 2)

        # update 2D plots
        self.update_2D_plot(raw_data, self.img1, self.cbar1)
        graph_data = [[self.ax2, self.img2, self.cbar2], [self.ax3, self.img3, self.cbar3]]
        for i in range(len(self.datasets)):
            temp_raw, temp_recovered = self.match_traces(raw_data, self.datasets[i])
            self.update_2D_plot(temp_recovered, *graph_data[i][-2:])
            graph_data[i][0].set_title(f"{labels[i]} Error: {self.calc_matrix_err(temp_raw['trace'], temp_recovered['trace']):0.02f}%")

    def add_FROG_data(self, frog_data, label = "FROG"):
        if len(self.t_amp_ax.lines) < 3:
            # first add lines to time and spectrum plots
            self.t_amp_line.append(self.add_line(self.t_amp_ax, "", "black"))
            self.t_pha_line.append(self.add_line(self.t_pha_ax, "", "black", style = "--"))
            self.w_amp_line.append(self.add_line(self.w_amp_ax, "", "black"))
            self.w_pha_line.append(self.add_line(self.w_pha_ax, "", "black", style = "--"))

        # prep dataset as before
        self.frog_data = self.prep_data(frog_data)

        # # # now plot data on new lines
        self.update_1D_plot(self.frog_data['time'], self.frog_data['pulse'], self.t_amp_ax, self.t_amp_line[-1], self.t_pha_ax, self.t_pha_line[-1])
        self.update_1D_plot(self.frog_data['wvl'], self.frog_data['spectrum'], self.w_amp_ax, self.w_amp_line[-1], self.w_pha_ax, self.w_pha_line[-1])
        self.update_fwhm(self.frog_data['time'], self.frog_data['pulse'], self.t_amp_ax, self.t_amp_line[-1], label)
        self.add_legend(f"{label}", self.w_amp_ax, self.w_amp_line[-1], "upper left", 3)

        if len(self.datasets) < 2:
            self.calc_all_errs(*self.datasets, self.frog_data)
        else:
            # change plot titles
            self.t_amp_ax.set_title("Predicted Pulses") 
            self.w_amp_ax.set_title("Predicted Spectras")

if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    import numpy as np
    import os
    
    app = QApplication(sys.argv)
    
    path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\DFN.pkl"
    path2 = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\PCP.pkl"
    path3 = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\FROG.pkl"

    # frog_path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20250115-SHG-FROG_ProbeAttWedge_1030nm_10fsStepSize"
    data = np.loadtxt(r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\raw_data.txt")
    raw_delay, raw_wvl, raw_trace = data[0, 1:], data[1:, 0], data[1:, 1:] / np.max(data[1:, 1:])

    window = ComparisonPlot()
    window.resize(1800, 800)
    
    data1 = window.load_pkl(path)
    data2 = window.load_pkl(path2)
    data3 = window.load_pkl(path3)

    raw_data = {'delay': raw_delay, 'wvl': raw_wvl, 'trace': raw_trace}

    window.plot(raw_data, data1, None, labels = ["DFN", "PCP"])

    window.add_FROG_data(data3)
    window.show()

    print(data1.keys())
    print(data2.keys())

    sys.exit(app.exec_())
