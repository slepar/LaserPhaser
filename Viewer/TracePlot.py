try:
    from Software.Viewer.BasePlot import BasePlot
except:
    from BasePlot import BasePlot

class TracePlot(BasePlot):
    """ This widget plots a single trace (FROG or FROSt). """
    def __init__(self, parent=None, cmap="binary"):
        super().__init__(parent, cmap)
        self.init()

    def init(self, label="Experimental Trace"):
        # Set up the plot using a single grid space
        gs = self.fig.add_gridspec(1, 1)  # Single subplot layout
        self.ax, self.img, self.cbar = self.init_2D_plot(gs[0], label)
        self.canvas.draw_idle()

    def plot(self, data):
        # Update the plot with the provided data
        self.update_2D_plot(data, self.img, self.cbar)
        self.canvas.draw_idle()

if __name__ == "__main__":
    import sys
    from PyQt5.QtWidgets import QApplication
    import numpy as np

    app = QApplication(sys.argv)
    
    path2 = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20250115-SHG-FROG_ProbeAttWedge_1030nm_HigherRes_5fsStepSize_Meas3\raw_data.txt"

    data = np.loadtxt(path2)
    raw_delay, raw_wvl, raw_trace = data[0, 1:], data[1:, 0], data[1:, 1:] / np.max(np.abs(data[1:, 1:]))  
    raw_data = {'delay': raw_delay, 'wvl': raw_wvl, 'trace': raw_trace} #, 'time': raw_time, 'freq': raw_freq}

    # path2 = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt\PCP.pkl"

    window = TracePlot()
    window.resize(800, 600)
    # data2 = window.load_pkl(path2)
    # print(data2.keys())

    window.plot(raw_data)
    window.show()
    sys.exit(app.exec_())
