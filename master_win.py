from PyQt5.QtWidgets import (QWidget, QVBoxLayout, QMainWindow, QDialog, QStackedWidget, 
                             QAction, QFileDialog, QApplication, QActionGroup, QLabel)

from Software.Processor import *
from Software.Viewer import *
from Software.Viewer.CMap.cmap_options import *
from Software.Predictor import *

from matplotlib.colors import LinearSegmentedColormap
import sys
import json
import scipy
import shutil
import pickle 
import os 
import numpy as np

class MasterWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.spectrum_op = "slice"
        self.cmap_name = "binary"

        self.cmap_options = {}
        self.parent_save = None
        self.raw_data, self.pty_data, self.ai_data = None, None, None
        self.init_ui()

    """ Initialization UI functions. """
    def init_ui(self):       
        self.setWindowTitle("Phase Recovery")
        self.resize(800, 600)
        self.showMaximized()
        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QVBoxLayout(main_widget)

        # label for trace name
        self.trace_label = QLabel("No data loaded")
        self.trace_label.setStyleSheet("color: gray; font: 16px")
        main_layout.addWidget(self.trace_label)

        # plot widget
        self.plot_widget = QStackedWidget(self)
        self.trace_plot = TracePlot(self, self.cmap_name)
        self.pred_plot = PredictionPlot(self, self.cmap_name)
        self.comp_plot = ComparisonPlot(self, self.cmap_name)
        self.plot_widget.addWidget(self.trace_plot)
        self.plot_widget.addWidget(self.pred_plot)
        self.plot_widget.addWidget(self.comp_plot)
        main_layout.addWidget(self.plot_widget)

        # menu bar with load button 
        self.menu_bar = self.menuBar()
        file_menu = self.menu_bar.addMenu("Load Data Files")
        self.add_action(file_menu, "Load XP Data", self.select_datafile)
        self.add_action(file_menu, "View Saved Data", self.open_dataviewer)   
        self.menu_bar.addMenu(file_menu)

        # predictions menu
        pred_menu = self.menu_bar.addMenu("Predictions")
        self.add_action(pred_menu, "FROG", self.init_FROG_pred)
        dfn_menu = pred_menu.addMenu("DFROStNET")
        self.add_action(dfn_menu, "Load Model", self.load_ai_model)
        self.add_action(dfn_menu, "Get Prediction", self.init_DFROStNET_pred)
        spectrum_submenu = dfn_menu.addMenu("Select Spectrum")
        spectrum_group = self.add_action_group()
        actions = {"mean": "Use mean", "slice": "Use slice", "upload": "Upload spectrum"}
        for name, label in actions.items():
            self.add_exclusive_action(label, spectrum_group, spectrum_submenu, self.handle_spectrum_selection, name)
        spectrum_options = {name: action for name, action in zip(actions.keys(), spectrum_group.actions())}
        spectrum_options.get(self.spectrum_op, list(actions.keys())[0]).setChecked(True)
        pred_menu.addMenu(dfn_menu)
        self.add_action(pred_menu, "PtyChoPy", self.init_PtyChoPy_pred)
        self.menu_bar.addMenu(pred_menu)

        # analysis menu
        anys_menu = self.menu_bar.addMenu("Analysis")
        self.add_action(anys_menu, "STFT", self.open_stft)
        self.add_action(anys_menu, "Compare Predictions", self.compare_live_predictions)
        self.menu_bar.addMenu(anys_menu)

        # color map menu
        cmap_submenu = self.menuBar().addMenu("CMap")
        cmap_group = self.add_action_group()
        for key, label in cmap_dict.items():
            self.cmap_options[key.lower()] = self.add_exclusive_action(label, cmap_group, cmap_submenu, self.handle_cmap_change)
        self.cmap_options.get(self.cmap_name, self.cmap_options[self.cmap_name]).setChecked(True)

    def add_action_group(self):
        group = QActionGroup(self)
        group.setExclusive(True) 
        return group
    
    def add_exclusive_action(self, label, group, submenu, function, name = None):
        action = QAction(label, self, checkable = True)
        group.addAction(action)
        submenu.addAction(action)
        name = label if name is None else name
        action.toggled.connect(lambda: function(name))
        return action
    
    def add_action(self, menu, action_name, function):
        dummy = QAction(action_name, self)
        dummy.triggered.connect(function)
        menu.addAction(dummy)

    """ Loading / Saving Functions """
    def select_datafile(self):
        self.file_path, _ = QFileDialog.getOpenFileName(self, "Open File", os.path.join(os.getcwd(), "XP_FROSt_traces"), 
                                                        "All Files (*)", 
                                                        options = QFileDialog.Options())

        if os.path.isfile(self.file_path):
            self.load_from_xp_path(self.file_path)

    def load_raw_data(self, path):
        data = np.loadtxt(path)
        raw_delay, raw_wvl, raw_trace = data[0, 1:], data[1:, 0], data[1:, 1:] / np.max(np.abs(data[1:, 1:]))  
        raw_time, raw_freq = self.gen_f_t_vector(np.flip(raw_wvl) / 1e3)
        raw_data = {'delay': raw_delay, 'wvl': raw_wvl, 'trace': raw_trace, 'time': raw_time, 'freq': raw_freq}
        return raw_data
    
    def load_from_xp_path(self, path):
        # Load data
        self.file_path = path
        self.raw_data = self.load_raw_data(self.file_path)
        
        # show in ui that trace is laoded. 
        self.switch_plot(0)
        self.trace_name = os.path.basename(path)
        self.trace_label.setText(f"'{self.trace_name}' loaded.")         
    
        # prep folders for outputs
        self.parent_save = os.path.join(os.getcwd(), "Recovered_Data", f"{self.trace_name}")
        os.makedirs(self.parent_save, exist_ok=True)
        destination_path = os.path.join(self.parent_save, "raw_data.txt")
        if not os.path.exists(destination_path):
            shutil.copy(self.file_path, destination_path)

    def upload_spectrum(self):
        # self.spectrum_file_path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\dummy_spectrum.txt"
        self.spectrum_file_path, _ = QFileDialog.getOpenFileName(self, "Select a Text File", "", "Text Files (*.txt);;All Files (*)")
        if self.spectrum_file_path:
            data = np.loadtxt(self.spectrum_file_path)  # Using numpy to read the file
            wvl = data[:, 0]
            est_spectrum = data[:, 1]
        return est_spectrum
    
    def load_ai_model(self):
        options = QFileDialog.Options()
        model_dir = os.path.join(os.getcwd(), "Software", "Predictor", "DFROStNET", "Models")
        self.model_path, _ = QFileDialog.getOpenFileName(self, "Open File", model_dir, "Weights Files (*.weights.h5)",  options=options)

    def dump_pkl_txt(self, file_name, data):
        self.save_pkl(data, file_name)          # save to pkl for simplicity
    
        # also save to txt in subfolder so we can see it
        file_name = file_name.split(".pkl")[0]
        os.makedirs(file_name, exist_ok = True) 
        for key, value in data.items():
            arr = np.array(value)
            txt_file_path = os.path.join(file_name, f"{key}.txt")
            np.savetxt(txt_file_path, arr, fmt = '%.6e')  # scientific notation, 6 decimals

    """ Prediction functions. """
    def init_DFROStNET_pred(self):
        try:
            self.dfn = DFROStNET(model_path = self.model_path, callback = self.handle_pred_data)
        except:
            self.dfn = DFROStNET(callback = self.handle_pred_data)
        self.handle_spectrum_selection(self.spectrum_op)
        self.dfn.run(self.raw_data, self.est_spectrum)

    def init_PtyChoPy_pred(self):
        self.pcp_win = PtyChoPy(self.cmap)
        self.pcp_win.iterations_complete.connect(self.handle_pred_data)
        self.pcp_win.init_vars(self.raw_data)

    def init_FROG_pred(self):
        # subfunction does processing and launches frog app
        self.frog_win = FROG(raw_data_path = self.file_path)
        self.frog_win.recovery_complete.connect(self.handle_pred_data)

    def handle_pred_data(self, data, source):
        self.switch_plot(1)
        self.pred_plot.plot(self.raw_data, data)

        if source == "PtyChoPy":
            self.pty_data = data           
            self.dump_pkl_txt(os.path.join(self.parent_save, 'PCP.pkl'), self.pty_data)

        elif source == "DFROStNET":
            self.ai_data = data           
            self.dump_pkl_txt(os.path.join(self.parent_save, 'DFN.pkl'), self.ai_data)

        elif source == "FROG":
            self.frog_data = data
            self.dump_pkl_txt(os.path.join(self.parent_save, 'FROG.pkl'), self.frog_data)

        QApplication.processEvents()    # Process pending UI events

    def handle_spectrum_selection(self, option):
        self.spectrum_op = option
        if self.raw_data is not None:
            trace = self.raw_data['trace']
            if option == "mean":
                est_spectrum = np.sqrt(np.mean(np.abs(trace), axis = 1))
            elif option == "slice":
                slice1 = np.mean(trace[:, :10], axis = 1)
                slice2 = np.mean(trace[:, -10:], axis = 1)
                est_spectrum = np.sqrt(np.abs(slice1 - slice2)) / np.max(np.abs(slice1 - slice2))
            elif option == "upload":
                est_spectrum = self.upload_spectrum()
            est_spectrum = scipy.ndimage.gaussian_filter1d(est_spectrum, sigma = 0.5)                  
            self.est_spectrum = est_spectrum / np.max(np.abs(est_spectrum))

    def compare_live_predictions(self):
        try:
            self.ai_data = self.load_pkl(os.path.join(self.parent_save, 'DFN.pkl')) if self.ai_data is None else self.ai_data
            self.pty_data = self.load_pkl(os.path.join(self.parent_save, 'PCP.pkl')) if self.pty_data is None else self.pty_data
            self.switch_plot(2)
        except:
            print("Make prediction prior to viewing.")

    """ CMap Stuff. """
    def handle_cmap_change(self, cmap_name):
        """Change the colormap of all plots."""
        self.cmap_name = cmap_name
        if self.cmap_name == "custom":                  # making a custom color now
            cmap_popup = CMapPopup()
            cmap_popup.show()
            if cmap_popup.exec_() == QDialog.Accepted:  
                self.cmap = cmap_popup.generate_cmap()           
        elif self.cmap_name.endswith(".json"):          # loading a previous custom cmap
            self.cmap = self.load_colormap(self.cmap_name)
        else:                                           # using a standard cmap
            self.cmap = self.cmap_name 

        for plot in [self.trace_plot, self.pred_plot, self.comp_plot]:
            plot.change_cmap(self.cmap)

    @staticmethod
    def load_colormap(cmap_name):
        path = os.path.join(custom_cmap_dir, cmap_name)
        with open(path, 'r') as f:
            cmap_data = json.load(f)
        
        # Reconstruct the colormap
        colors = [(entry["r"], entry["g"], entry["b"]) for entry in cmap_data]
        return LinearSegmentedColormap.from_list("dummy", colors)
    
    """ Misc. """
    def open_stft(self):
        if (self.pty_data or self.ai_data) is not None:
            pulse, time = (self.pty_data['pulse'], self.pty_data['time']) if self.pty_data is not None else (self.ai_data['pulse'], self.ai_data['time'])
            self.stft_widget = STFT(pulse, time, cmap = self.cmap_name)
            self.stft_widget.show()

    def switch_plot(self, new_index):
        if new_index == 0:              # plot raw data 
            self.trace_plot.clear_plot()
            self.plot_widget.setCurrentIndex(new_index)
            self.trace_plot.plot(self.raw_data)

        if new_index == 1:              # plot raw data in predplot and switch canvas
            self.pred_plot.clear_plot()
            self.plot_widget.setCurrentIndex(new_index)
            # plots when data is received from class

        elif new_index == 2:        # plot comparison
            self.comp_plot.clear_plot()
            self.plot_widget.setCurrentIndex(new_index)
            self.comp_plot.plot(self.raw_data, self.ai_data, self.pty_data, ["DFROStNET", "PtyChoPy"])      

    def open_dataviewer(self):
        self.viewer = DataViewer()

    @staticmethod
    def gen_f_t_vector(wavelength):
        freq = (3e8 / wavelength * 1e-9)        # wvl should be in nm from the spectrometer
        sorted_indices = np.argsort(freq)
        temp_freq = freq[sorted_indices]
        freq = np.linspace(temp_freq[0], temp_freq[-1], len(temp_freq))
        df = (freq[9] - freq[8])
        Dt = 1 / (len(freq)*df) 
        time = np.arange(-len(freq)//2, len(freq)//2) * Dt
    
        # freq = (3e8 / wavelength * 1e-9)                        # wvl in nm from the spectrometer
        # freq = np.linspace(freq[0], freq[-1], len(freq))        # constant step
        # df = (freq[9] - freq[8])

        # dt_delay = delay[9] - delay[8]
        # Dt_delay = (delay[-1] - delay[0])*4                     # *4 from og MATLAB code

        # # time = np.linspace(-Dt_delay/2, Dt_delay/2, int(dt_delay))
        # time = np.arange(start = -Dt_delay//2, stop = Dt_delay//2, step = int(dt_delay))

        # Df1 = (1 / dt_delay)*1e-15
        # freq = (((np.arange(len(time)) / (len(time) - 1)) - 0.5) * Df1) #+ f0

        # # Dt = 1 / (len(freq)*df) 
        # # time = np.arange(-len(freq)//2, len(freq)//2) * Dt
        return time, freq
        
    @staticmethod
    def save_pkl(data, path):
        with open(path, 'wb') as f:
            pickle.dump(data, f)

    @staticmethod
    def load_pkl(path):
        with open(path, 'rb') as f:
            data = pickle.load(f)
        return data

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MasterWindow()
    window.show()
    sys.exit(app.exec_())
