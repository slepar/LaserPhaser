import subprocess
import os 
import shutil
import pyautogui as ag
import time as clock
import ctypes
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
import threading
import numpy as np
import pickle 
from scipy.interpolate import interp1d
import math
from PyQt5.QtWidgets import QApplication, QWidget
import sys
from PyQt5.QtCore import pyqtSignal
import pyperclip
import matplotlib.pyplot as plt


class FileCreationHandler(FileSystemEventHandler):
    def __init__(self, extension = ".dat", callback = None):
        self.extension = extension
        self.callback = callback 
        self.timer = None  # Timer for debounce
        self.debounce_interval = 1.0  # Seconds to wait

    def on_created(self, event):
        if event.src_path.endswith(self.extension):
            self.reset_timer()

    def on_modified(self, event):
        if event.src_path.endswith(self.extension):
            self.reset_timer()

    def reset_timer(self):
        if self.timer:
            self.timer.cancel()  # Cancel any existing timer
        self.timer = threading.Timer(self.debounce_interval, self.callback)
        self.timer.start()
    
""" Function that handles the formatting and interfaces with the FROG software. 
    It must be a QWidget otherwise the FROG software shuts down. """
class FROG(QWidget):
    recovery_complete = pyqtSignal(object, str)

    def __init__(self, raw_data_path):
        super().__init__()
        self.frog_app_path = os.path.join(os.getcwd(), "Software\Predictor\FROG\FROG_soft\Frog3.exe")
        self.raw_file_path = raw_data_path

        self.prep_save_dir()
        self.processed_file_path = self.reformat_FROG(self.local_raw_path)
        self.start_monitoring_for_files(self.frog_save_subdir)

        # # launch software with mouse lock: ~ 4-5 seconds to execute
        self.launch_frog()

    def launch_frog(self):
        ctypes.windll.user32.BlockInput(True)
        subprocess.Popen(self.frog_app_path, shell = True)
        # clock.sleep(0.5)
        # self.close_popup()
        # clock.sleep(0.5)

        # load data
        ag.hotkey("ctrl", "f")          # opens the load data window
        pyperclip.copy(self.processed_file_path)            # copies reshaped file path to clipboard for user to paste
        ctypes.windll.user32.BlockInput(False)

    @staticmethod
    def close_popup():
        # auto detects and closes the popup
        popup_img_path = os.path.join(os.getcwd(), r"Software\Predictor\FROG\FROG_soft\popup_img.png") 
        try:
            location = ag.locateOnScreen(popup_img_path, confidence = 0.5)
            click_on_x = (location.left + location.width - 5, location.top + 15)
            ag.click(click_on_x)           # closes the invalid file popup
        except:
            pass    

    def prep_save_dir(self):
        parent_savedir = os.path.join(os.getcwd(), "Recovered_Data")
        self.file_name = os.path.basename(self.raw_file_path)

        # make a subfolder for this trace
        self.frog_save_subdir = os.path.join(parent_savedir, self.file_name, "Raw_FROG_Data")
        os.makedirs(self.frog_save_subdir, exist_ok = True)
        self.local_raw_path = os.path.join(self.frog_save_subdir, self.file_name)
        if not os.path.exists(self.local_raw_path):
            shutil.copy(self.raw_file_path, self.local_raw_path)           # copies raw data to subfolder

    @staticmethod
    def reformat_FROG(file_path):
        # extract raw data
        Image = np.loadtxt(file_path)
        delays = Image[0, 1:] 
        nb_delay_pts = np.size(delays)
        delay_inc = delays[1] - delays[0]
        wavelengths = Image[1:, 0]                             #900 to 2500nm
        spectrogram = Image[1:, 1:]  

        # make new wvl
        wvlmin = math.ceil(wavelengths[0]) 
        wvlmax = math.floor(wavelengths[-1])
        wvl_inc = round(wavelengths[9] - wavelengths[8], 1)     #nm 
        new_wvl = np.arange(wvlmin, wvlmax, wvl_inc)            #equally spaced wavelengths for interpolation
        nb_wvl_pts = np.size(new_wvl)
        wvl_central = new_wvl[nb_wvl_pts//2]

        # reshape spectrogram
        newspectro = np.zeros((np.shape(new_wvl)[0], np.shape(delays)[0]))
        for j, _ in enumerate(spectrogram[0, :]):
            f = interp1d(wavelengths, spectrogram[:, j])
            newspectro[:, j] = f(new_wvl)

        # make header
        Header = [nb_delay_pts, nb_wvl_pts, delay_inc, wvl_inc, wvl_central] #header for FROG retrieval

        # dump to new file
        processed_file_path = file_path + "_Reshaped.txt"
        with open(processed_file_path, 'w') as f: 
            for line in Header:
                    f.write(str(line) + '\n')
            for i, _ in enumerate (newspectro[:, 0]):
                for j, _ in enumerate (newspectro[0, :]):               
                    f.write(str(newspectro[i,j]) + '\n')
            f.close()
        return processed_file_path

    def start_monitoring_for_files(self, directory, extension = ".dat"):
        event_handler = FileCreationHandler(extension = extension, callback = self.handle_FROG_recovery)
        observer = Observer()
        observer.schedule(event_handler, directory, recursive=False)

        # Start the observer in a separate thread
        observer_thread = threading.Thread(target = observer.start, daemon=True)
        observer_thread.start()

    def handle_FROG_recovery(self):
        self.raw_wvl, self.raw_delay, self.raw_trace = self.read_trace_dat_file(os.path.join(self.frog_save_subdir, 'a.dat'))
        self.trace_wvl, self.delay, self.recovered_trace = self.read_trace_dat_file(os.path.join(self.frog_save_subdir, 'arecon.dat'))
        self.time, self.amp_t, self.phase_t, self.pulse = self.read_field_dat_file(os.path.join(self.frog_save_subdir, 'ek.dat'))
        self.wvl, self.amp_w, self.phase_w, self.spectrum = self.read_field_dat_file(os.path.join(self.frog_save_subdir, 'speck.dat'))
        self.recovered_trace /= np.max(np.abs(self.recovered_trace))
        self.raw_trace /= np.max(np.abs(self.raw_trace))

        data = {'pro_trace': self.raw_trace,
                'trace': self.recovered_trace,
                'sh_wvl': self.raw_wvl,
                'delay': self.delay,
                'time': self.time,
                'pulse': self.pulse,
                'wvl': self.wvl,
                'spectrum': self.spectrum}
        self.recovery_complete.emit(data, "FROG")

    @staticmethod
    def save_pkl(data, path):
        with open(path, 'wb') as f:
            pickle.dump(data, f)

    @staticmethod
    def load_pkl(path):
        with open(path, 'rb') as f:
            data = pickle.load(f)
        return data

    @staticmethod
    def read_trace_dat_file(filename):
        with open(filename, 'r') as file:
            num_points = int(file.readline().strip().split()[0])
            min_val, max_val = map(float, file.readline().strip().split())
            
            # Read wavelengths
            wavelengths = []
            for _ in range(num_points):
                wavelengths.append(float(file.readline().strip()))
            
            # Read delay values
            delays = []
            for _ in range(num_points):
                delays.append(float(file.readline().strip()))
            
            # Read the 2D data array
            data = []
            for _ in range(num_points):
                row = []
                for _ in range(num_points):
                    row.append(float(file.readline().strip()))
                data.append(row)       
        return np.array(wavelengths), np.array(delays), np.array(data)

    @staticmethod
    def read_field_dat_file(filename):
        data = np.loadtxt(filename)
        
        # Assign each column to a variable
        axis = data[:, 0]
        intensity = data[:, 1]
        phase = data[:, 2]
        pulse = data[:, 3] + 1j*data[:, 4]
        return axis, intensity, phase, pulse

if __name__ == "__main__":
    dummy_file = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Software\Predictor\FROG\Dummy_FROG.txt"
    dummy_file = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\XP_FROG_traces\20250115-SHG-FROG_ProbeAttWedge_1030nm_10fsStepSize"


    app = QApplication(sys.argv)
    window = FROG(dummy_file)
    sys.exit(app.exec_())