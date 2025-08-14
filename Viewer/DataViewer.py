import sys
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QHBoxLayout,
    QWidget, QStackedWidget, QMenu, QAction, QFileDialog,
    QTreeWidget, QTreeWidgetItem
)
from PyQt5.QtCore import Qt
import numpy as np

try:
    from Software.Viewer.PredictionPlot import PredictionPlot
    from Software.Viewer.ComparisonPlot import ComparisonPlot
    from Software.Viewer.TracePlot import TracePlot
except:
    from PredictionPlot import PredictionPlot
    from ComparisonPlot import ComparisonPlot
    from TracePlot import TracePlot
import os 
import pickle

class DataViewer(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Data Viewer")
        self.init_ui()
        self.showMaximized()

    def init_ui(self):
        # Main layout
        main_layout = QHBoxLayout()
        main_container = QWidget()
        main_container.setLayout(main_layout)
        self.setCentralWidget(main_container)

        # Sidebar for file management
        self.sidebar = QTreeWidget()
        self.sidebar.setHeaderLabels(["Data Available For Viewing"])
        self.sidebar.setFixedWidth(250)  # Set the sidebar width in pixels
        self.sidebar.setContextMenuPolicy(Qt.CustomContextMenu)
        self.sidebar.customContextMenuRequested.connect(self.show_context_menu)
        self.sidebar.itemChanged.connect(self.update_display_plot)
        main_layout.addWidget(self.sidebar)

        # Stacked widget for visualization
        self.stacked_widget = QStackedWidget()
        self.trace_widget = TracePlot()
        self.pred_widget = PredictionPlot()
        self.comp_widget = ComparisonPlot()
        self.stacked_widget.addWidget(self.trace_widget)
        self.stacked_widget.addWidget(self.pred_widget)
        self.stacked_widget.addWidget(self.comp_widget)
        main_layout.addWidget(self.stacked_widget)

        self.data = {}  # Dictionary to hold data

    def show_context_menu(self, position):
        # Get the item at the position where the context menu was requested
        item = self.sidebar.itemAt(position)

        # Create the context menu
        context_menu = QMenu()

        # Option for all items
        add_action = QAction("Add Data File", self)
        add_action.triggered.connect(self.add_data_file)
        context_menu.addAction(add_action)

        # Specific options for top-level items
        if item and not item.parent():
            delete_action = QAction("Delete", self)
            delete_action.triggered.connect(lambda: self.delete_item(item))
            context_menu.addAction(delete_action)

        # Show the context menu
        context_menu.exec_(self.sidebar.mapToGlobal(position))

    def delete_item(self, item):
        # Delete the specified top-level item
        index = self.sidebar.indexOfTopLevelItem(item)
        if index != -1:
            self.sidebar.takeTopLevelItem(index)

    def add_file_to_tree(self, file_ID, subitems):
        # add file_ID as parent
        root_item = QTreeWidgetItem([file_ID])
        root_item.setFlags(root_item.flags())
        self.sidebar.addTopLevelItem(root_item)

        # adds children for each recovery file
        for subitem_name in subitems:
            sub_item = QTreeWidgetItem([subitem_name])
            sub_item.setFlags(sub_item.flags() | Qt.ItemIsUserCheckable)
            sub_item.setCheckState(0, Qt.Unchecked)
            root_item.addChild(sub_item)
        root_item.setExpanded(True)

    def add_data_file(self):
        dir_path = QFileDialog.getExistingDirectory(self, "Open Data Directory", os.path.join(os.getcwd(), "Recovered_Data"))
        # dir_path = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt"

        # set up tree with files
        file_ID = os.path.basename(dir_path)
        subitems = [file for file in os.listdir(dir_path) if (file.endswith('.pkl') or file.endswith('.txt'))]
        self.add_file_to_tree(file_ID, subitems)
    
        # load data into self.data
        self.load_data_from_dir(dir_path, file_ID, subitems)

    def load_data_from_dir(self, dir_path, file_ID, subitems):
        # load all data into data dict with single file ID
        self.data[file_ID] = {}
        for item in subitems:
            if item.endswith(".pkl"):
                data = self.load_pkl(os.path.join(dir_path, item))
                self.data[file_ID][item] = data

        # load raw data
        raw_data = np.loadtxt(os.path.join(dir_path, 'raw_data.txt'))
        self.data[file_ID]['raw_data.txt'] = {'delay': raw_data[0, 1:], 'wvl': raw_data[1:, 0], 
                                              'trace': raw_data[1:, 1:] / np.max(raw_data[1:, 1:])}

    def update_display_plot(self):
        checked_items = self.count_checked_files(self.sidebar)

        if len(checked_items) == 1:      # one file is checked
            file_ID, subfile = checked_items[0]
            if subfile.endswith('.txt'):           # if that data is raw
                self.trace_widget.plot(self.data[file_ID][subfile])
                self.stacked_widget.setCurrentIndex(0)  # switch to trace plot

            elif subfile.endswith("FROG.pkl"):      # for FROG predictions
                self.pred_widget.clear_plot()
                self.pred_widget.plot(self.data[file_ID]['raw_data.txt'], self.data[file_ID][subfile])
                self.stacked_widget.setCurrentIndex(1)  

            elif subfile.endswith(".pkl") and "FROG" not in subfile:      # for FROSt predictions
                self.pred_widget.clear_plot()
                self.pred_widget.plot(self.data[file_ID]['raw_data.txt'], self.data[file_ID][subfile])
                self.stacked_widget.setCurrentIndex(1)  

        elif len(checked_items) >= 2:
            file_ID = {file_ID for file_ID, _ in checked_items}
            pkl_files = [item for item in checked_items if item[1].endswith('.pkl')]

            # plots are for the same raw trace
            if len(file_ID) == 1:
                file_ID = file_ID.pop()
                subfiles = [subfile for _, subfile in checked_items]               

                self.stacked_widget.setCurrentIndex(2)  # switch to comparison plot
                self.comp_widget.clear_plot()

                # for 2 FROSt plots
                if (all("FROG" not in item for item in subfiles)) and len(pkl_files) == 2:
                    self.comp_widget.plot(self.data[file_ID]['raw_data.txt'],
                                          self.data[file_ID][subfiles[0]], self.data[file_ID][subfiles[1]], 
                                          labels = [filename.split(".")[0] for filename in subfiles])

                # FROG and FROSt
                elif (any("FROG" in item for item in subfiles)) and len(pkl_files) == 2:
                    frog_files = [subfile for subfile in subfiles if "FROG" in subfile]
                    frost_files = [subfile for subfile in subfiles if subfile not in frog_files]
                    self.comp_widget.plot(raw_data = self.data[file_ID]['raw_data.txt'],
                                          data1 = self.data[file_ID][frost_files[0]],
                                          labels = [filename.split(".")[0] for filename in frost_files])
                    self.comp_widget.add_FROG_data(self.data[file_ID][frog_files[-1]])

                # handles 2 FROSt plots and a FROG
                elif len(pkl_files) == 3 and (any("FROG" in item for item in subfiles)):
                    frog_files = [subfile for subfile in subfiles if "FROG" in subfile]
                    frost_files = [subfile for subfile in subfiles if subfile not in frog_files]
                    self.comp_widget.plot(self.data[file_ID]['raw_data.txt'],
                                          self.data[file_ID][frost_files[0]], self.data[file_ID][frost_files[1]], 
                                          labels = [filename.split(".")[0] for filename in frost_files])
                    self.comp_widget.add_FROG_data(self.data[file_ID][frog_files[-1]])
    
    @staticmethod
    def count_checked_files(tree_widget):
        """Counts the number of checked items in a QTreeWidget ."""
        checked_items = []

        def check_item_recursive(item):
            if item.checkState(0) == Qt.Checked:
                parent = item.parent()
                parent_name = parent.text(0) if parent is not None else None
                checked_items.append((parent_name, item.text(0)))
            for i in range(item.childCount()):
                check_item_recursive(item.child(i))

        for i in range(tree_widget.topLevelItemCount()):
            check_item_recursive(tree_widget.topLevelItem(i))

        return checked_items

    @staticmethod
    def load_pkl(path):
        with open(path, 'rb') as f:
            data = pickle.load(f)
        return data

if __name__ == "__main__":
    import os
    import pickle
    
    app = QApplication(sys.argv)
    dir = r"C:\Users\lepar\Documents\MSc\Phase_Recovery_UIs\Recovered_Data\20201201_1erFROst.txt"
    data = np.loadtxt(os.path.join(dir, 'raw_data.txt'))
    raw_delay, raw_wvl, raw_trace = data[0, 1:], data[1:, 0], data[1:, 1:] / np.max(data[1:, 1:])

    viewer = DataViewer()

    viewer.add_data_file()

    sys.exit(app.exec_())
