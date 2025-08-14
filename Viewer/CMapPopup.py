
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from PyQt5.QtWidgets import QApplication, QDialog, QWidget, QVBoxLayout, QPushButton, QColorDialog, QLabel, QHBoxLayout, QScrollArea, QFrame, QMenu, QLineEdit, QFileDialog
import os
import json

class CMapPopup(QDialog):
    def __init__(self, parent = None):
        super().__init__(parent)

        self.setWindowTitle('Custom Cmap')
        self.setGeometry(100, 100, 400, 400)
        self.colors = []
        self.init_ui()
        # self.show()

    def init_ui(self):
        layout = QVBoxLayout()

        # Create a scroll area to hold the color squares
        self.color_scroll_area = QScrollArea(self)
        self.color_scroll_area.setWidgetResizable(True)
        self.color_container = QWidget()  # This will contain all the color squares
        self.color_container_layout = QVBoxLayout(self.color_container)
        self.color_scroll_area.setWidget(self.color_container)
        layout.addWidget(self.color_scroll_area)

        # Button to add a color
        self.add_color_button = QPushButton('Add Color', self)
        self.add_color_button.clicked.connect(self.add_color)
        layout.addWidget(self.add_color_button)

        # Button to create and display colormap
        self.create_button = QPushButton('Save Cmap', self)
        self.create_button.clicked.connect(self.save_colormap)
        layout.addWidget(self.create_button)

        # Set the layout
        self.setLayout(layout)

    def add_color(self):
        color = QColorDialog.getColor()
        if color.isValid():
            self.colors.append(color.name())
            self.update_color_squares()

    def update_color_squares(self):
        # Clear the existing color squares
        for i in reversed(range(self.color_container_layout.count())):
            widget = self.color_container_layout.itemAt(i).widget()
            if widget is not None:
                widget.deleteLater()

        # Add a new square for each selected color
        for color in self.colors:
            color_square = QLabel(self)
            color_square.setFixedSize(30, 30)  # Size of the color square
            color_square.setStyleSheet(f"background-color: {color};")
            color_square.setContextMenuPolicy(3)  # Enable right-click context menu
            color_square.customContextMenuRequested.connect(lambda pos, color = color: self.show_right_click_menu(pos, color_square, color))
            self.color_container_layout.addWidget(color_square)

    def show_right_click_menu(self, pos, color_square, color):
        # Create the context menu (right-click menu)
        context_menu = QMenu(self)
        delete_action = context_menu.addAction("Delete Color")
        flip_action = context_menu.addAction("Flip Order")
        demo_action = context_menu.addAction("Demo Cmap")

        action_dummy = context_menu.exec_(self.mapToGlobal(pos))

        if action_dummy == delete_action:
            self.delete_color(color_square, color)
        elif action_dummy == flip_action:
            self.flip_color_order()
        elif action_dummy == demo_action:
            self.show_colormap()

    def delete_color(self, color_square, color):
        self.colors.remove(color)
        self.update_color_squares()

    def flip_color_order(self):
        self.colors.reverse()
        self.update_color_squares()

    def show_colormap(self):
        plt.imshow(np.random.rand(10, 10), cmap = self.generate_cmap())
        plt.colorbar()
        plt.title("Custom Colormap Example")
        plt.show()

    def save_colormap(self):
        initial_directory = os.path.join(os.getcwd(), "Software/Viewer/CMap")
        cmap_path, _ = QFileDialog.getSaveFileName(self, "Save Colormap", initial_directory, "")
        
        custom_cmap = self.generate_cmap()
        samples = np.linspace(0, 1, 256)
        cmap_colors = [custom_cmap(x)[:3] for x in samples]  # Ignore the alpha channel
        cmap_data = [{"r": r, "g": g, "b": b} for r, g, b in cmap_colors]
        if not cmap_path.endswith(".json"):
            cmap_path += ".json"
        with open(cmap_path, 'w') as f:
            json.dump(cmap_data, f)

        print(f"Colormap saved to {cmap_path}.")
        
        plt.close()
        self.accept()  # Close the dialog and return success

    def generate_cmap(self):
        if len(self.colors) < 2:
            return  # Ensure at least two colors are selected
        return LinearSegmentedColormap.from_list("dummy", self.colors)


if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = CMapPopup()
    window.show()
    sys.exit(app.exec_())
