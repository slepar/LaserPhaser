from PyQt5.QtWidgets import (QVBoxLayout, QDialog, QDialogButtonBox, QLabel)

class ErrorPopup(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)

        # Set up the dialog layout
        self.setWindowTitle("Error: Divergence Detected")
        self.setGeometry(500, 500, 300, 150)
        
        layout = QVBoxLayout()

        # Add a label to show the error message
        label = QLabel("The parameters selected during pre-processing have lead to divergence.")
        layout.addWidget(label)

        # Add buttons
        self.button_box = QDialogButtonBox()
        self.close_button = self.button_box.addButton("Cancel", QDialogButtonBox.RejectRole)
        self.switch_button = self.button_box.addButton("Retry processing", QDialogButtonBox.AcceptRole)

        layout.addWidget(self.button_box)

        # Connect buttons to their respective actions
        self.close_button.clicked.connect(self.reject)  # Close popup
        self.switch_button.clicked.connect(self.accept)  # Signal to switch widget

        self.setLayout(layout)

        
    def return_to_screen(self):
        # Signal that this button was pressed (AcceptRole used as a placeholder)
        self.done(QDialog.Accepted)


# # Run the application
# if __name__ == "__main__":
#     app = QApplication(sys.argv)
#     # main_window = MainWindow()
#     main_window.show()
#     sys.exit(app.exec_())
