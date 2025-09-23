from PySide6.QtWidgets import QGroupBox, QFormLayout, QLineEdit, QListView



class PropertyForm(QGroupBox):

    def __init__(self, form=None):
        super().__init__(form.name)
        self.form = form
        layout = QFormLayout()
        self.setLayout(layout)

        self.setup_ui()

    def setup_ui(self):
        form_items = self.form.form()
        for name, val in form_items.items():
            self.layout().addRow(name, QLineEdit(str(val)))


class FluidList(QListView):
    
    def __intit__(self, parent=None, fluid_list=None):
        super().__init__(parent)
        self.setModel(fluid_list)
        
