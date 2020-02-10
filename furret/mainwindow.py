import os
import pickle
import shutil

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
import furret.config as config
from furret.preferences import load_settings
from furret.query import Query
from furret.preferences import PreferencesDialog


class MainWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        load_settings()
        self.table_widget = QTableWidget()
        self.init_ui()

    def update_table(self):
        self.table_widget.clearContents()

        candidates = [f.path for f in os.scandir(config.working_directory) if f.is_dir()]
        queries = []
        for c in candidates:
            if not os.path.isdir(c):
                continue
            if not os.path.isfile(os.path.join(c, 'query.pickle')):
                continue
            query_file = os.path.join(c, 'query.txt')
            if not os.path.isfile(query_file):
                self.statusBar().showMessage('query.txt file missing')
                continue
            with open(query_file) as text:
                lines = text.readlines()
            if not lines or len(lines) < 2 or len(lines[1]) < 19:
                self.statusBar().showMessage('query.txt file missing')
                continue
            queries.append((lines[0], lines[1][:19], c))
        content = sorted(queries, key=lambda xx: xx[1])

        self.table_widget.setRowCount(len(content))

        for i, line in enumerate(content):
            self.table_widget.setItem(i, 0, QTableWidgetItem(line[0]))
            self.table_widget.setItem(i, 1, QTableWidgetItem(line[1]))
            self.table_widget.setItem(i, 2, QTableWidgetItem(line[2]))
            self.table_widget.move(0, 0)

        self.table_widget.resizeColumnsToContents()
        x = self.table_widget.verticalHeader().size().width()
        for i in range(self.table_widget.columnCount()):
            x += self.table_widget.columnWidth(i)

        y = self.table_widget.horizontalHeader().size().height()
        for i in range(self.table_widget.rowCount()):
            y += self.table_widget.rowHeight(i)

        y += self.menuBar().height()
        y += self.statusBar().height()

        self.setFixedSize(max(x + 2, 200), y + 2)

    # noinspection PyUnresolvedReferences
    def init_ui(self):

        # Actions

        # New Query
        new_query_action = QAction("New...", self)
        new_query_action.setShortcut("Ctrl+N")
        new_query_action.triggered.connect(self.new_query)
        # Delete Query
        delete_query_action = QAction("Delete...", self)
        delete_query_action.triggered.connect(self.delete_query)
        # Exit
        exit_action = QAction('&Quit', self)
        exit_action.setShortcut('Ctrl+Q')
        exit_action.triggered.connect(self.quit)
        # Preferences
        preferences_action = QAction("&Preferences...", self)
        preferences_action.triggered.connect(self.preferences)
        # About
        about_action = QAction("&About...", self)
        about_action.triggered.connect(self.about)
        # Make tables
        make_tables_action = QAction("&Make Tables", self)
        make_tables_action.triggered.connect(self.make_tables)
        # Download structures
        download_structures_action = QAction("&Download PDB", self)
        download_structures_action.triggered.connect(self.download_structures)
        # Families Sequences Action
        families_sequences_action = QAction("&Generate Sequences", self)
        families_sequences_action.triggered.connect(self.families_sequences)
        # Families Structures Action
        families_structures_action = QAction("&Generate Structures", self)
        families_structures_action.triggered.connect(self.families_structures)
        # Families Motives Action
        families_motives_action = QAction("&Generate Motives", self)
        families_motives_action.triggered.connect(self.families_motives)
        summary_report_action = QAction("&Summary", self)
        summary_report_action.triggered.connect(self.report_summary)

        self.statusBar()

        # self.toolbar = self.addToolBar('Bar')

        menubar = self.menuBar()
        menubar.setNativeMenuBar(False)
        app_menu = menubar.addMenu(config.APPLICATION_NAME)
        app_menu.addAction(about_action)
        app_menu.addAction(preferences_action)
        app_menu.addAction(exit_action)
        fam_menu = QMenu("Families", self)
        fam_menu.addAction(families_sequences_action)
        fam_menu.addAction(families_structures_action)
        fam_menu.addAction(families_motives_action)
        query_menu = menubar.addMenu('&Query')
        query_menu.addAction(new_query_action)
        query_menu.addAction(delete_query_action)
        query_menu.addAction(download_structures_action)
        query_menu.addMenu(fam_menu)
        query_menu.addAction(summary_report_action)

        self.table_widget.setColumnCount(3)
        self.table_widget.setColumnHidden(2, True)
        self.table_widget.setRowCount(0)
        self.table_widget.setHorizontalHeaderLabels(['Query', 'Date'])
        self.table_widget.horizontalHeaderItem(0).setTextAlignment(Qt.AlignHCenter)
        self.table_widget.horizontalHeaderItem(1).setTextAlignment(Qt.AlignHCenter)
        self.table_widget.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
        self.table_widget.setAlternatingRowColors(True)
        self.table_widget.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.table_widget.setSelectionBehavior(QAbstractItemView.SelectRows)
        # self.table_widget.cellDoubleClicked.connect(self.got_double_click)
        self.table_widget.setSelectionMode(QAbstractItemView.SingleSelection)
        self.update_table()

        self.setWindowTitle(f'{config.APPLICATION_NAME} Main Window')

        self.setCentralWidget(self.table_widget)

        self.show()
        self.update_table()

    def new_query(self):
        dlg = QInputDialog(self)
        dlg.setInputMode(QInputDialog.TextInput)
        dlg.setLabelText("Query:")
        dlg.resize(500, 100)
        ok = dlg.exec_()
        query_text = dlg.textValue()

        if ok:
            _ = Query(query_text, self.statusBar())
            self.update_table()

    def delete_query(self):
        a = self.table_widget.selectedIndexes()
        if len(a) > 0:
            row = a[0].row()
            should_delete = QMessageBox.question(self,
                                                 "Deleting Query",
                                                 "Are You sure you want to delete the selected query?",
                                                 QMessageBox.Yes, QMessageBox.No)
            if should_delete == QMessageBox.Yes:
                the_dir = self.table_widget.item(row, 2).text()
                shutil.rmtree(the_dir)
                self.table_widget.removeRow(row)
                self.update_table()

        else:
            self.statusBar().showMessage('Select a query to be removed!')

    def make_tables(self):
        self.statusBar().showMessage('Make Tables Action triggered!')

    def download_structures(self):
        self.statusBar().showMessage('Download Structures Action triggered!')
        the_dir = self.get_selection_directory()
        if the_dir:
            with open(os.path.join(the_dir, 'query.pickle'), 'rb') as query_pickle:
                the_query = pickle.load(query_pickle)
                the_query.download_structures(self.statusBar())
        else:
            self.statusBar().showMessage('Select a query to download structures!')

    def families_sequences(self):
        the_dir = self.get_selection_directory()
        if the_dir:
            with open(os.path.join(the_dir, 'query.pickle'), 'rb') as query_pickle:
                the_query = pickle.load(query_pickle)
                the_query.gen_fam_seq(self.statusBar())
        else:
            self.statusBar().showMessage('Select a query')

    def families_structures(self):
        the_dir = self.get_selection_directory()
        if the_dir:
            with open(os.path.join(the_dir, 'query.pickle'), 'rb') as query_pickle:
                the_query = pickle.load(query_pickle)
                the_query.gen_fam_struct(self.statusBar())
        else:
            self.statusBar().showMessage('Select a query')

    def families_motives(self):
        the_dir = self.get_selection_directory()
        if the_dir:
            with open(os.path.join(the_dir, 'query.pickle'), 'rb') as query_pickle:
                the_query = pickle.load(query_pickle)
                the_query.gen_meme(self.statusBar())
        else:
            self.statusBar().showMessage('Select a query')

    def report_summary(self):
        the_dir = self.get_selection_directory()
        if the_dir:
            with open(os.path.join(the_dir, 'query.pickle'), 'rb') as query_pickle:
                the_query = pickle.load(query_pickle)
                the_query.gen_summary(self.statusBar())
        else:
            self.statusBar().showMessage('Select a query')

    def renew_query(self):
        self.statusBar().showMessage(f"To do: renew_query")
        return
        # the_dir = self.get_selection_directory()
        # the_file = ''
        # if the_dir:
        #     try:
        #         the_file = os.path.join(the_dir, 'query.txt')
        #         with open(the_file, 'rt') as query_file:
        #             lines = query_file.readlines()
        #             query_text = lines[0].strip()
        #             if query_text:
        #                 _ = Query(query_text, self.statusBar(), import_dir=os.path.join(the_dir, 'Imported'))
        #                 self.update_table()
        #             else:
        #                 self.statusBar().showMessage(f"No query found in {the_file}")
        #     except OSError:
        #         self.statusBar().showMessage(f"Can't read {the_file}")
        # else:
        #     self.statusBar().showMessage(f"No query selected")

    def get_selection_directory(self):
        result = ''
        a = self.table_widget.selectedIndexes()
        if len(a) > 0:
            n = a[0].row()
            result = self.table_widget.item(n, 2).text()
        return result

    @staticmethod
    def quit():
        the_app = QApplication.instance()
        the_app.quit()

    def preferences(self):
        dialog = PreferencesDialog(self)
        dialog.exec_()

    def about(self):
        QMessageBox.about(self, "", f"{config.APPLICATION_NAME}\n An application to help classify toxic proteins")
