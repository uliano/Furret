from PyQt5.QtWidgets import *
from pathlib import Path
from PyQt5.QtCore import QSettings
import furret.config as config
from typing import Optional
import os


class PreferencesDialog(QDialog):
    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        dir_label = QLabel("Working Directory:")
        self.working_directory = QLineEdit()
        meme_exe_label = QLabel("Meme executable:")
        self.meme_executable = QLineEdit()
        meme_options_label = QLabel("Meme options:")
        self.meme_options = QLineEdit()
        # moe_exe_label = QLabel("MOEbatch executable")
        # self.moe_executable = QLineEdit()
        entrez_email_label = QLabel("Email (for Entrez):")
        self.entrez_email = QLineEdit()

        ok_button = QPushButton("OK")
        cancel_button = QPushButton("Cancel")
        button_line = QHBoxLayout()
        button_line.addStretch()
        button_line.addWidget(ok_button)
        button_line.addWidget(cancel_button)
        grid = QGridLayout()
        grid.addWidget(dir_label, 0, 0)
        grid.addWidget(self.working_directory, 0, 1)
        grid.addWidget(meme_exe_label, 1, 0)
        grid.addWidget(self.meme_executable, 1, 1)
        grid.addWidget(meme_options_label, 2, 0)
        grid.addWidget(self.meme_options, 2, 1)
        # grid.addWidget(moe_exe_label, 3, 0)
        # grid.addWidget(self.moe_executable, 3, 1)
        grid.addWidget(entrez_email_label, 3, 0)
        grid.addWidget(self.entrez_email, 3, 1)

        main_layout = QVBoxLayout()
        main_layout.addLayout(grid)
        main_layout.addStretch()
        main_layout.addLayout(button_line)
        self.setLayout(main_layout)

        self.setWindowTitle("Preferences")
        self.setMinimumWidth(900)

        ok_button.clicked.connect(self.accept)
        cancel_button.clicked.connect(self.reject)

        self.working_directory.setText(config.working_directory)
        self.meme_executable.setText(config.meme_executable)
        self.meme_options.setText(config.meme_options)
        # self.moe_executable.setText(config.moe_executable)
        self.entrez_email.setText(config.entrez_email)

    def accept(self) -> None:
        settings = QSettings(config.APPLICATION_NAME, config.COMPANY_NAME)
        settings.setValue('workingDirectory', self.working_directory.text())
        settings.setValue('memeExecutable', self.meme_executable.text())
        settings.setValue('memeOptions', self.meme_options.text())
        # settings.setValue('moeExcecutable', self.moe_executable.text())
        settings.setValue('entrezEmail', self.entrez_email.text())
        load_settings()
        os.makedirs(self.working_directory.text(), exist_ok=True)
        super().accept()


def load_settings() -> None:
    """loads saved preferences into config module variables"""
    settings = QSettings(config.APPLICATION_NAME, config.COMPANY_NAME)

    config.working_directory = settings.value('workingDirectory', str(Path.home()))
    config.meme_executable = settings.value('memeExecutable', '')
    config.meme_options = settings.value('memeOptions',
                                         '-protein -oc . -nostatus -time 18000 -mod zoops -nmotifs 100'
                                         ' -minw 6 -maxw 50 -objfun classic -markov_order 0')
    config.entrez_email = settings.value('entrezEmail', 'my.name@my.domain')
    # config.moe_executable = settings.value('moeExcecutable', '')
