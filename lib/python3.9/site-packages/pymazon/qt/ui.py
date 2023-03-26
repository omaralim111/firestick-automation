"""
Pymazon - A Python based downloader for the Amazon.com MP3 store
Copyright (c) 2010 Steven C. Colbert

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
import os
import webbrowser

from PyQt4.QtCore import Qt, pyqtSignal, pyqtSlot   
from PyQt4.QtGui import QDialog, QStyledItemDelegate, QMainWindow, QFileDialog,\
                        QMessageBox, QIcon, QPixmap, QApplication, QStyle,\
                        QStyleOptionProgressBar                      

from pymazon.core.downloader import Downloader
from pymazon.core.item_model import Album
from pymazon.core.parser import AmzParser
from pymazon.core.settings import settings
from pymazon.qt import _ui, _nameformatdialog, _settingsdialog
from pymazon.qt.tree_model import TreeModel
from pymazon.resource import load_icon_path, download_icon_path,\
                             exit_icon_path, settings_icon_path,\
                             python_icon_path, show_icon_path


class NameFormatDialog(QDialog, _nameformatdialog.Ui_NameFormatDialog):
    def __init__(self, parent):
        super(QDialog, self).__init__(parent)
        self.setupUi(self)
        
    def launch(self):
        self.nameFormatLineEdit.setText(settings.name_template)
        res = self.exec_()
        if res:
            return str(self.nameFormatLineEdit.text())
        
        
class SettingsDialog(QDialog, _settingsdialog.Ui_SettingsDialog):
    def __init__(self, parent):
        super(QDialog, self).__init__(parent)
        self.setupUi(self)
        self.nameFormatDialog = NameFormatDialog(self)        
        
    @pyqtSlot()
    def on_nameFormatButton_clicked(self):
        new_format = self.nameFormatDialog.launch()
        if new_format:
            self.nameFormatLineEdit.setText(new_format)
            
    @pyqtSlot()
    def on_saveDirButton_clicked(self):
        res = str(QFileDialog.getExistingDirectory(self, 
                                                   'Choose Save Directory',
                                                   settings.save_dir))
        if res:
            if not os.access(res, os.W_OK):
                msg = 'No write access to specified directory. '
                msg += 'Try again.'
                QMessageBox.information(self, 'No Access!', msg)
                self.on_saveDirButton_clicked()
                return
            self.saveDirLineEdit.setText(res)
            
    def launch(self):
        self.saveDirLineEdit.setText(settings.save_dir)
        self.nameFormatLineEdit.setText(settings.name_template)
        self.numThreadsSpinBox.setValue(settings.num_threads)
        idx = self.toolkitComboBox.findText(settings.toolkit)
        self.toolkitComboBox.setCurrentIndex(idx)
        ok = self.exec_()
        if ok:
            settings.save_dir = str(self.saveDirLineEdit.text())
            settings.name_template = str(self.nameFormatLineEdit.text())
            settings.num_threads = str(self.numThreadsSpinBox.value())
            settings.toolkit = str(self.toolkitComboBox.currentText())
            settings.write_config_file()
        
        
class ProgressBarDelegate(QStyledItemDelegate):       
    def paint(self, painter, options, index):        
        opts = QStyleOptionProgressBar()
        opts.rect = options.rect        
        opts.minimum = 1
        opts.maximum = 100
        opts.textVisible = True        
        data = index.data().toPyObject()        
        percent, txt = data        
        if txt == 'Error!' or txt == 'Ready' or txt == 'Complete!':
            QApplication.style().drawItemText(painter, opts.rect, 4, 
                                              opts.palette, True, txt)
        else:
            opts.progress = percent
            opts.text = txt
            QApplication.style().drawControl(QStyle.CE_ProgressBar, opts, 
                                             painter)      
            

class MainWindow(QMainWindow, _ui.Ui_MainWindow):
    
    updateInfo = pyqtSignal()
    resetInfo = pyqtSignal()
    
    def __init__(self, amz_files, parent=None):
        QMainWindow.__init__(self, parent)
        self.setupUi(self)
        
        loadFileIcon = QIcon.fromTheme('document-new', QIcon(QPixmap(load_icon_path)))
        downloadIcon = QIcon.fromTheme('go-down', QIcon(QPixmap(download_icon_path)))
        settingsIcon = QIcon.fromTheme('preferences-other', QIcon(QPixmap(settings_icon_path)))
        showIcon = QIcon.fromTheme('emblem-downloads', QIcon(QPixmap(show_icon_path)))
        exitIcon = QIcon.fromTheme('window-close', QIcon(QPixmap(exit_icon_path)))
        self.pythonPixmap = QPixmap(python_icon_path)
        self.pymazon_text = self.nowDownloadingLabel.text()
        
        self.actionLoadFiles.setIcon(loadFileIcon)
        self.actionDownload.setIcon(downloadIcon)
        self.actionSettings.setIcon(settingsIcon)
        self.actionShowDownloads.setIcon(showIcon)
        self.actionQuit.setIcon(exitIcon)
        self.albumArtLabel.setPixmap(self.pythonPixmap)
        
        self.settingsDialog = SettingsDialog(self)
                
        self.tree_model = None        
        self.pbardelegate = ProgressBarDelegate()
        self.treeView.setItemDelegateForColumn(1, self.pbardelegate)
        
        self.downloader = None
        self.current_album = None
        self.old_albums = []
        self.updateInfo.connect(self.update_album_info)
        self.resetInfo.connect(self.reset_displaybar_info) 
        
        if amz_files:
            self.load_new_amz_files(amz_files)               
        
    @pyqtSlot()
    def on_actionLoadFiles_triggered(self):
        self.new_amz()
        
    @pyqtSlot()
    def on_actionSettings_triggered(self):        
        self.settingsDialog.launch()        
        
    @pyqtSlot()
    def on_actionDownload_triggered(self):        
        self.download_tracks()
        
    @pyqtSlot()
    def on_actionShowDownloads_triggered(self):
        webbrowser.open(settings.save_dir)
        
    @pyqtSlot()
    def on_actionQuit_triggered(self):
        self.close()
        
    @pyqtSlot()
    def on_actionAbout_triggered(self):
        webbrowser.open('http://code.google.com/p/pymazon')
    
    def closeEvent(self, evt):
        if self.downloader:
            self.downloader.kill()
        evt.accept()

    def new_amz(self):
        caption = 'Choose AMZ File'
        filefilter = "Amazon MP3 Download (*.amz)"
        files = [str(f) for f in 
                 QFileDialog.getOpenFileNames(self, caption, '', filefilter)]
        if files:
            self.load_new_amz_files(files)
            
    def load_new_amz_files(self, amz_files):
        parser = AmzParser()
        for f in amz_files:
            parser.parse(f)
        objects = parser.get_parsed_objects()        
        self.tree_model = TreeModel(objects)        
        self.treeView.setModel(self.tree_model) 
        self.setup_treeview()      
    
    def setup_treeview(self):       
        self.treeView.expandAll()
        self.treeView.resizeColumnToContents(0)
        width = self.treeView.columnWidth(0)
        self.treeView.setColumnWidth(1, 100)
        self.treeView.setColumnWidth(0, width + 50)
        
    def download_tracks(self):
        if not self.tree_model:
            return
        save_dir = settings.save_dir        
        if not os.access(save_dir, os.W_OK):
            msg = 'No write access to save directory. '
            msg += 'Choose a new directory with write privileges.'
            QMessageBox.information(self, 'No Access!', msg)
            return        
        self.downloader = Downloader(self.tree_model, self.update_cb, 
                                     self.finished_cb)
        self.actionDownload.setEnabled(False)
        self.actionSettings.setEnabled(False)
        self.actionLoadFiles.setEnabled(False)        
        self.downloader.start()     

    def update_cb(self, node):
        if isinstance(node.elem.obj, Album):
            if not node.elem.obj is self.current_album:
                if node.elem.obj not in self.old_albums:
                    self.old_albums.append(node.elem.obj)
                    self.current_album = node.elem.obj
                    self.updateInfo.emit()  
                else:
                    pass                          
        self.tree_model.update_node(node)
        
    def finished_cb(self):        
        self.downloader = None
        self.resetInfo.emit()
        self.actionDownload.setEnabled(True)
        self.actionSettings.setEnabled(True)
        self.actionLoadFiles.setEnabled(True)
        
    def reset_displaybar_info(self):
        self.current_album = None
        self.old_albums = []
        self.nowDownloadingLabel.setText(self.pymazon_text)
        self.albumArtLabel.setPixmap(self.pythonPixmap)
        
    def update_album_info(self):
        self.update_album_art()
        self.update_downloading_text()
        
    def update_album_art(self):
        img = self.current_album.image
        pm = QPixmap()
        pm.loadFromData(img, None, Qt.AutoColor)
        self.albumArtLabel.setPixmap(pm)        
        
    def update_downloading_text(self):
        name = self.current_album.artist
        title = self.current_album.title
        txt = 'Now dowloading <b>%s</b> by <b>%s</b>' % (title, name)
        self.nowDownloadingLabel.setText(txt)
        
def main(amzs):
    app = QApplication([])
    win = MainWindow(amzs)
    win.show()                
    app.exec_()
