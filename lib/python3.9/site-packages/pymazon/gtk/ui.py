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

import pygtk
pygtk.require('2.0')
import gobject
import gtk

from pymazon.core.downloader import Downloader
from pymazon.core.item_model import Album
from pymazon.core.parser import AmzParser
from pymazon.core.settings import settings
from pymazon.gtk.tree_model import TreeModel
from pymazon.resource import python_icon_path


class ProgressRenderer(gtk.CellRendererProgress):
    def do_render(self, window, widget, background_area, area, expose_area, flags):       
        gtk.CellRendererProgress.do_render(self, window, widget, 
                                       background_area, area,
                                       expose_area, flags)     


gobject.type_register(ProgressRenderer)


class NameFormatDialog(gobject.GObject):
    def __init__(self, builder):
        self.nameFormatDialog = builder.get_object('nameFormatDialog')
        self.nameFormatDialog.set_modal(True)
        self.nameFormatText = builder.get_object('nameFormatText')
        self.connect_signals(builder)
        
    def connect_signals(self, builder):
        nameFormatOkButton = builder.get_object('nameFormatOkButton')
        nameFormatOkButton.connect('clicked', self.on_nameFormatOkButton_clicked)
        nameFormatCancelButton = builder.get_object('nameFormatCancelButton')
        nameFormatCancelButton.connect('clicked', self.on_nameFormatCancelButton_clicked)
    
    def on_nameFormatOkButton_clicked(self, *args):
        self.nameFormatDialog.emit('response', 1)     
    
    def on_nameFormatCancelButton_clicked(self, *args):
        self.nameFormatDialog.emit('response', 0)

    def launch(self):
        self.nameFormatText.set_text(settings.name_template)
        if self.nameFormatDialog.run():
            txt = self.nameFormatText.get_text()
        else:
            txt = None
        self.nameFormatDialog.hide()
        return txt
        
    
class SettingsDialog(gobject.GObject):
    def __init__(self, builder):
        self.settingsDialog = builder.get_object('settingsDialog')        
        self.settingsDialog.set_modal(True)        
        self.nameFormatDialog = NameFormatDialog(builder)
        self.nameFormatLineEdit = builder.get_object('nameFormatLineEdit')
        self.saveDirLineEdit = builder.get_object('saveDirLineEdit')
        self.toolkitComboBox = builder.get_object('toolkitComboBox')
        self.numThreadsSpinBox = builder.get_object('numThreadsSpinBox')
        self.connect_signals(builder)
        
    def connect_signals(self, builder):
        settingsOkButton = builder.get_object('settingsOkButton')
        settingsOkButton.connect('clicked', self.on_settingsOkButton_clicked)
        settingsCancelButton = builder.get_object('settingsCancelButton')
        settingsCancelButton.connect('clicked', self.on_settingsCancelButton_clicked)
        saveDirButton = builder.get_object('saveDirButton')
        saveDirButton.connect('clicked', self.on_saveDirButton_clicked)
        nameFormatButton = builder.get_object('nameFormatButton')
        nameFormatButton.connect('clicked', self.on_nameFormatButton_clicked)        
        
    def on_settingsOkButton_clicked(self, *args):
        self.settingsDialog.emit('response', 1)   
        
    def on_settingsCancelButton_clicked(self, *args):
        self.settingsDialog.emit('response', 0)
    
    def on_saveDirButton_clicked(self, *args):
        caption = 'Choose Save Directory'              
        dialog = gtk.FileChooserDialog(caption, self.settingsDialog,                                       
                                       buttons=(gtk.STOCK_CANCEL,
                                                gtk.RESPONSE_REJECT,
                                                gtk.STOCK_OK,
                                                gtk.RESPONSE_ACCEPT))
        dialog.set_action(gtk.FILE_CHOOSER_ACTION_SELECT_FOLDER)
        dialog.set_current_folder(settings.save_dir)
        res = dialog.run()
        new_dir = dialog.get_current_folder()
        dialog.destroy()
        if res == gtk.RESPONSE_REJECT:
            return
        if not os.access(new_dir, os.W_OK):
            msg = 'No write access to specified directory. '
            msg += 'Try again.'
            dialog = gtk.MessageDialog(self.settingsDialog, 
                                       buttons=gtk.BUTTONS_OK)
            dialog.set_markup(msg)
            dialog.run()
            dialog.destroy()
            self.on_saveDirButton_clicked()
            return
        self.saveDirLineEdit.set_text(new_dir)
        
    def on_nameFormatButton_clicked(self, *args):
        new_format = self.nameFormatDialog.launch()        
        if new_format:            
            self.nameFormatLineEdit.set_text(new_format)
        
    def set_combo_box(self):
        # this is currently extremely hackish
        toolkits = {'qt4': 0, 'gtk': 1}
        idx = toolkits[settings.toolkit]
        self.toolkitComboBox.set_active(idx)
        
    def launch(self):
        self.saveDirLineEdit.set_text(settings.save_dir)
        self.nameFormatLineEdit.set_text(settings.name_template)
        self.set_combo_box()
        self.numThreadsSpinBox.set_value(settings.num_threads)
        if self.settingsDialog.run():
            settings.save_dir = self.saveDirLineEdit.get_text()
            settings.name_template = self.nameFormatLineEdit.get_text()
            settings.toolkit = self.toolkitComboBox.get_active_text()
            settings.num_threads = self.numThreadsSpinBox.get_value()
            settings.write_config_file()
        self.settingsDialog.hide()            
        
        
class MainWindow(gobject.GObject):    
    def __init__(self, amzs):
        glade_fpath = os.path.join(os.path.split(__file__)[0], '_ui.glade')
        builder = gtk.Builder()
        builder.add_from_file(glade_fpath)            
        
        self.mainWindow = builder.get_object('MainWindow')
        self.mainWindow.show()
          
        self.albumArtImage = builder.get_object('albumArtImage')
        self.nowDownloadingLabel = builder.get_object('nowDownloadingLabel')
        self.pythonPixbuf = gtk.gdk.pixbuf_new_from_file(python_icon_path)
        self.pymazon_text = '<span weight="bold" size="xx-large">Pymazon! w00t!</span>'
        self.reset_displaybar_info()        
        
        self.settings_dialog = SettingsDialog(builder)             

        self.treeView = builder.get_object('treeView')
        self.tree_model = None 
        
        self.pbarRenderer = ProgressRenderer()       
        self.colStatus = builder.get_object('colStatus')
        self.colStatus.clear()
        self.colStatus.pack_start(self.pbarRenderer, True)
        self.colStatus.add_attribute(self.pbarRenderer, 'value', 1)
        self.colStatus.add_attribute(self.pbarRenderer, 'text', 2)
        
        self.connect_signals(builder)
        
        self.downloader = None
        self.current_album = None
        self.old_albums = []
        if amzs:
            self.load_new_amz_files(amzs)
        
    def connect_signals(self, builder):
        builder.connect_signals({'on_MainWindow_destroy': self.on_destroy})
        self.actionAbout = builder.get_object('actionAbout')
        self.actionAbout.connect('activate', self.on_actionAbout_activate)
        self.actionLoadFiles = builder.get_object('actionLoadFiles')
        self.actionLoadFiles.connect('activate', self.on_actionLoadFiles_activate)
        self.actionSettings = builder.get_object('actionSettings')
        self.actionSettings.connect('activate', self.on_actionSettings_activate)
        self.actionShowDownloads = builder.get_object('actionShowDownloads')
        self.actionShowDownloads.connect('activate', self.on_actionShowDownloads_activate)
        self.actionDownload = builder.get_object('actionDownload')
        self.actionDownload.connect('activate', self.on_actionDownload_activate)      
        self.actionQuit = builder.get_object('actionQuit')
        self.actionQuit.connect('activate', self.on_destroy)        
        
    def on_actionAbout_activate(self, *args):
        webbrowser.open('http://code.google.com/p/pymazon')

    def on_actionLoadFiles_activate(self, *args):
        self.new_amz()

    def on_actionSettings_activate(self, *args):
        self.settings_dialog.launch()

    def on_actionShowDownloads_activate(self, *args):
        webbrowser.open(settings.save_dir)

    def on_actionDownload_activate(self, *args):
        self.download_tracks()

    def on_destroy(self, *arg):
        if self.downloader:
            self.downloader.kill()
        gtk.main_quit()    

    def new_amz(self):
        caption = 'Choose AMZ File'
        filefilter = gtk.FileFilter()        
        filefilter.add_pattern('*.amz')        
        dialog = gtk.FileChooserDialog(caption, self.mainWindow,                                       
                                       buttons=(gtk.STOCK_CANCEL,
                                                gtk.RESPONSE_REJECT,
                                                gtk.STOCK_OK,
                                                gtk.RESPONSE_ACCEPT))
        dialog.set_select_multiple(True)
        dialog.set_filter(filefilter)
        res = dialog.run()
        files = dialog.get_filenames()
        dialog.destroy()
        if res == gtk.RESPONSE_REJECT:
            return
        if files:
            self.load_new_amz_files(files)

    def load_new_amz_files(self, amz_files):
        parser = AmzParser()
        for f in amz_files:
            parser.parse(f)
        objects = parser.get_parsed_objects()        
        self.tree_model = TreeModel(objects)        
        self.treeView.set_model(self.tree_model)
        self.treeView.expand_all()

    def download_tracks(self):
        self.num_threads = 1        
        if not self.tree_model:
            return
        save_dir = settings.save_dir               
        if not os.access(save_dir, os.W_OK):
            print 'here'
            msg = 'No write access to save directory. '
            msg += 'Choose a new directory with write privileges.'
            dialog = gtk.MessageDialog(self.mainWindow, 
                                       buttons=gtk.BUTTONS_OK)
            dialog.set_markup(msg)
            dialog.run()
            dialog.destroy()
            return        
        self.downloader = Downloader(self.tree_model, self.update_cb, 
                                     self.finished_cb)
        self.actionDownload.set_sensitive(False)
        self.actionSettings.set_sensitive(False)
        self.actionLoadFiles.set_sensitive(False)        
        self.downloader.start()            

    def update_cb(self, node):
        if isinstance(node.elem.obj, Album):
            if not node.elem.obj is self.current_album:
                if node.elem.obj not in self.old_albums:
                    self.old_albums.append(node.elem.obj)
                    self.current_album = node.elem.obj
                    self.update_album_info()
                else:
                    pass     
        self.tree_model.update_node(node)
    
    def finished_cb(self):
        self.downloader = None
        self.reset_displaybar_info()
        self.actionDownload.set_sensitive(True)
        self.actionSettings.set_sensitive(True)
        self.actionLoadFiles.set_sensitive(True)
        
    def reset_displaybar_info(self):
        self.current_album = None
        self.old_albums = []        
        self.albumArtImage.set_from_pixbuf(self.pythonPixbuf)
        self.nowDownloadingLabel.set_markup(self.pymazon_text)
        
    def update_album_info(self):
        self.update_album_art()
        self.update_downloading_text()
        
    def update_album_art(self):
        img = self.current_album.image
        loader = gtk.gdk.PixbufLoader()
        loader.write(img)
        loader.close()
        pb = loader.get_pixbuf()
        self.albumArtImage.set_from_pixbuf(pb)
        self.update_downloading_text()
        
    def update_downloading_text(self):
        name = self.current_album.artist
        title = self.current_album.title
        txt = 'Now dowloading <b>%s</b> by <b>%s</b>' % (title, name)
        self.nowDownloadingLabel.set_markup(txt)
        
    
def main(amzs):
    gtk.gdk.threads_init()
    main_window = MainWindow(amzs)
    gtk.gdk.threads_enter()    
    gtk.main()
    gtk.gdk.threads_leave()
