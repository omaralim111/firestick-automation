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

from ConfigParser import ConfigParser, NoOptionError, NoSectionError
import os
import sys
import warnings


def pymazon_showwarning(message, category, *args):
    print 'Warning:'
    print message
warnings.showwarning = pymazon_showwarning


# so they can spell them however they want
_toolkits = {'qt':'qt4', 'Qt':'qt4', 'qt4':'qt4', 'Qt4':'qt4', 'pyqt':'qt4',
             'PyQt':'qt4', 'Pyqt':'qt4', 'pyqt4':'qt4', 'PyQt4':'qt4',
             'Pyqt4':'qt4', 'gtk':'gtk', 'Gtk':'gtk', 'pygtk':'gtk',
             'PyGtk':'gtk', 'Pygtk':'gtk'}             


def _get_pymazon_dir():
    try:
        pymazon_dir = os.path.join(os.path.expanduser('~'), '.pymazon')
        if not os.path.exists(pymazon_dir):
            os.makedirs(pymazon_dir)
        if not os.access(pymazon_dir, os.W_OK):
            raise
        return pymazon_dir
    except:
        print ('Pymazon is unable to write to the Pymazon config '
               'directory: %s. This directory is required for proper '
               'Pymazon operation. Please make sure it exists and write '
               'access is granted. Then, restart Pymazon.' % pymazon_dir)
        sys.exit(0)
        

def _get_default_toolkit():
    # try qt first, then gtk, if neither, use qt
    try:
        import PyQt4
        toolkit = 'qt4'
    except ImportError:
        try:
            import pygtk
            toolkit = 'gtk'
        except ImportError:
            toolkit = 'qt4'
    return toolkit

        
#------------------------------------------------------------------------------
# Program-wide settings Class 
#------------------------------------------------------------------------------

class PymazonSettingsError(Exception):
    pass


class _PymazonSettings(object):   
    def __init__(self):
        # some sensible defaults
        self._pymazon_dir = _get_pymazon_dir()
        self._config_file =  os.path.join(self._pymazon_dir, 'pymazonrc')
        self._save_dir = os.environ.get('XDG_MUSIC_DIR', os.getcwd())
        self._name_template = '${tracknum} - ${title}'
        self._toolkit = None # will be figured out in read/write_config_file
        self._num_threads = 1 
        self._read_config_file()
        
    def _get_pymazon_dir(self):
        return self._pymazon_dir
    
    pymazon_dir = property(_get_pymazon_dir)
    
    def _get_save_dir(self):
        return self._save_dir
    
    def _set_save_dir(self, value):
        self._save_dir = os.path.expanduser(value)
        
    save_dir = property(_get_save_dir, _set_save_dir)
    
    def _get_name_template(self):
        return self._name_template
    
    def _set_name_template(self, value):
        template_string = str(value).lstrip(os.path.sep)
        self._name_template = template_string
        
    name_template = property(_get_name_template, _set_name_template)
    
    def _get_toolkit(self):
        return self._toolkit
    
    def _set_toolkit(self, value):
        if value not in _toolkits:
            msg = 'Invalid toolkit specified. '
            msg += 'Valid toolkit identifiers are: %s ' % _toolkits.keys()
            msg += 'Reverting to toolkit: %s' % self.toolkit
            warnings.warn(msg)
            return
        self._toolkit = _toolkits[value]
        
    toolkit = property(_get_toolkit, _set_toolkit)
    
    def _get_num_threads(self):
        return self._num_threads
    
    def _set_num_threads(self, val):
        val = int(val)
        val = max(val, 1)
        val = min(val, 5)
        self._num_threads = val
            
    num_threads = property(_get_num_threads, _set_num_threads)
    
    def write_config_file(self):
        # don't do unnecessary imports if the user already has a config file.
        if self._toolkit is None:
            self._toolkit = _get_default_toolkit()
        cp = '''\
[config]
save_dir = %s
name_template = %s
num_threads = %s
toolkit = %s''' % (self.save_dir, self.name_template, 
                   self.num_threads, self.toolkit)
        f = open(self._config_file, 'w')
        f.write(cp)
        f.close()
    
    def _read_config_file(self):        
        if os.path.exists(self._config_file):
            cp = ConfigParser()
            cp.read(self._config_file)
            try:
                templ = cp.get('config', 'name_template')
                self.name_template = templ
            except (NoOptionError, NoSectionError):
                pass
            try:
                toolkit = cp.get('config', 'toolkit')
                self.toolkit = toolkit
            except (NoOptionError, NoSectionError):
                pass
            try:
                save_dir = cp.get('config', 'save_dir')
                self.save_dir = save_dir
            except (NoOptionError, NoSectionError):
                pass
            try:
                amz_dir = cp.get('config', 'amz_dir')
                settings.amz_dir = amz_dir
            except (NoOptionError, NoSectionError):
                pass
            try:
                num_threads = cp.get('config', 'num_threads')
                self.num_threads = num_threads
            except (NoOptionError, NoSectionError):
                pass
        else:
            self.write_config_file()     
        
settings = _PymazonSettings()


