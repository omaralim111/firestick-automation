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
import string

from pymazon.core.settings import settings
from pymazon.util.image import url_image_cache


class Container(object):
    def __init__(self):
        self.children = []
        
    def __iter__(self):
        return iter(self.children)
    

class HasStatus(Container):
    def __init__(self):
        super(HasStatus, self).__init__()
        self._progress = 0
        self._message = 'Ready'
        
    def _get_status(self):
        # if we have children, the status depends on them        
        if self.children:            
            progress = [child.status[0] for child in self.children]
            total_progress = sum([prog for prog in progress if prog != -1])
            n = len(self.children)
            perc = total_progress / n
            if perc == 100:
                return (perc, 'Complete!')
            return (perc, '%s%%' % perc)
        else:
            return (self._progress, self._message)
        
    def _set_status(self, value):
        self._progress, self._message = value
        
    status = property(_get_status, _set_status)
    
    
class Downloadable(HasStatus):
    def __init__(self, url=None):
        super(Downloadable, self).__init__()
        self.url = url       
        
    def _safe_save_name(self, fname):        
        # ensure we don't overwrite any files
        i = 1
        nfname = fname
        root, ext = os.path.splitext(fname)
        while True:
            if not os.path.isfile(nfname):
                break            
            nfname = ''.join([root, '(', str(i), ')', ext])
            i += 1
        return nfname
    
    def save(self , fname, data=None):
        fname = self._safe_save_name(fname)
        if not data:
            raise IOError('No data to save.')       
        # make any intermediate directories
        head, tail = os.path.split(fname)
        if not os.path.exists(head):
            os.makedirs(head)
        f = open(fname, 'wb')
        f.write(data)
        f.close()   
    
    
class Album(HasStatus):
    def __init__(self, title=None, artist=None, image_url=None, tracks=None):
        super(Album, self).__init__()
        self.title = title
        self.artist = artist        
        self.tracks = tracks or [] 
        self.image_url = image_url
        self.children = self.tracks
        
    @property
    def image(self):
        return url_image_cache.get(self.image_url)
    
    def __unicode__(self):
        return self.artist + ' - ' + self.title
    

class Track(Downloadable):
    def __init__(self, album, title=None, number=None, url=None, filesize=None, 
                 extension=None):
        super(Track, self).__init__(url=url)
        self.album = album
        self.title = title
        self.number = number        
        self.filesize = filesize
        self.extension = extension
        
    def save(self, data):        
        if not os.access(settings.save_dir, os.W_OK):
            raise IOError('No write access to save dir.')      
        
        template = string.Template(settings.name_template)
        sn = template.safe_substitute(artist=self.album.artist, 
                                      title=self.title,
                                      tracknum=self.number,
                                      album=self.album.title)        
        
        save_path = os.path.join(settings.save_dir, sn + '.' + self.extension)        
        super(Track, self).save(save_path, data)        
    
    def __unicode__(self):
        # for display purposes
        return self.number + '. ' + self.title  
    
    
class OtherMedia(HasStatus):
    def __init__(self, tracks=[]):
        self.tracks = tracks
        self.children = self.tracks
    
    def __unicode__(self):
        return 'Other Media'
    
