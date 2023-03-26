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

import urllib2


class UrlImageCache:
    '''A simple image cache for images residing at a url. 
    Extremely simple to use. 
    
    cache = ImageCache()
    image = cache.get(url) # url is valid text url (http://www.foo.com/foo.jpg)
    
    The cache has no size limit and image live in the cache forever. Thus, 
    this shouldn't be used for huge images if you value your memory. It is 
    intended to cache small album art thumbnails.
    
    '''    
    def __init__(self):
        self.cache = {}

    def _download(self, url):        
        handle = urllib2.urlopen(url)
        data = handle.read()
        handle.close()
        return data

    def _factory(self, url):
        try:
            pixbuf = self._download(url)
        except urllib2.URLError:
            pixbuf = ''
        return pixbuf
            
    def get(self, url):
        return self.cache.setdefault(url, self._factory(url))      
        
url_image_cache = UrlImageCache()