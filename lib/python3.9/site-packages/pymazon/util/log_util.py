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

import datetime
import logging
import os

from pymazon.core.settings import settings
    

class SimpleCache(object):    
    def __init__(self, missing_func):
        self.missing = missing_func
        self.cache = {}
        
    def __getitem__(self, item):
        if item not in self.cache:
            self.cache[item] = self.missing(item)
        return self.cache[item]
         
     
class PymazonLogger(object):
    '''A class interface for logging.'''
    cache = SimpleCache(logging.getLogger)
    inited = False    
    
    def __init__(self, who):
        PymazonLogger.setup_logging()
        self.log_handle = self.cache[who]    
     
    def __getattr__(self, attr):
        # on new style classes, this will only be called if the attribute
        # lookup raises an AttributeError exception, meaning the attribute 
        # was not found. In this case, we delegate the work to the 
        # logging handle, making the loggers methods available 
        # to the user
        return getattr(self.log_handle, attr)
    
    @classmethod   
    def setup_logging(cls):
        # only setup once
        if cls.inited:
            return
         
        # setup logging
        cls.log_dir = os.path.join(settings.pymazon_dir, 'logs')
        if not os.path.exists(cls.log_dir):
            os.makedirs(cls.log_dir)
        
        # this name filename should be unique unless you start
        # multiple instances of pymazon at exactly the same second
        # which is highly unlikely.
        today = datetime.datetime.today().strftime('pymazon_error_log_%d_%m_%y_%H_%M_%S')
        cls.log_filename = os.path.join(cls.log_dir, today + '.txt')
                        
        fmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        logging.basicConfig(filename=cls.log_filename, format=fmt)
        
        cls.inited = True
    
    @classmethod
    def shutdown(cls):
        # dont shutdown if not setup
        if not cls.inited:
            return
            
        # remove a log with no info
        logging.shutdown()
        f = open(cls.log_filename, 'r')
        n = len(f.read())
        f.close()
        if n == 0:
            os.remove(cls.log_filename)
            
        cls.inited = False
