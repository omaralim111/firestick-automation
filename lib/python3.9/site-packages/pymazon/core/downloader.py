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

from collections import deque
from cStringIO import StringIO
import threading
import urllib2

from pymazon.core.item_model import Downloadable
from pymazon.core.settings import settings
from pymazon.core.tree_model import TreeModel
from pymazon.util.log_util import PymazonLogger


logger = PymazonLogger('downloader')


class _DownloadWorker(threading.Thread):    
    def __init__(self, parent):
        super(_DownloadWorker, self).__init__()
        self.parent = parent        
        self._abort = False
        self._pause = False
        self.node = None
        
        # just in case the amazon servers want to get sneaky and redirect 
        self.redirects = urllib2.HTTPRedirectHandler()
        self.cookies =  urllib2.HTTPCookieProcessor()
        self.opener = urllib2.build_opener(self.redirects, self.cookies)

    def kill(self):        
        self._abort = True   

    def run(self):
        while not self._abort:
            self.node = self.parent.next_node()
            if not self.node:
                break
            self.do_work()
            
    def do_work(self):
        obj = self.node.elem.obj        
        
        handle = self._connect(obj)
        if not handle:
            return 
        
        data = self._download(obj, handle)        
        if not data:
            return
        
        try:
            obj.save(data)
        except IOError, e:
            logger.error('Saving file %s \n %s \n' % (dir(obj), e))
            obj.status = (0, 'Error!')
            self.parent.update(self.node)
            return
        
        obj.status = (100, 'Complete!')
        self.parent.update(self.node)        
        
    def _connect(self, obj):
        request = urllib2.Request(obj.url)        
        obj.status = (0, 'Connecting...')
        self.parent.update(self.node)
        try:
            handle = self.opener.open(request)
        except urllib2.URLError, e:
            logger.error('Opening request at url: %s \n %s \n' % (obj.url, e))
            obj.status = (0, 'Error!')
            self.parent.update(self.node)
            return
        return handle   

    def _download(self, obj, handle):
        fs = int(obj.filesize)
        chunk_size = fs / 100 #1% of file per chunk              
        buf = StringIO()
        perc_complete = 0  
        obj.status = (perc_complete, '%s%%' % perc_complete)
        self.parent.update(self.node)
        try:
            while True:
                if self._abort:
                    buf.close()
                    return None                
                chunk = handle.read(chunk_size)
                if not chunk: # we got the whole file
                    break
                buf.write(chunk)
                perc_complete += 1
                perc_complete = min(100, perc_complete)
                obj.status = (perc_complete, '%s%%' % perc_complete)
                self.parent.update(self.node)                
        # purposely swallow any and all exceptions during downloading so 
        # other tracks can continue
        except Exception, e:
            logger.error('Reading from opened url: %s \n %s\n' % (obj.url, e))
            obj.status = (0, 'Error!')
            self.parent.update(self.node)
            buf.close()
            return           
        data = buf.getvalue()
        buf.close()        
        return data 


class Downloader(threading.Thread):
    def __init__(self, tree, update_cb=lambda a: None, 
                 finished_cb=lambda: None):
        super(Downloader, self).__init__()
        
        if not isinstance(tree, TreeModel):
            raise ValueError('tree must be an instance of TreeModel')
        
        self.tree = tree
        self.update_cb = update_cb
        self.finished_cb = finished_cb   
        self.workers = []
        self.queue = deque(self.get_download_nodes())
                        
        
        # don't launch more threads than we can use
        self.num_threads = min(settings.num_threads, len(self.queue))        
        self._abort = False
        
    def get_download_nodes(self):
        def filter_func(node):
            return isinstance(node.elem.obj, Downloadable)
        return self.tree.filter_nodes(filter_func)
        
    def update(self, node):
        self.update_cb(node)
        if node.parent:
            self.update(node.parent)

    def kill(self):
        if self.workers:
            for worker in self.workers:
                worker.kill()
                worker.join()
        self._abort = True

    def next_node(self):
        try:
            node = self.queue.popleft()
            return node
        except IndexError:
            return None

    def run(self):
        self.workers = [_DownloadWorker(self) for i in range(self.num_threads)]
        for worker in self.workers:
            worker.start()
        for worker in self.workers:
            worker.join()
        if not self._abort:
            self.finished_cb()
