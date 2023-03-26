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
import sys

from pymazon.core.downloader import Downloader
from pymazon.core.item_model import Track
from pymazon.core.parser import AmzParser
from pymazon.core.settings import settings
from pymazon.core.tree_model import TreeModel


pymazon_splash = """
##############################################
#                                            #
#                PYMAZON!                    # 
#                                            #
##############################################
"""


pymazon_finished = """
##############################################

Files downloaded to:
{0}

Name format used:
{1}

##############################################
#                                            #
#               Finished!                    #
#                                            #
##############################################
"""


class CommandLineRunner(object):
    def __init__(self, amzs):
        self.tree_model = None
        self.downloader = None        
        self.current_node = None
        self.tab = '    '
        self.load_new_amz_files(amzs)        
        
    def run(self):
        for line in self.get_printable_download_list():
            print line
        print
        print 'Download these items? (Yes/No)'
        res = raw_input()
        if res in ['Yes', 'Y', 'y', 'yes', 'YES']:
            self.download_tracks()
        
    def load_new_amz_files(self, amz_files):
        parser = AmzParser()
        for f in amz_files:
            parser.parse(f)
        objects = parser.get_parsed_objects()        
        self.tree_model = TreeModel()
        self.tree_model.create_tree(objects)
        
    def finished_cb(self):
        finished = pymazon_finished.format(settings.save_dir, 
                                           settings.name_template)
        print
        print finished
    
    def update_cb(self, node):
        if isinstance(node.elem.obj, Track):
            artist = node.parent.elem.data()
            title = node.elem.data()
            status = node.elem.status()[1]
            if not node is self.current_node:
                self.current_node = node                
                print                
            else:
                txt = u'\r{0} {1} - {2}'.format(artist, title, status).encode('UTF-8', 'ignore')                
                sys.stdout.write(txt)
                sys.stdout.flush()                
            
    def download_tracks(self):
        settings.num_threads = 1 # until this printer becomes more sophisticated
        if not self.tree_model:
            return
        save_dir = settings.save_dir        
        if not os.access(save_dir, os.W_OK):
            msg = 'No write access to save directory. '
            msg += 'Choose a new directory with write privileges.'
            raise IOError(msg)        
        
        self.downloader = Downloader(self.tree_model, self.update_cb, 
                                     self.finished_cb)
        print 
        print '##############################################'
        self.downloader.start()
        self.downloader.join()     

    def get_printable_download_list(self):
        if not self.tree_model:
            return
        root_nodes = self.tree_model.get_root_nodes()
        lines = []
        def add_lines(node, ntabs):
            txt = node.elem.data()
            lines.append(ntabs * self.tab + txt)
            for snode in node.subnodes:
                add_lines(snode, ntabs+1)
        for node in root_nodes:
            add_lines(node, 0)
        return lines
    
    
def main(amzs):
    print pymazon_splash
    runner = CommandLineRunner(amzs)
    runner.run()
    