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


class TreeNode(object):
    def __init__(self, elem, parent, row):
        self.elem = elem
        self.parent = parent
        self.row = row
        self.subnodes = self._get_children()

    def _get_children(self):
        return [TreeNode(elem, self, idx) for idx, elem 
                in enumerate(self.elem.subelements)]

    def __getitem__(self, n):
        if not type(n) == int:
            raise ValueError('Can only index with integer values')
        if n < 0 or n >= len(self.subnodes):
            raise ValueError('Index out of Range')
        return self.subnodes[n] 
    

class TreeElement(object):
    def __init__(self, obj, parent):
        self.obj = obj
        self.parent = parent
        self.subelements = self._get_children()    
        
    def _get_children(self):
        return [TreeElement(child, self) for child in self.obj]
            
    def data(self):
        return unicode(self.obj)
    
    def status(self):
        return self.obj.status
    
    
class TreeModel(object):       
    def create_tree(self, parser_list):
        self.root_elements = self._get_root_elements(parser_list)       
        self.root_nodes = self._get_root_nodes()
        
    def _get_root_elements(self, parser_list):
        return [TreeElement(obj, None) for obj in parser_list]
        
    def _get_root_nodes(self):
        return [TreeNode(elem, None, idx) 
                for idx, elem in enumerate(self.root_elements)]    
    
    def get_root_elements(self):
        return self.root_elements
    
    def get_root_nodes(self):
        return self.root_nodes
    
    def flatten(self, node):
        if node.subnodes:
            for snode in node.subnodes:
                for ssnode in self.flatten(snode):
                    yield ssnode
        yield node                   
    
    def filter_nodes(self, filter_func):
        filtered_nodes = []
        for root_node in self.root_nodes:            
            passing = [node for node in self.flatten(root_node) 
                       if filter_func(node)]
            filtered_nodes.extend(passing)        
        return filtered_nodes
        
        
        
