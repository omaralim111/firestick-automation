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

from PyQt4.QtCore import QAbstractItemModel, QModelIndex, Qt

from pymazon.core.tree_model import TreeModel as PymazonTreeModel
    
    
class TreeModel(QAbstractItemModel, PymazonTreeModel):
    
    headers = ['Title', 'Status']
    
    def __init__(self, parser_list):        
        QAbstractItemModel.__init__(self)
        self.create_tree(parser_list)              
        
    def index(self, row, column, parent):
        if not parent.isValid():
            return self.createIndex(row, column, self.root_nodes[row])
        parent_node = parent.internalPointer()
        return self.createIndex(row, column, parent_node.subnodes[row])

    def parent(self, index):
        if not index.isValid():
            return QModelIndex()
        node = index.internalPointer()
        if node.parent is None:
            return QModelIndex()
        else:
            return self.createIndex(node.parent.row, 0, node.parent)

    def reset(self):
        self.root_nodes = self._get_root_nodes()
        QAbstractItemModel.reset(self)

    def rowCount(self, parent):
        if not parent.isValid():
            return len(self.root_nodes)
        node = parent.internalPointer()
        return len(node.subnodes) 
    
    def columnCount(self, parent):
        return 2    
    
    def data(self, index, role):
        if not index.isValid():
            return None
        node = index.internalPointer()
        if role == Qt.DisplayRole:
            if index.column() == 0:            
                return node.elem.data()
            if index.column() == 1:
                return node.elem.status()        
        return None   
    
    def headerData(self, section, orientation, role):
        if orientation == Qt.Horizontal and role == Qt.DisplayRole:
            return self.headers[section]
        return None
    
    def update_node(self, node):        
        idx = self.createIndex(node.row, 1, node)
        self.dataChanged.emit(idx, idx)        
        