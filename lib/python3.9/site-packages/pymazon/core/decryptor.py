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

import base64
from cStringIO import StringIO
import sys

try:
    from Crypto.Cipher import DES
except ImportError:
    msg = 'Unable to import PyCrypto. Make sure it is installed '
    msg += 'and available on the Python path.'
    print msg
    sys.exit()
    
from pymazon.util.log_util import PymazonLogger


logger = PymazonLogger('decryptor')


class DecryptException(Exception):
    pass


class AmzDecryptor(object):
    ''' Decrypt a base64 encoded & DES-CBC encrypted Amazon .amz file.
    
    Usage:
    
    >>> encrypted_txt = open('./amzfile.amz').read()
    >>> decryptor = AmzDecryptor()
    >>> valid_xml = decryptor.decrypt(encrypted_txt)
    
    '''
    
    # the key and initial value for the .amz DES CBC encryption
    KEY = '\x29\xAB\x9D\x18\xB2\x44\x9E\x31'
    IV = '\x5E\x72\xD7\x9A\x11\xB3\x4F\xEE'    
        
    def strip_trailing_bytes(self, amz_xml):
        # strips the octet even block padding to make
        # a valid xml file. Amazon doesn't follow the padding 
        # convention specified by the DES standard, 
        # rather, they add 8 -'\x08' bytes, and enough 
        # '\x00' bytes to make the file mod(8) even. 
        # The problem is that we don't know how many '\x00' bytes 
        # there are. So, we strip bytes until we reach the first
        # xml closing tag.
        buf = StringIO(amz_xml)
        n = -1
        buf.seek(n, 2)
        while buf.read(1) != '>':
            n -= 1
            buf.seek(n, 2)    
        buf.truncate()    
        return buf.getvalue()

    def decrypt(self, amz_data):
        try:
            d_obj = DES.new(self.KEY, DES.MODE_CBC, self.IV)
            cipher = base64.b64decode(amz_data)
            xml = d_obj.decrypt(cipher)
            valid_xml = self.strip_trailing_bytes(xml)
            return valid_xml
        except Exception, e:
            de = DecryptException(e)
            logger.error(de)
            raise de
    
