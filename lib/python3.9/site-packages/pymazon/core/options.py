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
from optparse import OptionParser

from pymazon.core.settings import settings


def validate_args(args):
    # validates that the only positional args are valid .amz files
    for arg in args:
        if not arg.endswith('.amz'):
            raise ValueError('Positional arguments must be .amz files.')        
        if not os.path.isfile(arg):
            raise ValueError('AMZ file %s does not exist. Please check the path '
                             'and try again.' % arg)       
    return args


def parse_options():    
    parser = OptionParser()
    parser.add_option('-d', '--dir', dest='save_dir',
                      help='Directory in which to save the downloads. '
                           'Defaults to current working directory.')
    parser.add_option('-c', '--command-line', dest='command_line',
                      help='Run Pymazon as a command line script. '
                           'If this option is used, it must be accompanied '
                           'by a path to one or more .amz files as a positional '
                           'argument.',
                      action='store_true', default=False)
    parser.add_option('-s', '--save-name-format', dest='fmt', help=\
'''A format string to tell Pymazon how to format the mp3 file names.
The default format is '${tracknum} - ${title}' which will generate a file
of the name '1 - foosong.mp3' for example. Valid tags are:

    ${tracknum}
    ${title}
    ${artist}
    ${album}

and can appear in any order.

The format given on the command line MUST be surrounded in single quotes:
$ pymazon --save-name-format='${tracknum} - ${artist} - ${title}'

or

$ pymazon -s '${tracknum} - ${artist} - ${title}'

This formatting command is flexible. For example, say you have a music folder
at ~/Music and you want to download the new Ludovico Einaudi album Nightbook.
And, let's say that you want your directory structure to look like this at the
end: ~/Music/Ludovico Einaudi/Nightbook/1 - song.mp3
Then, assuming ~/foo.amz is the path to the AMZ file, use the following
command:

As a command line script:
$ pymazon -c -d ~/Music -s '${artist}/${album}/${tracknum} - ${title}' ~/foo.amz

Or launch the gui ready-to-go with:
$ pymazon -d ~/Music -s '${artist}/${album}/${tracknum} - ${title}' ~/foo.amz

''', default='')
    parser.add_option('-g', '--gui-toolkit', dest='toolkit',
                      help='Specify which gui toolkit to use. '
                           'Valid options are: "qt4" and "gtk". '
                           'Default is qt4.')
    parser.add_option('-n', '--num-threads', dest='num_threads',
                      help='Specify the number of threads for use in the '
                           'downloader. i.e. the number of simultaneous '
                           'downloads. When using the command line interfaces, '
                           'the number of threads is always 1.')
   
    options, args = parser.parse_args()    
    
    if options.save_dir:
        settings.save_dir = options.save_dir        
    if options.fmt:
        settings.name_template = options.fmt
    if options.toolkit:
        settings.toolkit = options.toolkit
    if options.num_threads:
        settings.num_threads = options.num_threads

    amzs = validate_args(args)
    
    if options.command_line:
        entry = 'cmd'
    else:
        entry = settings.toolkit
        
    return (entry, amzs)
