"""
Pymazon - A Python based downloader for the Amazon.com MP3 store
Copyright (c) 2009 Steven C. Colbert

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

load_icon_path = os.path.join(*(os.path.split(__file__)[:-1] + ('icons', 'amz-load.png')))
download_icon_path = os.path.join(*(os.path.split(__file__)[:-1] + ('icons', 'download.png')))
exit_icon_path = os.path.join(*(os.path.split(__file__)[:-1] + ('icons', 'pymazon-exit.png')))
settings_icon_path = os.path.join(*(os.path.split(__file__)[:-1] + ('icons', 'pymazon-settings.png')))
python_icon_path = os.path.join(*(os.path.split(__file__)[:-1] + ('icons', 'python-powered.png')))
show_icon_path = os.path.join(*(os.path.split(__file__)[:-1] + ('icons', 'show-downloads.png')))
