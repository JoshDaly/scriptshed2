#!/usr/bin/env python
###############################################################################
#
# __simpleCheckmOutputParser__.py 
#
###############################################################################
# #
# This program is free software: you can redistribute it and/or modify #
# it under the terms of the GNU General Public License as published by #
# the Free Software Foundation, either version 3 of the License, or #
# (at your option) any later version. #
# #
# This program is distributed in the hope that it will be useful, #
# but WITHOUT ANY WARRANTY; without even the implied warranty of #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the #
# GNU General Public License for more details. #
# #
# You should have received a copy of the GNU General Public License #
# along with this program. If not, see <http://www.gnu.org/licenses/>. #
# #
###############################################################################

__author__ = "Joshua Daly"
__copyright__ = "Copyright 2017"
__credits__ = ["Joshua Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Joshua Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Development"

###############################################################################
# system imports

import argparse
import math
from multiprocessing import Pool
from subprocess import Popen, PIPE
import subprocess
#from Bio import SeqIO
import re
import os
import fnmatch
import glob
from datetime import datetime
import shutil

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class CheckmParser(object):
    def __init__(self,
                 checkm_file):
        self.parseCheckmFile(checkm_file)
        
    def parseCheckmFile(self,
                        checkm_file):
        with open(checkm_file) as fh:
            line_number = 0
            for l in fh:
                tabs = l.rstrip().split("\t")
                if line_number == 0:
                    print l.rstrip()
                else:
                    col_1  = tabs[0]
                    col_2  = tabs[1]
                    col_3  = tabs[2]
                    col_4  = tabs[3]
                    col_5  = tabs[4]
                    col_6  = tabs[5]
                    col_7  = tabs[6]
                    col_8  = tabs[7]
                    col_9  = tabs[8]
                    col_10 = tabs[9]
                    col_11 = tabs[10]
                    col_12 = tabs[11] # completeness
                    col_13 = tabs[12] # contamination
                    col_14 = tabs[13]
                    
                    if float(col_12) - (4 * float(col_13)) >=60:
                        print "\t".join(tabs)
                line_number += 1
                
###############################################################################
###############################################################################
###############################################################################
###############################################################################
def runCommand(cmd):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    args = cmd
    p = subprocess.Popen(args, shell=True, stdout=PIPE, stderr=PIPE) # shell=bash is not recommended. Only use when '>' must be in cmd.
    print cmd
    return p.communicate()
 
def doWork( args ):
    """ Main wrapper"""
    # run scripts
    CheckmParser(args.checkm_file)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    #------------------------------
    
    parser.add_argument('checkm_file',help='Checkm output file in tab format')
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################