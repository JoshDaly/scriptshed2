#!/usr/bin/env python
###############################################################################
#
# __checkSequencingIndex__.py 
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
import sys
import os
import fnmatch
import json

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class IndexChecker(object):
    def __init__(self,
                 index_file):
        self.indices = []
        self.bad_indexes = {}
        
        self.parseIndexFile(index_file)
        self.compareIndexes()
        self.printPassedIndexes()
        
    
    def parseIndexFile(self,
                       index_file):
        with open(index_file) as fh:
            for l in fh:
                index = l.rstrip()
                self.indices.append(index)
                
    def compareIndexes(self,):
        for i in range(len(self.indices)):
            for j in range(i+1, len(self.indices)):
                index1 = self.indices[i]
                index2 = self.indices[j]
                mismatches = self.calculateMismatches(index1,index2)

                if mismatches < 3:
                    for i in range(1000):
                        if i not in self.bad_indexes:
                            self.bad_indexes[i] = [index1, index2]
                            break
            
    def calculateMismatches(self,
                            index1,
                            index2):
        mismatches = 0
        for i,nucleotide in enumerate(index1):
            if nucleotide != index2[i]:
                mismatches += 1
        return mismatches

    def printPassedIndexes(self,):
        for i in self.bad_indexes.keys():
            index1 = self.bad_indexes[i][0]
            index2 = self.bad_indexes[i][1]
            
            print "\t".join([index1,index2])

###############################################################################
###############################################################################
###############################################################################
###############################################################################
   
def doWork( args ):
    """ Main wrapper"""
    
    IndexChecker(args.index_file)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    
    parser.add_argument('index_file',help='File containing list of indexes to compare')
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################