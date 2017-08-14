#!/usr/bin/env python
###############################################################################
#
# __calculateGenomeAbundanceInMetagenome__.py 
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
from multiprocessing import Pool
from subprocess import Popen, PIPE
import datatime

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class GenomeAbundance(object):
    
    def __init__(self,
                 bin_dir,
                 bamm_dir,
                 raw_data_dir):
        # data
        
        # functions
        self.Find = Find()
        
        self.parsePopulationBins(bin_dir)
        self.grabReadMappingForEachGenome(bamm_dir)
        self.calculateMetagenomeFraction(raw_data_dir)

    def parsePopulationBins(self,
                            bin_dir,
                            bin_extension):
        population_bins = self.Find("*%s" % bin_extension,bin_dir)
    
    def grabReadMappingForEachGenome(self,
                                     bamm_dir):
        pass
    
    def calculateMetagenomeFraction(self,
                                    raw_data_dir):
        pass
    
    
class Find(object):
    def find(self,
             pattern,
             path):
        result = []
        
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root, name))
        return result
###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand( cmd ):
    """Run a command and take care of stdout

expects 'cmd' to be a string like "foo -b ar"

returns (stdout, stderr)
"""
    print cmd
    args = shlex.split(cmd)
    p = subprocess.Popen(args) # shell=bash is not recommended. Only use when '>' must be in cmd. 
    return p.communicate()
    #p = Popen(cmd.split(' '), stdout=PIPE)
    #return p.communicate()

def doWork( args ):
    """ Main wrapper"""
    GenomeAbundance()

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    #------------------------------
    
    parser.add_argument('bin_dir',help='Directory containing population bins')
    parser.add_argument('bamm_dir',help='Directory containing bamm parse counts coverage files')
    parser.add_argument('raw_data_dir',help='File containing raw data in fastq.gz format')
    parser.add_argument('bin_extension',default="fa",help='Population bin fasta file extension')
    parser.add_argument('bamm_extension',default=,help='Bamm counts coverage file extension')

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
