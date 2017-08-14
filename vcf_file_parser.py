#!/usr/bin/env python
###############################################################################
#
# __bcf_file_parser__.py 
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
from pysam import VariantFile

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class BCF_PARSER(object):
    def __init__(self,
                 bcf_directory,
                 contig_id,
                 contig_start,
                 contig_end,
                 file_prefix):
        self.contig_id = contig_id
        self.contig_start = contig_start
        self.contig_end = contig_end
        self.file_prefix = file_prefix
        
        self.parse_bcf_files(bcf_directory)  
 
    def parse_bcf_files(self,
                        bcf_directory):
        bcf_files = Find().find("%s*.vcf" % self.file_prefix, bcf_directory)
        
        for bcf_file in bcf_files:
            self.filter_bcf_file(bcf_file)
            print "finished file %s" % bcf_file
            
    def filter_bcf_file(self,
                        bcf_file):
        bcf_in = VariantFile(bcf_file,'rb')
        bcf_out = VariantFile("%s.target.vcf" % bcf_file[:-4],'w',header=bcf_in.header)
        for rec in bcf_in.fetch():
            if rec.contig == self.contig_id:
                if self.contig_start == False and self.contig_end == False:
                    pass
                else:
                    if rec.pos >= self.contig_start and rec.pos <= self.contig_end:
                        bcf_out.write(rec)
        bcf_in.close()
        bcf_out.close()
        
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
   
def doWork( args ):
    """ Main wrapper"""
    
    # run scripts
    BCF_PARSER(args.bcf_directory,
               args.contig_id,
               args.contig_start,
               args.contig_end,
               args.file_prefix)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('bcf_directory',help='Directory containing bcf files to parse')
    parser.add_argument('contig_id',type=str,help='Target contig')
    parser.add_argument('-contig_start','--contig_start',type=int,default=False,help='Set region to call SNPs from (Start). Default is whole contig')
    parser.add_argument('-contig_end','--contig_end',type=int,default=False,help='Set region to call SNPs from (End). Default is whole contig')
    parser.add_argument('-fp','--file_prefix',default=False,help='BCF filename prefix')
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################