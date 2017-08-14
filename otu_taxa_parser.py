#!/usr/bin/env python
###############################################################################
#
# __otu_taxa_parser__.py - parse otu taxonomy to required level
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

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def parseOTUTable(otu_table,
                  column_header):
    taxa_column = 0
    
    with open(otu_table) as fh:
        line_number = 0 
        for l in fh:
            tabs = l.rstrip().split("\t")
            if line_number == 0:
                print l.rstrip()
                for index,column in enumerate(tabs):
                    if column == column_header:
                        taxa_column = index
                        
            else:
                array_to_print = []
                for index,value in enumerate(tabs):
                    if index == taxa_column:
                        value = parseTaxonomy(value)
                    else:
                        pass
                    array_to_print.append(value)
                print "\t".join(array_to_print)
            line_number += 1

def parseTaxonomy(taxonomy_str):
    taxonomy_ranks = taxonomy_str.split(";")
    taxa_to_return = taxonomy_ranks[-1]
    for index, value in enumerate(taxonomy_ranks):
        if len(value) <4:
            taxa_to_return = taxonomy_ranks[index-1]
            break
    return taxa_to_return

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
    parseOTUTable(args.otu_table,
                  args.column_header)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('otu_table',help='OTU table in tab-delimited format')
    parser.add_argument('-column_header','--column_header',default='Taxon',help='Column name containing taxonomy information')
    
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################