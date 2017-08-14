#!/usr/bin/env python
###############################################################################
#
# __make_otu_table__.py - Collate OTU output for GraftM or CommunityM
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

__author__ = "Josh Daly"
__copyright__ = "Copyright 2015"
__credits__ = ["Josh Daly"]
__license__ = "GPL3"
__version__ = "0.0.1"
__maintainer__ = "Josh Daly"
__email__ = "joshua.daly@uqconnect.edu.au"
__status__ = "Development"

###############################################################################

# system imports
import argparse
import sys
from multiprocessing import Pool
from subprocess import Popen, PIPE
import glob
import os
import datetime
import fnmatch

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class GraftMOTU(object):
    def __init__(self,
                 graftmDir,
                 outFile):
        self.graftmData         = {}
        self.sampleTotalCounts  = {}
        self.samples            = []
    
        # run scripts
        self.wrapper(graftmDir, outFile)
    
    def wrapper(self,
                graftmDir,
                outFile):
        
        # loop througt graftm dir, grabbing counts files
        self.grabGraftMData(graftmDir)
        
        # turn counts into relative abundance
        self.transformData()
        
        # print out GraftM data in OTU format
        self.writeDataToFile(outFile)
    
    def grabGraftMData(self,
                             graftmDir):
        countsFileSearch    = os.path.join(graftmDir, '*/*/*_count_table.txt')
        countsFiles         = glob.glob(countsFileSearch)
 
        for f in countsFiles:
            
            sample = f.split("/")[-2].split("_")[0]
            
            # add to samples dict
            self.samples.append(sample)
            
            # grab sample data
            with open(f) as fh:
                for l in fh:
                    if '#' not in l:
                        
                        tabs        = l.rstrip().split("\t")
                        ID          = tabs[0]
                        count       = int(tabs[1])
                        taxString   = tabs[2]
                        
                        self.addData(sample,
                                     taxString,
                                     count)
                        
    def addData(self,
                sample,
                taxString,
                count):
        try:
            self.graftmData[taxString][sample] = count
        except KeyError:
            self.graftmData[taxString] = {sample:count}
            
    def transformData(self, ):
        # transform counts data into relative abundance
        # loop through taxa, samples
        for taxa in self.graftmData.keys():
            for sample in self.graftmData[taxa]:
                count = int(self.graftmData[taxa][sample])
                self.addCounts(sample,
                               count
                               )
            
    def addCounts(self,
                  sample,
                  count):
        try:
            self.sampleTotalCounts[sample] += count
        except KeyError:
            self.sampleTotalCounts[sample] = count
    
    def writeDataToFile(self,
                        outFile):
        
        now = datetime.datetime.now().strftime("%Y_%m_%d")
        
        of = open(outFile, 'w')
        
        # write header
        str_to_print = 'taxa'
        for sample in self.samples:
            str_to_print += "\t%s" % sample
        of.write("%s\n" % str_to_print)
        
        for tax in self.graftmData.keys():
            
            str_to_print = ''
            
            str_to_print += tax
            
            for sample in self.samples:
                
                try:
                    count = self.graftmData[tax][sample]
                    #str_to_print += "\t%f" % (count/float(self.sampleTotalCounts[sample]))
                    str_to_print += "\t%d" % (count)
                
                except KeyError:
                    str_to_print += "\t0"
            
            of.write("%s\n" % str_to_print)
    
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

class CommunityMOTU(object):
    def __init__(self,
                 communitym_dir,
                 out_file,
                 extension
                 ):
        
        self.communitymData     = {}
        self.sampleTotalCounts  = {}
        self.samples            = []
        
        # run scripts
        self.wrapper(communitym_dir,out_file, extension)
    
    def wrapper(self,
                communitym_dir,
                out_file,
                extension):
        
        # loop through communitym directory, grabbing sample + counts
        self.parseCommunityMFiles(communitym_dir, extension)
        
        # combine communitym data by taxonomy, collating counts
        self.collateCounts()
    
        # write combine communitym output to file
        self.writeDataToFile(out_file)
    
    def parseCommunityMFiles(self,
                             communitym_dir,
                             extension):
        
        # change this depending on directory structure
        communitym_otu_tables = self.find('*%s' % extension,communitym_dir)
        
        for otu_table in communitym_otu_tables:
            
            sample = "_".join(otu_table.split("/")[-1].split("_")[0:2])
            
            # add to samples dict
            self.samples.append(sample)
            
            # grab sample data
            with open(otu_table) as fh:
                for l in fh:
                    if '#' not in l:
                        
                        tabs        = l.rstrip().split("\t")
                        ID          = tabs[0]
                        count       = int(tabs[1])
                        taxString   = tabs[2].replace("unresolved_by_lca","")
                        
                        self.addData(sample,
                                     taxString,
                                     count)

    def find(self,
             pattern,
             path):
        result = []
        
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root, name))
        return result
    
    def addData(self,
                sample,
                taxString,
                count):
        try:
            self.communitymData[taxString][sample] = count
        except KeyError:
            self.communitymData[taxString] = {sample:count}
    
    def collateCounts(self,):
        
        for taxa in self.communitymData.keys():
            
            for sample in self.communitymData[taxa]:
                
                count = int(self.communitymData[taxa][sample])
                
                self.addCounts(sample, count)
    
    def addCounts(self,
                  sample,
                  count):
        try:
            self.sampleTotalCounts[sample] += count
        except KeyError:
            self.sampleTotalCounts[sample] = count       

    def writeDataToFile(self,
                        outFile):
        
        of = open(outFile, 'w')
        
        # write header
        str_to_print = 'taxa'
        for sample in self.samples:
            str_to_print += "\t%s" % sample
        of.write("%s\n" % str_to_print)
        
        for tax in self.communitymData.keys():
            
            str_to_print = ''
            
            str_to_print += tax
            
            for sample in self.samples:
                
                try:
                    count = self.communitymData[tax][sample]
                    str_to_print += "\t%d" % (count)
                
                except KeyError:
                    str_to_print += "\t0"
            
            of.write("%s\n" % str_to_print)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def doWork( args ):
    """ Main wrapper"""
    
    if args.subparser_name == "graftm":
        GraftMOTU(args.graftmDir,
                    args.outFile)
    elif args.subparser_name == "communitym":
        CommunityMOTU(args.communitym_dir,
                      args.out_file,
                      args.extension)
                

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')
    
    # place holder for subparsing
    graftm_parser = subparsers.add_parser('graftm',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        help='Combine GraftM OTUs',
                                        description='Combine GraftM OTUs')
    
    graftm_parser.add_argument('graftmDir', help="Directory containig GraftM output files")
    graftm_parser.add_argument('-o','--outFile', default= 'graftmOTU.tsv',help="GraftM OTU table. default= graftmOTU.tsv")

    # place holder for subparsing
    communitym_parser = subparsers.add_parser('communitym',
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                        help='Combine CommunityM OTUs',
                                        description='Combine CommunityM OTUs')
    
    communitym_parser.add_argument('communitym_dir', help="Directory containig CommunityM output files")
    communitym_parser.add_argument('-o','--out_file', default= 'communitym_otu_table.tsv',help="CommunityM OTU table. default= communitym_otu_table.tsv")
    communitym_parser.add_argument('-x','--extension', default= 'communitym_otu_table.tab',help="CommunityM OTU table extension. default= communitym_otu_table.tab")

    # parse the arguments
    args = parser.parse_args()

    # do what we came here to do
    doWork(args)

###############################################################################
###############################################################################
###############################################################################
###############################################################################
