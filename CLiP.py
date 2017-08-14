#!/usr/bin/env python
###############################################################################
#
# __contigLinkagePlot__.py 
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
__copyright__ = "Copyright 2016"
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
import numpy as np
np.seterr(all='raise')
import math as math
import matplotlib.pyplot as plt
import os
import fnmatch
import shlex
from multiprocessing import Pool
from subprocess import Popen, PIPE
import subprocess
import colorsys

# local imports

###############################################################################
###############################################################################
###############################################################################
###############################################################################

class ContigLinkagePlot(object):
    def __init__(self,
                 bamm_dir,
                 assembly_dir,
                 outdir,
                 force,
                 extension,
                 highlight_contigs
                 ):
        # local objects
        self.Find       = FIND()
        self.BamMParser = BamMParser()
    
        # data storage
        self.bamm_files         = {}
        self.plot_data          = {}
        self.contig_len         = {}
        self.coords             = {}
        self.contig_linkage     = {}
        self.gc_content         = {}
        self.highlight_contigs  = {} 
    
        # run scripts
        if highlight_contigs:
            self.parseHighlightContigs(highlight_contigs)
        else:
            self.calculateGC(assembly_dir, extension)
        self.parseBamFiles(bamm_dir)
        self.createLinkagePlots(outdir, force, highlight_contigs)

    def calculateGC(self,
                    assembly_dir,
                    extension):

        self.createGCFiles(assembly_dir, extension)
        self.parseGCFiles(assembly_dir)
        self.removeGCTempFiles(assembly_dir)

    def parseHighlightContigs(self,
                              highlight_contigs):
        with open(highlight_contigs) as fh:
            for l in fh:
                contig_id = l.rstrip()
                self.highlight_contigs[contig_id] = 1

    def createGCFiles(self,
                      assembly_dir,
                      extension):
        assemblies = self.Find.find('*%s' % extension, assembly_dir)
        cmds = []
        
        # create GC files
        for assembly in assemblies:
            cmds.append('geecee -auto -sequence %s -outfile %s.gc.txt' % (assembly, assembly))
        
        for cmd in cmds:
            runCommand(cmd)
    
    def parseGCFiles(self,
                     assembly_dir):
        gc_files = self.Find.find('*.gc.txt', assembly_dir)
            
        for gc_file in gc_files:
            genome_id = gc_file.split("/")[-1].replace(".fa","")
            with open(gc_file) as fh:
                for l in fh:
                    if '#' not in l:
                        tabs = l.rstrip().split()
                        contig_id  = tabs[0]
                        gc_content = float(tabs[1])
                        self.addGC(genome_id,contig_id,gc_content)
    
    def pseudocolor(self,
                    val,
                    minval,
                    maxval):
        # convert val in range minval..maxval to the range 0..120 degrees which
        # correspond to the colors red..green in the HSV colorspace
        h = (float(val-minval) / (maxval-minval)) * 120
        # convert hsv color (h,1,1) to its rgb equivalent
        # note: the hsv_to_rgb() function expects h to be in the range 0..1 not 0..360
        r, g, b = colorsys.hsv_to_rgb(h/360, 1., 1.)
        return tuple([r, g, b])
    
    def addGC(self,
              genome_id,
              contig_id,
              gc_content):
        try:
            self.gc_content[genome_id][contig_id] = gc_content
        except KeyError:
            self.gc_content[genome_id] = {contig_id:gc_content}
    
    def removeGCTempFiles(self,
                          assembly_dir):
        gc_files = self.Find.find("*.gc.txt", assembly_dir)
        cmds = []
        for gc_file in gc_files:
            cmds.append("rm %s" % gc_file)
        for cmd in cmds:
            runCommand(cmd)
       
    def parseBamFiles(self,
                      bamm_dir):
        
        bamm_coverage_files = self.Find.find('*cov*', bamm_dir)
        bamm_links_files    = self.Find.find('*link*', bamm_dir)
        
        # loop through bamm coverage files
        for coverage_file in bamm_coverage_files:
            if "_" in coverage_file:
                genome_id = coverage_file.split("/")[-1].split('_')[0]
            else:
                genome_id = coverage_file.split("/")[-1].split('.')[0]
            try:
                self.bamm_files[genome_id]['coverages'] = coverage_file
            except KeyError:
                self.bamm_files[genome_id] = {'coverages':coverage_file}
        
        # loop through bamm links files
        for links_file in bamm_links_files:
            if "_" in coverage_file:
                genome_id = links_file.split("/")[-1].split('_')[0]
            else:
                genome_id = links_file.split("/")[-1].split('.')[0]
            try:
                self.bamm_files[genome_id]['links'] = links_file
            except KeyError:
                self.bamm_files[genome_id] = {'links':links_file}

    def createLinkagePlots(self,
                          outdir,
                          force,
                          highlight_contigs
                          ):
        """Create individual linkage plots for each population genome"""
        
        # loop through population genomes
        for genome_id in self.gc_content.keys():
            
            sample_id = genome_id.split("_")[0]
            
            coverage_file = self.bamm_files[sample_id]['coverages']
            links_file    = self.bamm_files[sample_id]['links']
        
            self.BamMParser.bammParser(links_file, coverage_file)
            
            # get plot data
            self.buildPlotData(self.BamMParser, genome_id)
            
            # try to run xyscatter
            self.checkXYScatter(outdir, genome_id, self.BamMParser, force, highlight_contigs)
            
    def buildPlotData(self,
                      bamm_data,
                      genome_id):
        """Build plot data object from BamMParser.bamm_links object"""
        # loop through contigs
        for contig1 in bamm_data.bamm_links.keys():
            
            for contig2 in bamm_data.bamm_links[contig1]:
                
                if self.checkContigsPresentInPopulationBin(contig1, contig2, genome_id):
                
                    for link in bamm_data.bamm_links[contig1][contig2]:
                        
                        try:
                            pos = link[2]
                            
                            try:
                                self.plot_data[contig1][pos] += 1 
                            
                            except KeyError:
                                
                                try:
                                    self.plot_data[contig1][pos] = 1
                                
                                except KeyError:
                                    self.plot_data[contig1]= {pos:1}
                        
                        except KeyError:
                            pass
    
    def checkContigsPresentInPopulationBin(self,
                                           contig1,
                                           contig2,
                                           genome_id):
        try:
            gc_1 = self.gc_content[genome_id][contig1]
            gc_2 = self.gc_content[genome_id][contig2]
            return True
        except KeyError:
            return False
    
    def checkXYScatter(self, outdir, gid, bamm_data, force, highlight_contigs):
        try:
            self.xyScatter(outdir, gid, bamm_data, force, highlight_contigs)
        except ZeroDivisionError:
            print 'No coverage data available for %s' % gid
            pass
    
    def xyScatter(self, outdir, gid, bamm_data, force, highlight_contigs):
        
        sample_id = gid.split("_")[0]
        
        # get the longest contig
        longest_contig_len  = 0
        longest_contig      = ''
        most_coverage       = 0
        chromosome_links    = 0
        total_links         = 0
        
        # initialise figure
        fig = plt.figure(figsize=(20,10),dpi=300)
        ax = fig.add_subplot(111)
        
        # initialise coords
        x = [] # contigsize
        y = [] # coverage
        
        # initialise colour array
        colours = []
        
        # contig size vs coverage + linkage
        for contig in self.gc_content[gid].keys():
            
            # colour by GC content
            if highlight_contigs:
                if contig in self.highlight_contigs:
                    colours.append('red')
                else:
                    colours.append('black')
            else:
                gc_content_rgb = self.pseudocolor(self.gc_content[gid][contig], 0.0, 1.0)
                colours.append(gc_content_rgb)
            
            if int(bamm_data.contig_len[contig]) > longest_contig_len:
                longest_contig_len  = int(bamm_data.contig_len[contig])
                longest_contig      = contig
            x.append(int(bamm_data.contig_len[contig]))
            y.append(self.returnAverage(bamm_data.bamm_cov_data[contig]))
            self.coords[contig] = [int(bamm_data.contig_len[contig]), self.returnAverage(bamm_data.bamm_cov_data[contig])]
            
        for i,v in enumerate(y):
            chromosome_coverage = self.returnAverage(bamm_data.bamm_cov_data[longest_contig])
            y[i] = float(v)/float(chromosome_coverage)
            if y[i] > most_coverage:
                most_coverage = y[i]
        
        # Grab longest contigs link count
        #for contig in bamm_data.bamm_links_data[longest_contig]:
        #    links = bamm_data.bamm_links_data[longest_contig][contig]
        #    chromosome_links += links
        
        # loop through contigs, and link with lines weighted by linkage
        for contig1 in bamm_data.bamm_links_data.keys():
            for contig2 in bamm_data.bamm_links_data[contig1]:
                if self.checkContigsPresentInPopulationBin(contig1, contig2, gid):
                    chromosome_coverage = self.returnAverage(bamm_data.bamm_cov_data[longest_contig])
                    total_links += bamm_data.bamm_links_data[contig1][contig2]
                    linkage = math.sqrt(bamm_data.bamm_links_data[contig1][contig2])/25
                    coords1 = self.coords[contig1]
                    coords2 = self.coords[contig2]
                    x2 = []
                    y2 = []
                    x2.append(coords1[0])
                    x2.append(coords2[0])
                    y2.append(coords1[1]/float(chromosome_coverage))
                    y2.append(coords2[1]/float(chromosome_coverage))
                    plt.plot(x2,y2,'k-',lw=linkage, alpha=0.5)
        
        # print out data
        try:
            self.writeStatsDataToFileWrapper(bamm_data, gid, force, outdir, longest_contig, total_links)
        except ZeroDivisionError:
            stats_file = os.path.join(outdir, '%s_cov_links_stats.csv' % (gid))
            print 'No coverage data available for %s' % stats_file
            pass
        
        plt.scatter(x,y, s=200, alpha=0.5, color=colours)
            
        # set axis limits
        #plt.ylim(-0.2, most_coverage+1)
        #plt.xlim(-(longest_contig_len/10), longest_contig_len+(longest_contig_len/10))
        
        # set axis labels
        plt.xlabel('Contig Size (bp)')
        plt.ylabel('Coverage Relative to Largest Contig (%)')
        
        # set axis to logarithmic
        ax.set_xscale('log')
        ax.set_yscale('log')
        
        # output directory/file
        outfile = os.path.join(outdir,"%s.scatter.ylog.png" % gid)
        plt.savefig(outfile,format='png')
        plt.close()
    
    def writeStatsDataToFileWrapper(self, bamm_data, gid, force, outdir, longest_contig, total_links):
        stats_file = os.path.join(outdir, '%s_cov_links_stats.csv' % gid)
    
        if force:
            # print header 
            self.writeStatsDataToFile(gid ,bamm_data, stats_file, longest_contig, total_links)
            
        else:
            if os.path.isfile(stats_file):
                print "File %s already exists, set force=True to overwrite" % (stats_file)
                    
            else:
                self.writeStatsDataToFile(gid, bamm_data, stats_file, longest_contig, total_links)
    
    def writeStatsDataToFile(self,
                             gid,
                             bamm_data,
                             stats_file,
                             longest_contig,
                             total_links):
        chromosome_coverage = self.returnAverage(bamm_data.bamm_cov_data[longest_contig])
        try:
            chromosome_linkage  = self.returnAverage(bamm_data.bamm_links_total[longest_contig])
            # print header 
            f = open(stats_file, 'w')
            f.write("\t".join(['contig',
                               'contig_len',
                               'contig_cov',
                               'contig_links',
                               'contig_links_div_chromo',
                               'contig_links_div_total',
                               'gc_content\n']))
            
            for contig in self.gc_content[gid].keys():
                contig_cov              = str(self.returnAverage(bamm_data.bamm_cov_data[contig])/float(chromosome_coverage))
                contig_len              = str(bamm_data.contig_len[contig])
                contig_links            = '0'
                contig_links_div_total  = '0'
                contig_links_div_chromo = '0'
                gc_content              = self.gc_content[gid][contig]
                try:
                    contig_links            = str(bamm_data.bamm_links_total[contig])
                    contig_links_div_total  = str(bamm_data.bamm_links_total[contig]/float(total_links))
                    contig_links_div_chromo = str(bamm_data.bamm_links_total[contig]/float(chromosome_linkage))
                except KeyError:
                    pass
                f.write("\t".join([contig,
                                   contig_len,
                                   contig_cov,
                                   contig_links,
                                   contig_links_div_chromo,
                                   contig_links_div_total,
                                   "%f\n" % gc_content]))
            
        except KeyError:
            # print header 
            f = open(stats_file, 'w')
            f.write("\t".join(['contig',
                               'contig_len',
                               'contig_cov\n']))
            for contig in bamm_data.bamm_cov_data.keys():
                contig_cov              = str(self.returnAverage(bamm_data.bamm_cov_data[contig])/float(chromosome_coverage))
                contig_len              = str(bamm_data.contig_len[contig])
                f.write("\t".join([contig,
                                   contig_len,
                                   "%s\n"  % contig_cov]))
    
    def returnAverage(self,array):
        a = np.array(array)
        return np.mean(a)
    
class FIND(object):
    def find(self,
             pattern,
             path):
        result = []
        
        for root, dirs, files in os.walk(path):
            for name in files:
                if fnmatch.fnmatch(name, pattern):
                    result.append(os.path.join(root, name))
        return result
    
class BamMParser(object):
    def __init__(self,):
        self.bamm_links         = {}
        self.bamm_links_data    = {}
        self.bamm_links_total   = {}
        self.bamm_cov_data      = {}
        self.contig_len         = {}
        
    def bammParser(self, bamm_links_file, bamm_cov_file):
        with open(bamm_links_file) as fh:
            for l in fh:
                if l[0] != '#':
                    BLFP = BamMLinksFileParser(l)
                    
                    self.buildLinksDataAll(BLFP.cid_1,
                                           BLFP.cid_2,
                                           BLFP.len_1,
                                           BLFP.len_2,
                                           BLFP.pos_1,
                                           BLFP.pos_2
                                           )
                    
                    # build links data
                    self.buildLinksData(BLFP.cid_1, BLFP.cid_2)
                    self.buildLinksData(BLFP.cid_2, BLFP.cid_1)
        
        with open(bamm_cov_file) as fh:
            for l in fh:
                if l[0] != '#':
                    BCFP = BamMCoverageFileParser(l)
                    
                    self.contig_len[BCFP.cid] = BCFP.cid_len
                    
                    # build coverage data
                    for cov in BCFP.covs:
                        self.buildCovData(BCFP.cid, cov)
                    
                    
    def buildLinksDataAll(self, cid1, cid2, len1, len2, pos1, pos2):
        try:
            self.bamm_links[cid1][cid2].append([len1, len2, pos1, pos2])
        except KeyError:
            try:
                self.bamm_links[cid1][cid2] = [[len1, len2, pos1, pos2]]
            except KeyError:
                self.bamm_links[cid1] = {cid2:[[len1, len2, pos1, pos2]]}
        try:
            self.bamm_links[cid1][cid1].append(len1, len2, pos1, pos2)
        except KeyError:
            try:
                self.bamm_links[cid2][cid1] = [[len2, len1, pos2, pos1]]
            except KeyError:
                self.bamm_links[cid2] = {cid1:[[len2, len1, pos2, pos1]]}
            
    
    def buildCovData(self, cid, cov):
        try:
            self.bamm_cov_data[cid].append(cov)
        except KeyError:
            self.bamm_cov_data[cid] = [cov]
                    
    def buildLinksData(self, cid1, cid2):
        try:
            self.bamm_links_data[cid1][cid2] += 1 
        except KeyError:
            try:
                self.bamm_links_data[cid1][cid2] = 1
            except KeyError:
                self.bamm_links_data[cid1] = {cid2:1}
        try:
            self.bamm_links_total[cid1] += 1 
        except KeyError:
            self.bamm_links_total[cid1] = 1 
            
        try:
            self.bamm_links_total[cid2] += 1
        except KeyError:
            self.bamm_links_total[cid2] = 1

class BamMLinksFileParser(object):
    def __init__(self, l):
        self.readBamMLinksFile(l)
        
    def readBamMLinksFile(self,l):
        tabs = l.rstrip().split("\t")
        #cid_1  cid_2   len_1   pos_1   rev_1   len_2   pos_2   rev_2   file
        self.cid_1 = tabs[0] 
        self.cid_2 = tabs[1]
        self.len_1 = tabs[2]
        self.pos_1 = tabs[3]
        self.rev_1 = tabs[4]
        self.len_2 = tabs[5]
        self.pos_2 = tabs[6]
        self.rev_2 = tabs[7]
        self.file  = tabs[8]

class BamMCoverageFileParser(object):
    def __init__(self, l):
        self.readBamMCoverageFile(l)
        
    def readBamMCoverageFile(self,l):
        tabs = l.rstrip().split("\t")
        #contig Length  649989915.SRR067990.bam 649989915.SRR072965_1.bam
        self.cid                = tabs[0]
        self.cid_len            = tabs[1]
        self.num_coverage_files = len(tabs)
        self.covs               = []
        for i in range(2, self.num_coverage_files):
            self.covs.append(float(tabs[i]))

###############################################################################
###############################################################################
###############################################################################
###############################################################################

def runCommand(cmd):
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
    ContigLinkagePlot(args.bamm_dir,
                      args.assembly_dir,
                      args.outdir,
                      args.force,
                      args.extension,
                      args.highlight_contigs)

###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    #subparsers = parser.add_subparsers(help="--", dest='subparser_name')
    
    #------------------------------
    
    parser.add_argument('bamm_dir', help="Directory containing bamm_cov and bamm_links files")
    parser.add_argument('outdir', help="Output directory")
    parser.add_argument('-assembly_dir','--assembly_dir', help="Directory containing genome fasta files")
    parser.add_argument('-f','--force', action='store_true', help="Force overwrite out file")
    parser.add_argument('-x','--extension', default='.fa', help="Fasta file suffix")
    parser.add_argument('-hc','--highlight_contigs', help="Line-separated list of contigs to highlight")
    #parser.add_argument('-t','--type', default='scatter',help="Set display type: bar or scatter")
    #parser.add_argument('-bc','--bamm_cov', default=False,help="Path to bamm coverage file")
    #parser.add_argument('-bl','--bamm_links', default=False,help="Path to bamm links file")
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################
