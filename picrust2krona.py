#!/usr/bin/env python
###############################################################################
#
# __picrust2krona__.py - parse picrust KO output and visualise using krona
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

methane_ko = {
    'K00399':1,
    'K00401':1,
    'K00402':1,
}

nitrate_reduction_ko = {
    'K00370':1,
    'K00371':1,
    'K00374':1,
    'K00373':1
}

sulfate_reduction_ko = {
    'K00958':1,
    'K00394':1,
    'K00395':1,
    'K11180':1,
    'K11181':1,
    'K11179':1,
    'K06039':1
}

sulfate_oxidation_ko = {
    'K17224':1,
    'K17222':1,
    'K17223':1
}
formate_ko = {
    'K00122':1,
    'K00123':1,
    'K00125':1,
    'K00126':1,
}
acetate_ko = {
    'K01905':1
}
propionate_ko  = {
    'K00932':1,
    'K01908':1    
}
butyrate_ko = {
    'K00929':1
}


class Picrust2KronaContributions(object):
    def __init__(self,
                 otu_table,
                 picrust_ko_file,
                 taxonomy_file,
                 ):
        self.gg_functions       = {} 
        self.taxonomy           = {}
        self.samples            = []
        self.otu_table          = {}
        self.krona_data         = {}
        self.otu_table_by_taxa  = {}
        
        self.parseTaxonomyFile(taxonomy_file)
        self.parsePicurstKO(picrust_ko_file)
        self.parseOTUtable(otu_table)
        self.makekronaplots()
    
    def parsePicurstKO(self,
                       picrust_ko_file):
        with open(picrust_ko_file) as fh:
            line_number = 0
            for l in fh:
                if line_number>0:
                    tabs   = l.rstrip().split("\t")
                    ko_id  = tabs[0]
                    gg_id  = tabs[2]
                    
                    self.initialseGG(gg_id)
                    domain = self.findDomain(gg_id)# bacteria or archaea
                    
                    # check for methane
                    if ko_id in methane_ko:
                         self.gg_functions[gg_id]['methanogens'] = 1
                            
                    # check for sulfate reduction
                    if domain == 'Bacteria':
                        if ko_id in sulfate_reduction_ko:
                            self.gg_functions[gg_id]['sulfate_reducing_bacteria'] = 1
                    else:
                        if ko_id in sulfate_reduction_ko:
                            self.gg_functions[gg_id]['sulfate_reducing_archaea'] = 1
                    
                    # check for nitrate reduction
                    if ko_id in nitrate_reduction_ko:
                        self.gg_functions[gg_id]['nitrate_reducing_bacteria'] = 1
                            
                    # check for sulfate oxidation
                    if ko_id in sulfate_oxidation_ko:
                        self.gg_functions[gg_id]['sulfate_oxidating_bacteria'] = 1
                    
                    # checkm for formate
                    if ko_id in formate_ko:
                        self.gg_functions[gg_id]['formate'] = 1
                    
                    # checkm for acetate    
                    if ko_id in acetate_ko:
                        self.gg_functions[gg_id]['acetate'] = 1
                        
                    # checkm for propionate
                    if ko_id in propionate_ko:
                        self.gg_functions[gg_id]['propionate'] = 1
                    
                    # checkm for butyrate
                    if ko_id in butyrate_ko:
                        self.gg_functions[gg_id]['butyrate'] = 1
                    
                line_number+=1
    
    def findDomain(self,
                   gg_id):
        if 'Bacteria' in self.taxonomy[gg_id]:
            return 'Bacteria'
        else:
            return 'Archaea'
    
    def initialseGG(self,
                    gg_id):
        if gg_id not in self.gg_functions:
            self.gg_functions[gg_id] = {}
            self.gg_functions[gg_id]['methanogens'] = 0
            self.gg_functions[gg_id]['sulfate_reducing_bacteria'] = 0
            self.gg_functions[gg_id]['sulfate_reducing_archaea'] = 0
            self.gg_functions[gg_id]['nitrate_reducing_bacteria'] = 0
            self.gg_functions[gg_id]['sulfate_oxidating_bacteria'] = 0
            self.gg_functions[gg_id]['formate'] = 0
            self.gg_functions[gg_id]['acetate'] = 0
            self.gg_functions[gg_id]['propionate'] = 0
            self.gg_functions[gg_id]['butyrate'] = 0
    
    def parseTaxonomyFile(self,
                          taxonomy_file):
        with open(taxonomy_file) as fh:
            for l in fh:
                tabs  = l.rstrip().split("\t")
                gg_id = tabs[0]
                taxa  = tabs[1]
                self.taxonomy[gg_id] = taxa

    def parseOTUtable(self,
                      otu_table):
        with open(otu_table) as fh:
            for l in fh:
                tabs = l.rstrip().split("\t")
                if "#OTU" in l:
                    for sample_id in tabs[1:]:
                        self.samples.append(sample_id)
                elif "#" not in l:
                    gg_id = tabs[0]
                    taxonomy = self.taxonomy[gg_id]
                    for i,abundance in enumerate(tabs[1:]):
                        sample_id = self.samples[i]
                        # for functional
                        try:
                            self.otu_table[sample_id][gg_id] = float(abundance)
                        except KeyError:
                            self.otu_table[sample_id]={gg_id:float(abundance)}
                            
                        # for community
                        try:
                            self.otu_table_by_taxa[sample_id][taxonomy] += float(abundance)
                        except KeyError:
                            try:
                                self.otu_table_by_taxa[sample_id][taxonomy] = float(abundance)
                            except KeyError:
                                self.otu_table_by_taxa[sample_id]={taxonomy:float(abundance)}

    def makekronaplots(self, ):
        
        for sample in self.samples:
            self.initialiseKronaData()
            
            of_presence_absence = open("%s.presence_absence.txt" % sample, 'w')
            
            # print header for file
            header_to_print = "\t".join(["OTU ID",
                                         "Taxonomy",
                                         "Abundance",
                                         "Methanogenesis",
                                         "Sulfate Reduction (Bacteria)",
                                         "Sulfate Reduction (Archaea)",
                                         "Sulfate Oxidation",
                                         "Nitrate Reduction (Bacteria)",
                                         "Formate",
                                         "Acetate",
                                         "Propionate",
                                         "Butyrate"])
            of_presence_absence.write("%s\n" % header_to_print)
            
            for gg_id in self.otu_table[sample]:
                otu_abundance = self.otu_table[sample][gg_id]
                taxonomy = self.taxonomy[gg_id]
                
                str_to_print = "\t".join([gg_id,
                                          taxonomy,
                                          str(otu_abundance)
                                          ])

                #### write to file ####
                # OTUID Abundance Functions....Presence/Absence
                
                # methane
                if self.gg_functions[gg_id]['methanogens'] == 1:
                    presence = "Y"
                    self.krona_data['methanogens_pos'] += otu_abundance
                    
                else:
                    presence = "N"
                    self.krona_data['methanogens_neg'] += otu_abundance
                
                str_to_print = "%s\t%s" % (str_to_print, presence)
                
                # sulfate reduction - bacteria
                if self.gg_functions[gg_id]['sulfate_reducing_bacteria'] == 1:
                    presence = "Y"
                    self.krona_data['sulfate_reducing_bacteria_pos'] += otu_abundance
                else:
                    presence = "N"
                    self.krona_data['sulfate_reducing_bacteria_neg'] += otu_abundance
                
                str_to_print = "%s\t%s" % (str_to_print, presence)
                
                # sulfate reduction - archaea
                if self.gg_functions[gg_id]['sulfate_reducing_archaea'] == 1:
                    presence = "Y"
                    self.krona_data['sulfate_reducing_archaea_pos'] += otu_abundance
                else:
                    presence = "N"
                    self.krona_data['sulfate_reducing_archaea_neg'] += otu_abundance
                
                str_to_print = "%s\t%s" % (str_to_print, presence)
                
                # sulfate oxidation
                if self.gg_functions[gg_id]['sulfate_oxidating_bacteria'] == 1:
                    presence = "Y"
                    self.krona_data['sulfate_oxidating_bacteria_pos'] += otu_abundance
                else:
                    presence = "N"
                    self.krona_data['sulfate_oxidating_bacteria_neg'] += otu_abundance
                
                str_to_print = "%s\t%s" % (str_to_print, presence)
                
                # nitrate reduction
                if self.gg_functions[gg_id]['nitrate_reducing_bacteria'] == 1:
                    presence = "Y"
                    self.krona_data['nitrate_reducing_bacteria_pos'] += otu_abundance
                else:
                    presence = "N"
                    self.krona_data['nitrate_reducing_bacteria_neg'] += otu_abundance
                    
                str_to_print = "%s\t%s" % (str_to_print, presence)
                    
                # formate
                if self.gg_functions[gg_id]['formate'] == 1:
                    presence = "Y"
                    self.krona_data['formate_pos'] += otu_abundance
                else:
                    presence = "N"
                    self.krona_data['formate_neg'] += otu_abundance
                
                str_to_print = "%s\t%s" % (str_to_print, presence)
                    
                # acetate
                if self.gg_functions[gg_id]['acetate'] == 1:
                    presence = "Y"
                    self.krona_data['acetate_pos'] += otu_abundance
                else:
                    presence = "N"
                    self.krona_data['acetate_neg'] += otu_abundance
                
                str_to_print = "%s\t%s" % (str_to_print, presence)
                    
                # propionate
                if self.gg_functions[gg_id]['propionate'] == 1:
                    presence = "Y"
                    self.krona_data['propionate_pos'] += otu_abundance
                else:
                    presence = "N"
                    self.krona_data['propionate_neg'] += otu_abundance
                
                str_to_print = "%s\t%s" % (str_to_print, presence)
                
                # butyrate
                if self.gg_functions[gg_id]['butyrate'] == 1:
                    presence = "Y"
                    self.krona_data['butyrate_pos'] += otu_abundance
                else:
                    presence = "N"
                    self.krona_data['butyrate_neg'] += otu_abundance
            
                str_to_print = "%s\t%s" % (str_to_print, presence)
                
                if otu_abundance > 0:
                    of_presence_absence.write("%s\n" % str_to_print)
            
            # output data to tsv file
            of = open("%s.functional.txt" % sample, 'w')
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # methane
            str_to_write = "\t".join(["Methanogens",
                                      "%f" % self.krona_data['methanogens_pos'],
                                      "%f" % self.krona_data['methanogens_neg']])
            of.write("%s\n" % str_to_write)
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # sulfate reduction - bacteria
            str_to_write = "\t".join(["Sulfate reducing bacteria",
                                    "%f" % self.krona_data['sulfate_reducing_bacteria_pos'],
                                    "%f" % self.krona_data['sulfate_reducing_bacteria_neg']])
            of.write("%s\n" % str_to_write)
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # sulfate reduction - archaea
            str_to_write = "\t".join(["Sulfate reducing archaea",
                                    "%f" % self.krona_data['sulfate_reducing_archaea_pos'],
                                    "%f" % self.krona_data['sulfate_reducing_archaea_neg']])
            of.write("%s\n" % str_to_write)
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # nitrate reduction
            str_to_write = "\t".join(["Nitrate reducing bacteria",
                                    "%f" % self.krona_data['nitrate_reducing_bacteria_pos'],
                                    "%f" % self.krona_data['nitrate_reducing_bacteria_neg']])            
            of.write("%s\n" % str_to_write)
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # sulfate oxidation
            str_to_write = "\t".join(["Sulfate oxidising bacteria",
                                    "%f" % self.krona_data['sulfate_oxidating_bacteria_pos'],
                                    "%f" % self.krona_data['sulfate_oxidating_bacteria_neg']])
            of.write("%s\n" % str_to_write)
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # formate
            str_to_write = "\t".join(["Formate",
                                    "%f" % self.krona_data['formate_pos'],
                                    "%f" % self.krona_data['formate_neg']])
            of.write("%s\n" % str_to_write)
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # acetate
            str_to_write = "\t".join(["Acetate",
                                    "%f" % self.krona_data['acetate_pos'],
                                    "%f" % self.krona_data['acetate_neg']])
            of.write("%s\n" % str_to_write)
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # propionate
            str_to_write = "\t".join(["Propionate",
                                    "%f" % self.krona_data['propionate_pos'],
                                    "%f" % self.krona_data['propionate_neg']])
            of.write("%s\n" % str_to_write)
            
            ################ header ###################
            header_to_write = "\t".join(["functional_class","Present","Absent"])
            of.write("%s\n" % header_to_write)
            
            # butyrate
            str_to_write = "\t".join(["Butyrate",
                                    "%f" % self.krona_data['butyrate_pos'],
                                    "%f" % self.krona_data['butyrate_neg']])
            of.write("%s\n" % str_to_write)
            
            of.close()
            
            
        for sample in self.samples:
            of = open("%s.krona.community.txt" % sample, 'w')
            for taxonomy in self.otu_table_by_taxa[sample]:
                abundance = self.otu_table_by_taxa[sample][taxonomy]
                str_to_write = "%f" % abundance
                for rank in taxonomy.split("; "):
                    if len(rank) > 4:
                        str_to_write = "%s\t%s" % (str_to_write, rank)
                of.write("%s\n" % str_to_write)
        
    def initialiseKronaData(self):
        self.krona_data['methanogens_pos'] = 0
        self.krona_data['methanogens_neg'] = 0
        self.krona_data['sulfate_reducing_bacteria_pos'] = 0
        self.krona_data['sulfate_reducing_bacteria_neg'] = 0
        self.krona_data['sulfate_reducing_archaea_pos'] = 0
        self.krona_data['sulfate_reducing_archaea_neg'] = 0
        self.krona_data['nitrate_reducing_bacteria_pos'] = 0
        self.krona_data['nitrate_reducing_bacteria_neg'] = 0
        self.krona_data['sulfate_oxidating_bacteria_pos'] = 0
        self.krona_data['sulfate_oxidating_bacteria_neg'] = 0
        self.krona_data['formate_pos'] = 0
        self.krona_data['formate_neg'] = 0
        self.krona_data['acetate_pos'] = 0
        self.krona_data['acetate_neg'] = 0
        self.krona_data['propionate_pos'] = 0
        self.krona_data['propionate_neg'] = 0
        self.krona_data['butyrate_pos'] = 0
        self.krona_data['butyrate_neg'] = 0
        
class Picrust2KronaPathway(object):
    def __init__(self,
                 otu_table):
        self.samples   = []
        self.otu_table = {}
    
        self.parseOtuTable(otu_table)
        self.outputKronaFiles()
    
    def parseOtuTable(self,
                      otu_table):
        with open(otu_table) as fh:
            for l in fh:
                tabs = l.rstrip().split("\t")
                if "#" in l:
                    if "#OTU" in l:
                        # header line
                        for sample_id in tabs[1:-1]:
                            self.samples.append(sample_id)    
                            # intialise otu_table
                            self.otu_table[sample_id] = {}
                else:
                    # data
                    kegg_pathway = tabs[0]
                    for index,abundance in enumerate(tabs[1:-1]):
                        sample_id = self.samples[index]
                        if float(abundance) > 0:
                            self.otu_table[sample_id][kegg_pathway] = abundance
    
    def outputKronaFiles(self, ):
        for sample_id in self.otu_table:
            of = open("%s.krona.pathway.txt" % sample_id, 'w')
            for kegg_pathway in self.otu_table[sample_id]:
                abundance = self.otu_table[sample_id][kegg_pathway]
                str_to_write = "\t".join([abundance,kegg_pathway])
                of.write("%s\n" % str_to_write)
            of.close()

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
    if args.subparser_name == 'contributions':
        Picrust2KronaContributions(args.otu_table,
                                    args.ko_table,
                                    args.taxonomy_file)
    elif args.subparser_name == 'pathway':
        Picrust2KronaPathway(args.otu_table)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(help="--", dest='subparser_name')
    
    contributions_parser = subparsers.add_parser('contributions',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    help='Create krona formatted text file from metagenome_contributions.py output',
                                    description='Create krona formatted text file from metagenome_contributions.py output')
    
    contributions_parser.add_argument('otu_table',help='File containing otu table from picrust: normalize_by_copy_number.py ')
    contributions_parser.add_argument('ko_table',help='File containing ko data from picrust: metagenome_contributions.py')
    contributions_parser.add_argument('-t','--taxonomy_file',default="/srv/db/gg/2013_05/gg_13_5_taxonomy.txt",help='File containing gg id -> taxonomy')
   
    pathway_parser = subparsers.add_parser('pathway',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    help='Create krona formatted text file from categorize_by_function.py output',
                                    description='Create krona formatted text file from categorize_by_function.py output')
    pathway_parser.add_argument('otu_table',help='File containing otu table from picrust: normalize_by_copy_number.py ')
    
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################