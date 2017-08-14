#!/usr/bin/env python
###############################################################################
#
# __mashderep__.py 
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

class Data(object):
    def __init__(self, name):
        self.__name  = name
        self.__links = set()

    @property
    def name(self):
        return self.__name

    @property
    def links(self):
        return set(self.__links)

    def add_link(self, other):
        self.__links.add(other)
        other.__links.add(self)

class MASH(object):
    def __init__(self,
                 mash_file,
                 genome_metadata_file,
                 outfile,
                 checkm,
                 gtdb,
                 file_delimeter,
                 threshold,
                 trusted_genomes
                 ):
        # data
        self.mash_data   = {}
        self.singletons  = {}
        self.genome_stats = {}        
        self.dereplicated_groups = {}
        self.trusted_genomes_list= {}
        
        # run scripts
        if trusted_genomes:
            self.parseTrustedGenomeList(trusted_genomes)
        if checkm:
           self.parseCheckmFile(genome_metadata_file,file_delimeter) 
        elif gtdb:
            self.parseGTDBmetadata(genome_metadata_file,file_delimeter)
        self.parseMashFile(mash_file,threshold,trusted_genomes)
        self.dereplicateGenomes(outfile)
    
    def parseTrustedGenomeList(self,
                               trusted_genomes):
        with open(trusted_genomes) as fh:
            for l in fh:
                genome = l.rstrip()
                self.trusted_genomes_list[genome] = 1
    
    def parseCheckmFile(self,
                        genome_metadata_file,
                        file_delimeter):
        with open(genome_metadata_file) as fh:
            line_number = 0
            for l in fh:
                if file_delimeter == "tab":
                    line_split = l.rstrip().split("\t")
                elif file_delimeter == "comma":
                    line_split = l.rstrip().split(",")
                if line_number > 0:
                    genome          = self.processGenomeString(line_split[0])
                    completion      = float(line_split[11])
                    contamination   = float(line_split[12])
                    
                    # add to checkm data dictionary
                    self.genome_stats[genome] = [completion, contamination]
                line_number += 1
        
    def parseGTDBmetadata(self,
                        genome_metadata_file,
                        file_delimeter):
        with open(genome_metadata_file) as fh:
            line_number = 0
            for l in fh:
                if file_delimeter == "tab":
                    line_split = l.rstrip().split("\t")
                elif file_delimeter == "comma":
                    line_split = l.rstrip().split(",")
                if line_number > 0:
                    genome          = self.processGenomeString(line_split[0])
                    completion      = float(line_split[30])
                    contamination   = float(line_split[31])
                    
                    # add to checkm data dictionary
                    self.genome_stats[genome] = [completion, contamination]
                line_number += 1
        
    def parseMashFile(self,
                      mash_file,
                      threshold,
                      trusted_genomes):
        
        with open(mash_file) as fh:
            line_number = 0
            for l in fh:
                tabs = l.rstrip().split("\t")
                if line_number > 0:
                    genome_a            = self.processGenomeString(tabs[0].split("/")[-1])
                    genome_b            = self.processGenomeString(tabs[1].split("/")[-1])
                    mash_dist           = (1 - float(tabs[2]))*100
                    p_value             = float(tabs[3])
                    matching_hashes     = tabs[4]
                    
                    if trusted_genomes:
                        if genome_a in self.trusted_genomes_list and genome_b in self.trusted_genomes_list:
                            if mash_dist >= threshold:
                                self.addConnection(genome_a, genome_b)
                                # add genomes to non-singelton list
                                self.addGenome(genome_a)
                                self.addGenome(genome_b)
                            else:
                                self.removeGenome(genome_a)
                                self.removeGenome(genome_b)
                    else:
                        if mash_dist >= threshold:
                                self.addConnection(genome_a, genome_b)
                                # add genomes to non-singelton list
                                self.addGenome(genome_a)
                                self.addGenome(genome_b)
                        else:
                            self.removeGenome(genome_a)
                            self.removeGenome(genome_b)
                line_number += 1
      
    def processGenomeString(self,
                            genome_string):
        if "genomic.fna" in genome_string:
            genome_name = genome_string.replace('_genomic.fna','')
        elif genome_string[-3:] == ".fa":
            genome_name = genome_string[:-3]
        elif genome_string[-4:] == ".fna":
            genome_name = genome_string[:-4]
        return genome_name
        
    def addGenome(self,
                  genome):
        self.singletons[genome] = 1
    
    def removeGenome(self,
                     genome):
        try:
            if self.singletons[genome] == 1:
                pass
        except KeyError:
            self.singletons[genome] = 0
    
    def dereplicateGenomes(self,
                           outfile):
        node_dict = {}

        nodes = set()

        # build trees from blast data

        regions_to_nodes= {}
        
        for genome_a in self.mash_data.keys():
            
            if genome_a not in regions_to_nodes:
                D = Data(genome_a)
                regions_to_nodes[genome_a] = D
                nodes |= {D}
                
                for genome_b in self.mash_data[genome_a].keys():
                    
                    if genome_b not in regions_to_nodes:
                        D = Data(genome_b)
                        regions_to_nodes[genome_b] = D
                        nodes |= {D}

        for genome_a in self.mash_data.keys():
            
            for genome_b in self.mash_data[genome_a].keys():
                regions_to_nodes[genome_a].add_link(regions_to_nodes[genome_b])

        # Find all the connected components
        number = 1

        # initialise outfile
        outfile = open(outfile, 'w')
        
        # print groups
        for components in self.connected_components(nodes):
            names = sorted(node.name for node in components)
            names_count= len(names)
            names = ",".join(names)
            genomes = names
            
            rep_genome = self.pickRepresentativeGenome(genomes, number)
            
            genomes_split = names.split(",")
            
            genomes_concatenated = []
            for genome in genomes_split:
                if genome != rep_genome:
                    genomes_concatenated.append(genome)
            outfile.write("%s\t%s\n" % (rep_genome, ",".join(genomes_concatenated)))
            number += 1        
    
        # print singletons
        for genome in self.singletons.keys():
            if self.singletons[genome] == 0:
                outfile.write("%s\t\n" % (genome))
    
    def pickRepresentativeGenome(self,
                                 genomes,
                                 group):
        genomes = genomes.split(",")
        
        for genome in genomes:
            if group not in self.dereplicated_groups:
                
                # add to dictionary
                self.dereplicated_groups[group] = genome
            
            else:
                # incumbent genome for group #
                genome_incumbent               = self.dereplicated_groups[group]
                genome_incumbent_contamination = self.genome_stats[genome_incumbent][1]
                genome_incumbent_completion    = self.genome_stats[genome_incumbent][0]  
                genome_incumbent_relative_completion = genome_incumbent_completion - genome_incumbent_contamination
                
                # new genome to compare to incumbent
                genome_contamination           = self.genome_stats[genome][1]
                genome_completion              = self.genome_stats[genome][0] 
                genome_relative_completion = genome_completion - genome_contamination
                
                if genome_relative_completion > genome_incumbent_relative_completion:
                    # replace incumbent with new genome
                    self.dereplicated_groups[group] = genome
        
        return self.dereplicated_groups[group]
        
    def addConnection(self, id1, id2):
        try:
            self.mash_data[id1][id2] = 1 
        except KeyError:
            self.mash_data[id1] = {id2:1}
        try:
            self.mash_data[id2][id1] = 1 
        except KeyError:
            self.mash_data[id2] = {id1:1}
    
    def connected_components(self,nodes):
        # List of connected components found. The order is random.
        result = []

        # Make a copy of the set, so we can modify it.
        nodes = set(nodes)

        nn = float(len(nodes))
        remains = nn

        # Iterate while we still have nodes to process.
        while nodes:

            # Get a random node and remove it from the global set.
            n = nodes.pop()

            # This set will contain the next group of nodes connected to each other.
            group = {n}

            # Build a queue with this node in it.
            queue = [n]

            # Iterate the queue.
            # When it's empty, we finished visiting a group of connected nodes.
            while queue:

                # Consume the next item from the queue.
                n = queue.pop(0)

                # Fetch the neighbors.
                neighbors = n.links

                # Remove the neighbors we already visited.
                neighbors.difference_update(group)

                # Remove the remaining nodes from the global set.
                remains -= float(len(neighbors))
                #print "%f completed" % (1 - remains/nn)
                nodes.difference_update(neighbors)

                # Add them to the group of connected nodes.
                group.update(neighbors)

                # Add them to the queue, so we visit them in the next iterations.
                queue.extend(neighbors)

            # Add the group to the list of groups.
            result.append(group)

        # Return the list of groups.
        return result

###############################################################################
###############################################################################
###############################################################################
###############################################################################
   
def doWork( args ):
    """ Main wrapper"""
    
    # run scripts
    MASH(args.mash_file,
         args.genome_metadata_file,
         args.outfile,
         args.checkm,
         args.gtdb,
         args.file_delimeter,
         args.threshold,
         args.trusted_genomes)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################

if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog='PROG',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument('mash_file',help='File containing all pairwise mash results ')
    parser.add_argument('genome_metadata_file',help='File containing genome completion and contamination stats. Set File type to either gtdb or checkm')
    parser.add_argument('outfile',help='Output filename')
    parser.add_argument('-checkm','--checkm',action="store_true",help='Set genome metadata file type to checkm')
    parser.add_argument('-gtdb','--gtdb',action="store_true",help='Set genome metadata file type to gtdb')
    parser.add_argument('-d','--file_delimeter',choices=['tab','comma'],help='Delimiter for gtdb metadata or checkm file.')
    parser.add_argument('-threshold','--threshold',type=int,default=95,help='Set mash percentage threshold for clustering')
    parser.add_argument('-tg','--trusted_genomes',help='File containing line-separated list of trusted genomes to use.')
    
    # parse the arguments
    args = parser.parse_args()
    
    # do what we came here to do
    doWork(args)
    
###############################################################################
###############################################################################
###############################################################################
###############################################################################
