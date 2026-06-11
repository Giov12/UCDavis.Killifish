#!/bin/env python3

from __future__ import annotations
import argparse
import os
import gzip
import copy
from collections import defaultdict

agp           = ''
seqs_file     = ''
contig_map    = dict() # contig -> Contig object
agp_map       = dict() # contig -> (chr, idx) # contig to the anchor chr & its index
contig_names  = dict() # contig -> contig name [OPTIONAL]
reverse_names = dict() # contig name -> contig [OPTIONAL]
valid_pairs   = set() 
names_file    = ''
nodes         = dict() # contig -> Node

class Contig:
    def __init__(self, name: str):
        self.name   = name
        self.dir_5p = defaultdict(int) # how many reads support
        self.dir_3p = defaultdict(int) # either the 5' or 3' directions

class Node:
    def __init__(self, name: str, placement: str, index: int):
        self.name       = name
        self.ref_seq    = placement
        self.idx        = index
        self.right      = list()
        self.left       = list()
        self.right_node = None
        self.left_node  = None

    def add_edge(self, node: Node, support: int, side: int) -> int:
        if (side == 0):
            self.left.append((node, support))
        else:
            self.right.append((node, support))

        return 0  

    def assign_most_supported(self) -> int:

        # func() to assign the most supported nodes to
        # as the left & right neighbors

        cur   = None
        score = -float("inf")
        for i in range(len(self.left)):
            entry = self.left[i]
            node  = entry[0]
            val   = entry[1]
            if (val > score):
                cur   = node
                score = val

        self.left_node = cur # if none, it will remain as none

        # repeat for right side
        cur   = None
        score = -float("inf")
        for i in range(len(self.right)):
            entry = self.right[i]
            node  = entry[0]
            val   = entry[1]
            if (val > score):
                cur   = node
                score = val

        self.right_node = cur

        return 0
    
    def __eq__(self, node: Node) -> bool:
        if (isinstance(node, Node) == False):
            return False
        return self.name == node.name
    
    def __hash__(self) -> int:
        return hash(self.name)
    
    def __repr__(self) -> str:
        return f"{self.name}_{self.ref_seq}_{self.idx}"

def set_arguments() -> int:
    """get & set the arguments"""

    global agp, seqs_file, names_file

    d = "Find putative neighoring contigs that were seperated"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-e", "--edge-seqs", help="A list of aligned sequences that are clipped and at the ends of the reference sequences", required=True, type=str, default='')
    parser.add_argument("-a", "--agp",       help="AGP file for the current reference assembly", default='', type=str, required=True)
    parser.add_argument("-n", "--names",     help="AGP file containing contig names to use instead [Optional]", default='', type=str)

    args = parser.parse_args()

    agp        = args.agp
    seqs_file  = args.edge_seqs
    names_file = args.names
  
    assert os.path.isfile(agp), f"Could not locate file {agp}"
    assert os.path.isfile(seqs_file), f"Could not locate file {seqs_file}"

    if (names_file != ''):
        assert os.path.isfile(names_file), f"Could not locate file {names_file}"
    
    return 0

def parse_edges() -> int:
    """fill the contig mapping"""

    global seqs_file, contig_map

    #
    # structure
    # read<tab>primary.conitg<tab>secondary.contigs\n
    #

    fh = open(seqs_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.strip().split('\t')
        pcntg  = fields[1] # primary contig/alignment
        if (pcntg not in contig_map):
            contig_map[pcntg] = Contig(pcntg)
        for scntg in fields[2].split(','):
            cntg = scntg[:-2]
            side = scntg[-1]
            if (side == '5'):
                contig_map[pcntg].dir_5p[cntg] += 1
            else:
                contig_map[pcntg].dir_3p[cntg] += 1

    fh.close()

    return 0

def read_names() -> int:
    """helper function to read in alternative contig names"""

    global names_file, contig_names, reverse_names

    if (names_file == ''):
        return 1 # no names to read in
    
    fh = gzip.open(names_file, "rt") if names_file.endswith(".gz") else open(names_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        if (fields[4] != 'W' and fields[4] != 'w'):
            continue
        contig    = fields[0]
        cntg_name = fields[5]
        contig_names[cntg_name] = contig
        reverse_names[contig]   = cntg_name # for reverse search
        
    fh.close()

    return 0

def load_agp() -> int:
    """read in the structure of the assembly"""

    global agp, agp_map, contig_names

    fh          = gzip.open(agp, "rt") if agp.endswith(".gz") else open(agp, 'r')
    check_names = len(contig_names) > 0
    idx         = 0
    prev        = ''

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields = line.split('\t')
        if (fields[4] != 'W' and fields[4] != 'w'):
            continue
        source    = fields[0] # e.g., chrom 1
        component = fields[5] # e.g., contig 1
        if (check_names and component in contig_names):
            component = contig_names[component] # swap the contig names
        if (prev != source):
            idx  = 1
            prev = source
        agp_map[component] = (source, idx)
        idx += 1

    fh.close()

    return 0

def find_agreements() -> int:
    """find neighboring contigs that show reciprocity"""

    global contig_map, valid_pairs

    for contig_1 in contig_map.values():
        # this contig is on the 5', so contig_1 should be on 3'
        for contig2_name in contig_1.dir_5p.keys():
            if (contig2_name in contig_map):
                contig_2 = contig_map[contig2_name]
                if (contig_1.name in contig_2.dir_3p):
                    tmp  = [contig_1.name, contig2_name]
                    tmp.sort()
                    pair = tuple(tmp)
                    if (pair not in valid_pairs):
                        valid_pairs.add(pair)

        # this contig is on the 3', so contig_1 should be on 5'
        for contig2_name in contig_1.dir_3p.keys():
            if (contig2_name in contig_map):
                contig_2 = contig_map[contig2_name]
                if (contig_1.name in contig_2.dir_5p):
                    tmp  = [contig_1.name, contig2_name]
                    tmp.sort()
                    pair = tuple(tmp)
                    if (pair not in valid_pairs):
                        valid_pairs.add(pair)

    return 0

def create_components() -> int:
    """function to iterate through all the nodes and construct components out of them"""

    global nodes

    # 
    # step 1: get nodes with no left neighbors
    #
    no_lefts = list()
    for node in nodes.values():
        node.assign_most_supported() # establish any neighbors
        if (node.left_node == None):
            no_lefts.append(node)

    if (len(no_lefts) == 0):
        print("Found no usuable contigs to start bridging contigs")
        return 1

    #
    # step 2: go through each node from
    # left to right & verify that they
    # are the most supported matches
    #
    components = list()
    for node in no_lefts:
        component = node
        cur       = component
        next      = component.right_node
        size      = 1
        components.append(cur) # this is the head node
        while (next != None):
            if (next.left_node == cur):
                cur = next
                next = cur.right_node
                size += 1
            else:
                cur.right_node = None # break the connection
                break
        # print("Added a component of size", size)

    fh = open("Connected.contigs.tsv", 'w')
    fh.write("#Component.ID\tNum.Contigs\tContig.IDs\n")

    # step 3, write and print to the console  the components
    for i, component in enumerate(components):
        temp = list()
        
        while (component != None):
            temp.append(repr(component))
            component = component.right_node
        size = len(temp)
        pstr = ' '.join(temp) # print string
        wstr = pstr.replace(' ', ',') # write string
        print("Component", i, f"Size: {size}")
        print(pstr, '\n')
        fh.write(f"{i}\t{size}\t{wstr}\n")

    fh.close()
    
    return 0

def identify_separations() -> int:
    """now check the pairings that are separated"""

    global agp_map, contig_map, valid_pairs, reverse_names, nodes

    total = 0
    ofh   = open("Separated_contigs.tsv", 'w')
    check = len(reverse_names) > 0
    ofh.write("#Contig.A\tChr.A\tIdx.A\tConitg.B\tChr.B\tIdx.B\tEdge\tRead.Support\n")

    for contig in contig_map.values():
        source, idx  = agp_map[contig.name]
        five_primes  = [(cntg, count) for cntg, count in contig.dir_5p.items()]
        three_primes = [(cntg, count) for cntg, count in contig.dir_3p.items()]
        five_primes.sort(key = lambda x : x[1], reverse=True)
        three_primes.sort(key = lambda x : x[1], reverse=True)
        d = "5'" # direction
        for pair in five_primes:
            left_cntg = pair[0]
            cntg_cnt  = pair[1] # number of reads supporting this match
            tmp       = [contig.name, left_cntg]
            tmp.sort()
            combo     = tuple(tmp)
            if (combo not in valid_pairs):
                continue
            cntg_src, cntg_idx = agp_map[left_cntg]
            if (cntg_src != source or abs(idx - cntg_idx) > 1):
                name1 = contig.name
                name2 = left_cntg
                # check if we need to update the contig names
                if (check):
                    if (name1 in reverse_names):
                        name1 = reverse_names[name1]
                    if (name2 in reverse_names):
                        name2 = reverse_names[name2]
                # create nodes to construct components later
                if (name1 not in nodes):
                    nodes[name1] = Node(name1, source, idx)
                if (name2 not in nodes):
                    nodes[name2] = Node(name2, cntg_src, cntg_idx)

                node1 = nodes[name1]
                node2 = nodes[name2]
                node1.add_edge(node2, cntg_cnt, 0) # add to left
                outline = f"{name1}\t{source}\t{idx}\t{name2}\t{cntg_src}\t{cntg_idx}\t{d}\t{cntg_cnt}\n"
                ofh.write(outline)
                total += 1
        d = "3'" # direction
        for pair in three_primes:
            right_cntg = pair[0]
            cntg_cnt   = pair[1] # number of reads supporting this match
            tmp        = [contig.name, right_cntg]
            tmp.sort()
            combo      = tuple(tmp)
            if (combo not in valid_pairs):
                continue
            cntg_src, cntg_idx = agp_map[right_cntg]
            if (cntg_src != source or abs(idx - cntg_idx) > 1):
                name1 = contig.name
                name2 = right_cntg
                if (check):
                    if (name1 in reverse_names):
                        name1 = reverse_names[name1]
                    if (name2 in reverse_names):
                        name2 = reverse_names[name2]
                if (name1 not in nodes):
                    nodes[name1] = Node(name1, source, idx)
                if (name2 not in nodes):
                    nodes[name2] = Node(name2, cntg_src, cntg_idx)
                node1 = nodes[name1]
                node2 = nodes[name2]
                node1.add_edge(node2, cntg_cnt, 1) # add to right
                outline = f"{name1}\t{source}\t{idx}\t{name2}\t{cntg_src}\t{cntg_idx}\t{d}\t{cntg_cnt}\n"
                ofh.write(outline)
                total += 1

    ofh.close()
    print(f"Found a total of {total} separations")

    return 0

def main() -> int:
    """Entry point to this small program"""

    # get arguments
    set_arguments()

    # read in the mapping info
    parse_edges()

    # read in any alternative names if present
    read_names()

    # load the agp file
    load_agp()

    # find reciprocal matches
    find_agreements()

    # now sieve through the pairings and placements
    identify_separations()

    # now connect the seperations
    create_components()

    return 0

if __name__ == "__main__":
    main()
