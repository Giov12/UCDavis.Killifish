#!/bin/env python3

from __future__ import annotations
import argparse
import os
import gzip
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
        self.dir_5p = {} # how many reads support
        self.dir_3p = {} # either the 5' or 3' directions

        # struct
        # direction map = {3: {cntg: count, cntg: count}, 5: {cntg: count, cntg: count}}
        #
        for d in ['3', '5']:
            self.dir_3p[d] = defaultdict(int)
            self.dir_5p[d] = defaultdict(int)

class Node:
    def __init__(self, name: str, placement: str, index: int):
        self.name       = name
        self.ref_seq    = placement
        self.idx        = index
        self.right      = list()
        self.left       = list()
        self._lefts     = list() # to hold onto equally supported nodes
        self._rights    = list() # hold onto equally supported nodes
        self.right_node = None
        self.right_orie = '' # orientation of right node
        self.left_node  = None
        self.left_orie  = '' # orientation of left node

    def add_edge(self, node: Node, support: int, self_side: int, other_side: str) -> int:
        #
        # other side -> 5' or 3' side of the other contig
        #
        
        if (self_side == 0):
            self.left.append((node, support, other_side))
        else:
            self.right.append((node, support, other_side))

        return 0  

    def find_most_supported(self) -> int:

        # func() to assign the most supported nodes
        # as the left & right neighbors

        score = -float("inf")
        lefts = list()
        for i in range(len(self.left)):
            entry = self.left[i]
            val   = entry[1]
            if (i == 0):
                score = val
            if (val == score):
                lefts.append(entry)
            elif (val > score):
                lefts.clear()
                lefts.append(entry)

        self._lefts.extend(lefts)

        # repeat for right side
        score  = -float("inf")
        rights = list()
        for i in range(len(self.right)):
            entry = self.right[i]
            val   = entry[1]
            if (i == 0):
                score = val
            if (val == score):
                rights.append(entry)
            elif (val > score):
                rights.clear()
                rights.append(entry)

        self._rights.extend(rights)

        return 0
    
    def no_left_neighbors(self) -> bool:
        return len(self._lefts) == 0
    
    def get_lefts(self) -> list[tuple[Node, int, str]]:
        return self._lefts
    
    def get_rights(self) -> list[tuple, None, int, str]:
        return self._rights
    
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
        fields    = line.strip().split('\t')
        pcntg_aln = fields[1] # primary contig/alignment
        pcntg, pside = pcntg_aln.split('_')
        if (pcntg not in contig_map):
            contig_map[pcntg] = Contig(pcntg)
        for scntg in fields[2].split(','):
            cntg = scntg[:-2]
            side = scntg[-1]
            if (pside == '5'):
                contig_map[pcntg].dir_5p[side][cntg] += 1
            else:
                contig_map[pcntg].dir_3p[side][cntg] += 1

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

def find_best_left_neighbor(node: Node) -> tuple[Node, str]:
    """find the most supported node on the 5' region"""

    global nodes

    left_nodes = node.get_lefts()

    for entry in left_nodes:
        left_node  = entry[0]
        support    = entry[1]
        left_orien = entry[2]
        left_node  = nodes[left_node.name]
        matched    = False
        for other_entries in left_node.get_lefts():
            if (other_entries[0].name == node.name):
                matched = True
                break
        if (matched):
            return (left_node, left_orien)
        if (matched == False):
            for other_entries in left_node.get_rights():
                if (other_entries[0].name == node.name):
                    matched = True
                    break
        if (matched):
            return (left_node, left_orien)

    return None  

def find_best_right_neighbor(node: Node) -> tuple[Node, str]:
    """find the most supported node on the 5' region"""

    global nodes

    right_nodes = node.get_rights()

    for entry in right_nodes:
        right_node = entry[0]
        support    = entry[1]
        right_orien = entry[2]
        right_node  = nodes[right_node.name]
        matched    = False
        for other_entries in right_node.get_lefts():
            if (other_entries[0].name == node.name):
                matched = True
                break
        if (matched):
            return (right_node, right_orien)
        if (matched == False):
            for other_entries in right_node.get_rights():
                if (other_entries[0].name == node.name):
                    matched = True
                    break
        if (matched):
            return (right_node, right_orien)

    return None

def create_components() -> int:
    """function to iterate through all the nodes and construct components out of them"""

    global nodes

    # 
    # step 1: get nodes with no left neighbors
    #
    no_lefts = list()
    for node in nodes.values():
        node.find_most_supported() # establish any neighbors
        if (node.no_left_neighbors()):
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
        component = [(node, '')]
        cur       = node
        visited    = set()
        visited.add(node.name)
        while (True):
            # now to get best match to the 3' regions
            next_neighbor = find_best_right_neighbor(cur)
            if (next_neighbor == None):
                break
            next_node  = next_neighbor[0]
            next_orien = next_neighbor[1]
            if (next_node.name in visited):
                break
            component.append((next_node, next_orien))
            visited.add(next_node.name)
            cur        = next_node
        components.append(component)

    fh = open("Connected.contigs.tsv", 'w')
    fh.write("#Component.ID\tNum.Contigs\tContig.IDs\n")

    # step 3, write and print to the console  the components
    for i, component in enumerate(components):
        size = len(component)
        if (size == 1):
            continue
        pstr = ' '.join(f"{repr(node)}({orien if orien else 'start'})" for node, orien in component)
        wstr = ','.join(f"{repr(node)}:{orien if orien else 'start'}" for node, orien in component)
        print("Component", i, f"Size: {size}")
        print(pstr, '\n')
        fh.write(f"{i}\t{size}\t{wstr}\n")

    fh.close()
    
    return 0

def create_edges_and_write_support() -> int:
    """a function to write the number of briding alignments"""

    global agp_map, contig_map, valid_pairs, reverse_names, nodes

    ofh   = open("Separated_contigs.tsv", 'w')
    check = len(reverse_names) > 0
    dirs  = ['5', '3']
    ofh.write("#Contig.A\tChr.A\tIdx.A\tSide.A\tConitg.B\tChr.B\tIdx.B\tSide.B\tRead.Support\n")

    for contig in contig_map.values():
        source, idx  = agp_map[contig.name]
        for d in dirs:
            five_primes  = [(cntg, count) for cntg, count in contig.dir_5p[d].items()]
            three_primes = [(cntg, count) for cntg, count in contig.dir_3p[d].items()]
            five_primes.sort(key = lambda x : x[1], reverse=True)
            three_primes.sort(key = lambda x : x[1], reverse=True)
        
            for pair in five_primes:
                left_cntg = pair[0]
                cntg_cnt  = pair[1] # number of reads supporting this match
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
                    outline = f"{name1}\t{source}\t{idx}\t5\t{name2}\t{cntg_src}\t{cntg_idx}\t{d}\t{cntg_cnt}\n"
                    ofh.write(outline)

                    node1 = nodes[name1]
                    node2 = nodes[name2]
                    node1.add_edge(node2, cntg_cnt, 0, d) # add to left
            
            for pair in three_primes:
                right_cntg = pair[0]
                cntg_cnt   = pair[1] # number of reads supporting this match
                cntg_src, cntg_idx = agp_map[right_cntg]
                if (cntg_src != source or abs(idx - cntg_idx) > 1):
                    name1 = contig.name
                    name2 = right_cntg
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
                    outline = f"{name1}\t{source}\t{idx}\t3\t{name2}\t{cntg_src}\t{cntg_idx}\t{d}\t{cntg_cnt}\n"
                    ofh.write(outline)

                    node1 = nodes[name1]
                    node2 = nodes[name2]
                    node1.add_edge(node2, cntg_cnt, 1, d) # add to right

    ofh.close()

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

    # create the edges and write it out for inspection later
    create_edges_and_write_support()

    # now connect edges
    create_components()

    return 0

if __name__ == "__main__":
    main()
