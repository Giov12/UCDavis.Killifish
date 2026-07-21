#!/bin/env python3

from __future__ import annotations
import argparse
import os
import gzip
from collections import defaultdict
from enum import Enum

agp           = ''
seqs_file     = ''
names_file    = ''
edge_weights  = defaultdict(int) # (end a, end b) -> read count
agp_map       = dict() # contig -> (chr, idx) # contig to the anchor chr & its index
contig_names  = dict() # contig -> contig name [OPTIONAL]
reverse_names = dict() # contig name -> contig [OPTIONAL]
nodes         = dict() # contig -> Node

class END(Enum):
    five_prime   = '5' 
    three_primer = '3' # how they will be represented

# below is a node struct that will represent an end (i.e., 5' or 3')
# of a contig. As such, a contig will actually be represented by 2 nodes
#
class Node:
    def __init__(self, contig: str, end: END, placement: str, index: int):
        self.contig  = contig
        self.end     = end
        self.ref_seq = placement
        self.idx     = index
        self.sibling = None
        self.left    = list()
        self._lefts  = list() # to hold onto equally supported nodes

    def add_edge(self, node: Node, support: int) -> int:
        self.left.append((node, support))
        return 0  

    def find_most_supported(self) -> int:

        # func() to assign the most supported nodes
        # on the left side

        score = -float("inf")
        lefts = list()

        for entry in self.left:
            val = entry[1]
            if (val > score):
                score = val
                lefts = [entry]
            elif (val == score):
                lefts.append(entry)
        self._lefts = lefts

        return 0
    
    def no_left_neighbors(self) -> bool:
        return len(self._lefts) == 0
    
    def get_lefts(self) -> list[tuple[Node, int]]:
        return self._lefts
    
    def __eq__(self, node: Node) -> bool:
        if (isinstance(node, Node) == False):
            return False
        return self.contig == node.contig
    
    def __hash__(self) -> int:
        return hash((self.contig, self.end))
    
    def __repr__(self) -> str:
        return f"{self.contig}_{self.ref_seq}_{self.idx}"

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

    global seqs_file, edge_weights

    #
    # structure
    # read<tab>primary.conitg<tab>secondary.contigs_end\n
    #

    fh = open(seqs_file, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        fields    = line.strip().split('\t')
        pcntg_aln = fields[1] # primary contig/alignment
        for scntg in fields[2].split(','):
            key = tuple(sorted((pcntg_aln, scntg)))
            edge_weights[key] += 1

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

def get_node(contig: str, end: END, src: str, idx: int) -> Node:
    """return the end node and/or create both nodes"""

    global nodes

    key = f"{contig}_{end.value}"

    # do we need to create it?
    if (key not in nodes):
        other      = END.three_primer if (end == END.five_prime) else END.five_prime
        okey       = f"{contig}_{other.value}" # key to the other side (5' vs 3')
        nodes[key] = Node(contig, end, src, idx)
        if (okey not in nodes): # creating both sides
            nodes[okey] = Node(contig, other, src, idx)
        nodes[key].sibling  = nodes[okey]
        nodes[okey].sibling = nodes[key]

    return nodes[key]

def find_best_left_neighbor(node: Node) -> Node | None:
    """find the best & mutally supported node of this end"""

    global nodes

    left_nodes = node.get_lefts()

    for entry in left_nodes:
        other_node = entry[0]
        for other in other_node.get_lefts():
            if (other[0].contig == node.contig and other[0].end == node.end):
                return other_node

    return None

def create_components() -> int:
    """function to iterate through all the nodes and construct components out of them"""

    global nodes

    # 
    # step 1: establish the best partner to every node/end
    #
    for node in nodes.values():
        node.find_most_supported()

    #
    # step 2, establish and get the starts
    # starts: ends with a terminal end
    #
    starts = list()
    for node in nodes.values():
        if (node.no_left_neighbors()):
            starts.append(node)

    if (len(starts) == 0):
        print("Found no usuable contigs to start bridging contigs")
        return 1

    #
    # step 3: establish the component
    # starting from a start node (i.e., terminal node with no left)
    # cross the internal edge to the sibling node (node from other end of the contig)
    # and then bridge from sibling to the next node. Orientation is '+'
    # when we enter at the 5' end and '-' when we enter a 3' end
    # left to right & verify that they are the most supported matches
    #
    components = list()
    visited    = set()
    for node in starts:
        if (node.contig in visited):
            continue
        component = list()
        cur       = node
        while (cur != None and cur.contig not in visited):
            sib    = cur.sibling # internal edge to other end of contig
            orient = '+' if (cur.end == END.five_prime) else '-'
            component.append((cur, orient))
            visited.add(cur.contig)
            visited.add(sib.contig)
            next_node = find_best_left_neighbor(sib)
            if (next_node == None or next_node.contig in visited):
                cur = None
            else:
                cur = next_node
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

    global agp_map, edge_weights, reverse_names, nodes

    ofh = open("Separated_contigs.tsv", 'w')
    ofh.write("#Contig.A\tChr.A\tIdx.A\tSide.A\tConitg.B\tChr.B\tIdx.B\tSide.B\tRead.Support\n")

    for (key1, key2), weight in edge_weights.items():
        contig1, edge1 = key1.rsplit('_', 1)
        contig2, edge2 = key2.rsplit('_', 1)
        agpKey1        = agp_map.get(contig1) # tuple (ref seq, index)
        agpKey2        = agp_map.get(contig2)

        # skip contigs already neighboring in the currently assembly
        if (agpKey1 != None and agpKey2 != None):
            if (agpKey1[0] == agpKey2[0] and abs(agpKey1[1] - agpKey2[1]) <= 1):
                continue

        # write out this information
        if (agpKey1 != None):
            src1 = agpKey1[0]
            idx1 = agpKey1[1]
        else:
            src1 = "NA"
            idx1 = -1
        
        if (agpKey2 != None):
            src2 = agpKey2[0]
            idx2 = agpKey2[1]
        else:
            src2 = "NA"
            idx2 = -1

        # write out the record
        ofh.write(f"{contig1}\t{src1}\t{idx1}\t{edge1}\t{contig2}\t{src2}\t{idx2}\t{edge2}\t{weight}\n")

        node1 = get_node(contig1, END(edge1), src1, idx1)
        node2 = get_node(contig2, END(edge2), src2, idx2)

        # create the edges
        node1.add_edge(node2, weight)
        node2.add_edge(node1, weight)

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
