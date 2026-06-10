#!/bin/env python3

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

class Contig:
    def __init__(self, name: str):
        self.name   = name
        self.dir_5p = defaultdict(int) # how many reads support
        self.dir_3p = defaultdict(int) # either the 5' or 3' directions

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

def identify_separations() -> int:
    """now check the pairings that are separated"""

    global agp_map, contig_map, valid_pairs, reverse_names

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
                if (check):
                    name1 = contig.name
                    name2 = left_cntg
                    if (name1 in reverse_names):
                        name1 = reverse_names[name1]
                    if (name2 in reverse_names):
                        name2 = reverse_names[name2]
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
                if (check):
                    name1 = contig.name
                    name2 = right_cntg
                    if (name1 in reverse_names):
                        name1 = reverse_names[name1]
                    if (name2 in reverse_names):
                        name2 = reverse_names[name2]
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

    return 0

if __name__ == "__main__":
    main()
