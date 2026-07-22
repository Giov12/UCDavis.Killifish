#!/bin/env python3

import argparse
import os
import gzip
import sys
from enum import Enum

fasta      = ''
cfile      = ''
components = list()
records    = dict()

# create once
_DNA_COMPLEMENT = str.maketrans("ATCGatcgNn", "TAGCtagcNn")
class Orientation(Enum):
    reverse = 0
    forward = 1

class Record:
    __slots__ = ("id", "seq")
    def __init__(self, id_: str) -> None:
        self.id  = id_
        self.seq = ''
    
class Component:
    def __init__(self, id_: str) -> None: 
        self.id          = id_
        self.components = list()

    def add_component(self, seq_id: str, orientation: Orientation) -> int:
        self.components.append((seq_id, orientation))
        return 0
    
    def get_components(self) -> list[tuple[str, Orientation]]:
        return self.components

def load_components() -> int:
    """load the components"""

    global cfile, records, components

    fh = open(cfile, 'r')

    for line in fh:
        if (len(line) == 0 or line[0] == '#'):
            continue
        
        fields = line.strip().split()
        comp   = Component(fields[0])
        for entry in fields[2].split(','):
            idx  = entry.find('_') # find first occurance
            seq  = entry[:idx]
            oren = Orientation.forward if (entry[-1] == '+') else Orientation.reverse
            comp.add_component(seq, oren)
            records[seq] = None # place holder

        components.append(comp)

    fh.close()

    print(f"Loaded {len(records)} target IDs for {len(components)} components")

    return 0

def set_arguments() -> int:
    """get the fasta file and optional targets"""

    global fasta, cfile

    d = "Construct a fasta file based on the contigs in each component"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-f", "--fasta", help="fasta file containing the contigs", required=True, type=str, default='')
    parser.add_argument("-c", "--components", help="A file containing the components to create", type=str, default='', required=True)

    args  = parser.parse_args()
    fasta = args.fasta
    cfile = args.components
  
    assert os.path.isfile(fasta), f"Could not locate file {fasta}"
    assert os.path.isfile(cfile), f"Could not locate file {cfile}"
    
    return 0

def get_header_id(header: str) -> str:
    """return the first characters before any spacing"""

    id_ = ''
    idx = header.find(' ')

    if (idx == -1):
        id_ = header[1:]
    else:
        id_ = header[1:idx]

    return id_

def load_fasta() -> int:
    """load the sequence records into memory"""

    global fasta, records

    fh     = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, 'r')
    curRec = ''
    seq    = list()
    count  = 0 

    for line in fh:

        if (len(line) == 0 or line[0] == '#'):
            continue
        line = line.strip()
        if (line[0] == '>'):
            if (curRec == ''):
                curRec = get_header_id(line)
                if (curRec not in records):
                    curRec = ''
            else:
                rec             = Record(curRec)
                rec.seq         = ''.join(seq)
                records[curRec] = rec # replace
                count          += 1
                seq.clear()
                curRec = get_header_id(line)
                if (curRec not in records):
                    curRec = ''

        elif (curRec != ''):
            seq.append(line)

    fh.close()

    if (curRec != ''):
        rec             = Record(curRec)
        rec.seq         = ''.join(seq)
        records[curRec] = rec
        count          += 1
        seq.clear()

    fh.close()

    if (len(records) != count):
        msg = f"Failed to load all target sequences needed to create components"
        sys.exit(msg)

    print(f"Loaded {len(records)} records")

    return 0

def write_components() -> int:
    """write each component as a sequence record"""

    global records, components, _DNA_COMPLEMENT

    fh  = open("components.fa", 'w')
    gap = 'N' * 100

    for component in components:
        seq   = list()
        _id   = ">component_" + component.id + '\n'
        comps = component.get_components()
        for comp in comps:
            target = comp[0]
            oren   = comp[1]
            sequen = records[target].seq
            sequen = sequen if (oren == Orientation.forward) else sequen.translate(_DNA_COMPLEMENT)[::-1]
            seq.append(sequen)
        sequence = gap.join(seq)
        record   = '\n'.join(sequence[i:i+60] for i in range(0, len(sequence), 60))
        fh.write(_id)
        fh.write(record)
        fh.write('\n')

    fh.close()

    print(f"Wrote a total of {len(components)} records")

    return 0

def main() -> int:
    """Entry point to this small program"""

    # get arguments
    set_arguments()

    # get the components and sequence record ids
    load_components()

    # read in the sequence records
    load_fasta()

    # construct the components fasta file
    write_components()

    return 0

if __name__ == "__main__":
    main()
