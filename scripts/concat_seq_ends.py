#!/bin/env python3

import argparse
import os
import gzip
import textwrap

fasta = ''
seqs  = list()

# create once
rc_dict = {'A':'T',
           'T':'A',
           'C':'G',
           'G':'C',
           'N':'N',
           '-':'-',
           'K':'M', # G||T : C||A
           'M':'K', # C||A : G||T
           'W':'W', # A||T : T||A
           'S':'S', # G||C : C||G
           'R':'Y', # A||G : T||C
           'Y':'R', # T||C : A||G
           'B':'V', # C||G||T : A||C||G
           'V':'B', # A||C||G : C||G||T
           'D':'H', # A||G||T : T||C||A
           'H':'D'}

class Record:
    __slots__ = ("id", "seq")
    def __init__(self, id_: str):
        self.id  = id_
        self.seq = ''

    def get_ends(self, size: int) -> tuple[str, str]: 
        """get both ends of the sequence"""
        if (size > len(self.seq)):
            length = len(self.seq) // 2
        else:
            length = size
        first  = self.seq[:length]
        second = self.seq[-length:]
        return (first, second)

def set_arguments() -> int:
    """get the fasta file"""

    global fasta

    d = "Construct a fasta file by concatenating the ends of each pairwise combination of sequences"

    parser = argparse.ArgumentParser(description = d)
    parser.add_argument("-f", "--fasta", help="fasta file", required=True, type=str, default='')

    args  = parser.parse_args()
    fasta = args.fasta
  
    assert os.path.isfile(fasta), f"Could not locate file {fasta}"

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

    global fasta, seqs

    fh     = gzip.open(fasta, "rt") if fasta.endswith(".gz") else open(fasta, 'r')
    curRec = ''
    seq    = list()

    for line in fh:

        if (len(line) == 0 or line[0] == '#'):
            continue
        line = line.strip()
        if (line[0] == '>'):
            if (curRec == ''):
                curRec = get_header_id(line)
            else:
                rec = Record(curRec)
                rec.seq = ''.join(seq)
                seqs.append(rec)
                seq.clear()
                curRec = get_header_id(line)

        elif (curRec != ''):
            seq.append(line)

    fh.close()

    if (curRec != ''):
        rec = Record(curRec)
        rec.seq = ''.join(seq)
        seqs.append(rec)
        seq.clear()

    fh.close()

    print(f"Loaded {len(seqs)} records")

    return 0

def reverse_seq(seq: str) -> str:
    """return the reverse complement of the sequence"""
    
    global rc_dict

    tmp = list()

    for nuc in seq[::-1]:
        tmp.append(rc_dict[nuc.upper()])

    return ''.join(tmp)

def write_seqs() -> int:
    """write each pairwise combination in every orientation"""

    global seqs

    fh     = open("concatenated.ends.fa", 'w')
    size   = 2_000_000
    labels = ["5'", "3'", "5'_r", "3'_r"]
    total  = 0 

    for i in range(len(seqs) - 1):
        rec  = seqs[i]
        ends = rec.get_ends(size)
        rev1 = reverse_seq(ends[0])
        rev2 = reverse_seq(ends[1])
        entries = [ends[0], ends[1], rev1, rev2]
        for j in range(i + 1, len(seqs)):
            rec2     = seqs[j]
            ends2    = rec2.get_ends(size)
            rev1_2   = reverse_seq(ends2[0])
            rev2_2   = reverse_seq(ends2[1])
            entries2 = [ends2[0], ends2[1], rev1_2, rev2_2]
            for k, entry in enumerate(entries):
                for m, entry2 in enumerate(entries2):
                    new_entry = entry + entry2 + '\n'
                    lab1      = labels[k]
                    lab2      = labels[m]
                    header    = f">{rec.id}_{lab1}_{rec2.id}_{lab2}\n"
                    fh.write(header)
                    fh.write(textwrap.fill(new_entry, width=60))
                    total    += 1

    fh.close()

    print(f"Wrote a total of {total} records")

    return 0


def main() -> int:
    """Entry point to this small program"""

    # get arguments
    set_arguments()

    # read in the sequence records
    load_fasta()

    # write out the ends
    write_seqs()

    return 0

if __name__ == "__main__":
    main()
